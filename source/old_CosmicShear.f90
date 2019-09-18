! Likelihood code for weak lensing including astrophysical systematics
! Dataset: Kilo Degree Survey (KiDS)
! Written by Shahab Joudaki (2015), adapted for KiDS (2016)

    module CosmicShear
    use CosmologyTypes
    use CAMB, only : ComovingRadialDistance, AngularDiameterDistance, AngularDiameterDistance2, f_K, Hofz  !distance also in Mpc no h units
    use constants
    use Precision
    use likelihood
    use settings
    use Interpolation
    use MatrixUtils
    use omp_lib
    use Likelihood_Cosmology
    use Calculator_Cosmology
    use CosmoTheory
    use ModelParams
    implicit none
!    private

    TYPE, extends(TCosmoCalcLikelihood) :: CosmicShearLikelihood
        real(mcp), dimension(2,70) :: arraysjorig1,arraysjorig2,arraysjorig3,arraysjorig4
        real(mcp), dimension(2,70,4) :: arraysjfull
        real(mcp), allocatable, dimension(:) :: xipm !7 tom-bins, 7 ang-bins, 2 for +/- gives 28*7*2 = 392
        real(mcp), allocatable, dimension(:,:) :: covxipm, invcovxipm,covxipminv !nzbins*(nzbins+1)*nangbins
        real(mcp), dimension(9) :: thetacfhtini, thetaradcfhtini !nangbins
        real(mcp), dimension(58999) :: ellgentestarrini
        real(mcp), allocatable,dimension(:) :: maskelements
        real(mcp), dimension(9,58999) :: bes0arr,bes4arr,bes2arr
        integer :: size_cov,size_covmask,size_covmaskplanck,sizcov,sizcovpremask,klinesum,set_scenario
        logical :: use_morell, use_rombint
        logical :: use_nl
    contains

    procedure :: LogLike => CosmicShear_LnLike
    procedure :: ReadIni => CosmicShear_ReadIni

    END TYPE CosmicShearLikelihood

    logical :: use_CosmicShear = .false.
!    public CosmicShearLikelihood, CosmicShearLikelihood_Add

    contains

    subroutine CosmicShearLikelihood_Add(LikeList, Ini)
    use settings
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: Ini
    Type(CosmicShearLikelihood), pointer :: this
    Type(TSettingIni) :: DataSets

	if(.not. Ini%Read_Logical('use_CosmicShear',.false.)) return
    call Ini%TagValuesForName('CosmicShear_dataset', DataSets)
    allocate(this)
    this%LikelihoodType = 'CosmicShear'
    this%needs_nonlinear_pk = .true.
    this%needs_background_functions = .true.
    this%needs_powerspectra = .true.
    this%needs_exact_z = .true.
    this%speed = -1
    this%num_z = 37
    this%size_cov = 180 !nzbins*(nzbins+1)*nangbins --- 4*5*(9+9)
    if (this%use_nl == .false.) then
    !  this%needs_nonlinear_pk = .false.
     this%needs_weylpower = .true. !NHmod
    endif
    call this%ReadDatasetFile(DataSets%Value(1))
!    this%LikelihoodType = 'CosmicShear'
    this%tag = DataSets%Name(1)
    call this%loadParamNames(trim(DataDir)//'CosmicShear.paramnames')
    call LikeList%Add(this)
    if (Feedback>1) write(*,*) 'Imported CosmicShear info'

    end subroutine CosmicShearLikelihood_Add

    subroutine CosmicShear_ReadIni(this, Ini)
    use MatrixUtils
    use settings
    use omp_lib
    external DGETRF, DGETRI
    class(CosmicShearLikelihood) this
    class(TSettingIni) :: Ini
    integer, allocatable, dimension(:) :: ipiv
    real(mcp), allocatable, dimension(:) :: work_inverse
    real(mcp), allocatable, dimension(:,:) :: xipmtemp,covxipmtemp,covxipmtempinv,masktemp
    integer :: info,lwork,zinit,hay,iii,jjj,it,jt,kt,setscenario,sizcovish,sizcovishsq,sizcovishpremask

    this%set_scenario = Ini%Read_Int('set_scenario')
    this%use_morell = Ini%Read_Logical('use_morell',.false.)
    this%use_nl = Ini%Read_Logical('use_nl', .false.)
    this%use_rombint = Ini%Read_Logical('use_rombint',.false.)
    print *, 'this%use_morell, this%use_rombint', this%use_morell, this%use_rombint
    print *, 'this%use_nl', this%use_nl !NHmod
!    setscenario = 1
    setscenario = this%set_scenario
    print *, 'setscenario', setscenario

    if (this%use_nl) then
      this%size_covmask = 130
    else
      this%size_covmask = 30
    end if
    this%size_covmaskplanck = 56

    if(setscenario == 0) then !no masking
        sizcovish = this%size_cov !180 elements post-masking
        sizcovishpremask = 180 !elements pre-masking, same number as post-masking, as this is the non-masking xi_+/- scenario
    end if
    if(setscenario == 1) then !fiducial masking advocated in Hildebrandt et al (2016) and Joudaki et al (2016)
        sizcovish = this%size_covmask !130 elements post-masking
        sizcovishpremask = 180 !elements pre-masking
    end if

    sizcovishsq = sizcovish**2
    this%sizcov = sizcovish
    this%sizcovpremask = sizcovishpremask

    allocate(this%xipm(sizcovish))
    allocate(this%covxipm(sizcovish,sizcovish))
    allocate(this%invcovxipm(sizcovish,sizcovish))
    allocate(this%covxipminv(sizcovish,sizcovish))
    allocate(this%maskelements(sizcovishpremask))
    allocate(masktemp(2,sizcovishpremask))
    allocate(xipmtemp(4,sizcovish))
    allocate(covxipmtemp(3,sizcovishsq))

    !!!Reading in source distributions
    open(7,file='data/KIDS/nz_z1_kids_binned_hendrik.dat') !bootstrap z
    READ (7,*) this%arraysjorig1
    close(7)
    open(7,file='data/KIDS/nz_z2_kids_binned_hendrik.dat') !bootstrap z
    READ (7,*) this%arraysjorig2
    close(7)
    open(7,file='data/KIDS/nz_z3_kids_binned_hendrik.dat') !bootstrap z
    READ (7,*) this%arraysjorig3
    close(7)
    open(7,file='data/KIDS/nz_z4_kids_binned_hendrik.dat') !bootstrap z
    READ (7,*) this%arraysjorig4
    close(7)

    this%arraysjorig1(2,:) = this%arraysjorig1(2,:)/(sum(this%arraysjorig1(2,:))*0.05d0)
    this%arraysjorig2(2,:) = this%arraysjorig2(2,:)/(sum(this%arraysjorig2(2,:))*0.05d0)
    this%arraysjorig3(2,:) = this%arraysjorig3(2,:)/(sum(this%arraysjorig3(2,:))*0.05d0)
    this%arraysjorig4(2,:) = this%arraysjorig4(2,:)/(sum(this%arraysjorig4(2,:))*0.05d0)
    this%arraysjfull(:,:,1) = this%arraysjorig1
    this%arraysjfull(:,:,2) = this%arraysjorig2
    this%arraysjfull(:,:,3) = this%arraysjorig3
    this%arraysjfull(:,:,4) = this%arraysjorig4

    !!!Reading in measurements
!    if(setscenario == 0) then !do not use this option
!        open(7,file='data/KIDS/xipm_filename.dat')
!    end if
    if(setscenario == 1) then
      if (this%use_nl) then
        open(7,file='data/KIDS/xipmcut_kids_blind1.dat') !USED
      else
        open(7,file='data/KIDS/xipmcut_kids_blind1_largecut.dat') !USED
      end if
!        open(7,file='data/KIDS/xipmcut_kids_blind2.dat')
!        open(7,file='data/KIDS/xipmcut_kids_blind3.dat')
    end if
    READ (7,*) xipmtemp
    close(7)
    this%xipm = xipmtemp(2,:)

    !!!Masking file
    if(setscenario == 0) then
        open(7,file='data/KIDS/xipm_kids4tom_selectsjnocut.dat') !USED
    end if
    if(setscenario == 1) then
      if (this%use_nl) then
        open(7,file='data/KIDS/xipm_kids4tom_selectsj.dat') !USED
      else
        open(7,file='data/KIDS/xipm_kids4tom_selectsj_largecut.dat') ! SP: USED
      end if
    end if
    READ (7,*) masktemp
    close(7)
    this%maskelements = masktemp(2,:)

    !!!Reading in covariance
    !if(setscenario == 0) then !do not use this option
     !   open(7,file='xipmcov_filename.dat')
    !end if
    if(setscenario == 1) then
      if (this%use_nl) then
        open(7,file='data/KIDS/xipmcutcov_kids_analytic_inc_m_blind1.dat') !USED
      else
        open(7,file='data/KIDS/xipmcutcov_kids_analytic_inc_m_blind1_largecut.dat') ! SP: USED
      end if
 !       open(7,file='xipmcutcov_kids_analytic_inc_m_blind2.dat')
 !       open(7,file='xipmcutcov_kids_analytic_inc_m_blind3.dat')
 !       open(7,file='xipmcutcov_kids_regcomb_blind1_nbodysj.dat')
 !       open(7,file='xipmcutcov_kids_regcomb_blind2_nbodysj.dat')
 !       open(7,file='xipmcutcov_kids_regcomb_blind3_nbodysj.dat')
    end if
    READ (7,*) covxipmtemp
    close(7)

    kt = 0
    do it=1,sizcovish
        do jt=1,sizcovish
            kt = kt+1
            this%covxipm(it,jt) = covxipmtemp(3,kt)
        end do
    end do

    !converted from arcmin to degrees
    this%thetacfhtini = (/ 0.71336d0, 1.45210d0, 2.95582d0, 6.01675d0, 12.24745d0, 24.93039d0, 50.74726d0, 103.29898d0, 210.27107d0 /)/60.0d0 !athena angular scales (deg) --- 9 ang bins
    this%thetaradcfhtini = this%thetacfhtini*3.14159265359d0/180.0d0 !converted to rad
    this%ellgentestarrini = (/ (hay,hay=2,59000,1) /)

print *, 'this%thetacfhtini',this%thetacfhtini

    !compute Bessel functions
    do iii=1,9
        do jjj=1,58999
            this%bes0arr(iii,jjj) = bessel_j0(this%ellgentestarrini(jjj)*this%thetaradcfhtini(iii))
            this%bes4arr(iii,jjj) = bessel_jn(4,this%ellgentestarrini(jjj)*this%thetaradcfhtini(iii))
            this%bes2arr(iii,jjj) = bessel_jn(2,this%ellgentestarrini(jjj)*this%thetaradcfhtini(iii))
        end do
    end do

    allocate(this%exact_z(this%num_z)) !allocate array for matter power spectrum redshifts
    !assign redshifts for matter power spectrum
    this%exact_z(1) = 0.0d0
    this%exact_z(2) = 0.025d0
    do zinit=3,37
        this%exact_z(zinit) = this%exact_z(zinit-1) + 0.154d0*this%exact_z(zinit-1)
    end do
    this%exact_z(37) = 3.474999d0

    print *, 'shahab this%exact_z', this%exact_z

    !calculate inverse of covariance
    allocate(work_inverse(sizcovish*sizcovish))         !for inverse calculation
    allocate(ipiv(sizcovish))                           !for inverse calculation
    this%invcovxipm = this%covxipm
    LWORK = sizcovish**2
    call DGETRF(sizcovish,sizcovish,this%invcovxipm,sizcovish,ipiv,info)         !LU decomposition
    call DGETRI(sizcovish,this%invcovxipm,sizcovish,ipiv,work_inverse,lwork,info)    !inverse from LU decompostion
    ipiv(:) =  0
    work_inverse(:) = 0
    if (info .ne. 0) stop 'Problem with the matrix inverse calculation'

!    if(setscenario == 0) then
!        this%invcovxipm = (this%invcovxipm)*(930.0-180.0-2.0)/(930.0-1.0) !Hartlap correction factor
!        this%invcovxipm = this%invcovxipm !Sellentin-Heavens or Analytic Covariance instead
!    end if
    if(setscenario == 1) then
!        this%invcovxipm = (this%invcovxipm)*(930.0-130.0-2.0)/(930.0-1.0) !Hartlap correction factor
        this%invcovxipm = this%invcovxipm !Sellentin-Heavens or Analytic Covariance instead
    end if

    DEALLOCATE(work_inverse)
    DEALLOCATE(ipiv)
    DEALLOCATE(xipmtemp)
    DEALLOCATE(covxipmtemp)

    this%klinesum = 0

    end subroutine CosmicShear_ReadIni

    function CosmicShear_LnLike(this,CMB,Theory,DataParams)
    implicit none
    class(CosmicShearLikelihood) :: this
    class(CMBParams) :: CMB
    class(TCosmoTheoryPredictions), target :: Theory
    type my_type
        type(TCosmoTheoryPredictions), allocatable :: theoryomp
        real(mcp) :: pisjomp,h0omp,homp,omdmomp,ombomp,omkomp,distlensflatomp,mumin2vomp,mumaxomp,tol_erroromp,ampiayesomp,redziayesomp,lumiayesomp,nellbinsomp
        real(mcp), allocatable, dimension(:) :: ellarromp, aphotarromp,exact_zmodomp,lumarromp
        integer, allocatable, dimension(:) :: exact_z_index,momp,m1arromp,m2arromp
        integer :: wchooseomp,ellgenomp,bnumomp,kline,m1omp,m2omp,binwomp,wtrapmaxomp
        real(mcp) :: mumax1omp,mumax2omp,mumaxlensomp,mumaxwomp
        real(mcp), allocatable, dimension(:,:) :: weightarromp,weightpsarromp,weightnoarromp
        type(SCubicSpline), allocatable, dimension(:) :: psourcetypeomp !NH compiler having trouble here
        type(SCubicSpline), allocatable :: gcubicomp
        Type(TCubicSpline), allocatable :: clinterptypeomp,clinterptypeiiomp,clinterptypegiomp !keep this as tcubicspline
    end type my_type
    real(mcp) rombint
    real(mcp) rombint_obj
    external rombint
    external rombint_obj
    type(my_type) obj
    real(mcp) :: CosmicShear_LnLike,ampiayes,redziayes,lumiayes
    real(mcp) :: tol_error
    real(mcp) :: a1phot,a2phot,a3phot,a4phot
    real(mcp) :: mumax,mumin,mumin2v, pisj,likesj,chisqsj,chisqsjsellentin,start,finish,mumax1,mumax2,mumaxw,morellfac,ckmz
    real(mcp), allocatable, dimension(:,:,:) :: nonlinspec
    real(mcp), allocatable, dimension(:) :: xiplus,xiplusii,xiplusgi
    real(mcp), allocatable, dimension(:) :: ximinus,ximinusii,ximinusgi
    real(mcp), allocatable, dimension(:) :: xiplusminus, xiplusminusgg,xiplusminusii,xiplusminusgi,dataminustheory,dataminustheoryhihi,dataminustheoryhoho,finaltheory
    real(mcp), allocatable, dimension(:,:) :: weightarr,weightpsarr,weightnoarr
    real(mcp), allocatable, dimension(:) :: m1arr,m2arr,binwarr
    real(mcp), allocatable, dimension(:) :: ellarr,trapzarr,trapzarrii,trapzarrgi
    real(mcp), allocatable, dimension(:) :: mumax1arr,mumax2arr, aphotarr,garr,rcslum,cfhtlum,lumarr,mnomp
    real(mcp), allocatable, dimension(:,:) :: clarr,clarrii,clarrgi
    real(mcp), allocatable, dimension(:) :: clarrbig,ellgentestarr,clarrbigii,clarrbiggi
    integer :: hoj,bnum,gigi,how,howin
    integer :: loopee,yeye,ellgen,intnum,bin,biin,wchoose, ho, m1, m2, mint,ttt,wtrapmax,wttt,binw
    integer :: nangbins, nzbins, nellbins,  nell, nzp
    real(mcp) :: ellgentest
    logical :: use_cfht = .true.
    logical :: use_rombintz = .true. !if true use Romberg integration instead of trapezoidal for outer integral
!    logical :: use_rombintz = .false.
!    logical :: use_morellz = .true. !if true use 101 ell values, otherwise 31 ell values
     logical :: use_morellz = .false. !if true use 101 ell values, otherwise 31 ell values
    real(mcp) :: DataParams(:)
    ampiayes = DataParams(1)
    redziayes = DataParams(2)
    lumiayes = DataParams(3) !not used for KiDS
    a1phot = DataParams(4) !not used for KiDS
    a2phot = DataParams(5) !not used for KiDS
    a3phot = DataParams(6) !not	used for KiDS
    a4phot = DataParams(7) !not	used for KiDS

    use_morellz = this%use_morell
    use_rombintz = this%use_rombint
    nangbins = 9
    nzbins = 4
    if(use_morellz == .false.) nellbins = 31
    if(use_morellz == .true.) nellbins = 101
    obj%nellbinsomp = nellbins
    nell = 58999
    nzp = 70
    wtrapmax = 37
    obj%wtrapmaxomp = wtrapmax

    allocate(obj%theoryomp)
    allocate(obj%ellarromp(nellbins))
    allocate(obj%exact_z_index(this%num_z))
    allocate(obj%momp(2))
    allocate(obj%psourcetypeomp(nzbins))
    allocate(obj%aphotarromp(nzbins))
    allocate(obj%lumarromp(nzbins))
    allocate(obj%exact_zmodomp(wtrapmax))
    allocate(obj%weightarromp(nzbins,wtrapmax-1))
    allocate(obj%weightpsarromp(nzbins,wtrapmax-1))
    allocate(obj%m1arromp(nzbins*(nzbins+1)/2))
    allocate(obj%m2arromp(nzbins*(nzbins+1)/2))
    allocate(obj%weightnoarromp(nellbins,wtrapmax-1))
    allocate(obj%clinterptypeomp)
    allocate(obj%clinterptypeiiomp)
    allocate(obj%clinterptypegiomp)
    allocate(obj%gcubicomp)

    allocate(clarrbig(nell))
    allocate(clarrbiggi(nell))
    allocate(clarrbigii(nell))
    allocate(ellgentestarr(nell))
    allocate(xiplusminus(nzbins*(nzbins+1)/2*nangbins*2))
    allocate(xiplusminusgg(nzbins*(nzbins+1)/2*nangbins*2))
    allocate(xiplusminusii(nzbins*(nzbins+1)/2*nangbins*2))
    allocate(xiplusminusgi(nzbins*(nzbins+1)/2*nangbins*2))
    allocate(dataminustheory(this%sizcov))
    allocate(dataminustheoryhihi(this%sizcov))
    allocate(dataminustheoryhoho(this%sizcov))
    allocate(finaltheory(this%sizcovpremask))
    allocate(clarr(nellbins,nzbins*(nzbins+1)/2))
    allocate(clarrii(nellbins,nzbins*(nzbins+1)/2))
    allocate(clarrgi(nellbins,nzbins*(nzbins+1)/2))
    allocate(ellarr(nellbins))
    allocate(xiplus(nangbins))
    allocate(xiplusii(nangbins))
    allocate(xiplusgi(nangbins))
    allocate(ximinus(nangbins))
    allocate(ximinusii(nangbins))
    allocate(ximinusgi(nangbins))
    allocate(mumax1arr(nzbins*(nzbins+1)/2))
    allocate(mumax2arr(nzbins*(nzbins+1)/2))
    allocate(aphotarr(nzbins))
    allocate(rcslum(nzbins))
    allocate(cfhtlum(nzbins))
    allocate(lumarr(nzbins))
    allocate(m1arr(nzbins*(nzbins+1)/2))
    allocate(m2arr(nzbins*(nzbins+1)/2))
    allocate(mnomp(2))
    allocate(weightarr(nzbins,wtrapmax-1))
    allocate(weightpsarr(nzbins,wtrapmax-1))
    allocate(binwarr(nzbins))
    allocate(trapzarr(wtrapmax-2))
    allocate(trapzarrii(wtrapmax-2))
    allocate(trapzarrgi(wtrapmax-2))
    allocate(weightnoarr(nellbins,wtrapmax-1))
    allocate(garr(Theory%MPK%ny))

    aphotarr(1) = a1phot
    aphotarr(2) = a2phot
    aphotarr(3) = a3phot
    aphotarr(4) = a4phot
    obj%aphotarromp(1) = a1phot
    obj%aphotarromp(2) = a2phot
    obj%aphotarromp(3) = a3phot
    obj%aphotarromp(4) = a4phot

    start  = OMP_get_wtime()
    cfhtlum = (/ 0.0174d0, 0.0694d0, 0.1518d0, 0.2182d0 /) !not used for KiDS
    if(use_cfht == .true.) lumarr = cfhtlum
    obj%lumarromp = lumarr

    obj%exact_zmodomp = this%exact_z
    tol_error = 0.005d0 !Numerical integration relative error, for Romberg
    pisj = 3.14159265359d0
    obj%pisjomp = 3.14159265359d0
    obj%h0omp = CMB%H0
    obj%homp = CMB%h
    obj%omdmomp = CMB%omdm
    obj%ombomp = CMB%omb
    obj%omkomp = CMB%omk
    obj%theoryomp = Theory


    mumin=0.0d0
    mumin2v=0.025d0 !min integration redshift
    mumax=3.474999d0 !max integration redshift

    !Generate array containing l-values where C(l) is evaluated
    if(use_morellz == .false.) morellfac = 0.5
    if(use_morellz == .true.) morellfac = 0.1
    do ellgen=1,nellbins
        if(ellgen < 10) ellarr(ellgen) = ellgen+1
        if(ellgen > 10 .or. ellgen == 10) ellarr(ellgen) = ellarr(ellgen-1)+morellfac*ellarr(ellgen-1)
    end do
    ellarr = int(ellarr)
    obj%ellarromp = ellarr
    obj%tol_erroromp = tol_error
    obj%mumin2vomp = mumin2v
    obj%mumaxomp = mumax

    obj%ampiayesomp = ampiayes
    obj%redziayesomp = redziayes
    obj%lumiayesomp = lumiayes
    obj%exact_z_index = this%exact_z_index
    wchoose = 1
    obj%wchooseomp = 1

    !interpolation for source redshift distributions
    do ho=1,nzbins
        call obj%psourcetypeomp(ho)%Init(this%arraysjfull(1,:,ho),this%arraysjfull(2,:,ho),nzp)
    end do

    obj%kline = 0

    !tomographic ordering
    m1arr = (/ 1, 1, 1, 1, 2, 2, 2, 3, 3, 4 /)
    m2arr = (/ 1, 2, 3, 4, 2, 3, 4, 3, 4, 4 /)

    !compute integration range accounting for photo-z error
    mint = 0
    do m1=1,nzbins
        do m2=1,nzbins
            if(m2 >= m1) then
                mint = mint + 1
                mumax1arr(mint) = mumax + aphotarr(m1arr(mint))
                mumax2arr(mint) = mumax + aphotarr(m2arr(mint))
            end if
        end do
    end do

    !requires the redshift in CosmologyTypes.f90 file to be set to z>3.5 in case of photo-z varying
    do gigi=1,Theory%MPK%ny
        garr(gigi) = exp(rombint(sjgrowtha,1.0d0/(1.0d0+Theory%MPK%y(gigi)),1.0d0,tol_error))
    end do
    call obj%gcubicomp%Init(Theory%MPK%y,garr,Theory%MPK%ny)


    if(use_rombintz == .true.) then !Romberg Integration for the outer integral

        !for nz tomographic bins, compute nz*(nz+1)/2 C(l)-combinations
        do bin = 1,nzbins*(nzbins+1)/2
            mint = 0
            do m1=1,nzbins
                do m2=1,nzbins
                    if(m2 >= m1) then
                        mint = mint + 1
                        if(mint == bin) then
                            obj%m1omp=m1arr(mint)
                            obj%m2omp=m2arr(mint)
                            obj%momp(1)=m1arr(mint)
                            obj%momp(2)=m2arr(mint)
                        end if
                    end if
                end do
            end do
            mnomp = obj%momp
            bnum = bin
            obj%bnumomp = bin
            mumax1 = mumax1arr(bin)
            mumax2 = mumax2arr(bin)
            obj%mumax1omp = mumax1arr(bin)
            obj%mumax2omp = mumax2arr(bin)
            obj%mumaxlensomp = min(mumax1,mumax2)

            !!!##############COMPUTE LENSING POWER SPECTRA (GG, GI, II, Gg) AT DISTINCT L-VALUES. PARALLELIZATION EMPLOYED############
            !$OMP PARALLEL DO FIRSTPRIVATE(ellgen,obj)
                do ellgen=1,nellbins
                    obj%ellgenomp = ellgen
                    clarr(ellgen,obj%bnumomp) = rombint_obj(obj,sjclsobjsf,1.0d0/(1.0d0+obj%mumaxlensomp),1.0d0/(1.0d0+obj%mumin2vomp),obj%tol_erroromp) !!good!!
                    clarrii(ellgen,obj%bnumomp) = rombint_obj(obj,sjclsiiobjsf,1.0d0/(1.0d0+obj%mumaxlensomp),1.0d0/(1.0d0+obj%mumin2vomp),obj%tol_erroromp) !good!!
                    clarrgi(ellgen,obj%bnumomp) = rombint_obj(obj,sjclsgiobjsf,1.0d0/(1.0d0+obj%mumaxlensomp),1.0d0/(1.0d0+obj%mumin2vomp),obj%tol_erroromp) !!good!!
                end do
            !$OMP END PARALLEL DO
        end do
    end if

    !Trapezoid Integration for the outer integral!!! for the case use_rombintz=F !!! ------ do not use this option when varying photo-z params, in that case set use_rombintz = T
    if(use_rombintz == .false.) then

        binwarr = (/ 1, 5, 8, 10 /) !4 tomographic bins, m1arr(i) = m2arr(i) for these indices --- need to change this array for different number of tomographic bins
        weightarr(:,:) = 0.0d0
        weightpsarr(:,:) = 0.0d0
        weightnoarr(:,:) = 0.0d0

        do binw = 1,nzbins
            mumaxw = mumax1arr(binwarr(binw))
            obj%mumaxwomp = mumaxw
            bnum = binwarr(binw)
            obj%bnumomp = bnum
            obj%m1omp=m1arr(binwarr(binw))
            obj%m2omp=m2arr(binwarr(binw))
            obj%momp(1)=m1arr(binwarr(binw))
            obj%momp(2)=m2arr(binwarr(binw))
            mnomp = obj%momp
            obj%binwomp = binw

            !$OMP PARALLEL DO FIRSTPRIVATE(wttt,obj)
                do wttt=2,obj%wtrapmaxomp
                    weightarr(obj%binwomp,wttt-1) = sjclsobjonlyweight(obj,obj%exact_zmodomp(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                    weightpsarr(obj%binwomp,wttt-1) = sjclsiiobjonlyweight(obj,obj%exact_zmodomp(wttt))   !have to make sure mumaxw > obj%exact_zmodomp(wttt)
                end do
            !$OMP END PARALLEL DO
        end do

        do ellgen=1,nellbins
            obj%ellgenomp = ellgen
            !$OMP PARALLEL DO FIRSTPRIVATE(wttt,obj)
                do wttt=2,obj%wtrapmaxomp
                    weightnoarr(obj%ellgenomp,wttt-1) = sjclsobjnoweight(obj,obj%exact_zmodomp(wttt))   !have to make sure mumaxw > this%exact_z(wttt)
                end do
            !$OMP END PARALLEL DO
        end do

        !for nz tomographic bins, compute nz*(nz+1)/2 C(l)-combinations
        do bin = 1,nzbins*(nzbins+1)/2
            mint = 0
            do m1=1,nzbins
                do m2=1,nzbins
                    if(m2 >= m1) then
                        mint = mint + 1
                        if(mint == bin) then
                            obj%m1omp=m1arr(mint)
                            obj%m2omp=m2arr(mint)
                            obj%momp(1)=m1arr(mint)
                            obj%momp(2)=m2arr(mint)
                        end if
                    end if
                end do
            end do
            mnomp = obj%momp
            bnum = bin
            obj%bnumomp = bin
            mumax1 = mumax1arr(bin)
            mumax2 = mumax2arr(bin)
            obj%mumax1omp = mumax1arr(bin)
            obj%mumax2omp = mumax2arr(bin)
            obj%mumaxlensomp = min(mumax1,mumax2)
            obj%weightarromp = weightarr
            obj%weightpsarromp = weightpsarr
            obj%weightnoarromp = weightnoarr
            obj%m1arromp = m1arr
            obj%m2arromp = m2arr

            !!!##############COMPUTE LENSING POWER SPECTRA (GG, GI, II, Gg) AT DISTINCT L-VALUES. PARALLELIZATION EMPLOYED############
            do ellgen=1,obj%nellbinsomp
                obj%ellgenomp = ellgen
                !$OMP PARALLEL DO FIRSTPRIVATE(ttt,obj)
                    do ttt=2,36
                        trapzarr(ttt-1) = 0.5d0*(obj%weightnoarromp(obj%ellgenomp,ttt-1)*obj%weightarromp(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightarromp(obj%m2arromp(obj%bnumomp),ttt-1) + obj%weightnoarromp(obj%ellgenomp,ttt)*obj%weightarromp(obj%m1arromp(obj%bnumomp),ttt)*obj%weightarromp(obj%m2arromp(obj%bnumomp),ttt))*(obj%exact_zmodomp(ttt+1)-obj%exact_zmodomp(ttt))
                        trapzarrii(ttt-1) = 0.5d0*(obj%weightnoarromp(obj%ellgenomp,ttt-1)*obj%weightpsarromp(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightpsarromp(obj%m2arromp(obj%bnumomp),ttt-1) + obj%weightnoarromp(obj%ellgenomp,ttt)*obj%weightpsarromp(obj%m1arromp(obj%bnumomp),ttt)*obj%weightpsarromp(obj%m2arromp(obj%bnumomp),ttt))*(obj%exact_zmodomp(ttt+1)-obj%exact_zmodomp(ttt))
                        trapzarrgi(ttt-1) = 0.5d0*(obj%weightnoarromp(obj%ellgenomp,ttt-1)*(obj%weightarromp(obj%m1arromp(obj%bnumomp),ttt-1)*obj%weightpsarromp(obj%m2arromp(obj%bnumomp),ttt-1)+obj%weightarromp(obj%m2arromp(obj%bnumomp),ttt-1)*obj%weightpsarromp(obj%m1arromp(obj%bnumomp),ttt-1)) + obj%weightnoarromp(obj%ellgenomp,ttt)*(obj%weightarromp(obj%m1arromp(obj%bnumomp),ttt)*obj%weightpsarromp(obj%m2arromp(obj%bnumomp),ttt)+obj%weightarromp(obj%m2arromp(obj%bnumomp),ttt)*obj%weightpsarromp(obj%m1arromp(obj%bnumomp),ttt)))*(obj%exact_zmodomp(ttt+1)-obj%exact_zmodomp(ttt))
                    end do
                !$OMP END PARALLEL DO
                clarr(ellgen,obj%bnumomp)=sum(trapzarr)
                clarrii(ellgen,obj%bnumomp)=sum(trapzarrii)
                clarrgi(ellgen,obj%bnumomp)=sum(trapzarrgi)
            end do
        end do
    end if

    !Interpolation for C(l)
    ellgentestarr = this%ellgentestarrini
    intnum=nellbins
    do biin = 1,nzbins*(nzbins+1)/2
        call obj%clinterptypeomp%Init(ellarr,clarr(:,biin),intnum)
        call obj%clinterptypeiiomp%Init(ellarr,clarrii(:,biin),intnum)
        call obj%clinterptypegiomp%Init(ellarr,clarrgi(:,biin),intnum)
        ellgentest = 1.0d0
        !$OMP PARALLEL DO FIRSTPRIVATE(obj,ellgen)
            do ellgen=1,58999
                clarrbig(ellgen) = obj%clinterptypeomp%Value(ellgen+1.0d0)
                clarrbigii(ellgen) = obj%clinterptypeiiomp%Value(ellgen+1.0d0)
                clarrbiggi(ellgen) = obj%clinterptypegiomp%Value(ellgen+1.0d0)
            end do
        !$OMP END PARALLEL DO

        xiplus(:) = 0.0d0
        ximinus(:) = 0.0d0
        xiplusii(:) = 0.0d0
        ximinusii(:) = 0.0d0
        xiplusgi(:) = 0.0d0
        ximinusgi(:) = 0.0d0
        hoj=0
        do yeye=1,nangbins
            do loopee=1,nell
                xiplus(yeye) = xiplus(yeye) + 1.0d0/pisj/2.0d0*clarrbig(loopee)*ellgentestarr(loopee)*this%bes0arr(yeye,loopee)
                ximinus(yeye) = ximinus(yeye) + 1.0d0/pisj/2.0d0*clarrbig(loopee)*ellgentestarr(loopee)*this%bes4arr(yeye,loopee)
                xiplusii(yeye) = xiplusii(yeye) + 1.0d0/pisj/2.0d0*clarrbigii(loopee)*ellgentestarr(loopee)*this%bes0arr(yeye,loopee)
                ximinusii(yeye) = ximinusii(yeye) + 1.0d0/pisj/2.0d0*clarrbigii(loopee)*ellgentestarr(loopee)*this%bes4arr(yeye,loopee)
                xiplusgi(yeye) = xiplusgi(yeye) + 1.0d0/pisj/2.0d0*clarrbiggi(loopee)*ellgentestarr(loopee)*this%bes0arr(yeye,loopee)
                ximinusgi(yeye) = ximinusgi(yeye) + 1.0d0/pisj/2.0d0*clarrbiggi(loopee)*ellgentestarr(loopee)*this%bes4arr(yeye,loopee)
            end do
        end do

        !old Catherine's scheme
        xiplusminus(1+(biin-1)*nangbins*2:nangbins+(biin-1)*nangbins*2) = xiplus + xiplusii + xiplusgi !full xi+ accounting for intrinsic alignments
        xiplusminus((nangbins+1)+(biin-1)*nangbins*2:nangbins*2+(biin-1)*nangbins*2) = ximinus + ximinusii + ximinusgi !full xi- accounting for intrinsic alignments
        xiplusminusgg(1+(biin-1)*nangbins*2:nangbins+(biin-1)*nangbins*2) = xiplus
        xiplusminusgg((nangbins+1)+(biin-1)*nangbins*2:nangbins*2+(biin-1)*nangbins*2) = ximinus
        xiplusminusii(1+(biin-1)*nangbins*2:nangbins+(biin-1)*nangbins*2) = xiplusii
        xiplusminusii((nangbins+1)+(biin-1)*nangbins*2:nangbins*2+(biin-1)*nangbins*2) = ximinusii
        xiplusminusgi(1+(biin-1)*nangbins*2:nangbins+(biin-1)*nangbins*2) = xiplusgi
        xiplusminusgi((nangbins+1)+(biin-1)*nangbins*2:nangbins*2+(biin-1)*nangbins*2) = ximinusgi

    end do

    this%klinesum = this%klinesum + obj%kline

    finaltheory(1:size(xiplusminus)) = xiplusminus
    howin = 1
    do how=1,this%sizcovpremask
        if(this%maskelements(how) == 1) then
            dataminustheory(howin) = this%xipm(howin) - finaltheory(how)
            dataminustheoryhoho(howin) = finaltheory(how)
            dataminustheoryhihi(howin) = this%xipm(howin)
            howin = howin + 1
        end if
    end do

    likesj = (Matrix_QuadForm(this%invcovxipm,dataminustheory))/2.0d0 !compute likelihood
    chisqsj = likesj*2.0
    chisqsjsellentin = 930.0d0*dlog(1.0d0 + likesj*2.0d0/(930.0d0-1.0d0)) !KiDS
    CosmicShear_LnLike = likesj !!Analytic or Hartlap
!!!    CosmicShear_LnLike = 930.0d0/2.0d0*dlog(1.0d0 + likesj*2.0d0/(930.0d0-1.0d0)) !!KiDS: Sellentin-Heavens (only with N-body covariance)
!!!    CosmicShear_LnLike = 0.0d0

    finish  = OMP_get_wtime()
    print *, 'finish-start, chisqsj,chisqsjsellentin,tol_error', finish-start, chisqsj,chisqsjsellentin,tol_error

    deallocate(obj%theoryomp)
    deallocate(obj%ellarromp)
    deallocate(obj%exact_z_index)
    deallocate(obj%psourcetypeomp)
    deallocate(obj%momp)
    deallocate(obj%aphotarromp)
    deallocate(obj%weightarromp)
    deallocate(obj%weightpsarromp)
    deallocate(obj%weightnoarromp)
    deallocate(obj%m1arromp)
    deallocate(obj%m2arromp)
    deallocate(obj%clinterptypeomp)
    deallocate(obj%clinterptypeiiomp)
    deallocate(obj%clinterptypegiomp)
    deallocate(obj%exact_zmodomp)
    deallocate(obj%gcubicomp)
    deallocate(obj%lumarromp)

    deallocate(clarrbig)
    deallocate(clarrbiggi)
    deallocate(clarrbigii)
    deallocate(ellgentestarr)
    deallocate(xiplusminus)
    deallocate(xiplusminusgg)
    deallocate(xiplusminusii)
    deallocate(xiplusminusgi)
    deallocate(dataminustheory)
    deallocate(dataminustheoryhihi)
    deallocate(dataminustheoryhoho)
    deallocate(finaltheory)
    deallocate(clarr)
    deallocate(clarrii)
    deallocate(clarrgi)
    deallocate(ellarr)
    deallocate(xiplus)
    deallocate(xiplusii)
    deallocate(xiplusgi)
    deallocate(ximinus)
    deallocate(ximinusii)
    deallocate(ximinusgi)
    deallocate(mumax1arr)
    deallocate(mumax2arr)
    deallocate(aphotarr)
    deallocate(m1arr)
    deallocate(m2arr)
    deallocate(mnomp)
    deallocate(rcslum)
    deallocate(cfhtlum)
    deallocate(lumarr)
    deallocate(weightarr)
    deallocate(weightpsarr)
    deallocate(weightnoarr)
    deallocate(binwarr)
    deallocate(trapzarr)
    deallocate(trapzarrii)
    deallocate(trapzarrgi)
    deallocate(garr)

    contains


    function sjgrowtha(reda) !integrate to obtain growth function D(0)/D(z)
        REAL(mcp) :: reda,sjgrowtha

        sjgrowtha = (1.0d0/reda)*Theory%growth_z%Value(1.0d0/reda-1.0d0)/Theory%sigma8_z%Value(1.0d0/reda-1.0d0)
    END function sjgrowtha

    !non-weight component for trapezoid integration
    function sjclsobjnoweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclsobjnoweight,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz
        distz = f_K(distz) !with curvature
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)

        kminsj = 1.0d-5
        kmaxsj = 100.0d0 !kmaxsjhoho
        deltafinal = 0.0d0
        if(kval < kminsj) obj%kline = obj%kline + 1
        if (this%use_nl== .true.) then
         if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK%PowerAt(kval,zlens)
        else
         if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK_WEYL%PowerAt(kval,zlens) !NHmod
        endif
        hubblez = Hofz(zlens)
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        if (this%use_nl ==.true.) then
         sjclsobjnoweight = 1.0d0/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0
        else
         sjclsobjnoweight = 1.0d0/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0/(9.0d0/4.0d0*(obj%h0omp/ckms)**4.0d0*(obj%omdmomp &
                            & +obj%ombomp)**2.0d0/(obj%homp)**3.0d0*(1.0d0+zlens)**2.0d0)
        endif
    END function sjclsobjnoweight

    !weight for trapezoid integration
    function sjclsobjonlyweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclsobjonlyweight,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz !assuming flatness
        distz = f_K(distz) !with curvature
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1 !doesn't matter here because imposed momp(1) = momp(2)
        sjclsobjonlyweight = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumaxwomp,obj%tol_erroromp) !note changed mumax to mumax1

    END function sjclsobjonlyweight

    !IA-weight for trapezoid integration
    function sjclsiiobjonlyweight(obj,zlens)
        type(my_type) :: obj
        REAL(mcp) :: zlens,sjclsiiobjonlyweight,kval,deltafinal,hubblez,growthsf,growthnorm
        real(mcp) :: kminsj,kmaxsj

        growthnorm = obj%gcubicomp%Value(zlens)	!this is D(0)/D(z)
	hubblez = Hofz(zlens)
        obj%wchooseomp = 1 !doesn't matter here
        sjclsiiobjonlyweight = psourceobjcubic(obj,zlens)*hubblez*(-obj%ampiayesomp*5.0d-14*2.77536627d11*(obj%omdmomp+obj%ombomp)*growthnorm*((1.0d0+zlens)/(1.0d0+0.3d0))**obj%redziayesomp)*(obj%lumarromp(obj%momp(1)))**obj%lumiayesomp !enforced momp(1) = momp(2), so either is fine

    END function sjclsiiobjonlyweight


    !WEAK LENSING INTEGRAND (GG) ---- wrt scale factor instead
    function sjclsobjsf(obj,alens)
        type(my_type) :: obj
        REAL(mcp) :: alens,zlens,sjclsobjsf,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms
        real(mcp) :: kminsj,kmaxsj

        zlens = 1.0d0/alens - 1.0d0
        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz !assuming flatness
        distz = f_K(distz) !with curvature
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)
        kminsj = 1.0d-5
        kmaxsj = 100.0d0 !kmaxsjhihi
        deltafinal = 0.0d0
        if(kval < kminsj) obj%kline = obj%kline + 1
        if (this%use_nl== .true.) then
         if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK%PowerAt(kval,zlens)
       	else
         if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK_WEYL%PowerAt(kval,zlens) !NHmod
        endif
	hubblez = Hofz(zlens)
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1
        weight1 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumax1omp,obj%tol_erroromp) !note changed mumax to mumax1
        if(obj%bnumomp == 1 .or. obj%bnumomp == 5 .or. obj%bnumomp == 8 .or. obj%bnumomp == 10) then !generalize later
            weight2 = weight1
        else
            obj%wchooseomp = 2
            weight2 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumax2omp,obj%tol_erroromp)
        end if
        if (this%use_nl ==.true.) then !NHmod
        sjclsobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0
        else
       sjclsobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0/(9.0d0/4.0d0*(obj%h0omp/ckms)**4.0d0*(obj%omdmomp &
                            & +obj%ombomp)**2.0d0/(obj%homp)**3.0d0*(1.0d0+zlens)**2.0d0)
       endif
 END function sjclsobjsf


    !INTRINSIC ALIGNMENT INTEGRAND (II) ---- wrt scale factor instead
    function sjclsiiobjsf(obj,alens)
        type(my_type) :: obj
        REAL(mcp) :: alens,zlens,sjclsiiobjsf,growthnorm,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms,growthsf
        real(mcp) :: kminsj,kmaxsj

        zlens = 1.0d0/alens - 1.0d0
        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz !assuming flatness ---- need to fix this
        distz = f_K(distz) !with curvature
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)
        kminsj = 1.0d-5
        kmaxsj = 100.0d0 !kmaxsjhihi
        deltafinal = 0.0d0
        if (this%use_nl== .true.) then
         if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK%PowerAt(kval,zlens)
        else
         if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK_WEYL%PowerAt(kval,zlens) !NHmod
        endif
        growthnorm = obj%gcubicomp%Value(zlens)
        deltafinal = deltafinal*((obj%lumarromp(obj%momp(1))*obj%lumarromp(obj%momp(2)))**obj%lumiayesomp)*(-obj%ampiayesomp*5.0d-14*2.77536627d11*(obj%omdmomp+obj%ombomp)*growthnorm*((1.0d0+zlens)/(1.0d0+0.3d0))**obj%redziayesomp)**2.0d0
	hubblez = Hofz(zlens)
        obj%wchooseomp = 1
        weight1 = psourceobjcubic(obj,zlens)*hubblez
        if(obj%bnumomp == 1 .or. obj%bnumomp == 5 .or. obj%bnumomp == 8 .or. obj%bnumomp == 10) then !generalize later
            weight2 = weight1
        else
            obj%wchooseomp = 2
            weight2 = psourceobjcubic(obj,zlens)*hubblez
        end if
        if (this%use_nl ==.true.) then
        sjclsiiobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0
        else
       sjclsiiobjsf = ((1.0d0+zlens)**2.0d0)*weight1*weight2/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0/(9.0d0/4.0d0*(obj%h0omp/ckms)**4.0d0*(obj%omdmomp+ &
                       & obj%ombomp)**2.0d0/(obj%homp)**3.0d0*(1.0d0+zlens)**2.0d0)
       end if !NHmod
    END function sjclsiiobjsf

    !LENSING - INTRINSIC ALIGNMENT INTEGRAND (GI) --- wrt scale factor instead
    function sjclsgiobjsf(obj,alens)
        type(my_type) :: obj
        REAL(mcp) :: alens,zlens,sjclsgiobjsf,weightpart,weight1,weight2,distz,kval,deltafinal,hubblez,ckms,weight3,weight4,growthsf,growthnorm
        real(mcp) :: kminsj,kmaxsj

        zlens = 1.0d0/alens - 1.0d0
        ckms = 299792.458d0
        distz = ComovingRadialDistance(zlens)
        obj%distlensflatomp = distz !assuming flatness
        distz = f_K(distz) !with curvature
        kval = (obj%ellarromp(obj%ellgenomp)+1.0d0/2.0d0)/((obj%h0omp/100.0d0)*distz)
        kminsj = 1.0d-5
        kmaxsj = 100.0d0 !kmaxsjhihi
        deltafinal = 0.0d0
        if (this%use_nl== .true.) then
         if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK%PowerAt(kval,zlens)
        else
         if((kval >= kminsj) .and. (kval <= kmaxsj)) deltafinal = obj%theoryomp%NL_MPK_WEYL%PowerAt(kval,zlens) !NHmod
        endif
        growthnorm = obj%gcubicomp%Value(zlens)	!this is D(0)/D(z)
        deltafinal = deltafinal*(-obj%ampiayesomp*5.0d-14*2.77536627d11*(obj%omdmomp+obj%ombomp)*growthnorm*((1.0d0+zlens)/(1.0d0+0.3d0))**obj%redziayesomp)
	hubblez = Hofz(zlens)
        weightpart = 3.0d0/2.0d0*(obj%omdmomp+obj%ombomp)*(obj%h0omp/ckms)**2.0d0*distz*(1.0d0+zlens)
        obj%wchooseomp = 1
        weight1 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumax1omp,obj%tol_erroromp) !note changed mumax to mumax1
        weight3 = psourceobjcubic(obj,zlens)*hubblez
        if(obj%bnumomp == 1 .or. obj%bnumomp == 5 .or. obj%bnumomp == 8 .or. obj%bnumomp == 10) then !generalize later
            weight2 = weight1
            weight4 = weight3
        else
            obj%wchooseomp = 2
            weight2 = weightpart*rombint_obj(obj,weightobjcubic,zlens,obj%mumax2omp,obj%tol_erroromp) !note changed mumax to mumax2
            weight4 = psourceobjcubic(obj,zlens)*hubblez
        end if
        if (this%use_nl== .true.) then !NHmod
        sjclsgiobjsf = ((1.0d0+zlens)**2.0d0)*(weight1*weight4*(obj%lumarromp(obj%momp(2)))**obj%lumiayesomp +&
                      & weight2*weight3*(obj%lumarromp(obj%momp(1)))**obj%lumiayesomp)/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0
        else
        sjclsgiobjsf = ((1.0d0+zlens)**2.0d0)*(weight1*weight4*(obj%lumarromp(obj%momp(2)))**obj%lumiayesomp +&
                      & weight2*weight3*(obj%lumarromp(obj%momp(1)))**obj%lumiayesomp)/(distz**2.0d0)/hubblez*deltafinal/(obj%homp)**3.0d0/(9.0d0/4.0d0*(obj%h0omp/ckms)**4.0d0*(obj%omdmomp+&
                      & obj%ombomp)**2.0d0/(obj%homp)**3.0d0*(1.0d0+zlens)**2.0d0)
        endif

 END function sjclsgiobjsf


    !lensing weight, cubic spline of source distribution
    function weightobjcubic(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,weightobjcubic,chis,diff,ckms

        chis = ComovingRadialDistance(zs)
        diff = f_K(chis-obj%distlensflatomp) !with curvature

        weightobjcubic=0.0d0
        if((zs-obj%aphotarromp(obj%momp(obj%wchooseomp))) >= obj%mumin2vomp)  weightobjcubic = obj%psourcetypeomp(obj%momp(obj%wchooseomp))%Value(zs-obj%aphotarromp(obj%momp(obj%wchooseomp)))
        weightobjcubic = weightobjcubic*diff/f_K(chis)

    END function weightobjcubic


    !cubic spline of source distribution
    function psourceobjcubic(obj,zs)
        type(my_type) :: obj
        REAL(mcp) :: zs,psourceobjcubic

        psourceobjcubic=0.0d0
        if((zs-obj%aphotarromp(obj%momp(obj%wchooseomp))) >= obj%mumin2vomp)  psourceobjcubic = obj%psourcetypeomp(obj%momp(obj%wchooseomp))%Value(zs-obj%aphotarromp(obj%momp(obj%wchooseomp)))

    END function psourceobjcubic


    END function CosmicShear_LnLike

    end module CosmicShear
