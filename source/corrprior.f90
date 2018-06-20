module corrprior

    use settings
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology
    use Likelihood_Cosmology
    use MatrixUtils
    implicit none
    private

    type, extends(TCosmoCalcLikelihood) :: CPLikelihood
        integer :: nbins
        real(mcp), allocatable, dimension(:,:) :: inv_cp_cov,cp_cov
        real(mcp), allocatable, dimension(:) :: theta_bins
        real(mcp), allocatable, dimension(:) :: z_bins
        real(mcp), allocatable, dimension(:) :: xi_obs ! Observed correlation functions
        real(mcp) :: alpha_c, sigma_alpha
    contains
    procedure :: LogLike => CP_LnLike
    procedure :: rombint_corr
    procedure :: rombint_corr_2d
    !procedure :: ReadIni => CP_ReadIni
    end type CPLikelihood

    logical, parameter :: CPdebugging = .false.

    public CPLikelihood, CPLikelihood_Add
    contains

    function csi(a,void_redshift, nbins,alpha_c,sigma_alpha)
      real (mcp) :: a
      real (mcp) csi
      integer :: nbins
      real (mcp) csi_0,a_max, a_min
      real(mcp) alpha_c,sigma_alpha
      real (mcp), dimension(nbins+1), intent(in) :: void_redshift
      !
      a_max = 1./(1.+void_redshift(1))
      a_min = 1./(1.+void_redshift(nbins+1))
      ! ! if (CPdebugging) then
      ! !   print*, a_max, a_min
      ! ! end if
      csi_0 = sigma_alpha**2./(alpha_c*PI)*(a_max - a_min)
      csi = csi_0/(1. +(abs(a)/alpha_c)**2.)

    end function csi

    function rombint_corr(this,f,a,b,x,tol)
    use Precision
    !  rombint_corr returns the integral from a to b of using Romberg integration.
    !  The method converges provided that f(x) is continuous in (a,b).
    !  f must be real(dl) and must be declared external in the calling
    !  routine.  tol indicates the desired relative accuracy in the integral.
    !
    implicit none
    integer, parameter :: MAXITER=20
    integer, parameter :: MAXJ=5
    real(mcp), dimension(MAXJ+1):: g
    real(mcp) f
    external f
    Class(CPLikelihood) :: this

    real(mcp) :: rombint_corr
    real(mcp), intent(in) :: a,b,tol, x
    integer :: nint, i, k, jmax, j
    real(mcp) :: h, gmax, error,  g0, g1, fourj
    !
    h=0.5d0*(b-a)
    gmax=h*(f(a-x,this%z_bins, this%nbins,this%alpha_c,this%sigma_alpha)+f(b-x,this%z_bins, this%nbins,this%alpha_c,this%sigma_alpha))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
        go to 40
    !  Calculate next trapezoidal rule approximation to integral.
    g0=0._dl
    do 20 k=1,nint
        g0=g0+f(a-x+(k+k-1)*h,this%z_bins, this%nbins,this%alpha_c,this%sigma_alpha)
20  continue
    g0=0.5d0*g(1)+h*g0
    h=0.5d0*h
    nint=nint+nint
    jmax=min(i,MAXJ)
    fourj=1._dl
    do 30 j=1,jmax
        !  Use Richardson extrapolation.
        fourj=4._dl*fourj
        g1=g0+(g0-g(j))/(fourj-1._dl)
        g(j)=g0
        g0=g1
30  continue
    if (abs(g0).gt.tol) then
        error=1._dl-gmax/g0
    else
        error=gmax
    end if
    gmax=g0
    g(jmax+1)=g0
    go to 10
40  rombint_corr=g0
    if (i.gt.MAXITER.and.abs(error).gt.tol)  then
        write(*,*) 'Warning: rombint_corr failed to converge; '
        write (*,*)'integral, error, tol:', rombint_corr,error, tol
    end if

    end function rombint_corr

    function rombint_corr_2d(this,xi,a,b,c,d,tol)
    use Precision
    !  rombint_corr returns the integral from a to b of using Romberg integration.
    !  The method converges provided that f(x) is continuous in (a,b).
    !  f must be real(dl) and must be declared external in the calling
    !  routine.  tol indicates the desired relative accuracy in the integral.
    !
    implicit none
    integer, parameter :: MAXITER=20
    integer, parameter :: MAXJ=5
    real(mcp), dimension(MAXJ+1):: g
    real(mcp) xi
    external xi
    Class(CPLikelihood) :: this

    real(mcp) :: rombint_corr_2d
    real(mcp), intent(in) :: a,b,tol,c,d
    integer :: nint, i, k, jmax, j
    real(mcp) :: h, gmax, error,  g0, g1, fourj
    !
    h=0.5d0*(b-a)
    gmax=h*(this%rombint_corr(xi,c ,d, a,tol)+this%rombint_corr(xi,c ,d, b,tol))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) &
        go to 40
    !  Calculate next trapezoidal rule approximation to integral.
    g0=0._dl
    do 20 k=1,nint
        g0=g0+this%rombint_corr(xi,c ,d, a+(k+k-1)*h,tol)
20  continue
    g0=0.5d0*g(1)+h*g0
    h=0.5d0*h
    nint=nint+nint
    jmax=min(i,MAXJ)
    fourj=1._dl
    do 30 j=1,jmax
        !  Use Richardson extrapolation.
        fourj=4._dl*fourj
        g1=g0+(g0-g(j))/(fourj-1._dl)
        g(j)=g0
        g0=g1
30  continue
    if (abs(g0).gt.tol) then
        error=1._dl-gmax/g0
    else
        error=gmax
    end if
    gmax=g0
    g(jmax+1)=g0
    go to 10
40  rombint_corr_2d=g0
    if (i.gt.MAXITER.and.abs(error).gt.tol)  then
        write(*,*) 'Warning: rombint_corr_2d failed to converge; '
        write (*,*)'integral, error, tol:', rombint_corr_2d,error, tol
    end if

  end function rombint_corr_2d

    subroutine CPLikelihood_Add(LikeList, voidbins, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(TTextFile) :: F
    Type(CMBParams) CMB
    Type(CPLikelihood), pointer :: this
    integer :: voidbins
    character(LEN=:), allocatable :: covariance_file
    integer i,j
    real (mcp), allocatable, dimension(:) :: Delta, scale_factor, void_redshift
    real(mcp), PARAMETER :: tol=1.0E-5

    if (Ini%Read_Logical('use_corrprior',.false.)) then
        allocate(this)

        this%nbins = voidbins
        !llocate the arrays
        allocate(this%cp_cov(this%nbins,this%nbins))
        allocate(this%inv_cp_cov(this%nbins,this%nbins))
        allocate(Delta(this%nbins))
        allocate(scale_factor(this%nbins+1))
        allocate(this%z_bins(this%nbins+1))
        !SPmod: insert here the computation of the covariance matrix
        scale_factor(1)= 1.
        this%z_bins(1) = 0.
        this%z_bins(2) = Ini%Read_Double('param[void_z0]')
        this%z_bins(3) = Ini%Read_Double('param[void_z1]')
        this%z_bins(4) = Ini%Read_Double('param[void_z2]')
        this%z_bins(5) = Ini%Read_Double('param[void_z3]')
        this%alpha_c = Ini%Read_Double('alpha_c')
        this%sigma_alpha = Ini%Read_Double('sigma_alpha')
        do i = 1, this%nbins
          Delta(i) = 1./(1.+this%z_bins(i)) - 1./(1.+this%z_bins(i+1))
          scale_factor(i+1) = 1./(1.+this%z_bins(i+1))
        end do

        if (CPdebugging) then
          do i=1,this%nbins
             write(*,"(100E15.8)") Delta(i),scale_factor(i),this%z_bins(i)
          end do
          write(*,*) '-----------------------------------------'
          do i=1,this%nbins+1
            write(*,"(100E15.8)") scale_factor(i), csi(scale_factor(i),this%z_bins, this%nbins,this%alpha_c,this%sigma_alpha)
          end do
          write(*,*) '-----------------------------------------'
          do i=1,this%nbins
            do j=1,this%nbins
              write(*,*) scale_factor(i),1./(Delta(i)*Delta(j))*this%rombint_corr_2d( csi ,scale_factor(i+1) ,scale_factor(i),scale_factor(j+1) ,scale_factor(j), tol )
            end do
          end do
        end if

        do i = 1, this%nbins
          do j = 1, this%nbins

            this%cp_cov(i,j) =1./(Delta(i)*Delta(j))*this%rombint_corr_2d( csi ,scale_factor(i+1) ,scale_factor(i),scale_factor(j+1) ,scale_factor(j), tol )

          end do
        end do

        this%inv_cp_cov = this%cp_cov
        call Matrix_Inverse(this%inv_cp_cov)

        if (CPdebugging) then
           do i=1,this%nbins
              write(*,"(100E15.8)") (this%inv_cp_cov(i,j),j=1,this%nbins)
           end do
        end if

        this%LikelihoodType = 'CP'
        call LikeList%Add(this)
        if (Feedback>1) write(*,*) 'read Correlation Prior dataset'
    end if

    end subroutine CPLikelihood_Add

    function CP_LnLike(this, CMB, Theory, DataParams)
    use MatrixUtils
    Class(CPLikelihood) :: this
    Class(CMBParams) CMB
    Class(TCosmoTheoryPredictions), target :: Theory
    real(mcp) :: DataParams(:)
    real(mcp) CP_LnLike
    real(mcp) vec(this%nbins)
    integer i,j
    real(mcp) void_mean

    if ( CMB%void_mean_fiducial == 2 ) then !SPmod: fiducial as mean of neighbours bins

      i=1
      void_mean = (CMB%void_qV(i) + CMB%void_qV(i+1))/2._mcp
      vec(i) = CMB%void_qV(i)-void_mean

      do i=2, this%nbins-1
        void_mean = (CMB%void_qV(i-1) + CMB%void_qV(i) + CMB%void_qV(i+1))/3._mcp
        vec(i) = CMB%void_qV(i)-void_mean
      end do

      i= this%nbins
      void_mean = (2.*CMB%void_qV(i) + CMB%void_qV(i-1))/3._mcp
      vec(i) = CMB%void_qV(i)-void_mean

    else if ( CMB%void_mean_fiducial == 1 ) then !fixed fiducial

      do i=1, this%nbins
        vec(i) = CMB%void_qV(i)-CMB%void_fiducial
      end do

    else

      call MpiStop(' no choice for q fiducial in correlation prior computation ')

    end if

    if (CPdebugging) then
       do i=1,this%nbins
          write(*,*) CMB%void_qV(i),vec(i)
       end do
    end if
    CP_LnLike = dot_product( vec, MatMul(this%inv_cp_cov,vec))
    CP_LnLike = CP_LnLike/2._mcp
    if (Feedback>1) write(*,*) 'CPLnLike=',CP_LnLike

    if (CPdebugging) stop 'bye bye'

    end function CP_LnLike


end module corrprior
