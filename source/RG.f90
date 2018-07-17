    ! RG dataset
    !
    ! Radio Galaxies Dataset from
    ! "Daly et al. (SCP) 2011, arXiv:0710.5112v2".

    module RG
    use cmbtypes
    use MatrixUtils
    use likelihood
    implicit none

    integer, parameter :: RG_num = 30

    type, extends(CosmologyLikelihood) :: RGLikelihood
        double precision :: RG_z(RG_num), RG_moduli(RG_num), RG_modulierr(RG_num)
    contains
    procedure :: LogLikeTheory => RG_LnLike
    end type RGLikelihood

    contains

    subroutine RGLikelihood_Add(LikeList, Ini)
    use IniFile
    use settings
    class(LikelihoodList) :: LikeList
    Type(TIniFile) :: ini
    Type(RGLikelihood), pointer :: like
    integer i
    ! The following line selects which error estimate to use
    ! default .True. = with systematic errors

    if (.not. Ini_Read_Logical_File(Ini, 'use_RG',.false.)) return

    allocate(like)
    Like%LikelihoodType = 'RG'
    Like%name='Daly2009'
    like%needs_background_functions = .true.
    call LikeList%Add(like)

    if (Feedback > 0) write (*,*) 'Reading: RG data'
    call OpenTxtFile(trim(DataDir)//'radio_z_mu_dmu.dat',tmp_file_unit)
    do i=1,  RG_num
        read(tmp_file_unit, *) Like%RG_z(i),Like%RG_moduli(i), Like%RG_modulierr(i)
    end do
    close(tmp_file_unit)

    end subroutine RGLikelihood_Add

    function RG_LnLike(like, CMB)
    use camb
    !Assume this is called just after CAMB with the correct model  use camb
    Class(CMBParams) CMB
    Class(RGLikelihood) :: like
    real(mcp) RG_LnLike
    integer i
    double precision z
    real(mcp) diffs(RG_num), chisq

    chisq=0.
    do i=1, RG_num
        z= Like%RG_z(i)
        diffs(i) = 5*log10((1+z)**2*AngularDiameterDistance(z))+25 -Like%RG_moduli(i)
        chisq=chisq+(diffs(i)**2)/(Like%RG_modulierr(i)**2)
    end do

    !! H0 normalisation alla Bridle and co.

    if (Feedback > 1) write (*,*) 'RG chisq: ', chisq

    RG_LnLike = chisq/2


    end function RG_LnLike


    end module RG
