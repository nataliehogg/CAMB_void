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
        integer :: num_z_p ! Source galaxy distribution p(z,bin)
        real(mcp), allocatable, dimension(:) :: z_p
    contains
    procedure :: LogLike => CP_LnLike
    !procedure :: ReadIni => CP_ReadIni
    end type CPLikelihood

    logical, parameter :: CPdebugging = .true.

    public CPLikelihood, CPLikelihood_Add
    contains

    subroutine CPLikelihood_Add(LikeList, voidbins, Ini)
    class(TLikelihoodList) :: LikeList
    class(TSettingIni) :: ini
    Type(TTextFile) :: F
    Type(CMBParams) CMB
    Type(CPLikelihood), pointer :: this
    Type(TSettingIni) :: DataSets
    integer :: voidbins
    character(LEN=:), allocatable :: covariance_file
    integer i,j

    if (Ini%Read_Logical('use_corrprior',.false.)) then
        call Ini%TagValuesForName('corrprior_dataset', DataSets)
        allocate(this)
        covariance_file  = Ini%ReadFileName('covariance_file')

        this%nbins = voidbins 
        !Reading covariance matrix
        call F%Open(covariance_file)
        
        allocate(this%cp_cov(this%nbins,this%nbins))
        allocate(this%inv_cp_cov(this%nbins,this%nbins))
        do i=1,this%nbins
            read (F%unit,*) (this%cp_cov(i,j),j=1,this%nbins)
        end do
        call F%Close()

        if (CPdebugging) then
           do i=1,this%nbins
              write(*,"(100E15.8)") (this%cp_cov(i,j),j=1,this%nbins)
           end do
        end if

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


    do i=1, this%nbins
       vec(i) = CMB%void_qV(i)-CMB%void_fiducial
    end do
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
