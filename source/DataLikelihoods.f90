    module DataLikelihoodList
    use likelihood
    use settings
    use CosmologyTypes
    implicit none

    contains

    subroutine SetDataLikelihoods(Ini)
    use HST
    use snovae
    use CMBLikelihoods
    use bao
    use mpk
    !use CosmicShear !SJ   
    use wigglez
    use szcounts !Anna
    use wl
    use ElementAbundances
    use corrprior !MMmod: VOID
    class(TSettingIni), intent(in) :: Ini

    CosmoSettings%get_sigma8 = Ini%Read_Logical('get_sigma8',.false.)

    call CMBLikelihood_Add(DataLikelihoods, Ini)

    call AbundanceLikelihood_Add(DataLikelihoods, Ini)

    call HSTLikelihood_Add(DataLikelihoods, Ini)

    call SNLikelihood_Add(DataLikelihoods, Ini)

    call MPKLikelihood_Add(DataLikelihoods, Ini)

    if (use_mpk) call WiggleZLikelihood_Add(DataLikelihoods, Ini)

    call BAOLikelihood_Add(DataLikelihoods, Ini)
    !call CosmicShearLikelihood_Add(DataLikelihoods, Ini) !SJ
    call SZLikelihood_Add(DataLikelihoods, Ini) !Anna

    call WLLikelihood_Add(DataLikelihoods, Ini)

    CosmoSettings%void_n = Ini%Read_Int('number_of_bins')
    call CPLikelihood_Add(DataLikelihoods, CosmoSettings%void_n, Ini) !MMmod: void

    end subroutine SetDataLikelihoods


    end module DataLikelihoodList
