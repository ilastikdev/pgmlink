#include "pgmlink/inferencemodel/perturbation/divmbest_perturbation.h"

namespace pgmlink
{

DivMBestPerturbation::DivMBestPerturbation(const Parameter &perturbation_param,
                           const InferenceModel::Parameter &inf_param):
    Perturbation(perturbation_param, inf_param),
    deterministic_offset_(NULL)
{
    if(perturbation_param_.distributionId != DiverseMbest)
        throw std::runtime_error("Cannot construct DivMBestPerturbation "
                                 "when another distribution is specified");
}

void DivMBestPerturbation::perturb(Perturbation::DeterministicOffset *det_off)
{
    deterministic_offset_ = det_off;
}

double DivMBestPerturbation::generateRandomOffset(EnergyType energyIndex,
                                          double energy,
                                          Traxel tr,
                                          Traxel tr2,
                                          size_t state,
                                          boost::shared_ptr<InferenceModel::TransitionPredictionsMap> transition_predictions)
{

    LOG(logDEBUG4) << "generateRandomOffset()";

    LOG(logDEBUG4) << "DiverseMBest/MBestCPLEX: random offset 0";
    return 0;
}

size_t DivMBestPerturbation::add_div_m_best_perturbation(marray::Marray<double>& energies,
        EnergyType energy_type,
        size_t factorIndex)
{
    if (perturbation_param_.distributionId == DiverseMbest)
    {
        assert(deterministic_offset_ != NULL);
        std::vector<std::vector<size_t> >& indexlist = (*deterministic_offset_)[factorIndex];
        for (std::vector<std::vector<size_t> >::iterator index = indexlist.begin(); index != indexlist.end(); index++)
        {
            energies(index->begin()) += perturbation_param_.distributionParam[energy_type];
        }
        factorIndex++;
    }
    return factorIndex;
}

} // namespace pgmlink
