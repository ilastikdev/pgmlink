#include "pgmlink/inferencemodel/perturbation/divmbest_perturbation.h"

namespace pgmlink
{

DivMBestPerturbation::DivMBestPerturbation(const Parameter &perturbation_param,
                           const InferenceModel::Parameter &inf_param):
    Perturbation(perturbation_param, inf_param)
{
    if(perturbation_param_.distributionId != DiverseMbest)
        throw std::runtime_error("Cannot construct DivMBestPerturbation "
                                 "when another distribution is specified");
}

void DivMBestPerturbation::push_away_from_solution(const PertGmType &model, std::vector<size_t> solution)
{
    size_t num_factors = model.numberOfFactors();
    if(deterministic_offset_.size() != num_factors)
        deterministic_offset_.resize(num_factors);

    for (size_t factorId = 0; factorId < num_factors; ++factorId)
    {
        PertGmType::FactorType factor = model[factorId];
        vector<size_t> varIndices;
        for (PertGmType::FactorType::VariablesIteratorType ind = factor.variableIndicesBegin();
                ind != factor.variableIndicesEnd();
                ++ind)
        {
            varIndices.push_back(solution[*ind]);
        }
        deterministic_offset_[factorId].push_back(varIndices);
    }
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
        std::vector<std::vector<size_t> >& indexlist = deterministic_offset_[factorIndex];
        for (std::vector<std::vector<size_t> >::iterator index = indexlist.begin(); index != indexlist.end(); index++)
        {
            energies(index->begin()) += perturbation_param_.distributionParam[energy_type];
        }
        factorIndex++;
    }
    return factorIndex;
}

} // namespace pgmlink
