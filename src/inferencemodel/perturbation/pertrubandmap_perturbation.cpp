#include "pgmlink/inferencemodel/perturbation/perturbandmap_perturbation.h"

namespace pgmlink
{

PerturbAndMapPerturbation::PerturbAndMapPerturbation(const Perturbation::Parameter &perturbation_param,
                           const pgmlink::Parameter &inf_param):
    Perturbation(perturbation_param, inf_param)
{
    if(perturbation_param_.distributionId != PerturbAndMAP)
        throw std::runtime_error("Cannot construct PerturbAndMapPerturbation "
                                 "when another distribution is specified");
}

double PerturbAndMapPerturbation::generateRandomOffset(EnergyType energyIndex,
                                          double energy,
                                          Traxel tr,
                                          Traxel tr2,
                                          size_t state,
                                          boost::shared_ptr<InferenceModel::TransitionPredictionsMap> transition_predictions)
{
    double rand;

    LOG(logDEBUG4) << "PerturbAndMAP";
    //distribution parameter: beta
    if (energyIndex >= perturbation_param_.distributionParam.size())
    {
        throw std::runtime_error("sigma is not set correctly");
    }
    rand = random_uniform_();
    //throw std::runtime_error("I don't think this formula is correct; debug when needed; check whether rand>=0");

    return perturbation_param_.distributionParam[energyIndex] * log(-log(rand));
}

} // namespace pgmlink
