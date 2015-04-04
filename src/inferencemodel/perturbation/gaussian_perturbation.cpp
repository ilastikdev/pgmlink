#include "pgmlink/inferencemodel/perturbation/gaussian_perturbation.h"

namespace pgmlink
{

GaussianPerturbation::GaussianPerturbation(const Parameter &perturbation_param,
                           const InferenceModel::Parameter &inf_param):
    Perturbation(perturbation_param, inf_param)
{
    if(perturbation_param_.distributionId != Gaussian)
        throw std::runtime_error("Cannot construct GaussianPerturbation "
                                 "when another distribution is specified");
}

double GaussianPerturbation::generateRandomOffset(EnergyType energyIndex,
                                          double energy,
                                          Traxel tr,
                                          Traxel tr2,
                                          size_t state,
                                          boost::shared_ptr<InferenceModel::TransitionPredictionsMap> transition_predictions)
{

    LOG(logDEBUG4) << "generateRandomOffset()";

    double rand;

    LOG(logDEBUG4) << "GaussianPerturbation";
    if (energyIndex >= perturbation_param_.distributionParam.size())
    {
        throw std::runtime_error("sigma is not set correctly");
    }
    return random_normal_() * perturbation_param_.distributionParam[energyIndex];
}

} // namespace pgmlink
