#include "pgmlink/inferencemodel/perturbedinferencemodel.h"

namespace pgmlink
{

PerturbedInferenceModel::PerturbedInferenceModel(const ConsTrackingInferenceModel::Parameter& param,
                                                 boost::shared_ptr<Perturbation> perturbation,
                                                 double ep_gap,
                                                 double cplex_timeout,
                                                 unsigned int num_threads):
    ConsTrackingInferenceModel(param, ep_gap, cplex_timeout, num_threads),
    perturbation_(perturbation)
{
}

double PerturbedInferenceModel::generateRandomOffset(EnergyType parameterIndex,
                                                     double energy,
                                                     Traxel tr,
                                                     Traxel tr2,
                                                     size_t state)
{
    return perturbation_->generateRandomOffset(parameterIndex, energy, tr, tr2, state, transition_predictions_);
}

size_t PerturbedInferenceModel::add_div_m_best_perturbation(marray::Marray<double> &energies,
                                                            EnergyType energy_type,
                                                            size_t factorIndex)
{
    return perturbation_->add_div_m_best_perturbation(energies, energy_type, factorIndex);
}

} // namespace pgmlink
