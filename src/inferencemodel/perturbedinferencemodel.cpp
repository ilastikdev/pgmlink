#include "pgmlink/inferencemodel/perturbedinferencemodel.h"

namespace pgmlink
{

PerturbedInferenceModel::PerturbedInferenceModel(
    const ConsTrackingInferenceModel::Parameter& param,
    const Perturbation::Parameter &perturbation_param,
    double ep_gap,
    double cplex_timeout):
    ConsTrackingInferenceModel(param, ep_gap, cplex_timeout),
    perturbation_(new Perturbation(perturbation_param, param))
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

void PerturbedInferenceModel::perturb(Perturbation::DeterministicOffset *det_off)
{
    perturbation_->perturb(det_off);
}

} // namespace pgmlink
