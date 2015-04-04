#include "pgmlink/inferencemodel/perturbedinferencemodel.h"
#include "pgmlink/inferencemodel/perturbation/gaussian_perturbation.h"
#include "pgmlink/inferencemodel/perturbation/divmbest_perturbation.h"
#include "pgmlink/inferencemodel/perturbation/classifier_uncertainty_perturbation.h"
#include "pgmlink/inferencemodel/perturbation/perturbandmap_perturbation.h"

namespace pgmlink
{

PerturbedInferenceModel::PerturbedInferenceModel(const ConsTrackingInferenceModel::Parameter& param,
                                                 const Perturbation::Parameter &perturbation_param,
                                                 double ep_gap,
                                                 double cplex_timeout):
    ConsTrackingInferenceModel(param, ep_gap, cplex_timeout)
{
    // instanciate perturbation depending on distribution type
    switch(perturbation_param.distributionId)
    {
        case Gaussian:
            perturbation_ = boost::make_shared<GaussianPerturbation>(perturbation_param, param);
            break;
        case PerturbAndMAP:
            perturbation_ = boost::make_shared<PerturbAndMapPerturbation>(perturbation_param, param);
            break;
        case DiverseMbest:
            perturbation_ = boost::make_shared<DivMBestPerturbation>(perturbation_param, param);
            break;
        case ClassifierUncertainty:
            perturbation_ = boost::make_shared<ClassifierUncertaintyPerturbation>(perturbation_param, param);
            break;
        default:
            throw std::runtime_error("The chosen perturbation distribution is not available "
                                     "with the current inference model");
    }
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
