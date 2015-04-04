#include "pgmlink/inferencemodel/dynprog_perturbedinferencemodel.h"
#include "pgmlink/inferencemodel/perturbation/gaussian_perturbation.h"
#include "pgmlink/inferencemodel/perturbation/divmbest_perturbation.h"
#include "pgmlink/inferencemodel/perturbation/classifier_uncertainty_perturbation.h"
#include "pgmlink/inferencemodel/perturbation/perturbandmap_perturbation.h"

#ifdef WITH_DPCT

namespace pgmlink
{

DynProgPerturbedInferenceModel::DynProgPerturbedInferenceModel(const InferenceModel::Parameter& param,
                                                               const Perturbation::Parameter &perturbation_param):
    DynProgConsTrackInferenceModel(param)
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
            throw std::runtime_error("Diverse M Best perturbation does not work with Magnusson yet");
//            perturbation_ = boost::make_shared<DivMBestPerturbation>(perturbation_param, param);
            break;
        case ClassifierUncertainty:
            perturbation_ = boost::make_shared<ClassifierUncertaintyPerturbation>(perturbation_param, param);
            break;
        default:
            throw std::runtime_error("The chosen perturbation distribution is not available "
                                     "with the current inference model");
    }
}

double DynProgPerturbedInferenceModel::generateRandomOffset(EnergyType parameterIndex,
                                                     double energy,
                                                     Traxel tr,
                                                     Traxel tr2,
                                                     size_t state)
{
    return perturbation_->generateRandomOffset(parameterIndex, energy, tr, tr2, state, transition_predictions_);
}

} // namespace pgmlink

#endif // WITH_DPCT
