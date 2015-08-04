#include "pgmlink/inferencemodel/dynprog_perturbedinferencemodel.h"
#include "pgmlink/inferencemodel/perturbation/gaussian_perturbation.h"
#include "pgmlink/inferencemodel/perturbation/divmbest_perturbation.h"
#include "pgmlink/inferencemodel/perturbation/classifier_uncertainty_perturbation.h"
#include "pgmlink/inferencemodel/perturbation/perturbandmap_perturbation.h"

#ifdef WITH_DPCT

namespace pgmlink
{

DynProgPerturbedInferenceModel::DynProgPerturbedInferenceModel(const InferenceModel::Parameter& param,
                                                               boost::shared_ptr<Perturbation> perturbation):
    DynProgConsTrackInferenceModel(param),
    perturbation_(perturbation)
{
}

double DynProgPerturbedInferenceModel::generateRandomOffset(EnergyType parameterIndex,
                                                     double energy,
                                                     Traxel tr,
                                                     Traxel tr2,
                                                     size_t state)
{
    return perturbation_->generateRandomOffset(parameterIndex, energy, tr, tr2, state, transition_predictions_);
}

double DynProgPerturbedInferenceModel::getDivMBestOffset(EnergyType energy_type,
                                                         const HypothesesGraph& g,
                                                         HypothesesGraph::Node n,
                                                         HypothesesGraph::Arc a,
                                                         size_t state)
{
    return perturbation_->getDivMBestOffset(energy_type, g, n, a, state);
}

} // namespace pgmlink

#endif // WITH_DPCT
