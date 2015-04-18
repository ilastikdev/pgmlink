#ifndef DYNPROG_PERTURBEDINFERENCEMODEL_H
#define DYNPROG_PERTURBEDINFERENCEMODEL_H

#ifdef WITH_DPCT

#include <vector>
#include <iostream>
#include <opengm/datastructures/marray/marray.hxx>
#include "dynprog_constrackinginferencemodel.h"
#include "pgmlink/uncertaintyParameter.h"
#include "pgmlink/inferencemodel/perturbation/perturbation.h"

namespace pgmlink
{

/**
 * @brief The DynProg_PerturbedInferenceModel class specializes the DynProg_ConsTrackingInferenceModel
 * and adds perturbation to the unaries by providing non-zero perturbations in virtual methods.
 */
class DynProgPerturbedInferenceModel : public DynProgConsTrackInferenceModel
{
public: // API
    DynProgPerturbedInferenceModel(const InferenceModel::Parameter& param,
                                   boost::shared_ptr<Perturbation> perturbation);

    virtual double generateRandomOffset(EnergyType parameterIndex,
                                        double energy = 0,
                                        Traxel tr = 0,
                                        Traxel tr2 = 0,
                                        size_t state = 0);

    virtual double getDivMBestOffset(EnergyType energy_type,
                                     const HypothesesGraph& g,
                                     HypothesesGraph::Node n,
                                     HypothesesGraph::Arc a,
                                     size_t state);

protected:
    boost::shared_ptr<Perturbation> perturbation_;
};

} // namespace pgmlink

#endif // WITH_DPCT

#endif // PERTURBEDINFERENCEMODEL_H
