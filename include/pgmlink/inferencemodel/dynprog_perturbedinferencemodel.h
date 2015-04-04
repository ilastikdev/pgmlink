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
                                   const Perturbation::Parameter& perturbation_param);

    virtual double generateRandomOffset(EnergyType parameterIndex,
                                        double energy = 0,
                                        Traxel tr = 0,
                                        Traxel tr2 = 0,
                                        size_t state = 0);

//    virtual size_t add_div_m_best_perturbation(marray::Marray<double> &energies,
//            EnergyType energy_type,
//            size_t factorIndex);

//    void perturb(Perturbation::DeterministicOffset* det_off);

protected:
    boost::shared_ptr<Perturbation> perturbation_;
};

} // namespace pgmlink

#endif // WITH_DPCT

#endif // PERTURBEDINFERENCEMODEL_H
