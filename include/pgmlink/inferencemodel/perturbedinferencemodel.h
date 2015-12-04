#ifndef PERTURBEDINFERENCEMODEL_H
#define PERTURBEDINFERENCEMODEL_H

#include <vector>
#include <iostream>
#include <opengm/datastructures/marray/marray.hxx>
#include "constrackinginferencemodel.h"
#include "pgmlink/uncertaintyParameter.h"
#include "pgmlink/inferencemodel/perturbation/perturbation.h"

namespace pgmlink
{

/**
 * @brief The PerturbedInferenceModel class specializes the ConsTrackingInferenceModel
 * and adds perturbation to the unaries by providing non-zero perturbations in virtual methods.
 */
class PerturbedInferenceModel : public ConsTrackingInferenceModel
{
public: // API
//<<<<<<< HEAD
//    PerturbedInferenceModel(const ConsTrackingInferenceModel::Parameter& param,
//        boost::shared_ptr<Perturbation> perturbation,
//        double ep_gap,
//        double cplex_timeout,
//        unsigned int num_threads = 1);
//=======
    PerturbedInferenceModel(const Parameter& param,
        boost::shared_ptr<Perturbation> perturbation);
//>>>>>>> c0ae1ffa3bed35ac471972fc3c7c0dcd5a44ffe7

    virtual double generateRandomOffset(EnergyType parameterIndex,
                                        double energy = 0,
                                        Traxel tr = 0,
                                        Traxel tr2 = 0,
                                        size_t state = 0);
    virtual size_t add_div_m_best_perturbation(marray::Marray<double> &energies,
            EnergyType energy_type,
            size_t factorIndex);

protected:
    boost::shared_ptr<Perturbation> perturbation_;
};

} // namespace pgmlink

#endif // PERTURBEDINFERENCEMODEL_H
