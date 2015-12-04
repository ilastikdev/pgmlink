#ifndef PERTURBANDMAP_PERTURBATION_H
#define PERTURBANDMAP_PERTURBATION_H

//Random distribs
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/math/distributions/normal.hpp>

#include <vector>
#include <iostream>

#include <boost/function.hpp>

#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/inference/inference.hxx>

#include "pgmlink/traxels.h"
#include "pgmlink/uncertaintyParameter.h"
#include "pgmlink/inferencemodel/inferencemodel.h"
#include "pgmlink/inferencemodel/perturbation/perturbation.h"

namespace pgmlink
{

class PerturbAndMapPerturbation : public Perturbation
{
public: // API
    PerturbAndMapPerturbation(const Perturbation::Parameter& perturbation_param, const pgmlink::Parameter& inf_param);

    virtual double generateRandomOffset(EnergyType parameterIndex,
                                        double energy = 0,
                                        Traxel tr = 0,
                                        Traxel tr2 = 0,
                                        size_t state = 0,
                                        boost::shared_ptr<InferenceModel::TransitionPredictionsMap> transition_predictions
                                            = boost::shared_ptr<InferenceModel::TransitionPredictionsMap>());
};

} // namespace pgmlink

#endif // PERTURBATION_H
