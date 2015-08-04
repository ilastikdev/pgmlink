#ifndef CLASSIFIER_UNCERTAINTY_PERTURBATION_H
#define CLASSIFIER_UNCERTAINTY_PERTURBATION_H

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

class ClassifierUncertaintyPerturbation : public Perturbation
{
public: // API
    ClassifierUncertaintyPerturbation(const Parameter& perturbation_param, const InferenceModel::Parameter& inf_param);

    virtual double generateRandomOffset(EnergyType parameterIndex,
                                        double energy = 0,
                                        Traxel tr = 0,
                                        Traxel tr2 = 0,
                                        size_t state = 0,
                                        boost::shared_ptr<InferenceModel::TransitionPredictionsMap> transition_predictions
                                            = boost::shared_ptr<InferenceModel::TransitionPredictionsMap>());

protected: // Methods
    double weightedNegLog(double prob, double weight);
    double inverseWeightedNegLog(double energy, double weight);
    double sigmoid(double x);
    double inverse_sigmoid(double x);
    double get_transition_variance(Traxel &tr1,
                                   Traxel &tr2,
                                   boost::shared_ptr<InferenceModel::TransitionPredictionsMap> transition_predictions);
    double sample_with_classifier_variance(double mean, double variance);
};

} // namespace pgmlink

#endif // PERTURBATION_H
