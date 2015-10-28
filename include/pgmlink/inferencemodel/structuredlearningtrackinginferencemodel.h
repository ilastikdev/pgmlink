#ifndef STRUCTUREDLEARNINGTRACKINGINFERENCEMODEL_H
#define STRUCTUREDLEARNINGTRACKINGINFERENCEMODEL_H

#include <boost/function.hpp>

#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/inference/inference.hxx>
#include <opengm/inference/lpcplex2.hxx>

#include "pgmlink/hypotheses.h"
#include "pgmlink/pgm.h"
#include "pgmlink/inferencemodel/constraint_pool.hxx"
#include "pgmlink/inferencemodel/inferencemodel.h"
#include "pgmlink/inferencemodel/constrackingexplicitinferencemodel.h"
#include "pgmlink/reasoner_constracking.h"

namespace pgmlink
{

/**
 * @brief The StructuredLearningTrackingInferenceModel inherits from ConsTrackingExplicitInferenceModel class which builds the OpenGM model needed to run basic conservation tracking.
 * StructuredLearningTrackingInferenceModel uses learnable functions in overrides of add_*_factor methods.
 */
  class StructuredLearningTrackingInferenceModel : public ConsTrackingExplicitInferenceModel
{
public:
    StructuredLearningTrackingInferenceModel(
        const Parameter& inferenceParam,
        double ep_gap,
        double cplex_timeout,
        opengm::learning::Weights<double>& inferenceWeights,
        bool withNormalization):
        ConsTrackingExplicitInferenceModel(
            inferenceParam,
            ep_gap,
            cplex_timeout,
            inferenceWeights),
        withNormalization_(withNormalization)
    {}

    virtual size_t add_division_factors(const HypothesesGraph&, size_t);
    virtual size_t add_transition_factors(const HypothesesGraph&, size_t);
    virtual size_t add_detection_factors(const HypothesesGraph&, size_t);

    void setModelStartTime(size_t time){
        modelStartTime_ = time;
    }

    void setModelEndTime(size_t time){
        modelEndTime_ = time;
    }

    size_t getModelStartTime(){
        return modelStartTime_;
    }

    size_t getModelEndTime(){
        return modelEndTime_;
    }

protected:
    bool withNormalization_;
    size_t modelStartTime_;
    size_t modelEndTime_;
};
} // namespace pgmlink

#endif // STRUCTUREDLEARNINGTRACKINGINFERENCEMODEL_H
