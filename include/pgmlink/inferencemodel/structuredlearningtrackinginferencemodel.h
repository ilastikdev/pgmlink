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
        double cplex_timeout):
        ConsTrackingExplicitInferenceModel(inferenceParam, ep_gap, cplex_timeout)
    {}

    virtual size_t add_division_factors(const HypothesesGraph&, size_t);
    virtual size_t add_transition_factors(const HypothesesGraph&, size_t);
    virtual size_t add_detection_factors(const HypothesesGraph&, size_t);
};
} // namespace pgmlink

#endif // STRUCTUREDLEARNINGTRACKINGINFERENCEMODEL_H
