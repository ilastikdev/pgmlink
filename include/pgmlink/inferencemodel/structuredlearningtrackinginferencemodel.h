#ifndef STRUCTUREDLEARNINGTRACKINGINFERENCEMODEL_H
#define STRUCTUREDLEARNINGTRACKINGINFERENCEMODEL_H

#include <boost/function.hpp>

#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/inference/inference.hxx>
#include <opengm/inference/lpcplex.hxx>

#include "pgmlink/hypotheses.h"
#include "pgmlink/pgm.h"
#include "pgmlink/inferencemodel/constraint_pool.hxx"
#include "pgmlink/inferencemodel/constrackinginferencemodel.h"
#include "pgmlink/reasoner_constracking.h"

#include "inferencemodel.h"

namespace pgmlink
{

/**
 * @brief The StructuredLearningTrackingInferenceModel inherits from ConsTrackingInferenceModel class which builds the OpenGM model needed to run basic conservation tracking.
 * StructuredLearningTrackingInferenceModel uses learnable functions in overrides of add_*_factor methods.
 */
  class StructuredLearningTrackingInferenceModel : public ConsTrackingInferenceModel
{
public:
    StructuredLearningTrackingInferenceModel(
        const Parameter& inferenceParam,
        double ep_gap,
        double cplex_timeout):
        ConsTrackingInferenceModel(inferenceParam, ep_gap, cplex_timeout),
        weights_((size_t)5) // { [0]detection_weight,[1]division_weight,[2]transition_weight,[3]appearance_weight,[4]disappearance_weight}

    {
        std::cout << "Constructor StructuredLearningTrackingInferenceModel" << std::endl;
    }

    opengm::learning::Weights<double> weights_;
     void setWeight ( size_t, double );
     double weight ( size_t );

    virtual size_t add_division_factors(const HypothesesGraph&, size_t);
    virtual size_t add_transition_factors(const HypothesesGraph&, size_t);
    virtual size_t add_detection_factors(const HypothesesGraph&, size_t);
};
} // namespace pgmlink

#endif // STRUCTUREDLEARNINGTRACKINGINFERENCEMODEL_H
