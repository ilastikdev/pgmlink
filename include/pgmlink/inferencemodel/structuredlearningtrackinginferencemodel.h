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
#include "pgmlink/inferencemodel/constrackinginferencemodel.h"
#include "pgmlink/reasoner_constracking.h"

namespace pgmlink
{

/**
 * @brief The StructuredLearningTrackingInferenceModel inherits from ConsTrackingExplicitInferenceModel class which builds the OpenGM model needed to run basic conservation tracking.
 * StructuredLearningTrackingInferenceModel uses learnable functions in overrides of add_*_factor methods.
 */
class StructuredLearningTrackingInferenceModel : public ConsTrackingInferenceModel
{
public:
    StructuredLearningTrackingInferenceModel(
        const Parameter& inferenceParam,
        double ep_gap,
        double cplex_timeout,
        opengm::learning::Weights<double>& learningWeights,
        bool withNormalization,
        FieldOfView fov,
        double borderWidth,
        unsigned int num_threads):

        ConsTrackingInferenceModel(inferenceParam, ep_gap,cplex_timeout,num_threads),
        withNormalization_(withNormalization),
        fov_(fov),
        borderWidth_(borderWidth),
        learningWeights_(learningWeights)
    {
        cplex2_param_.verbose_ = true;
        cplex2_param_.integerConstraintNodeVar_ = true;
        cplex2_param_.epGap_ = ep_gap;
        cplex2_param_.timeLimit_ = cplex_timeout;
        cplex2_param_.numberOfThreads_ = num_threads;
    }


    virtual void fixFirstDisappearanceNodesToLabels(
            const HypothesesGraph& g,
            const HypothesesGraph &tracklet_graph,
            std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > &traxel2tracklet_map);

    virtual size_t add_division_factors(const HypothesesGraph&, size_t);
    virtual size_t add_transition_factors(const HypothesesGraph&, size_t);
    virtual size_t add_detection_factors(const HypothesesGraph&, size_t);

    virtual void add_constraints_to_pool(const HypothesesGraph& );

    template<class INF>
    void add_constraints(INF& optimizer);

    virtual IlpSolution infer();
    void set_inference_params(
            size_t numberOfSolutions,
            const std::string& feature_filename,
            const std::string& constraints_filename,
            const std::string& ground_truth_filename);

    IlpSolution extractSolution(size_t k, const std::string& ground_trugh_filename);

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
    opengm::learning::Weights<double>& getLearningWeights(){ return learningWeights_; }

protected:
    bool withNormalization_;
    size_t modelStartTime_;
    size_t modelEndTime_;
    FieldOfView fov_;
    double borderWidth_;
    opengm::learning::Weights<double>& learningWeights_;
};

template<class INF>
void StructuredLearningTrackingInferenceModel::add_constraints(INF &optimizer)
{
    linear_constraint_pool_.add_constraints_to_model(model_, optimizer);
    LOG(logDEBUG) << "StructuredLearningTrackingInferenceModel::add_constraints";
}

} // namespace pgmlink

#endif // STRUCTUREDLEARNINGTRACKINGINFERENCEMODEL_H
