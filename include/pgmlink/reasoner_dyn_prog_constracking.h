#ifndef REASONER_DYN_PROG_CONSTRACKING_H
#define REASONER_DYN_PROG_CONSTRACKING_H

#include "reasoner_constracking.h"
#include <dpct/graph.h>
#include <dpct/trackingalgorithm.h>

namespace pgmlink
{

template<class T>
class ConservationTrackingUserData : public dpct::UserData
{
public:
    ConservationTrackingUserData(T& t):
        ref_(t)
    {}

    virtual std::string toString() const { return std::string(); }
    T getRef() { return ref_; }

private:
    T ref_;
};

typedef ConservationTrackingUserData<HypothesesGraph::Node> ConservationTrackingNodeData;
typedef ConservationTrackingUserData<HypothesesGraph::Arc> ConservationTrackingArcData;

class DynProgConservationTracking : public ConservationTracking
{
public:
    DynProgConservationTracking(const Parameter& param);
    ~DynProgConservationTracking();

    virtual void infer();
    virtual void conclude(HypothesesGraph&g);
    virtual void formulate( const HypothesesGraph& );

    // methods derived from ConservationTracking which do not apply here
    virtual void conclude(HypothesesGraph&, boost::shared_ptr<ConsTrackingInferenceModel> inference_model);
    virtual void perturbedInference(HypothesesGraph&, bool with_inference = true);

    template<class ArcIterator>
    double getTransitionArcScore(const HypothesesGraph& g, ArcIterator a);
protected:
    dpct::Graph inference_graph_;
    std::vector<dpct::TrackingAlgorithm::Path> solution_paths_;

    double get_transition_probability(Traxel &tr1, Traxel &tr2, size_t state);
    double get_transition_prob(double distance, size_t state, double alpha);
};

} // namespace pgmlink

#endif // REASONER_DYN_PROG_CONSTRACKING_H

