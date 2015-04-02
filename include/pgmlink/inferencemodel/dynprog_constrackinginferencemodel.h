#ifndef REASONER_DYN_PROG_CONSTRACKING_H
#define REASONER_DYN_PROG_CONSTRACKING_H

#include "inferencemodel.h"
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

    virtual std::string toString() const
    {
        return std::string();
    }
    T getRef()
    {
        return ref_;
    }

private:
    T ref_;
};

typedef ConservationTrackingUserData<HypothesesGraph::Node> ConservationTrackingNodeData;
typedef ConservationTrackingUserData<HypothesesGraph::Arc> ConservationTrackingArcData;

class DynProgConsTrackInferenceModel : public InferenceModel
{
public:
    DynProgConsTrackInferenceModel(const Parameter& param);
    ~DynProgConsTrackInferenceModel();

    virtual std::vector<size_t> infer();
    virtual void conclude(HypothesesGraph&g,
                          HypothesesGraph &tracklet_graph,
                          std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > &tracklet2traxel_node_map,
                          std::vector<size_t>& solution);

    virtual void build_from_graph(const HypothesesGraph&);

    virtual void write_labeledgraph_to_file(const HypothesesGraph & g,
                            const std::string &ground_truth_filename){};

    template<class ArcIterator>
    double getTransitionArcScore(const HypothesesGraph& g, ArcIterator a);
protected:
    // dpct inference members
    dpct::Graph inference_graph_;
    std::vector<dpct::TrackingAlgorithm::Path> solution_paths_;
};

template<class ArcIterator>
double DynProgConsTrackInferenceModel::getTransitionArcScore(const HypothesesGraph& g, ArcIterator a)
{
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = g.get(node_tracklet());

    Traxel tr1, tr2;
    if (param_.with_tracklets)
    {
        tr1 = tracklet_map[g.source(a)].back();
        tr2 = tracklet_map[g.target(a)].front();
    }
    else
    {
        tr1 = traxel_map[g.source(a)];
        tr2 = traxel_map[g.target(a)];
    }

    return param_.transition(get_transition_probability(tr1, tr2, 0)) - param_.transition(get_transition_probability(tr1, tr2, 1));
}

} // namespace pgmlink

#endif // REASONER_DYN_PROG_CONSTRACKING_H

