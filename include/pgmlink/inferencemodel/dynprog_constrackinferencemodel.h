#ifndef REASONER_DYN_PROG_CONSTRACKING_H
#define REASONER_DYN_PROG_CONSTRACKING_H

#ifdef WITH_DPCT
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

class ConservationTrackingNodeData : public ConservationTrackingUserData<HypothesesGraph::Node>
{
public:
    ConservationTrackingNodeData(HypothesesGraph::Node& n, Traxel& t):
        ConservationTrackingUserData<HypothesesGraph::Node>(n),
        traxel_(t)
    {}

    const Traxel& getTraxel() const
    {
        return traxel_;
    }

private:
    Traxel traxel_;
};

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
    virtual void fixFirstDisappearanceNodesToLabels(
            const HypothesesGraph& g,
            const HypothesesGraph &tracklet_graph,
            std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > &tracklet2traxel_map);

    template<class ArcIterator>
    double getTransitionArcScore(const HypothesesGraph& g, ArcIterator a);

    virtual double generateRandomOffset(EnergyType parameterIndex,
                                        double energy = 0,
                                        Traxel tr = 0,
                                        Traxel tr2 = 0,
                                        size_t state = 0);

    virtual double getDivMBestOffset(EnergyType energy_type,
                                     const HypothesesGraph &g,
                                     HypothesesGraph::Node n,
                                     HypothesesGraph::Arc a,
                                     size_t state);

private:
    double evaluate_motion_model(dpct::Node* a,
                                dpct::Node* b, 
                                dpct::Node* c) const;

protected:
    // dpct inference members
    dpct::Graph inference_graph_;
    std::vector<dpct::TrackingAlgorithm::Path> solution_paths_;
    std::map<HypothesesGraph::Node, dpct::Graph::NodePtr> node_reference_map_;
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

    // transProb(state=1) - transProb(state=0) would be energy or cost. score is the negative of that
    double score = param_.transition(get_transition_probability(tr1, tr2, 0))
            - param_.transition(get_transition_probability(tr1, tr2, 1));
    return score - generateRandomOffset(Transition, -score, tr1, tr2);
}

} // namespace pgmlink

#endif // WITH_DPCT

#endif // REASONER_DYN_PROG_CONSTRACKING_H

