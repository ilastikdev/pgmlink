#ifndef REASONER_FLOW_CONSTRACKING_H
#define REASONER_FLOW_CONSTRACKING_H

#ifdef WITH_DPCT
#include "inferencemodel.h"
#include <dpct/flowgraph.h>

namespace pgmlink
{

class FlowConsTrackInferenceModel : public InferenceModel
{
public:
    FlowConsTrackInferenceModel(const Parameter& param);
    ~FlowConsTrackInferenceModel();

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
    double getTransitionArcCost(const HypothesesGraph& g, ArcIterator a);

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

protected:
    // dpct inference members
    dpct::FlowGraph inference_graph_;
    std::map<HypothesesGraph::Node, dpct::FlowGraph::Node> node_reference_map_;
    std::map<HypothesesGraph::Node, dpct::FlowGraph::Arc> app_reference_map_;
    std::map<HypothesesGraph::Node, dpct::FlowGraph::Arc> dis_reference_map_;
    std::map<HypothesesGraph::Node, dpct::FlowGraph::Arc> div_reference_map_;
    std::map<HypothesesGraph::Arc, dpct::FlowGraph::Arc> arc_reference_map_;
};

template<class ArcIterator>
double FlowConsTrackInferenceModel::getTransitionArcCost(const HypothesesGraph& g, ArcIterator a)
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

    double cost = param_.transition(get_transition_probability(tr1, tr2, 1))
            - param_.transition(get_transition_probability(tr1, tr2, 0));
    return cost - generateRandomOffset(Transition, cost, tr1, tr2);
}

} // namespace pgmlink

#endif // WITH_DPCT

#endif // REASONER_DYN_PROG_CONSTRACKING_H

