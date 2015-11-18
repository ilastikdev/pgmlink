#ifdef WITH_DPCT

#include "pgmlink/inferencemodel/flow_constrackinferencemodel.h"
#include <stdexcept>
#include <memory>
#include <dpct/flowgraph.h>

namespace pgmlink
{

FlowConsTrackInferenceModel::FlowConsTrackInferenceModel(const Parameter &param):
    InferenceModel(param)
{
}

FlowConsTrackInferenceModel::~FlowConsTrackInferenceModel()
{
}

std::vector<size_t> FlowConsTrackInferenceModel::infer()
{
    LOG(logINFO) << "Starting Tracking...";
    
    inference_graph_.maxFlowMinCostTracking();

    return std::vector<size_t>();
}


void FlowConsTrackInferenceModel::build_from_graph(const HypothesesGraph& g)
{
    const HypothesesGraph *graph = &g;

    HypothesesGraph::node_timestep_map& timestep_map = graph->get(node_timestep());
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = graph->get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = graph->get(node_tracklet());

    size_t first_timestep = graph->earliest_timestep();
    size_t last_timestep = graph->latest_timestep();

    LOG(logINFO) << "Creating Flow Graph nodes";

    // add all nodes
    for (HypothesesGraph::NodeIt n(*graph); n != lemon::INVALID; ++n)
    {
        LOG(logDEBUG3) << "Adding node in timestep " << timestep_map[n] << " to FlowGraph" << std::endl;
        std::vector<float> costs;

        if(param_.with_tracklets)
            costs = tracklet_map[n].front().features["detEnergy"];
        else
            costs = traxel_map[n].features["detEnergy"];

        // for merger resolving: paths can start and end at other timeframes!
        // if(!inference_graph_.getConfig().withAppearance)
        // {
        //     if(lemon::countInArcs(*graph, n) == 0)
        //     {
        //         in_first_frame = true;
        //     }
        // }

        // if(!inference_graph_.getConfig().withDisappearance)
        // {
        //     if(lemon::countOutArcs(*graph, n) == 0)
        //     {
        //         in_last_frame = true;
        //     }
        // }

        // appearance and disappearance cost
        Traxel tr;
        if(param_.with_tracklets)
        {
            tr = tracklet_map[n].front();
        }
        else
        {
            tr = traxel_map[n];
        }

        std::vector<float> app_costs = tr.features["appEnergy"];
        std::vector<float> dis_costs = tr.features["disEnergy"];
        
        std::vector<double> costDeltas;
        std::vector<double> appearanceCostDeltas;
        std::vector<double> disappearanceCostDeltas;
        for(size_t i = 1; i < costs.size(); i++)
        {
            costDeltas.push_back(costs[i] - costs[i-1]);
            appearanceCostDeltas.push_back(app_costs[i] - app_costs[i-1]);
            disappearanceCostDeltas.push_back(dis_costs[i] - dis_costs[i-1]);
        }
        assert(costDeltas.size() == param_.max_number_objects);

        dpct::FlowGraph::Node inf_node = inference_graph_.addNode(costDeltas);
        dpct::FlowGraph::Arc inf_app = inference_graph_.addArc(inference_graph_.getSource(), inf_node, appearanceCostDeltas);
        dpct::FlowGraph::Arc inf_dis = inference_graph_.addArc(inf_node, inference_graph_.getTarget(), disappearanceCostDeltas);

        node_reference_map_[n] = inf_node;
        app_reference_map_[n] = inf_app;
        dis_reference_map_[n] = inf_dis;
    }

    LOG(logINFO) << "Creating DPCT arcs";

    // add all transition arcs
    for (HypothesesGraph::ArcIt a(*graph); a != lemon::INVALID; ++a)
    {
        dpct::FlowGraph::Node source = node_reference_map_[graph->source(a)];
        dpct::FlowGraph::Node target = node_reference_map_[graph->target(a)];

        Traxel tr1, tr2;
        if (param_.with_tracklets)
        {
            tr1 = tracklet_map[graph->source(a)].back();
            tr2 = tracklet_map[graph->target(a)].front();
        }
        else
        {
            tr1 = traxel_map[graph->source(a)];
            tr2 = traxel_map[graph->target(a)];
        }

        std::vector<float> costs = tr1.get_feature_store()->get_traxel_features(tr1, tr2)["transEnergy"];
        std::vector<double> costDeltas;
        for(size_t i = 1; i < costs.size(); i++)
        {
            costDeltas.push_back(costs[i] - costs[i-1]);
        }

        dpct::FlowGraph::Arc inf_arc = inference_graph_.addArc(source,
                                                                  target,
                                                                  costDeltas);
        arc_reference_map_[a] = inf_arc;
    }

    if(param_.with_divisions)
    {
        LOG(logINFO) << "Preparing division arcs";
        // allow division where a node has more than one output
        for (HypothesesGraph::NodeIt n(*graph); n != lemon::INVALID; ++n)
        {
            size_t number_of_outarcs = 0;
            for (HypothesesGraph::OutArcIt a(*graph, n); a != lemon::INVALID; ++a)
            {
                ++number_of_outarcs;
            }

            if (number_of_outarcs > 1)
            {
                Traxel tr;
                if (param_.with_tracklets)
                {
                    tr = tracklet_map[n].back();
                }
                else
                {
                    tr = traxel_map[n];
                }
                std::vector<float> costs = tr.features["divEnergy"];
                assert(costs.size() == 2);
                double division_cost = costs[1] - costs[0];
                // double perturb_div = generateRandomOffset(Division, division_cost, tr);
                // perturb_div += getDivMBestOffset(Division, g, n, HypothesesGraph::Arc(), size_t(true));
                LOG(logDEBUG3) << "Adding possible division from " << tr << " with score: "
                               << division_cost << "(div)" << std::endl;

                dpct::FlowGraph::Arc inf_div = inference_graph_.allowMitosis(node_reference_map_[n],
                                                division_cost);
                                              // division_cost + perturb_div);

                div_reference_map_[n] = inf_div;
            }
        }
    }

    // LOG(logINFO) << "Constructed DPCT graph with " << inference_graph_.getNumNodes() << " nodes and " << inference_graph_.getNumArcs()
    //              << " arcs on " << inference_graph_.getNumTimesteps() << " timesteps" << std::endl;
}

void FlowConsTrackInferenceModel::fixFirstDisappearanceNodesToLabels(
        const HypothesesGraph& g,
        const HypothesesGraph &tracklet_graph,
        std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > &tracklet2traxel_map)
{
    throw std::runtime_error("Not implemented");
}

double FlowConsTrackInferenceModel::generateRandomOffset(EnergyType parameterIndex,
                                                            double energy,
                                                            Traxel tr,
                                                            Traxel tr2,
                                                            size_t state)
{
    return 0.0;
}

double FlowConsTrackInferenceModel::getDivMBestOffset(EnergyType energy_type,
                                                         const HypothesesGraph &g,
                                                         HypothesesGraph::Node n,
                                                         HypothesesGraph::Arc a,
                                                         size_t state)
{
    return 0.0;
}

void FlowConsTrackInferenceModel::conclude(HypothesesGraph& g,
        HypothesesGraph &tracklet_graph,
        std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > &tracklet2traxel_node_map,
        std::vector<size_t> &solution)
{
    g.add(node_active2()).add(arc_active()).add(division_active());
    property_map<node_active2, HypothesesGraph::base_graph>::type& active_nodes = g.get(node_active2());
    property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = g.get(arc_active());
    property_map<division_active, HypothesesGraph::base_graph>::type& active_divisions = g.get(division_active());

    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());

    // add counting properties for analysis of perturbed models
    g.add(arc_active_count()).add(node_active_count()).add(division_active_count()).add(arc_value_count());
    property_map<arc_active_count, HypothesesGraph::base_graph>::type& active_arcs_count =
        g.get(arc_active_count());
    property_map<arc_value_count, HypothesesGraph::base_graph>::type& arc_values =
        g.get(arc_value_count());
    property_map<node_active_count, HypothesesGraph::base_graph>::type& active_nodes_count =
        g.get(node_active_count());
    property_map<division_active_count, HypothesesGraph::base_graph>::type& active_divisions_count =
        g.get(division_active_count());

    int iterStep = active_nodes_count[g.nodeFromId(0)].size();
    if (iterStep == 0)
    {
        //initialize vectors for storing optimizer results
        for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
        {
            active_arcs_count.set(a, std::vector<bool>());
            arc_values.set(a, std::vector<size_t>());
        }
        for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
        {
            active_nodes_count.set(n, std::vector<long unsigned int>());
            active_divisions_count.set(n, std::vector<bool>());
        }
    }

    if (!param_.with_tracklets)
    {
        tracklet_graph.add(tracklet_intern_arc_ids()).add(traxel_arc_id());
    }

    property_map<tracklet_intern_arc_ids, HypothesesGraph::base_graph>::type& tracklet_arc_id_map
            = tracklet_graph.get(tracklet_intern_arc_ids());
    property_map<traxel_arc_id, HypothesesGraph::base_graph>::type& traxel_arc_id_map
            = tracklet_graph.get(traxel_arc_id());

    // initialize nodes and divisions to 0
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        active_nodes.set(n, 0);
        active_divisions.set(n, false);
        active_nodes_count.get_value(n).push_back(0);
        active_divisions_count.get_value(n).push_back(0);
    }

    //initialize arc counts by 0
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
        active_arcs.set(a, false);
        active_arcs_count.get_value(a).push_back(0);
        arc_values.get_value(a).push_back(0);
    }


    typedef std::map<HypothesesGraph::Node, dpct::FlowGraph::Arc>::iterator NodeToFlowArcMapIt;
    typedef std::map<HypothesesGraph::Arc, dpct::FlowGraph::Arc>::iterator ArcToFlowArcMapIt;

    // detections
    for(HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        int flow = 0;
        if(inference_graph_.getFlowMap()[app_reference_map_[n]] > 0)
        {
            // appearance
            flow = inference_graph_.getFlowMap()[app_reference_map_[n]];
        }
        else if(inference_graph_.getFlowMap()[dis_reference_map_[n]] > 0)
        {
            // disappearance
            inference_graph_.getFlowMap()[dis_reference_map_[n]];
        }
        else
        {
            // standard detection usage
            flow = inference_graph_.sumInFlow(node_reference_map_[n]);
            int outFlow = inference_graph_.sumOutFlow(node_reference_map_[n]);
            assert(outFlow == flow || outFlow == flow + 1);
        }

        if(flow > 0)
        {
            if (param_.with_tracklets)
            {
                // set state of tracklet nodes
                std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map[n];

                for (std::vector<HypothesesGraph::Node>::const_iterator tr_n_it = traxel_nodes.begin();
                     tr_n_it != traxel_nodes.end();
                     ++tr_n_it)
                {
                    HypothesesGraph::Node no = *tr_n_it;
                    active_nodes.set(no, flow);
                    active_nodes_count.get_value(no)[iterStep] = active_nodes[no];
                }

                // set state of tracklet internal arcs
                std::vector<int> arc_ids = tracklet_arc_id_map[n];
                for (std::vector<int>::const_iterator arc_id_it = arc_ids.begin();
                     arc_id_it != arc_ids.end();
                     ++arc_id_it)
                {
                    HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
                    active_arcs.set(a, true);
                    active_arcs_count.get_value(a)[iterStep] = true;
                    arc_values.get_value(a)[iterStep] = flow;
                }
            }
            else
            {
                active_nodes.set(n, flow);
                active_nodes_count.get_value(n)[iterStep] = active_nodes[n];
            }
        }
    }

    // transitions
    for(ArcToFlowArcMapIt transIt = arc_reference_map_.begin(); transIt != arc_reference_map_.end(); ++transIt)
    {
        HypothesesGraph::Arc a = transIt->first;
        dpct::FlowGraph::Arc fa = transIt->second;
        int flow = inference_graph_.getFlowMap()[fa];

        if(flow > 0)
        {
            if(param_.with_tracklets)
            {
                active_arcs.set(g.arcFromId((traxel_arc_id_map[a])), true);
                active_arcs_count.get_value(g.arcFromId((traxel_arc_id_map[a])))[iterStep] = true;
                arc_values.get_value(g.arcFromId((traxel_arc_id_map[a])))[iterStep] = flow;
            }
            else
            {
                active_arcs.set(a, true);
                active_arcs_count.get_value(a)[iterStep] = true;
                arc_values.get_value(a)[iterStep] = flow;
            }
        }
    }

    // divisions
    for(NodeToFlowArcMapIt divIt = div_reference_map_.begin(); divIt != div_reference_map_.end(); ++divIt)
    {
        HypothesesGraph::Node n = divIt->first;
        dpct::FlowGraph::Arc fa = divIt->second;
        int flow = inference_graph_.getFlowMap()[fa];

        if(flow > 0)
        {
            assert(flow == 1);

            if(param_.with_tracklets)
            {
                n = tracklet2traxel_node_map[n].back();
                active_divisions.set(n, true);
                active_divisions_count.get_value(n)[iterStep] = true;
            }
            else
            {
                active_divisions.set(n, true);
                active_divisions_count.get_value(n)[iterStep] = true;
            }
        }
    }
}

} // namespace pgmlink

#endif // WITH_DPCT
