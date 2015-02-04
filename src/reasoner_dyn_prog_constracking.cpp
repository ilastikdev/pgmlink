#include "pgmlink/reasoner_dyn_prog_constracking.h"
#include <stdexcept>
#include <memory>
#include <dpct/magnusson.h>

namespace pgmlink
{

DynProgConservationTracking::DynProgConservationTracking(const ConservationTracking::Parameter &param):
    ConservationTracking(param),
    inference_graph_(dpct::Graph::Configuration(param.with_appearance, param.with_disappearance, param.with_divisions))
{
}

DynProgConservationTracking::~DynProgConservationTracking()
{
}

void DynProgConservationTracking::infer()
{
    LOG(logINFO) << "Starting Tracking...";
    dpct::Magnusson tracker(&inference_graph_, true, true);
    double score = tracker.track(solution_paths_);
    LOG(logINFO) << "Done Tracking in " << tracker.getElapsedSeconds() << " secs with score " << score << " !";
}

// assumes detection, division, appearance and disappearance cost given as NegLnXXX functions
template<class ArcIterator>
double DynProgConservationTracking::getTransitionArcScore(const HypothesesGraph& g, ArcIterator a)
{
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = g.get(node_tracklet());

    Traxel tr1, tr2;

    if (inference_model_param_.with_tracklets)
    {
        tr1 = tracklet_map[g.source(a)].back();
        tr2 = tracklet_map[g.target(a)].front();
    }
    else
    {
        tr1 = traxel_map[g.source(a)];
        tr2 = traxel_map[g.target(a)];
    }

    return inference_model_param_.transition(get_transition_probability(tr1, tr2, 0)) - inference_model_param_.transition(get_transition_probability(tr1, tr2, 1));
}

void DynProgConservationTracking::formulate(const HypothesesGraph& g)
{
    const HypothesesGraph *graph;

    // for formulate, add_constraints, add_finite_factors: distinguish graph & tracklet_graph
    if (with_tracklets_)
    {
        LOG(logINFO) << "ConservationTracking::perturbedInference: generating tracklet graph";
        tracklet2traxel_node_map_ = generateTrackletGraph2(g, tracklet_graph_);
        graph = &tracklet_graph_;
    }
    else
    {
        graph = &g;
    }

    std::map<HypothesesGraph::Node, dpct::Graph::NodePtr> node_reference_map;
    HypothesesGraph::node_timestep_map& timestep_map = graph->get(node_timestep());
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = graph->get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map = graph->get(node_tracklet());

    size_t first_timestep = graph->earliest_timestep();
    size_t last_timestep = graph->latest_timestep();

    LOG(logINFO) << "Creating DPCT nodes";

    // add all nodes
    for (HypothesesGraph::NodeIt n(*graph); n != lemon::INVALID; ++n)
    {
        LOG(logDEBUG3) << "Adding node in timestep " << timestep_map[n] << " to DPCT" << std::endl;
        std::vector<double> scoreDeltas;

        for(size_t state = 0; state <= inference_model_param_.max_number_objects; state++)
        {
            double energy = 0.0;
            if (inference_model_param_.with_tracklets)
            {
                // add all detection factors of the internal nodes
                for (std::vector<Traxel>::const_iterator trax_it = tracklet_map[n].begin(); trax_it != tracklet_map[n].end(); ++trax_it)
                {
                    energy += inference_model_param_.detection(*trax_it, state);
                }

                // add all transition factors of the internal arcs
                Traxel tr_prev;
                bool first = true;
                for (std::vector<Traxel>::const_iterator trax_it = tracklet_map[n].begin(); trax_it != tracklet_map[n].end(); ++trax_it)
                {
                    LOG(logDEBUG4) << "internal arcs traxel " << *trax_it;
                    Traxel tr = *trax_it;
                    if (!first)
                    {
                        energy += inference_model_param_.transition( get_transition_probability(tr_prev, tr, state) );
                    }
                    else
                    {
                        first = false;
                    }
                    tr_prev = tr;
                }
            }
            else
            {
                // only look at this single traxel
                energy = inference_model_param_.detection(traxel_map[n], state);
            }

            scoreDeltas.push_back(-1.0 * energy);

            LOG(logDEBUG3) << "\tstate " << state << " has score " << scoreDeltas.back() << std::endl;
        }

        // add source and sink links
        int node_begin_time = -1;
        int node_end_time = -1;
        if (inference_model_param_.with_tracklets)
        {
            node_begin_time = tracklet_map[n].front().Timestep;
            node_end_time = tracklet_map[n].back().Timestep;
        }
        else
        {
            node_begin_time = traxel_map[n].Timestep;
            node_end_time = traxel_map[n].Timestep;
        }

        bool in_first_frame = node_begin_time == first_timestep;
        bool in_last_frame = node_end_time == last_timestep;

        // appearance and disappearance score
        Traxel tr;
        if(inference_model_param_.with_tracklets)
        {
            tr = tracklet_map[n].front();
        }
        else
        {
            tr = traxel_map[n];
        }

        double app_score = -1.0 * inference_model_param_.appearance_cost(tr);

        if(inference_model_param_.with_tracklets)
        {
            tr = tracklet_map[n].back();
        }
        double dis_score = -1.0 * inference_model_param_.disappearance_cost(tr);
        LOG(logDEBUG3) << "\tapp-score " << app_score << std::endl;
        LOG(logDEBUG3) << "\tdis-score " << dis_score << std::endl;

        dpct::Graph::NodePtr inf_node = inference_graph_.addNode(timestep_map[n] - first_timestep, scoreDeltas, app_score, dis_score, in_first_frame, in_last_frame, std::make_shared<ConservationTrackingNodeData>(n));
        node_reference_map[n] = inf_node;
    }

    LOG(logINFO) << "Creating DPCT arcs";

    // add all transition arcs
    for (HypothesesGraph::ArcIt a(*graph); a != lemon::INVALID; ++a)
    {
        dpct::Graph::NodePtr source = node_reference_map[graph->source(a)];
        dpct::Graph::NodePtr target = node_reference_map[graph->target(a)];

        double score = getTransitionArcScore(*graph, a);

        dpct::Graph::ArcPtr inf_arc = inference_graph_.addMoveArc(source, target, score, std::make_shared<ConservationTrackingArcData>(a));
    }

    if(inference_model_param_.with_divisions)
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
                if (inference_model_param_.with_tracklets)
                {
                    tr = tracklet_map[n].back();
                }
                else
                {
                    tr = traxel_map[n];
                }
                double division_score = inference_model_param_.division(tr, 0) - inference_model_param_.division(tr, 1);

                for (HypothesesGraph::OutArcIt a(*graph, n); a != lemon::INVALID; ++a)
                {
                    // division arc score = division score + move score
                    double move_score = getTransitionArcScore(*graph, a);
                    LOG(logDEBUG3) << "Adding possible division from " << tr << " with score: " << move_score << "(move) + " << division_score << "(div)" << std::endl;
                    inference_graph_.allowMitosis(node_reference_map[n], node_reference_map[graph->target(a)], move_score + division_score);
                }
            }
        }
    }

    LOG(logINFO) << "Constructed DPCT graph with " << inference_graph_.getNumNodes() << " nodes and " << inference_graph_.getNumArcs()
                 << " arcs on " << inference_graph_.getNumTimesteps() << " timesteps" << std::endl;
}

// copied from ConsTrackingInferenceModel, but without the python TransitionClassifier
double DynProgConservationTracking::get_transition_probability(Traxel& tr1, Traxel& tr2, size_t state) {
    LOG(logDEBUG4) << "get_transition_probability()";

    double prob;

    //read the FeatureMaps from Traxels
    double distance = 0;
    if (inference_model_param_.with_optical_correction)
    {
        distance = tr1.distance_to_corr(tr2);
    }
    else
    {
        distance = tr1.distance_to(tr2);
    }
    prob = get_transition_prob(distance, state, inference_model_param_.transition_parameter);
    LOG(logDEBUG4) << "get_transition_probability(): using deterministic function: " << tr1
                   << " " << tr2 << " [" << state << "] = " << prob << "; distance = " << distance;
    assert(prob >= 0 && prob <= 1);
    return prob;
}

// copied from ConsTrackingInferenceModel
double DynProgConservationTracking::get_transition_prob(double distance, size_t state, double alpha) {
    double prob = exp(-distance / alpha);
    if (state == 0)
    {
        return 1 - prob;
    }
    return prob;
}

void DynProgConservationTracking::conclude(HypothesesGraph& g)
{
    g.add(node_active2()).add(arc_active()).add(division_active());
    property_map<node_active2, HypothesesGraph::base_graph>::type& active_nodes = g.get(node_active2());
    property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = g.get(arc_active());
    property_map<division_active, HypothesesGraph::base_graph>::type& active_divisions = g.get(division_active());

    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());

    if (!with_tracklets_)
        tracklet_graph_.add(tracklet_intern_arc_ids()).add(traxel_arc_id());

    property_map<tracklet_intern_arc_ids, HypothesesGraph::base_graph>::type& tracklet_arc_id_map = tracklet_graph_.get(tracklet_intern_arc_ids());
    property_map<traxel_arc_id, HypothesesGraph::base_graph>::type& traxel_arc_id_map = tracklet_graph_.get(traxel_arc_id());

    // initialize nodes and divisions to 0
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        active_nodes.set(n, 0);
        active_divisions.set(n, false);
    }

    //initialize arc counts by 0
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
        active_arcs.set(a, false);
    }

    // small function used to increase number of objects per node
    std::function<void(dpct::Node*)> increase_object_count = [&](dpct::Node* node)
    {
        std::shared_ptr<ConservationTrackingNodeData> nd = std::static_pointer_cast<ConservationTrackingNodeData>(node->getUserData());
        HypothesesGraph::Node n = nd->getRef();
        LOG(logDEBUG3) << "increasing use count of " << traxel_map[n] << std::endl;

        if (inference_model_param_.with_tracklets)
        {
            // set state of tracklet nodes
            std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map_[n];

            for (std::vector<HypothesesGraph::Node>::const_iterator tr_n_it = traxel_nodes.begin(); tr_n_it != traxel_nodes.end(); ++tr_n_it)
            {
                HypothesesGraph::Node no = *tr_n_it;
                active_nodes.set(no, active_nodes[no] + 1);
            }

            // set state of tracklet internal arcs
            std::vector<int> arc_ids = tracklet_arc_id_map[n];
            for (std::vector<int>::const_iterator arc_id_it = arc_ids.begin(); arc_id_it != arc_ids.end(); ++arc_id_it)
            {
                HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
//                assert(active_arcs[a] == false);
                active_arcs.set(a, true);
            }
        }
        else
        {
            active_nodes.set(n, active_nodes[n] + 1);
        }
    };

    std::function<void(dpct::Arc*)> activate_arc = [&](dpct::Arc* arc)
    {
        std::shared_ptr<ConservationTrackingArcData> ad = std::static_pointer_cast<ConservationTrackingArcData>(arc->getUserData());
        HypothesesGraph::Arc a = ad->getRef();

        if(inference_model_param_.with_tracklets)
        {
            active_arcs.set(g.arcFromId((traxel_arc_id_map[a])), true);
        }
        else
        {
            active_arcs.set(a, true);
        }
    };

    size_t num_paths = 0;
    // for each path, increment the number of cells the nodes and arcs along the path
    for(dpct::TrackingAlgorithm::Path& p : solution_paths_)
    {
        std::cout << "\rLooking at path " << num_paths++  << " of length " << p.size() << std::flush;
        // a path starts at the dummy-source and goes to the dummy-sink. these arcs are of type dummy, and thus skipped
        bool first_arc_on_path = true;
        for(dpct::Arc* a: p)
        {
            assert(a != nullptr);
            assert(a->getType() != dpct::Arc::Swap);
            assert(a->getSourceNode() != nullptr);
            assert(a->getTargetNode() != nullptr);

            switch(a->getType())
            {
                case dpct::Arc::Move:
                {
                    // send one cell through the nodes
                    if(first_arc_on_path)
                    {
                        increase_object_count(a->getSourceNode());
                        first_arc_on_path = false;
                    }
                    increase_object_count(a->getTargetNode());

                    // set arc to active
                    activate_arc(a);
                } break;
                case dpct::Arc::Appearance:
                {
                    // the node that appeared is set active here, so detections without further path are active as well
                    increase_object_count(a->getTargetNode());
                    first_arc_on_path = false;
                } break;
                case dpct::Arc::Disappearance:
                {
                    if(first_arc_on_path)
                    {
                        increase_object_count(a->getSourceNode());
                    }
                    first_arc_on_path = false;
                    // nothing to do, last node on path was already set active by previous move or appearance
                }  break;
                case dpct::Arc::Division:
                {
                    // set as active division
                    std::shared_ptr<ConservationTrackingNodeData> nd = std::static_pointer_cast<ConservationTrackingNodeData>(a->getObservedNode()->getUserData());
                    HypothesesGraph::Node parent = nd->getRef();
                    nd = std::static_pointer_cast<ConservationTrackingNodeData>(a->getTargetNode()->getUserData());
                    HypothesesGraph::Node child = nd->getRef();

                    if(inference_model_param_.with_tracklets)
                    {
                        parent = tracklet2traxel_node_map_[parent].back();
                        child = tracklet2traxel_node_map_[child].front();
                    }

                    LOG(logINFO) << "activating division for " << traxel_map[parent] << std::endl;
                    active_divisions.set(parent, true);
                    increase_object_count(a->getTargetNode());
                    first_arc_on_path = false;

                    // activate corresponding arc, which is not active in dpct!
                    for (HypothesesGraph::OutArcIt oa(g, parent); oa != lemon::INVALID; ++oa)
                    {
                        if(g.target(oa) == child)
                        {
                            LOG(logINFO) << "Found arc to activate for division!" << std::endl;
                            active_arcs.set(oa, true);
                            break;
                        }
                    }
                } break;
                case dpct::Arc::Swap:
                {
                    throw std::runtime_error("Got a swap arc even though it should have been cleaned up!");
                } break;
                case dpct::Arc::Dummy:
                {
                    // do nothing
                } break;
                default:
                {
                    throw std::runtime_error("Unkown arc type");
                } break;
            }
        }
    }
}

void DynProgConservationTracking::conclude(HypothesesGraph &, boost::shared_ptr<ConsTrackingInferenceModel> inference_model)
{
    throw std::runtime_error("Not implemented - should not be used");
}

void DynProgConservationTracking::perturbedInference(HypothesesGraph &, bool with_inference)
{
    throw std::runtime_error("Not implemented - should not be used");
}


} // namespace pgmlink
