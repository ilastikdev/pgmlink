#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <string.h>
#include <memory.h>
#include <opengm/inference/lpcplex.hxx>
#include <opengm/datastructures/marray/marray.hxx>
#include <opengm/graphicalmodel/graphicalmodel_hdf5.hxx>

#include <boost/date_time/posix_time/posix_time.hpp>

#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/traxels.h"
#include "pgmlink/constraint_pool.hxx"

using namespace std;

namespace pgmlink {
ConservationTracking::~ConservationTracking() {
    //	if (pgm_ != NULL) {
    //		delete pgm_;
    //		pgm_ = NULL;
    //	}
    //   if (optimizer_ != NULL) {
    //      delete optimizer_;
    //      optimizer_ = NULL;
    //   }
}

double ConservationTracking::forbidden_cost() const {
    return forbidden_cost_;
}

void ConservationTracking::formulate(const HypothesesGraph& hypotheses) {
    LOG(logDEBUG) << "ConservationTracking::formulate: entered";
    reset();
    pgm_ = boost::shared_ptr < pgm::OpengmModelDeprecated > (new pgm::OpengmModelDeprecated());

    HypothesesGraph const *graph;
    if (with_tracklets_) {
        LOG(logINFO) << "ConservationTracking::formulate: generating tracklet graph";
        tracklet2traxel_node_map_ = generateTrackletGraph2(hypotheses, tracklet_graph_);
        graph = &tracklet_graph_;
    } else {
        graph = &hypotheses;
    }

    LOG(logDEBUG) << "ConservationTracking::formulate: add_transition_nodes";
    add_transition_nodes(*graph);
    LOG(logDEBUG) << "ConservationTracking::formulate: add_appearance_nodes";
    add_appearance_nodes(*graph);
    LOG(logDEBUG) << "ConservationTracking::formulate: add_disappearance_nodes";
    add_disappearance_nodes(*graph);

    LOG(logDEBUG) << "ConservationTracking::formulate: add_division_nodes";
    if (with_divisions_) {
        add_division_nodes(*graph);
    }
    pgm::OpengmModelDeprecated::ogmGraphicalModel* model = pgm_->Model();

    LOG(logDEBUG) << "ConservationTracking::formulate: add_finite_factors";
    add_finite_factors(*graph);
    LOG(logDEBUG) << "ConservationTracking::formulate: finished add_finite_factors";
    typedef opengm::LPCplex<pgm::OpengmModelDeprecated::ogmGraphicalModel,
            pgm::OpengmModelDeprecated::ogmAccumulator> cplex_optimizer;
    cplex_optimizer::Parameter param;
    param.verbose_ = true;
    param.integerConstraint_ = true;
    param.epGap_ = ep_gap_;
    param.timeLimit_ = cplex_timeout_;
    LOG(logDEBUG) << "ConservationTracking::formulate ep_gap = " << param.epGap_;

    optimizer_ = new cplex_optimizer(*model, param);

    LOG(logDEBUG) << "ConservationTracking::formulate: add_constraints";
    if (with_constraints_) {
        add_constraints(*graph);
    }

    LOG(logINFO) << "number_of_transition_nodes_ = " << number_of_transition_nodes_;
    LOG(logINFO) << "number_of_appearance_nodes_ = " << number_of_appearance_nodes_;
    LOG(logINFO) << "number_of_disappearance_nodes_ = " << number_of_disappearance_nodes_;
    LOG(logINFO) << "number_of_division_nodes_ = " << number_of_division_nodes_;

}

void ConservationTracking::infer() {
    if (!with_constraints_) {
        opengm::hdf5::save(optimizer_->graphicalModel(), "./conservationTracking.h5", "conservationTracking");
        throw std::runtime_error("GraphicalModel::infer(): inference with soft constraints is not implemented yet. The conservation tracking factor graph has been saved to file");
    }
    opengm::InferenceTermination status = optimizer_->infer();
    if (status != opengm::NORMAL) {
        throw std::runtime_error("GraphicalModel::infer(): optimizer terminated abnormally");
    }
}

void ConservationTracking::conclude(HypothesesGraph& g) {
    // extract solution from optimizer
    vector<pgm::OpengmModelDeprecated::ogmInference::LabelType> solution;
    opengm::InferenceTermination status = optimizer_->arg(solution);
    if (status != opengm::NORMAL) {
        throw runtime_error("GraphicalModel::infer(): solution extraction terminated abnormally");
    }

    // add 'active' properties to graph
    g.add(node_active2()).add(arc_active()).add(division_active());

    property_map<node_active2, HypothesesGraph::base_graph>::type& active_nodes =
            g.get(node_active2());
    property_map<arc_active, HypothesesGraph::base_graph>::type& active_arcs = g.get(arc_active());
    property_map<division_active, HypothesesGraph::base_graph>::type& division_nodes =
            g.get(division_active());
    if (!with_tracklets_) {
        tracklet_graph_.add(tracklet_intern_arc_ids()).add(traxel_arc_id());
    }
    property_map<tracklet_intern_arc_ids, HypothesesGraph::base_graph>::type& tracklet_arc_id_map =
            tracklet_graph_.get(tracklet_intern_arc_ids());
    property_map<traxel_arc_id, HypothesesGraph::base_graph>::type& traxel_arc_id_map =
            tracklet_graph_.get(traxel_arc_id());

    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
        active_arcs.set(a, false);
    }

    // write state after inference into 'active'-property maps
    // the node is also active if its appearance node is active
    for (std::map<HypothesesGraph::Node, size_t>::const_iterator it = app_node_map_.begin();
            it != app_node_map_.end(); ++it) {
        if (with_tracklets_) {
            // set state of tracklet nodes
            std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map_[it->first];
            for (std::vector<HypothesesGraph::Node>::const_iterator tr_n_it = traxel_nodes.begin();
                    tr_n_it != traxel_nodes.end(); ++tr_n_it) {
                HypothesesGraph::Node n = *tr_n_it;
                active_nodes.set(n, solution[it->second]);
            }
            // set state of tracklet internal arcs
            std::vector<int> arc_ids = tracklet_arc_id_map[it->first];
            for (std::vector<int>::const_iterator arc_id_it = arc_ids.begin();
                    arc_id_it != arc_ids.end(); ++arc_id_it) {
                HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
                assert(active_arcs[a] == false);
                if (solution[it->second] > 0) {
                    active_arcs.set(a, true);
                    assert(active_nodes[g.source(a)] == solution[it->second]
                            && "tracklet internal arcs must have the same flow as their connected nodes");
                    assert(active_nodes[g.target(a)] == solution[it->second]
                            && "tracklet internal arcs must have the same flow as their connected nodes");
                }
            }
        } else {
            active_nodes.set(it->first, solution[it->second]);
        }
    }

    // the node is also active if its disappearance node is active
    for (std::map<HypothesesGraph::Node, size_t>::const_iterator it = dis_node_map_.begin();
            it != dis_node_map_.end(); ++it) {

        if (solution[it->second] > 0) {
            if (with_tracklets_) {
                // set state of tracklet nodes
                std::vector<HypothesesGraph::Node> traxel_nodes = tracklet2traxel_node_map_[it->first];
                for (std::vector<HypothesesGraph::Node>::const_iterator tr_n_it =
                        traxel_nodes.begin(); tr_n_it != traxel_nodes.end(); ++tr_n_it) {
                    HypothesesGraph::Node n = *tr_n_it;
                    if (active_nodes[n] == 0) {
                        active_nodes.set(n, solution[it->second]);
                    } else {
                        assert(active_nodes[n] == solution[it->second]);
                    }
                }
                // set state of tracklet internal arcs
                std::vector<int> arc_ids = tracklet_arc_id_map[it->first];
                for (std::vector<int>::const_iterator arc_id_it = arc_ids.begin();
                        arc_id_it != arc_ids.end(); ++arc_id_it) {
                    HypothesesGraph::Arc a = g.arcFromId(*arc_id_it);
                    if (solution[it->second] > 0) {
                        active_arcs.set(a, true);
                        assert(active_nodes[g.source(a)] == solution[it->second]
                                && "tracklet internal arcs must have the same flow as their connected nodes");
                        assert(active_nodes[g.target(a)] == solution[it->second]
                                && "tracklet internal arcs must have the same flow as their connected nodes");
                    }
                }
            } else {
                if (active_nodes[it->first] == 0) {
                    active_nodes.set(it->first, solution[it->second]);
                } else {
                    assert(active_nodes[it->first] == solution[it->second]);
                }
            }
        }
    }

    for (std::map<HypothesesGraph::Arc, size_t>::const_iterator it = arc_map_.begin();
            it != arc_map_.end(); ++it) {
        if (solution[it->second] >= 1) {
            if (with_tracklets_) {
                active_arcs.set(g.arcFromId((traxel_arc_id_map[it->first])), true);
            } else {
                active_arcs.set(it->first, true);
            }
        }
    }
    // initialize division node map
    if (with_divisions_) {
        for (std::map<HypothesesGraph::Node, size_t>::const_iterator it = div_node_map_.begin();
                it != div_node_map_.end(); ++it) {
            division_nodes.set(it->first, false);
        }
        for (std::map<HypothesesGraph::Node, size_t>::const_iterator it = div_node_map_.begin();
                it != div_node_map_.end(); ++it) {
            if (solution[it->second] >= 1) {
                if (with_tracklets_) {
                    // set division property for the last node in the tracklet
                    division_nodes.set(tracklet2traxel_node_map_[it->first].back(), true);
                } else {
                    division_nodes.set(it->first, true);
                }
            }
        }
    }
}

const std::map<HypothesesGraph::Arc, size_t>& ConservationTracking::get_arc_map() const {
    return arc_map_;
}

void ConservationTracking::reset() {
    if (optimizer_ != NULL) {
        delete optimizer_;
        optimizer_ = NULL;
    }
    arc_map_.clear();
    div_node_map_.clear();
    app_node_map_.clear();
    dis_node_map_.clear();
}

void ConservationTracking::add_appearance_nodes(const HypothesesGraph& g) {
    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        pgm_->Model()->addVariable(max_number_objects_ + 1);
        app_node_map_[n] = pgm_->Model()->numberOfVariables() - 1;
        assert(pgm_->Model()->numberOfLabels(app_node_map_[n]) == max_number_objects_ + 1);
        ++count;
    }
    number_of_appearance_nodes_ = count;
}

void ConservationTracking::add_disappearance_nodes(const HypothesesGraph& g) {
    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        pgm_->Model()->addVariable(max_number_objects_ + 1);
        dis_node_map_[n] = pgm_->Model()->numberOfVariables() - 1;
        assert(pgm_->Model()->numberOfLabels(dis_node_map_[n]) == max_number_objects_ + 1);
        ++count;
    }
    number_of_disappearance_nodes_ = count;
}

void ConservationTracking::add_transition_nodes(const HypothesesGraph& g) {
    size_t count = 0;
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
        pgm_->Model()->addVariable(max_number_objects_ + 1);
        arc_map_[a] = pgm_->Model()->numberOfVariables() - 1;
        assert(pgm_->Model()->numberOfLabels(arc_map_[a]) == max_number_objects_ + 1);
        ++count;
    }
    number_of_transition_nodes_ = count;
}

void ConservationTracking::add_division_nodes(const HypothesesGraph& g) {
    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        size_t number_of_outarcs = 0;
        for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
            ++number_of_outarcs;
        }
        if (number_of_outarcs > 1) {
            pgm_->Model()->addVariable(2);
            div_node_map_[n] = pgm_->Model()->numberOfVariables() - 1;
            assert(pgm_->Model()->numberOfLabels(div_node_map_[n]) == 2);
            ++count;
        }
    }
    number_of_division_nodes_ = count;
}

namespace {
double get_transition_prob(double distance, size_t state, double alpha) {
    double prob = exp(-distance / alpha);
    if (state == 0) {
        return 1 - prob;
    }
    return prob;
}
}

void ConservationTracking::add_finite_factors(const HypothesesGraph& g) {
    LOG(logDEBUG) << "ConservationTracking::add_finite_factors: entered";
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    property_map<node_tracklet, HypothesesGraph::base_graph>::type& tracklet_map =
            g.get(node_tracklet());
    property_map<tracklet_intern_dist, HypothesesGraph::base_graph>::type& tracklet_intern_dist_map =
            g.get(tracklet_intern_dist());

    ////
    //// add detection factors
    ////
    LOG(logDEBUG) << "ConservationTracking::add_finite_factors: add detection factors";
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
        size_t num_vars = 0;
        vector<size_t> vi;
        vector<double> cost;

        int node_begin_time = -1;
        int node_end_time = -1;
        if (with_tracklets_) {
            node_begin_time = tracklet_map[n].front().Timestep;
            node_end_time = tracklet_map[n].back().Timestep;
        } else {
            node_begin_time = traxel_map[n].Timestep;
            node_end_time = traxel_map[n].Timestep;
        }

        if (app_node_map_.count(n) > 0) {
            vi.push_back(app_node_map_[n]);
            if (node_begin_time <= g.earliest_timestep()) {  // "<" holds if there are only tracklets in the first frame
                // pay no appearance costs in the first timestep
                cost.push_back(0.);
            } else {
                if (with_tracklets_) {
                    cost.push_back(appearance_cost_(tracklet_map[n].front()));
                    LOG(logDEBUG4) << "App-costs 1: " << appearance_cost_(tracklet_map[n].front()) << ", " << tracklet_map[n].front();
                } else {
                    cost.push_back(appearance_cost_(traxel_map[n]));
                    LOG(logDEBUG4) << "App-costs 2: " << appearance_cost_(traxel_map[n]) << ", " << traxel_map[n];
                }
            }
            ++num_vars;
        }
        if (dis_node_map_.count(n) > 0) {
            vi.push_back(dis_node_map_[n]);
            double c = 0;
            if (node_end_time < g.latest_timestep()) { // "<" holds if there are only tracklets in the last frame
                if (with_tracklets_) {
                    c += disappearance_cost_(tracklet_map[n].back());
                    LOG(logDEBUG4) << "Disapp-costs 1: " << disappearance_cost_(tracklet_map[n].back()) << ", " << tracklet_map[n].back();
                } else {
                    c += disappearance_cost_(traxel_map[n]);
                    LOG(logDEBUG4) << "Disapp-costs 2: " << disappearance_cost_(traxel_map[n]) << ", " << traxel_map[n];
                }
            }
            cost.push_back(c);
            ++num_vars;
        }

        // convert vector to array
        vector<size_t> coords(num_vars, 0); // number of variables
        // ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
        pgm::OpengmExplicitFactor<double> table(vi.begin(), vi.end(), forbidden_cost_, (max_number_objects_ + 1));
        for (size_t state = 0; state <= max_number_objects_; ++state) {
            double energy = 0;
            if (with_tracklets_) {
                // add all detection factors of the internal nodes
                for (std::vector<Traxel>::const_iterator trax_it = tracklet_map[n].begin();
                        trax_it != tracklet_map[n].end(); ++trax_it) {
                    energy += detection_(*trax_it, state);
                }
                // add all transition factors of the internal arcs
                for (std::vector<double>::const_iterator intern_dist_it =
                        tracklet_intern_dist_map[n].begin();
                        intern_dist_it != tracklet_intern_dist_map[n].end(); ++intern_dist_it) {
                    energy += transition_(
                            get_transition_prob(*intern_dist_it, state, transition_parameter_));
                }
            } else {
                energy = detection_(traxel_map[n], state);
            }
            LOG(logDEBUG2) << "ConservationTracking::add_finite_factors: detection[" << state
                    << "] = " << energy;
            for (size_t var_idx = 0; var_idx < num_vars; ++var_idx) {
                coords[var_idx] = state;
                // if only one of the variables is > 0, then it is an appearance in this time frame
                // or a disappearance in the next timeframe. Hence, add the cost of appearance/disappearance
                // to the detection cost
                table.set_value(coords, energy + state * cost[var_idx]);
                coords[var_idx] = 0;
                LOG(logDEBUG4) << "ConservationTracking::add_finite_factors: var_idx "
                        << var_idx << " = " << energy;
            }
            // also this energy if both variables have the same state
            if (num_vars == 2) {
                coords[0] = state;
                coords[1] = state;
                // only pay detection energy if both variables are on
                table.set_value(coords, energy);
                coords[0] = 0;
                coords[1] = 0;

                LOG(logDEBUG4) << "ConservationTracking::add_finite_factors: var_idxs 0 and var_idx 1 = "
                        << energy;
            }
        }

        LOG(logDEBUG3) << "ConservationTracking::add_finite_factors: adding table to pgm";
        table.add_to(*(pgm_->Model()));
    }

    ////
    //// add transition factors
    ////
    LOG(logDEBUG) << "ConservationTracking::add_finite_factors: add transition factors";
    property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_distances = g.get(
            arc_distance());
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a) {
        size_t vi[] = { arc_map_[a] };
        vector<size_t> coords(1, 0); // number of variables
        // ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
        pgm::OpengmExplicitFactor<double> table(vi, vi + 1, forbidden_cost_, (max_number_objects_ + 1));
        for (size_t state = 0; state <= max_number_objects_; ++state) {
            double energy = transition_(get_transition_prob(arc_distances[a], state, transition_parameter_));
            LOG(logDEBUG2) << "ConservationTracking::add_finite_factors: transition[" << state
                    << "] = " << energy;
            coords[0] = state;
            table.set_value(coords, energy);
            coords[0] = 0;
        }
        table.add_to(*(pgm_->Model()));
    }

    ////
    //// add division factors
    ////
    if (with_divisions_) {
        LOG(logDEBUG) << "ConservationTracking::add_finite_factors: add division factors";
        for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n) {
            if (div_node_map_.count(n) == 0) {
                continue;
            }
            size_t vi[] = { div_node_map_[n] };
            vector<size_t> coords(1, 0); // number of variables
            // ITER first_ogm_idx, ITER last_ogm_idx, VALUE init, size_t states_per_var
            pgm::OpengmExplicitFactor<double> table(vi, vi + 1, forbidden_cost_, 2);
            for (size_t state = 0; state <= 1; ++state) {
                double energy = 0;
                if (with_tracklets_) {
                    energy = division_(tracklet_map[n].back(), state);
                } else {
                    energy = division_(traxel_map[n], state);
                }
                LOG(logDEBUG2) << "ConservationTracking::add_finite_factors: division[" << state
                        << "] = " << energy;
                coords[0] = state;
                table.set_value(coords, energy);
                coords[0] = 0;
            }
            table.add_to(*(pgm_->Model()));
        }
    }
}

size_t ConservationTracking::cplex_id(size_t opengm_id, size_t state) {
    return optimizer_->lpNodeVi(opengm_id, state);
}

void ConservationTracking::add_constraints(const HypothesesGraph& g) {
    LOG(logDEBUG) << "ConservationTracking::add_constraints: entered";

    pgm::ConstraintPool constraint_pool(forbidden_cost_, with_divisions_, with_appearance_, with_disappearance_, with_misdetections_allowed_);

    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        ////
        //// outgoing transitions
        ////
        {
            std::vector<size_t> transition_nodes;
            for (HypothesesGraph::OutArcIt a(g, n); a != lemon::INVALID; ++a) {
                transition_nodes.push_back(arc_map_[a]);
            }

            int division_node = -1;
            if(div_node_map_.count(n) > 0)
            {
                division_node = div_node_map_[n];
            }
            size_t appearance_node = app_node_map_[n];

            constraint_pool.add_constraint(pgm::ConstraintPool::OutgoingConstraint(appearance_node, division_node, transition_nodes));
        }


        ////
        //// incoming transitions
        ////
        {
            std::vector<size_t> transition_nodes;
            for (HypothesesGraph::InArcIt a(g, n); a != lemon::INVALID; ++a) {
                transition_nodes.push_back(arc_map_[a]);
            }
            size_t disappearance_node = dis_node_map_[n];

            constraint_pool.add_constraint(pgm::ConstraintPool::IncomingConstraint(transition_nodes, disappearance_node));
        }

        ////
        //// disappearance/appearance coupling
        ////
        if (app_node_map_.count(n) > 0 && dis_node_map_.count(n) > 0)
        {
            constraint_pool.add_constraint(pgm::ConstraintPool::DetectionConstraint((size_t)dis_node_map_[n], (size_t)app_node_map_[n]));
        }
    }

    constraint_pool.force_softconstraint(!with_constraints_);

    // For testing purposes:
    boost::posix_time::time_facet *facet = new boost::posix_time::time_facet("%Y%m%d-%H%M%S");
    std::stringstream filename_model;
    filename_model.imbue(std::locale(std::cout.getloc(), facet));
    filename_model << "./constracking_model-" << boost::posix_time::second_clock::local_time() << ".h5";

    opengm::hdf5::save(*pgm_->Model(), filename_model.str(), "model");
    
    std::stringstream filename_constraints;
    filename_constraints.imbue(std::locale(std::cout.getloc(), facet));
    filename_constraints << "./constracking_constraints-" << boost::posix_time::second_clock::local_time() << ".cp";

    std::ofstream out_stream(filename_constraints.str().c_str());
    boost::archive::text_oarchive oa(out_stream);
    oa & constraint_pool;
    out_stream.close();

    constraint_pool.add_constraints_to_problem(*pgm_->Model(), *optimizer_);
}

} /* namespace pgmlink */
