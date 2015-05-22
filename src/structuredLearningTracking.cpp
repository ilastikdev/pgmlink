#include <cassert>
#include <memory>
#include <set>
#include <string>
#include <iostream>
#include <sstream>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "pgmlink/randomforest.h"
#include "pgmlink/features/feature.h"
#include "pgmlink/pgm.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_pgm.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/merger_resolving.h"
#include "pgmlink/structuredLearningTracking.h"
#include "pgmlink/tracking.h"
#include <boost/python.hpp>
#include "pgmlink/field_of_view.h"

#include <stdio.h>

using boost::shared_ptr;
using boost::shared_array;

namespace pgmlink
{

////
//// class StructuredLearningTracking
////
EventVectorVectorVector StructuredLearningTracking::operator()(
    TraxelStore& ts,
    double forbidden_cost,
    double ep_gap,
    bool with_tracklets,
    double division_weight,
    double transition_weight,
    double disappearance_cost,
    double appearance_cost,
    bool with_merger_resolution,
    int n_dim,
    double transition_parameter,
    double border_width,
    bool with_constraints,
    UncertaintyParameter uncertaintyParam,
    double cplex_timeout,
    TimestepIdCoordinateMapPtr coordinates,
    boost::python::object transition_classifier)
{

    //build_hypo_graph(ts);

    // TODO need solution without copying the event vector


    /*
    EventVectorVectorVector events = track(
                                         forbidden_cost,
                                         ep_gap,
                                         with_tracklets,
                                         10.,
                                         division_weight,
                                         transition_weight,
                                         disappearance_cost,
                                         appearance_cost,
                                         with_merger_resolution,
                                         n_dim,
                                         transition_parameter,
                                         border_width,
                                         with_constraints,
                                         uncertaintyParam,
                                         cplex_timeout,
                                         transition_classifier
                                     );


    if (with_merger_resolution)
    {
        EventVectorVectorVector merger_resolved_events;

        for(auto& event : events)
        {
            merger_resolved_events.push_back(resolve_mergers(
                                                 event,
                                                 coordinates,
                                                 ep_gap,
                                                 transition_weight,
                                                 with_tracklets,
                                                 n_dim,
                                                 transition_parameter,
                                                 with_constraints,
                                                 transition_classifier
                                             ));
        }

        return merger_resolved_events;
    }
    else
    {
        return events;
    }
    */
}

namespace
{
std::vector<double> computeDetProb(double vol, std::vector<double> means, std::vector<double> s2)
{
    std::vector<double> result;

    double sum = 0;
    for (size_t k = 0; k < means.size(); ++k)
    {
        double val = vol - means[k];
        val = exp(-(val * val) / s2[k]);
        result.push_back(val);
        sum += val;
    }

    // normalize
    for(std::vector<double>::iterator it = result.begin(); it != result.end(); ++it)
    {
        (*it) /= sum;
    }

    return result;
}
}

void StructuredLearningTracking::hypothesesGraphTest(const HypothesesGraph& g)
{
    //boost::shared_ptr<std::vector< std::map<unsigned int, bool> > > stateOfNodes = state_of_nodes(g);
    HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());

    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        //model_.addVariable(param_.max_number_objects + 1);
        //app_node_map_[n] = model_.numberOfVariables() - 1;

        timestep_map[n];

        //assert(model_.numberOfLabels(app_node_map_[n]) == param_.max_number_objects + 1);

        //std::cout << count << " : " << timestep_map[n] << std::endl;
        ++count;
    }
    std::cout << "TOTAL NUMBER OF NODES: " << count << std::endl;

    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    typedef property_map<disappearance_label, HypothesesGraph::base_graph>::type disappearance_label_map;
    node_traxel_map& traxel_map = g.get(node_traxel());
    //disappearance_label_map& disappearance_labels_map = g.get(disappearance_label());

    for(int t = g.earliest_timestep(); t <= g.latest_timestep(); ++t)
    {
        std::cout << "TIME: " << t << std::endl;

        count = 0;
        for(node_timestep_map_t::ItemIt node(timestep_map, t); node != lemon::INVALID; ++node)
        {
            std::cout << "   Traxel Id: " << traxel_map[node].Id << "   Center: (" << traxel_map[node].X() << "," << traxel_map[node].Y() << "," << traxel_map[node].Z() << ")" << std::endl;
            //std::cout << "   Dissappearance Label: " << disappearance_labels_map[node] << std::endl;
            ++count;
        }
        std::cout << "   Number of Nodes: " << count << std::endl;
    }

    count = 0;
    for(HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
        HypothesesGraph::Node from = (&g)->source(a);
        HypothesesGraph::Node to = (&g)->target(a);
//        Traxel from_tr = traxel_map[from];
//        Traxel to_tr = traxel_map[to];

        ++count;
    }
    std::cout << "TOTAL NUMBER OF ARCS: " << count << std::endl;

}

void StructuredLearningTracking::addLabels(HypothesesGraph& g)
{
    g.add(appearance_label());
    g.add(disappearance_label());
    g.add(division_label());
    g.add(arc_label());
}

void StructuredLearningTracking::addAppearanceLabel(HypothesesGraph& g, int time, int label, double detectionProbability)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = g.get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            std::cout << " APPEARANCE Label   : [" << time << "] : " << traxel_map[node].Id << ": "  << detectionProbability << std::endl;
            g.add_appearance_label(node, detectionProbability);
        }
}

void StructuredLearningTracking::addDisappearanceLabel(HypothesesGraph& g, int time, int label, double detectionProbability)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = g.get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            std::cout << " DISAPPEARANCE Label: [" << time << "] : " << traxel_map[node].Id << ": "  << detectionProbability << std::endl;
            g.add_disappearance_label(node, detectionProbability);
        }
}

void StructuredLearningTracking::addDivisionLabel(HypothesesGraph& g, int time, int label, double detectionProbability)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = g.get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            std::cout << " DIVISION Label     : [" << time << "] : " << traxel_map[node].Id << ": "  << detectionProbability << std::endl;
            g.add_division_label(node, detectionProbability);
        }
}

void StructuredLearningTracking::addArcLabel(HypothesesGraph& g, int time, int startLabel, int endLabel, double detectionProbability)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    //typedef property_map<traxel_arc_id, HypothesesGraph::base_graph>::type traxel_arc_id_map;
    node_traxel_map& traxel_map = g.get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());
    //traxel_arc_id_map& arc_id_map = g.get(traxel_arc_id());

    HypothesesGraph::Node to;
    for(node_timestep_map_t::ItemIt node(timestep_map, time-1); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == startLabel){
            for(HypothesesGraph::base_graph::OutArcIt arc(g, node); arc != lemon::INVALID; ++arc){
                to = (&g)->target(arc);
                if (traxel_map[to].Id == endLabel){
                    std::cout << " ARC Label          : [" << time-1 << "," << time << "] : (" << traxel_map[node].Id << " ---> " << traxel_map[to].Id << "): "  << detectionProbability << std::endl;
                    g.add_arc_label(arc, detectionProbability);
                }
            }
        }
}

void StructuredLearningTracking::addFirstLabels(HypothesesGraph& g, int time, int label, double detectionProbability)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = g.get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            std::cout << " APPEARANCE Label   : [" << time << "] : " << traxel_map[node].Id << ": "  << detectionProbability << std::endl;
            g.add_appearance_label(node, detectionProbability);
            std::cout << " DISAPPEARANCE Label: [" << time << "] : " << traxel_map[node].Id << ": "  << 0 << std::endl;
            g.add_disappearance_label(node,0);
        }
}

void StructuredLearningTracking::addLastLabels(HypothesesGraph& g, int time, int label, double detectionProbability)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = g.get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            std::cout << " APPEARANCE Label   : [" << time << "] : " << traxel_map[node].Id << ": "  << 0 << std::endl;
            g.add_appearance_label(node,0);
            std::cout << " DISAPPEARANCE Label: [" << time << "] : " << traxel_map[node].Id << ": "  << detectionProbability << std::endl;
            g.add_disappearance_label(node,detectionProbability);
        }
}

void StructuredLearningTracking::addSingletonLabels(HypothesesGraph& g, int time, int label, double detectionProbability)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = g.get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            std::cout << " APPEARANCE Label   : [" << time << "] : " << traxel_map[node].Id << ": "  << detectionProbability << std::endl;
            g.add_appearance_label(node,detectionProbability);
            std::cout << " DISAPPEARANCE Label: [" << time << "] : " << traxel_map[node].Id << ": "  << detectionProbability << std::endl;
            g.add_disappearance_label(node,detectionProbability);
        }
}

void StructuredLearningTracking::addIntermediateLabels(HypothesesGraph& g, int time, int label, double detectionProbability)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = g.get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            std::cout << " APPEARANCE Label   : [" << time << "] : " << traxel_map[node].Id << ": "  << detectionProbability << std::endl;
            g.add_appearance_label(node,detectionProbability);
            std::cout << " DISAPPEARANCE Label: [" << time << "] : " << traxel_map[node].Id << ": "  << detectionProbability << std::endl;
            g.add_disappearance_label(node,detectionProbability);
        }
}


bool StructuredLearningTracking::exportCrop(FieldOfView crop)//, const std::string& name)
{
    crops_.push_back(crop);
    cout << "C++   ====> Crop starts: " << crops_.back().lower_bound()[0] << " " << crops_.back().lower_bound()[1] << " " << crops_.back().lower_bound()[2] << " " << crops_.back().lower_bound()[3] << endl;
    cout << "                 stops : " << crops_.back().upper_bound()[0] << " " << crops_.back().upper_bound()[1] << " " << crops_.back().upper_bound()[2] << " " << crops_.back().upper_bound()[3] << endl;

    return true;
}

/*
boost::shared_ptr<HypothesesGraph> StructuredLearningTracking::build_hypo_graph(TraxelStore& ts)
{

    LOG(logDEBUG3) << "entering build_hypo_graph" << endl;;

    LOG(logDEBUG1) << "max_number_objects  \t" << max_number_objects_  ;
    LOG(logDEBUG1) << "size_dependent_detection_prob\t" <<  use_size_dependent_detection_ ;
    LOG(logDEBUG1) << "avg_obj_size\t" <<      avg_obj_size_;
    LOG(logDEBUG1) << "with_divisions\t" <<      with_divisions_;
    LOG(logDEBUG1) << "division_threshold\t" <<      division_threshold_;

    traxel_store_ = &ts;

    use_classifier_prior_ = false;
    Traxel trax = *(traxel_store_->begin());
    FeatureMap::const_iterator it = trax.features.find("detProb");
    LOG(logDEBUG4) << "available features of " << trax << ":";
    for(auto it : trax.features.get())
    {
        LOG(logDEBUG4) << "\t" << it.first << std::endl;
    }

    if(it != trax.features.end())
    {
        use_classifier_prior_ = true;
        LOG(logDEBUG3) << "COULD find detProb, using classifier prior";
    }

    if(not use_classifier_prior_ and use_size_dependent_detection_)
    {
        LOG(logDEBUG3) << "creating detProb feature in traxel store from "
                          "deterministic detection probability function";
        std::vector<double> means;
        if (means_.size() == 0 )
        {
            for(int i = 0; i < max_number_objects_ + 1; ++i)
            {
                means.push_back(i * avg_obj_size_);
                LOG(logDEBUG4) << "mean[" << i << "] = " << means[i];
            }
        }
        else
        {
            assert(sigmas_.size() != 0);
            for(int i = 0; i < max_number_objects_ + 1; ++i)
            {
                means.push_back(means_[i]);
                LOG(logDEBUG4) << "mean[" << i << "] = " << means[i];
            }
        }

        std::vector<double> sigma2;
        if (sigmas_.size() == 0)
        {
            double s2 = (avg_obj_size_ * avg_obj_size_) / 4.0;
            if (s2 < 0.0001)
            {
                s2 = 0.0001;
            }
            for(int i = 0; i < max_number_objects_ + 1; ++i)
            {
                sigma2.push_back(s2);
                LOG(logDEBUG4) << "sigma2[" << i << "] = "  << sigma2[i];
            }
        }
        else
        {
            for (int i = 0; i < max_number_objects_ + 1; ++i)
            {
                sigma2.push_back(sigmas_[i]);
                LOG(logDEBUG4) << "sigma2[" << i << "] = "  << sigma2[i];
            }
        }

        for(TraxelStore::iterator tr = traxel_store_->begin(); tr != traxel_store_->end(); ++tr)
        {
            Traxel trax = *tr;
            FeatureMap::const_iterator it = trax.features.find("count");
            if(it == trax.features.end())
            {
                throw runtime_error("get_detection_prob(): cellness feature not in traxel");
            }
            double vol = it->second[0];
            std::vector<double> detProb;
            detProb = computeDetProb(vol, means, sigma2);
            feature_array detProbFeat(feature_array::difference_type(max_number_objects_ + 1));
            for(int i = 0; i <= max_number_objects_; ++i)
            {
                double d = detProb[i];
                if (d < 0.01)
                {
                    d = 0.01;
                }
                else if (d > 0.99)
                {
                    d = 0.99;
                }
                LOG(logDEBUG2) << "detection probability for " << trax.Id << "[" << i << "] = " << d;
                detProbFeat[i] = d;
            }
            trax.features["detProb"] = detProbFeat;
            traxel_store_->replace(tr, trax);
        }
    }

    LOG(logDEBUG1) << "-> building hypotheses" << endl;
    SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(1, // max_nearest_neighbors
            max_dist_,
            true, // forward_backward
            with_divisions_, // consider_divisions
            division_threshold_
                                                                );
    SingleTimestepTraxel_HypothesesBuilder hyp_builder(traxel_store_, builder_opts);
    hypotheses_graph_ = boost::shared_ptr<HypothesesGraph>(hyp_builder.build());

    hypotheses_graph_->add(arc_distance()).add(tracklet_intern_dist()).add(node_tracklet())
            .add(tracklet_intern_arc_ids()).add(traxel_arc_id());

    property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_distances =
            (hypotheses_graph_)->get(arc_distance());
    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map =
            (hypotheses_graph_)->get(node_traxel());


    Traxel some_traxel = (*traxel_map.beginValue());
    if (some_traxel.features.find("com_corrected") != some_traxel.features.end())
    {
        LOG(logINFO) << "optical correction enabled";
        with_optical_correction_ = true;
    }

    for(HypothesesGraph::ArcIt a(*hypotheses_graph_); a != lemon::INVALID; ++a)
    {
        HypothesesGraph::Node from = (hypotheses_graph_)->source(a);
        HypothesesGraph::Node to = (hypotheses_graph_)->target(a);
        Traxel from_tr = traxel_map[from];
        Traxel to_tr = traxel_map[to];

        if (with_optical_correction_)
        {
            arc_distances.set(a, from_tr.distance_to_corr(to_tr));
        }
        else
        {
            arc_distances.set(a, from_tr.distance_to(to_tr));
        }
    }

    if(event_vector_dump_filename_ != "none")
    {
        // store the traxel store and the resulting event vector
        std::ofstream ofs(event_vector_dump_filename_.c_str());
        boost::archive::text_oarchive out_archive(ofs);
        out_archive << ts;
    }
    return hypotheses_graph_;
}
*/
/*
boost::shared_ptr<HypothesesGraph> StructuredLearningTracking::get_hypo_graph()
{
    return hypotheses_graph_;
}

boost::shared_ptr<HypothesesGraph> StructuredLearningTracking::get_resolved_hypotheses_graph()
{
    if(!resolved_graph_)
    {
        throw std::runtime_error("Merger Resolving was not run, cannot get resolved graph.");
    }
    return resolved_graph_;
}



*/


/*

EventVectorVectorVector StructuredLearningTracking::track(double forbidden_cost,
        double ep_gap,
        bool with_tracklets,
        double detection_weight,
        double division_weight,
        double transition_weight,
        double disappearance_cost,
        double appearance_cost,
        bool with_merger_resolution,
        int n_dim,
        double transition_parameter,
        double border_width,
        bool with_constraints,
        UncertaintyParameter uncertaintyParam,
        double cplex_timeout,
        boost::python::object transition_classifier)
{
    ConservationTracking::Parameter param = get_conservation_tracking_parameters(
            forbidden_cost,
            ep_gap,
            with_tracklets,
            detection_weight,
            division_weight,
            transition_weight,
            disappearance_cost,
            appearance_cost,
            with_merger_resolution,
            n_dim,
            transition_parameter,
            border_width,
            with_constraints,
            uncertaintyParam,
            cplex_timeout,
            transition_classifier,
            solver_);
    uncertainty_param_ = uncertaintyParam;

    return track_from_param(param);
}

EventVectorVectorVector StructuredLearningTracking::track_from_param(ConservationTracking::Parameter& param,
                                                       bool fixLabeledNodes)
{
    original_hypotheses_graph_ = boost::make_shared<HypothesesGraph>();
    HypothesesGraph::copy(*hypotheses_graph_, *original_hypotheses_graph_);

//	PyEval_InitThreads();
//	PyGILState_STATE gilstate = PyGILState_Ensure();


    ConservationTracking pgm(param);

    pgm.labels_export_file_name_ = tracking_labels_export_file_name_;
    if(fixLabeledNodes)
    {
        pgm.enableFixingLabeledAppearanceNodes();
    }
    pgm.perturbedInference(*hypotheses_graph_);

    size_t num_solutions = uncertainty_param_.numberOfIterations;
    if (num_solutions == 1)
    {
        cout << "-> storing state of detection vars" << endl;
        last_detections_ = state_of_nodes(*hypotheses_graph_);
    }

    //	PyGILState_Release(gilstate);


    //TODO: conceptual problem here:
    //revise prune_inactive//events

    cout << "-> constructing unresolved events" << endl;

    EventVectorVectorVector all_ev(num_solutions);
    for (size_t i = 0; i < num_solutions; ++i)
    {
        all_ev[i] = *events(*hypotheses_graph_, i);
    }

    if(event_vector_dump_filename_ != "none")
    {
        // store the traxel store and the resulting event vector
        std::ofstream ofs(event_vector_dump_filename_.c_str());
        boost::archive::text_oarchive out_archive(ofs);
        out_archive << all_ev[0];
    }

    return all_ev;

}

ConservationTracking::Parameter StructuredLearningTracking::get_conservation_tracking_parameters(double forbidden_cost,
        double ep_gap,
        bool with_tracklets,
        double detection_weight,
        double division_weight,
        double transition_weight,
        double disappearance_cost,
        double appearance_cost,
        bool with_merger_resolution,
        int n_dim,
        double transition_parameter,
        double border_width,
        bool with_constraints,
        UncertaintyParameter uncertaintyParam,
        double cplex_timeout,
        boost::python::api::object transition_classifier,
        ConservationTracking::SolverType solver)
{
    LOG(logDEBUG1) << "max_number_objects  \t" << max_number_objects_  ;
    LOG(logDEBUG1) << "size_dependent_detection_prob\t" <<  use_size_dependent_detection_ ;
    LOG(logDEBUG1) << "forbidden_cost\t" <<      forbidden_cost;
    LOG(logDEBUG1) << "ep_gap\t" <<      ep_gap;
    LOG(logDEBUG1) << "avg_obj_size\t" <<      avg_obj_size_;
    LOG(logDEBUG1) << "with_tracklets\t" <<      with_tracklets;
    LOG(logDEBUG1) << "detection_weight\t" <<      detection_weight;
    LOG(logDEBUG1) << "division_weight\t" <<      division_weight;
    LOG(logDEBUG1) << "transition_weight\t" <<      transition_weight;
    LOG(logDEBUG1) << "with_divisions\t" <<      with_divisions_;
    LOG(logDEBUG1) << "disappearance_cost\t" <<      disappearance_cost;
    LOG(logDEBUG1) << "appearance_cost\t" <<      appearance_cost;
    LOG(logDEBUG1) << "with_merger_resolution\t" <<      with_merger_resolution;
    LOG(logDEBUG1) << "n_dim\t" <<      n_dim;
    LOG(logDEBUG1) << "transition_parameter\t" <<      transition_parameter;
    LOG(logDEBUG1) << "border_width\t" <<      border_width;
    LOG(logDEBUG1) << "with_constraints\t" <<      with_constraints;
    LOG(logDEBUG1) << "cplex_timeout\t" <<      cplex_timeout;
    uncertaintyParam.print();

    Traxels empty;
    boost::function<double(const Traxel&, const size_t)> detection, division;
    boost::function<double(const double)> transition;
    boost::function<double(const Traxel&)> appearance_cost_fn, disappearance_cost_fn;
    
    LOG(logDEBUG1) << "division_weight = " << division_weight;
    LOG(logDEBUG1) << "transition_weight = " << transition_weight;
    //border_width_ is given in normalized scale, 1 corresponds to a maximal distance of dim_range/2
    LOG(logINFO) << "using border-aware appearance and disappearance costs, with absolute margin: " << border_width;


    ConservationTracking::Parameter param(
        max_number_objects_,
        detection,
        division,
        transition,
        forbidden_cost,
        ep_gap,
        with_tracklets,
        with_divisions_,
        disappearance_cost_fn,
        appearance_cost_fn,
        true, // with_misdetections_allowed
        true, // with_appearance
        true, // with_disappearance
        transition_parameter,
        with_constraints,
        uncertaintyParam,
        cplex_timeout,
        division_weight,
        detection_weight,
        transition_weight,
        border_width,
        transition_classifier,
        with_optical_correction_,
        solver
    );

    std::vector<double> model_weights = {detection_weight, division_weight, transition_weight, disappearance_cost, appearance_cost};

    setParameterWeights(param,model_weights); 
    return param;
}

void StructuredLearningTracking::setParameterWeights(ConservationTracking::Parameter& param,std::vector<double> weights)
{

    param.detection_weight  =weights[0];
    param.division_weight   =weights[1];
    param.transition_weight =weights[2];

    size_t tmin = hypotheses_graph_->earliest_timestep();
    size_t tmax = hypotheses_graph_->latest_timestep();

    if (use_classifier_prior_)
    {
        LOG(logINFO) << "Using classifier prior";
        param.detection = NegLnDetection(weights[0]);
    }
    else if (use_size_dependent_detection_)
    {
        LOG(logINFO) << "Using size dependent prior";
        param.detection = NegLnDetection(weights[0]); // weight
    }
    else
    {
        LOG(logINFO) << "Using hard prior";
        // assume a quasi geometric distribution
        std::vector<double> prob_vector;
        double p = 0.7; // e.g. for max_number_objects=3, p=0.7: P(X=(0,1,2,3)) = (0.027, 0.7, 0.21, 0.063)
        double sum = 0;
        for(double state = 0; state < max_number_objects_; ++state)
        {
            double prob = p * pow(1 - p, state);
            prob_vector.push_back(prob);
            sum += prob;
        }
        prob_vector.insert(prob_vector.begin(), 1 - sum);

        param.detection = boost::bind<double>(NegLnConstant(weights[0], prob_vector), _2);
    }

    param.division = NegLnDivision(weights[1]);
    param.transition = NegLnTransition(weights[2]);

    param.appearance_cost_fn = SpatialBorderAwareWeight(weights[4],
                             param.border_width,
                             false, // true if relative margin to border
                             fov_,
                             tmin);// set appearance cost to zero at t = tmin

    param.disappearance_cost_fn = SpatialBorderAwareWeight(weights[3],
                            param.border_width,
                            false, // true if relative margin to border
                            fov_,
                            tmax);// set disappearance cost to zero at t = tmax
}

////
//// helper function needed for boost::algorithm::all_of
////
template <typename T>
bool has_data(const std::vector<T>& vector)
{
    return vector.size() > 0;
}


////
//// helper function equivalent to all_of (c++11)
////
template<class InputIterator, class UnaryPredicate>
bool all_true (InputIterator first, InputIterator last, UnaryPredicate pred)
{
    while (first != last)
    {
        if (!pred(*first))
        {
            return false;
        }
        ++first;
    }
    return true;
}

EventVectorVector StructuredLearningTracking::resolve_mergers(
    EventVectorVector& events,
    TimestepIdCoordinateMapPtr coordinates,
    double ep_gap,
    double transition_weight,
    bool with_tracklets,
    int n_dim,
    double transition_parameter,
    bool with_constraints,
    boost::python::object transitionClassifier
)
{
    // TODO Redundancy to track(). -> Problem?
    boost::function<double(const double)> transition;
    transition = NegLnTransition(transition_weight);

    cout << "-> resolving mergers" << endl;
    // TODO why doesn't it check for empty vectors in the event vector from the
    // first element on?
    if ( not all_true(events.begin() + 1, events.end(), has_data<Event>))
    {
        LOG(logDEBUG) << "Nothing to be done in ConstTracking::resolve_mergers:";
        LOG(logDEBUG) << "Empty vector in event vector";
    }
    else if (max_number_objects_ == 1)
    {
        LOG(logDEBUG) << "Nothing to resolve in ConstTracking::resolve_mergers:";
        LOG(logDEBUG) << "max_number_objects = 1";
    }
    else
    {
        // create a copy of the hypotheses graph to perform merger resolution without destroying the old graph
        resolved_graph_ = boost::make_shared<HypothesesGraph>();
        HypothesesGraph::copy(*hypotheses_graph_, *resolved_graph_);

        MergerResolver m(resolved_graph_.get());
        FeatureExtractorBase* extractor;
        DistanceFromCOMs distance;
        if (coordinates)
        {
            extractor = new FeatureExtractorArmadillo(coordinates);
        }
        else
        {
            calculate_gmm_beforehand(*resolved_graph_, 1, n_dim);
            extractor = new FeatureExtractorMCOMsFromMCOMs;
        }
        FeatureHandlerFromTraxels handler(*extractor, distance, traxel_store_);

        m.resolve_mergers(handler);

        HypothesesGraph g_res;
        resolve_graph(*resolved_graph_,
                      g_res,
                      transition,
                      ep_gap,
                      with_tracklets,
                      transition_parameter,
                      with_constraints,
                      transitionClassifier,
                      solver_);
//            prune_inactive(resolved_graph);

        cout << "-> constructing resolved events" << endl;
        boost::shared_ptr<std::vector< std::vector<Event> > > multi_frame_moves = multi_frame_move_events(*resolved_graph_);
        boost::shared_ptr<std::vector< std::vector<Event> > > resolved_tos = resolved_to_events(*resolved_graph_);

        cout << "-> merging unresolved and resolved events" << endl;
        // delete extractor; // TO DELETE FIRST CREATE VIRTUAL DTORS
        boost::shared_ptr<EventVectorVector> events_tmp = merge_event_vectors(events, *multi_frame_moves);
        boost::shared_ptr<EventVectorVector> events_ptr = merge_event_vectors(*events_tmp, *resolved_tos);
        //      all_ev[0] = *merge_event_vectors(*ev, *multi_frame_moves);

        // TODO The in serialized event vector written in the track() function
        // will be overwritten. Is this the desired behaviour?
        if(event_vector_dump_filename_ != "none")
        {
            // store the traxel store and the resulting event vector
            std::ofstream ofs(event_vector_dump_filename_.c_str());
            boost::archive::text_oarchive out_archive(ofs);
            out_archive << *events_ptr;
        }

        // cleanup extractor
        delete extractor;

        return *events_ptr;
    }
    cout << "-> done resolving mergers" << endl;
    return events;
}

std::vector<std::map<unsigned int, bool> > StructuredLearningTracking::detections()
{
    std::vector<std::map<unsigned int, bool> > res;
    if (last_detections_)
    {
        return *last_detections_;
    }
    else
    {
        throw std::runtime_error(
            "MrfTracking::detections(): previous tracking result required");
    }
}

void StructuredLearningTracking::save_ilp_solutions(const std::string& filename)
{
    std::ofstream result_file(filename.c_str());

    if(!result_file.good())
    {
        throw std::runtime_error("Couldn't open file to save ILP solutions");
    }

    for(size_t variable_idx = 0; variable_idx < ilp_solutions_[0].size(); variable_idx++)
    {
        for(size_t result_idx = 0; result_idx < ilp_solutions_.size(); result_idx++)
        {
            result_file << ilp_solutions_[result_idx][variable_idx] << " ";
        }

        result_file << "\n";
    }
}

void StructuredLearningTracking::setTrackLabelingExportFile(std::string file_name)
{
    tracking_labels_export_file_name_ = file_name;
}

void StructuredLearningTracking::createStructuredLearningFiles(std::string feature_file_name,
                                      std::string constraints_file_name,
                                      std::string ground_truth_file_name)
{

    if(not  constraints_file_name.empty())
    {
        std::ofstream constraints_file;
        constraints_file.open (constraints_file_name);
        constraints_file.close();
    }

    if(not feature_file_name.empty())
    {
        std::ofstream feature_file;
        feature_file.open (feature_file_name);
        feature_file.close();
    }

    if(not ground_truth_file_name.empty())
    {
        std::ofstream labels_file;
        labels_file.open (ground_truth_file_name);
        labels_file.close();
    }
}
*/



/*
void StructuredLearningTracking::writeStructuredLearningFiles(std::string feature_file_name,
                                      std::string constraints_file_name,
                                      std::string ground_truth_file_name,
                                      ConservationTracking::Parameter param)
{

    //create empty files that opengm can append to
    createStructuredLearningFiles(feature_file_name,constraints_file_name,ground_truth_file_name);

    //every iteration creates one line in the feature matrix file
    // if no feature file is written we only need to create the modle once
    const size_t number_of_weights = 5;
    size_t num_writes = number_of_weights;
    if(feature_file_name.empty())
        num_writes = 1;

    //writing main loop
    for (int i = 0; i < num_writes; ++i)
    {
        std::vector<double> model_weights =  std::vector<double>(number_of_weights, 0.);
        model_weights[i] = 1.0;
        setParameterWeights(param,model_weights);

        ConservationTracking pgm(param);

        pgm.features_file_      = feature_file_name;
        pgm.constraints_file_   = constraints_file_name;

        //during the model creation the feature and constraint files are filled
        HypothesesGraph *graph = pgm.get_prepared_graph(*hypotheses_graph_);
        boost::shared_ptr<InferenceModel> inference_model = pgm.create_inference_model();
        inference_model->build_from_graph(*graph);

        boost::static_pointer_cast<ConsTrackingInferenceModel>(inference_model)->set_inference_params(
            0,
            feature_file_name,
            constraints_file_name,
            "");

        //write graph label to file (only once)
        if (!ground_truth_file_name.empty())
        {   
            inference_model->write_labeledgraph_to_file(*hypotheses_graph_, ground_truth_file_name);
            ground_truth_file_name.clear();
        }

        //only fill constraint file once
        if (!constraints_file_name.empty())
        {
            constraints_file_name.clear();
        }
    }
    
    // feature file has been created line by line. Transpose to get one column per feature
    if(not feature_file_name.empty())
    {
        transpose_matrix_in_file(feature_file_name);
    }
}
*/

/*
std::vector<double> StructuredLearningTracking::learnTrackingWeights(std::string feature_file_name,
                                      std::string constraints_file_name,
                                      std::string ground_truth_file_name,
                                      std::string lossweights,
                                      std::string options)
{

    std::vector<double> out;
    std::string command = std::string("/home/swolf/local/src/sbmrm/build/binaries/sbmrm") + " --featuresFile=" + feature_file_name + " --constraintsFile=" + constraints_file_name + " --labelsFile=" + ground_truth_file_name + " " +  options;
    if(not lossweights.empty())
    {
        command += " --weightCostsFile=" + ground_truth_file_name + " --weightCostString=" + "\"" + lossweights + "\" ";
    }

    LOG(logINFO) << "calling funkey with " << command;
    std::string shell_output =  exec(command.c_str());
    LOG(logINFO) << shell_output << endl;
    int start = shell_output.find("optimial w is [") + 15;
    int end = shell_output.find("]", start);
    std::string numlist = shell_output.substr(start, end - start);


    for(int i = 0; i != std::string::npos ; i = numlist.find(",", i))
    {
        if(i > 0)
            i += 2;
        out.push_back(string_to_double(numlist.substr(i, numlist.find(",", i) - i).c_str()));
    }
    return out;
}
*/


/*


double StructuredLearningTracking::hammingloss_of_files(std::string f1, std::string f2)
{
    ifstream in(f1);
    ifstream in2(f2);
    double loss = 0.;

    while ((!in.eof()) && (!in2.eof()))
    {
        string line, line2;
        getline(in, line);
        getline(in2, line2);
        LOG(logDEBUG4) << "StructuredLearningTracking::hammingloss_of_files: comparing  " << line[0] << " and " << line2[0] ;
        
        if(line[0] != line2[0])
            loss += 1;

    }
    return loss;
}
*/
} // namespace tracking
