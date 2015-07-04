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

#include <opengm/learning/loss/hammingloss.hxx>

#include "pgmlink/randomforest.h"
#include "pgmlink/features/feature.h"
#include "pgmlink/pgm.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_pgm.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/merger_resolving.h"
#include "pgmlink/structuredLearningTracking.h"
#include "pgmlink/inferencemodel/structuredlearningtrackinginferencemodel.h"
//#include "pgmlink/tracking.h"
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
{}

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
void StructuredLearningTracking::prepareTracking(ConservationTracking& pgm, ConservationTracking::Parameter& param)
{
    inference_model_param_.max_number_objects = param.max_number_objects;
    inference_model_param_.with_constraints = param.with_constraints;
    inference_model_param_.with_tracklets = param.with_tracklets;
    inference_model_param_.with_divisions = param.with_divisions;
    inference_model_param_.with_appearance = param.with_appearance;
    inference_model_param_.with_disappearance = param.with_disappearance;
    inference_model_param_.with_misdetections_allowed = param.with_misdetections_allowed;
    inference_model_param_.with_optical_correction = param.with_optical_correction;
    inference_model_param_.detection = param.detection;
    inference_model_param_.detectionNoWeight = param.detectionNoWeight;
    inference_model_param_.division = param.division;
    inference_model_param_.transition = param.transition;
    inference_model_param_.transition_parameter = param.transition_parameter;
    inference_model_param_.transition_classifier = param.transition_classifier;
    inference_model_param_.forbidden_cost = param.forbidden_cost;
    inference_model_param_.appearance_cost = param.appearance_cost_fn;
    inference_model_param_.disappearance_cost = param.disappearance_cost_fn;

    boost::shared_ptr<InferenceModel> inference_model =
        boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(create_inference_model(param));

    boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model)->setWeight((size_t)0,param.detection_weight);
    boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model)->setWeight((size_t)1,param.division_weight);
    boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model)->setWeight((size_t)2,param.transition_weight);
    boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model)->setWeight((size_t)3,param.appearance_weight);
    boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model)->setWeight((size_t)4,param.disappearance_weight);
    numWeights_ = boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model)->weights_.size();

    pgm.setInferenceModel(inference_model);
    std::cout << " AFTER  --------------------- setInferenceModel" << std::endl;
}

boost::shared_ptr<InferenceModel> StructuredLearningTracking::create_inference_model(ConservationTracking::Parameter& param)
{

    ep_gap_ = param.ep_gap;
    cplex_timeout_ = param.cplex_timeout;

        return boost::make_shared<StructuredLearningTrackingInferenceModel>(
            inference_model_param_,
            ep_gap_,
            cplex_timeout_);
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

void StructuredLearningTracking::addLabels()
{
    hypotheses_graph_->add(appearance_label());
    hypotheses_graph_->add(disappearance_label());
    hypotheses_graph_->add(division_label());
    hypotheses_graph_->add(arc_label());
}

void StructuredLearningTracking::addAppearanceLabel(int time, int label, double cellCount)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    HypothesesGraph::node_timestep_map& timestep_map = hypotheses_graph_->get(node_timestep());

    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = hypotheses_graph_->get(node_traxel());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            //std::cout << " APPEARANCE Label   : [" << time << "] : " << traxel_map[node].Id << ": "  << cellCount << std::endl;
            hypotheses_graph_->add_appearance_label(node, cellCount);
        }
}

void StructuredLearningTracking::addDisappearanceLabel(int time, int label, double cellCount)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = hypotheses_graph_->get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = hypotheses_graph_->get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            //std::cout << " DISAPPEARANCE Label: [" << time << "] : " << traxel_map[node].Id << ": "  << cellCount << std::endl;
            hypotheses_graph_->add_disappearance_label(node, cellCount);
        }
}

void StructuredLearningTracking::addDivisionLabel(int time, int label, double cellCount)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = hypotheses_graph_->get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = hypotheses_graph_->get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            //std::cout << " DIVISION Label     : [" << time << "] : " << traxel_map[node].Id << ": "  << cellCount << std::endl;
            hypotheses_graph_->add_division_label(node, cellCount);
        }
}

void StructuredLearningTracking::addArcLabel(int startTime, int startLabel, int endLabel, double cellCount)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    //typedef property_map<traxel_arc_id, HypothesesGraph::base_graph>::type traxel_arc_id_map;
    node_traxel_map& traxel_map = hypotheses_graph_->get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = hypotheses_graph_->get(node_timestep());
    //traxel_arc_id_map& arc_id_map = g.get(traxel_arc_id());

    HypothesesGraph::Node to;
    for(node_timestep_map_t::ItemIt node(timestep_map, startTime); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == startLabel){
            for(HypothesesGraph::base_graph::OutArcIt arc(*hypotheses_graph_, node); arc != lemon::INVALID; ++arc){
                to = hypotheses_graph_->target(arc);
                if (traxel_map[to].Id == endLabel){
                    //std::cout << " ARC Label          : [" << startTime << "," << startTime+1 << "] : (" << traxel_map[node].Id << " ---> " << traxel_map[to].Id << "): "  << cellCount << std::endl;
                    hypotheses_graph_->add_arc_label(arc, cellCount);
                }
            }
        }
}

void StructuredLearningTracking::addFirstLabels(int time, int label, double cellCount)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = hypotheses_graph_->get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = hypotheses_graph_->get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            //std::cout << " DISAPPEARANCE Label: [" << time << "] : " << traxel_map[node].Id << ": "  << 0 << std::endl;
            hypotheses_graph_->add_disappearance_label(node,0);
            //std::cout << " APPEARANCE Label   : [" << time << "] : " << traxel_map[node].Id << ": "  << cellCount << std::endl;
            hypotheses_graph_->add_appearance_label(node, cellCount);
        }
}

void StructuredLearningTracking::addLastLabels(int time, int label, double cellCount)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = hypotheses_graph_->get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = hypotheses_graph_->get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            //std::cout << " DISAPPEARANCE Label: [" << time << "] : " << traxel_map[node].Id << ": "  << cellCount << std::endl;
            hypotheses_graph_->add_disappearance_label(node,cellCount);
            //std::cout << " APPEARANCE Label   : [" << time << "] : " << traxel_map[node].Id << ": "  << 0 << std::endl;
            hypotheses_graph_->add_appearance_label(node,0);
        }
}

void StructuredLearningTracking::addIntermediateLabels(int time, int label, double cellCount)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = hypotheses_graph_->get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = hypotheses_graph_->get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            //std::cout << " DISAPPEARANCE Label: [" << time << "] : " << traxel_map[node].Id << ": "  << cellCount << std::endl;
            hypotheses_graph_->add_disappearance_label(node,cellCount);
            //std::cout << " APPEARANCE Label   : [" << time << "] : " << traxel_map[node].Id << ": "  << cellCount << std::endl;
            hypotheses_graph_->add_appearance_label(node,cellCount);
        }
}

bool StructuredLearningTracking::exportCrop(FieldOfView crop)//, const std::string& name)
{
    crops_.push_back(crop);
    std::cout << "C++   ====> Crop starts: " << crops_.back().lower_bound()[0] << " " << crops_.back().lower_bound()[1] << " " << crops_.back().lower_bound()[2] << " " << crops_.back().lower_bound()[3] << std::endl;
    std::cout << "                 stops : " << crops_.back().upper_bound()[0] << " " << crops_.back().upper_bound()[1] << " " << crops_.back().upper_bound()[2] << " " << crops_.back().upper_bound()[3] << std::endl;

    numCrops_ = crops_.size();

    return true;
}

void StructuredLearningTracking::makeStructuredLearningTrackingDataset(
        double forbidden_cost,
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
        boost::python::object transition_classifier
        )
{
    std::cout << "C++  makeStructuredLearningTrackingDataset IN" << std::endl;
    opengm::datasets::StructuredLearningTrackingDataset<StructuredLearningTrackingInferenceModel::GraphicalModelType,opengm::learning::HammingLoss> sltDataset(numCrops_, crops_, numWeights_, numLabels_, ndim_, hypotheses_graph_);

    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = hypotheses_graph_->get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = hypotheses_graph_->get(node_timestep());

    // set up graphical models
    for(size_t m=0; m<numCrops_; ++m){
        std::cout << "___________________________________________________MODEL/CROP: " << m << std::endl;

        boost::shared_ptr<HypothesesGraph> hypothesesSubGraph = boost::make_shared<HypothesesGraph>();
        // boost::shared_ptr<HypothesesGraph> other = boost::make_shared<HypothesesGraph>();

        HypothesesGraph::base_graph::NodeMap<bool> selected_nodes(*hypotheses_graph_);
        for (HypothesesGraph::NodeIt n(*hypotheses_graph_); n != lemon::INVALID; ++n)
            selected_nodes[n] = false;

        //for(int t = hypothesesGraph.earliest_timestep(); t <= hypothesesGraph.latest_timestep(); ++t)
        for(int t = crops_[m].lower_bound()[0]; t <= crops_[m].upper_bound()[0]; ++t)
        {
            std::cout << "TIME: " << t << std::endl;

            //count = 0;
            for(node_timestep_map_t::ItemIt node(timestep_map, t); node != lemon::INVALID; ++node){

                if(ndim_==2 and
                   crops_[m].lower_bound()[1] <= traxel_map[node].X() and traxel_map[node].X() <= crops_[m].upper_bound()[1] and
                   crops_[m].lower_bound()[2] <= traxel_map[node].Y() and traxel_map[node].Y() <= crops_[m].upper_bound()[2] or
                   ndim_==3 and
                   crops_[m].lower_bound()[1] <= traxel_map[node].X() and traxel_map[node].X() <= crops_[m].upper_bound()[1] and
                   crops_[m].lower_bound()[2] <= traxel_map[node].Y() and traxel_map[node].Y() <= crops_[m].upper_bound()[2] and
                   crops_[m].lower_bound()[3] <= traxel_map[node].Z() and traxel_map[node].Z() <= crops_[m].upper_bound()[3] ){

                    // selected_nodes
                    std::cout << " Node:" << traxel_map[node].Id << std::endl;
                    selected_nodes[node] = true;
                }
            }
        }

        std::cout << "______________" << std::endl;
        HypothesesGraph::base_graph::ArcMap<bool> selected_arcs(*hypotheses_graph_);
        for(HypothesesGraph::ArcIt a(*hypotheses_graph_); a != lemon::INVALID; ++a)
        {
            HypothesesGraph::Node from = hypotheses_graph_->source(a);
            HypothesesGraph::Node to = hypotheses_graph_->target(a);

            if(selected_nodes[from] and selected_nodes[to])
                selected_arcs[a] = true;
            else
                selected_arcs[a] = false;
        }

        HypothesesGraph::copy_subgraph(*hypotheses_graph_, *hypothesesSubGraph,selected_nodes,selected_arcs);

        // structured learning tracking inference model
        // from ConsTracking::track
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

        // from ConsTracking::track_from_param
        ConservationTracking pgm(param);



        // from ConservationTracking::perturbedInference
        HypothesesGraph *graph = pgm.get_prepared_graph(*hypothesesSubGraph);

        //traxel_map = graph->get(node_traxel());

        size_t numNodes = 0;
        for (HypothesesGraph::NodeIt n(*graph); n != lemon::INVALID; ++n)
        {
            std::cout << " Node:" << traxel_map[n].Id << std::endl;
            ++numNodes;
        }
        std::cout << "TOTAL NUMBER OF NODES prepared_graph: " << numNodes << std::endl;

        size_t numArcs = 0;
        for(HypothesesGraph::ArcIt a(*graph); a != lemon::INVALID; ++a)
        {
            std::cout << " Arc (" << traxel_map[graph->source(a)].Id << " ---> " << traxel_map[graph->target(a)].Id << " ) " << std::endl;
            ++numArcs;
        }
        std::cout << "TOTAL NUMBER OF ARCS prepared_graph: " << numArcs << std::endl;

        prepareTracking(pgm, param);
        boost::shared_ptr<InferenceModel> inference_model = pgm.getInferenceModel();
        inference_model->build_from_graph(*graph);
        boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model)->set_inference_params(
            1,//numberOfSolutions,
            "",//get_export_filename(0, features_file_),
            "",//constraints_file_,
            "");//get_export_filename(0, labels_export_file_name_));
        //inference_model->infer();

        sltDataset.setGraphicalModel(m, boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model)->model());



        //sltDataset.build_model_with_loss(m);  // <--- THIS NEEDS to be called, but: ARE my factors set up properly with LEARNABLE functions???

        //set up ground truth
        LabelType numberOfLabels = max_number_objects_;
        sltDataset.resizeGTS(numCrops_);




        typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
        node_traxel_map& traxel_map = graph->get(node_traxel());
        HypothesesGraph::node_timestep_map& timestep_map = graph->get(node_timestep());

        property_map< appearance_label, HypothesesGraph::base_graph>::type& appearance_labels = graph->get(appearance_label());
        property_map< disappearance_label, HypothesesGraph::base_graph>::type& disappearance_labels = graph->get(disappearance_label());
        property_map< division_label, HypothesesGraph::base_graph>::type& division_labels = graph->get(division_label());
        property_map< arc_label, HypothesesGraph::base_graph>::type& arc_labels = graph->get(arc_label());
        //property_map< appearance_label, HypothesesGraph::base_graph>::type& detection_labels = graph->get(detection_label());

        numNodes = 0;
        for (HypothesesGraph::NodeIt n(*graph); n != lemon::INVALID; ++n){
            LabelType appearance_label=-1, disappearance_label=-1, division_label=-1;

            appearance_label = appearance_labels[n];
            disappearance_label = disappearance_labels[n];
            division_label = division_labels[n];

            std::cout << "t= " << timestep_map[n] << "   " << traxel_map[n].Id <<" (appearance, disappearance, division) = (" << appearance_label << "," << disappearance_label << "," << division_label << ")" << std::endl;

            //sltDataset.setGTS(m,traxel_map[n].Id,label);
            ++numNodes;
        }
        std::cout << " Total number of nodes in the subgraph = " << numNodes << std::endl;

        numArcs = 0;
        for(HypothesesGraph::ArcIt a(*graph); a != lemon::INVALID; ++a)
        {
            LabelType arc_label=-1;

            arc_label = arc_labels[a];
            std::cout << "t= (" << timestep_map[graph->source(a)] << "," << timestep_map[graph->target(a)] << ")   Arc:  ( " << traxel_map[graph->source(a)].Id << " ---> " << traxel_map[graph->target(a)].Id << " ) " << arc_label << std::endl;
            ++numArcs;
        }
        std::cout << "TOTAL NUMBER OF ARCS in the subgraph = " << numArcs << std::endl;
    } // for model m


    std::cout << "C++  makeStructuredLearningTrackingDataset OUT" << std::endl;
}

} // namespace tracking
