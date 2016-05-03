#ifndef OPENGM_UNSIGNED_INTEGER_POW_HXX_
#define OPENGM_UNSIGNED_INTEGER_POW_HXX_
#endif

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

#ifdef WITH_GUROBI
#include <opengm/inference/lpgurobi.hxx>
#else
#include <opengm/inference/lpcplex2.hxx>
#endif

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include <opengm/learning/loss/hammingloss.hxx>
#include <opengm/learning/solver/CplexBackend.h>
#include <opengm/learning/solver/QuadraticSolverBackend.h>
#include <opengm/learning/struct-max-margin.hxx>

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
void StructuredLearningTracking::prepareTracking(
        ConservationTracking& pgm,
        Parameter& param,
        opengm::learning::Weights<double>& trackingWeights,
        bool withNormalization,
        bool withClassifierPrior)
{
    use_classifier_prior_ = withClassifierPrior;

    boost::shared_ptr<InferenceModel> inference_model =
        boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(create_inference_model(param,trackingWeights,withNormalization));

    numWeights_ = boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model)->getLearningWeights().size();
    pgm.setInferenceModel(inference_model);
}

boost::shared_ptr<InferenceModel> StructuredLearningTracking::create_inference_model(
        Parameter& param,
        opengm::learning::Weights<double>& trackingWeights,
        bool withNormalization)
{

    ep_gap_ = param.ep_gap;
    cplex_timeout_ = param.cplex_timeout;

        return boost::make_shared<StructuredLearningTrackingInferenceModel>(
            param,
            ep_gap_,
            cplex_timeout_,
            trackingWeights,
            withNormalization,
            fov_,
            param.border_width,
            numThreads_);
}


void StructuredLearningTracking::hypothesesGraphTest(const HypothesesGraph& g)
{
    HypothesesGraph::node_timestep_map& timestep_map = g.get(node_timestep());

    size_t count = 0;
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {

        timestep_map[n];

        ++count;
    }

    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    typedef property_map<disappearance_label, HypothesesGraph::base_graph>::type disappearance_label_map;
    node_traxel_map& traxel_map = g.get(node_traxel());

    for(int t = g.earliest_timestep(); t <= g.latest_timestep(); ++t)
    {

        count = 0;
        for(node_timestep_map_t::ItemIt node(timestep_map, t); node != lemon::INVALID; ++node)
        {
            ++count;
        }
    }

    count = 0;
    for(HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
        HypothesesGraph::Node from = (&g)->source(a);
        HypothesesGraph::Node to = (&g)->target(a);

        ++count;
    }
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
            hypotheses_graph_->add_division_label(node, cellCount);
        }
}

bool StructuredLearningTracking::addArcLabel(int startTime, int startLabel, int endLabel, double cellCount)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = hypotheses_graph_->get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = hypotheses_graph_->get(node_timestep());

    bool found = false;
    HypothesesGraph::Node to;
    for(node_timestep_map_t::ItemIt node(timestep_map, startTime); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == startLabel){
            for(HypothesesGraph::base_graph::OutArcIt arc(*hypotheses_graph_, node); arc != lemon::INVALID; ++arc){
                to = hypotheses_graph_->target(arc);
                if (traxel_map[to].Id == endLabel){
                    hypotheses_graph_->add_arc_label(arc, cellCount);
                    return true;
                }
            }
        }
//    if(not found){
//        std::cout << "[StructuredLearningTracking::addArcLabel] You have tried to set a label of an arc that does not exist!" << std::endl;
//    }
    return found;
}

void StructuredLearningTracking::addFirstLabels(int time, int label, double cellCount)
{
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = hypotheses_graph_->get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = hypotheses_graph_->get(node_timestep());

    for(node_timestep_map_t::ItemIt node(timestep_map, time); node != lemon::INVALID; ++node)
        if (traxel_map[node].Id == label){
            //hypotheses_graph_->add_disappearance_label(node,0);
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
            hypotheses_graph_->add_disappearance_label(node,cellCount);
            //hypotheses_graph_->add_appearance_label(node,0);
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
            hypotheses_graph_->add_disappearance_label(node,cellCount);
            hypotheses_graph_->add_appearance_label(node,cellCount);
        }
}

bool StructuredLearningTracking::exportCrop(FieldOfView crop)
{
    crops_.push_back(crop);

    numCrops_ = crops_.size();

    return true;
}

Parameter StructuredLearningTracking::get_structured_learning_tracking_parameters(
        double forbidden_cost,
        double ep_gap,
        bool with_tracklets,
        double detection_weight,
        double division_weight,
        double transition_weight,
        double disappearance_cost,
        double appearance_cost,
        bool with_merger_resolution,
        unsigned int n_dim,
        double transition_parameter,
        double border_width,
        bool with_constraints,
        UncertaintyParameter uncertaintyParam,
        double cplex_timeout,
        boost::python::api::object transition_classifier,
        SolverType solver,
        bool training_to_hard_constraints,
        unsigned int num_threads,
        bool withNormalization,
        bool withClassifierPrior,
        bool verbose)
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
    boost::function<double (const double)> transition;
//    boost::function<double(const Traxel&, const Traxel&, const size_t)> transition;
    boost::function<double(const Traxel&)> appearance_cost_fn, disappearance_cost_fn;

    LOG(logDEBUG1) << "division_weight = " << division_weight;
    LOG(logDEBUG1) << "transition_weight = " << transition_weight;
    //border_width_ is given in normalized scale, 1 corresponds to a maximal distance of dim_range/2

    LOG(logINFO) << "using border-aware appearance and disappearance costs, with absolute margin: " << border_width;

    Parameter param(
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
        with_merger_resolution,
        n_dim,
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
        solver,
        training_to_hard_constraints,
        num_threads,
        withNormalization,
        withClassifierPrior,
        verbose
    );
    use_classifier_prior_ = withClassifierPrior;

    std::vector<double> model_weights;
	model_weights.push_back(detection_weight);
	model_weights.push_back(division_weight);
	model_weights.push_back(transition_weight);
	model_weights.push_back(disappearance_cost);
	model_weights.push_back(appearance_cost);

    StructuredLearningTracking::setParameterWeights(
        param,
        model_weights);

    return param;
}

void StructuredLearningTracking::structuredLearningFromParam(Parameter& param)
{
    //transition_= param.transition;
    StructuredLearningTracking::structuredLearning(
        param.forbidden_cost,
        param.ep_gap,
        param.with_tracklets,
        param.detection_weight,
        param.division_weight,
        param.transition_weight,
        param.disappearance_weight,
        param.appearance_weight,
        param.with_merger_resolution,
        param.n_dim,
        param.transition_parameter,
        param.border_width,
        param.with_constraints,
        param.uncertainty_param,
        param.cplex_timeout,
        param.transition_classifier,
        param.num_threads,
        param.withNormalization,
        param.withClassifierPrior,
        param.verbose
    );
}

void StructuredLearningTracking::setParameterWeights(Parameter& param,std::vector<double> sltWeights)
{

    param.detection_weight  =sltWeights[0];
    param.division_weight   =sltWeights[1];
    param.transition_weight =sltWeights[2];
    param.appearance_weight = sltWeights[3];
    param.disappearance_weight = sltWeights[4];

    size_t tmin = hypotheses_graph_->earliest_timestep();
    size_t tmax = hypotheses_graph_->latest_timestep();

    if (use_classifier_prior_)
    {
        LOG(logINFO) << "[StructuredLearningTracking::setParameterWeights] Using classifier prior";
        param.detection = NegLnDetection(sltWeights[0]);
        param.detectionNoWeight = NegLnDetectionNoWeight(sltWeights[0]);
    }
    else if (use_size_dependent_detection_)
    {
        LOG(logINFO) << "Using size dependent prior";
        param.detection = NegLnDetection(sltWeights[0]);
        param.detectionNoWeight = NegLnDetectionNoWeight(sltWeights[0]);
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

        param.detection = boost::bind<double>(NegLnConstant(sltWeights[0], prob_vector), _2);
        param.detectionNoWeight = boost::bind<double>(NegLnConstantNoWeight(sltWeights[0], prob_vector), _2);
    }

    param.division = NegLnDivision(sltWeights[1]);
    param.divisionNoWeight = NegLnDivisionNoWeight(sltWeights[1]);

    param.transition = NegLnTransition(sltWeights[2]); // not needed if defined in python via register...
    param.transitionNoWeight = NegLnTransitionNoWeight(sltWeights[2]); // not needed if defined in python via register...

    param.appearance_cost_fn = SpatialBorderAwareWeight(sltWeights[3],
                             param.border_width,
                             false, // true if relative margin to border
                             fov_,
                             tmin);// set appearance cost to zero at t = tmin

    param.disappearance_cost_fn = SpatialBorderAwareWeight(sltWeights[4],
                            param.border_width,
                            false, // true if relative margin to border
                            fov_,
                            tmax);// set disappearance cost to zero at t = tmax
}

void StructuredLearningTracking::structuredLearning(
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
        boost::python::object transition_classifier,
        unsigned int num_threads,
        bool withNormalization,
        bool withClassifierPrior,
        bool verbose
        )
{

    typedef opengm::datasets::StructuredLearningTrackingDataset<
            StructuredLearningTrackingInferenceModel::GraphicalModelType,opengm::learning::HammingLoss> DSS;

    trackingWeights_.setWeight(0,detection_weight);
    trackingWeights_.setWeight(1,division_weight);
    trackingWeights_.setWeight(2,transition_weight);
    trackingWeights_.setWeight(3,appearance_cost);
    trackingWeights_.setWeight(4,disappearance_cost);

    DSS sltDataset(crops_, trackingWeights_);

    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map_t;
    typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_map;
    node_traxel_map& traxel_map = hypotheses_graph_->get(node_traxel());
    HypothesesGraph::node_timestep_map& timestep_map = hypotheses_graph_->get(node_timestep());

    property_map<traxel_arc_id, HypothesesGraph::base_graph>::type& arc_id_map = hypotheses_graph_->get(traxel_arc_id());

    std::vector<boost::shared_ptr<HypothesesGraph>> hypothesesSubGraph;
    std::vector<HypothesesGraph*> graph;
    std::vector<boost::shared_ptr<InferenceModel>> inference_model;
        
    // set up graphical models in parallel
    #ifdef WITH_OPENMP
    omp_lock_t modelLock;
    omp_init_lock(&modelLock);
    #pragma omp parallel for
    #endif
    for(size_t m=0; m<numCrops_; ++m){
        LOG(logINFO) << "\n\n GRAPHICAL MODEL.............. " << m << "\n";

        // lock the model
        #ifdef WITH_OPENMP
        omp_set_lock(&modelLock);
        sltDataset.lockModel(m);
        omp_unset_lock(&modelLock);
        #else
        sltDataset.lockModel(m);
        #endif

        HypothesesGraph::base_graph::NodeMap<bool> selected_nodes(*hypotheses_graph_);
        for (HypothesesGraph::NodeIt n(*hypotheses_graph_); n != lemon::INVALID; ++n)
            selected_nodes[n] = false;

        LOG(logINFO) << "[StructuredLearningTracking] Control print: select crop nodes";
        int count = 0;
        for(int t = crops_[m].lower_bound()[0]; t <= crops_[m].upper_bound()[0]; ++t)
        {
            for(node_timestep_map_t::ItemIt node(timestep_map, t); node != lemon::INVALID; ++node){

                if(ndim_==2 and
                   crops_[m].lower_bound()[1] <= traxel_map[node].X() and traxel_map[node].X() <= crops_[m].upper_bound()[1] and
                   crops_[m].lower_bound()[2] <= traxel_map[node].Y() and traxel_map[node].Y() <= crops_[m].upper_bound()[2] or
                   ndim_==3 and
                   crops_[m].lower_bound()[1] <= traxel_map[node].X() and traxel_map[node].X() <= crops_[m].upper_bound()[1] and
                   crops_[m].lower_bound()[2] <= traxel_map[node].Y() and traxel_map[node].Y() <= crops_[m].upper_bound()[2] and
                   crops_[m].lower_bound()[3] <= traxel_map[node].Z() and traxel_map[node].Z() <= crops_[m].upper_bound()[3] ){

                    selected_nodes[node] = true;
                    //LOG(logINFO) << "app. node: " << traxel_map[node].Id << "    t " << timestep_map[node] << "    node num " << count;

                    count++;
                }
            }
        }

        LOG(logINFO) << "[StructuredLearningTracking] Control print: select arcs";
        HypothesesGraph::base_graph::ArcMap<bool> selected_arcs(*hypotheses_graph_);
        for(HypothesesGraph::ArcIt a(*hypotheses_graph_); a != lemon::INVALID; ++a)
        {
            HypothesesGraph::Node from = hypotheses_graph_->source(a);
            HypothesesGraph::Node to = hypotheses_graph_->target(a);

            if(selected_nodes[from] and selected_nodes[to]){
                selected_arcs[a] = true;
            }
            else{
                selected_arcs[a] = false;
            }
        }

        LOG(logINFO) << "[StructuredLearningTracking] Control print: copy subgraph";
        hypothesesSubGraph.push_back(boost::make_shared<HypothesesGraph>());
        HypothesesGraph::copy_subgraph(*hypotheses_graph_, *(hypothesesSubGraph[m]),selected_nodes,selected_arcs);

        LOG(logINFO) << "[StructuredLearningTracking] Control print: parameters";
        Parameter param = get_structured_learning_tracking_parameters(
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
            solver_,
            false,
            num_threads,
            withNormalization,
            withClassifierPrior,
            verbose);
        uncertainty_param_ = uncertaintyParam;

        //param.transition = transition_;
        ConservationTracking pgm(param);

        LOG(logINFO) << "[StructuredLearningTracking] Control print: prepared graph";

        graph.push_back( pgm.get_prepared_graph(*(hypothesesSubGraph[m])) );

        size_t numNodes = 0;
        for (HypothesesGraph::NodeIt n(*(graph[m])); n != lemon::INVALID; ++n)
        {
            ++numNodes;
        }

        size_t numArcs = 0;
        for(HypothesesGraph::ArcIt a(*(graph[m])); a != lemon::INVALID; ++a)
        {
            ++numArcs;
        }

        LOG(logINFO) << "[StructuredLearningTracking] Control print: prepare tracking";

        prepareTracking(pgm, param, sltDataset.getWeights(),withNormalization,withClassifierPrior);
        LOG(logINFO) << "[StructuredLearningTracking] Control print: get inference model";
        inference_model.push_back(pgm.getInferenceModel());

        boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model[m])->setModelStartTime(crops_[m].lower_bound()[0]);
        boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model[m])->setModelEndTime(crops_[m].upper_bound()[0]);

        LOG(logINFO) << "[StructuredLearningTracking] Control print: build from graph";
        inference_model[m]->build_from_graph(*(graph[m]));

        LOG(logINFO) << "[StructuredLearningTracking] Control print: set inference parameters";
        boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model[m])->set_inference_params(
            1,//numberOfSolutions,
            "",//get_export_filename(0, features_file_),
            "",//constraints_file_,
            "");//get_export_filename(0, labels_export_file_name_));

        LOG(logINFO) << "[StructuredLearningTracking] Control print: set graphical model";
        sltDataset.setGraphicalModel(m, boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model[m])->model());
        sltDataset.resizeGTS(m);

        property_map< appearance_label, HypothesesGraph::base_graph>::type& appearance_labels = graph[m]->get(appearance_label());
        property_map< disappearance_label, HypothesesGraph::base_graph>::type& disappearance_labels = graph[m]->get(disappearance_label());
        property_map< division_label, HypothesesGraph::base_graph>::type& division_labels = graph[m]->get(division_label());
        property_map< arc_label, HypothesesGraph::base_graph>::type& arc_labels = graph[m]->get(arc_label());

        node_traxel_map& traxel_map_gm = graph[m]->get(node_traxel());
        HypothesesGraph::node_timestep_map& timestep_map_gm = graph[m]->get(node_timestep());

        LOG(logINFO) << "[StructuredLearningTracking] Control print: arcs";
        size_t number_of_transition_nodes = boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model[m])->get_number_of_transition_nodes();
        size_t indexArcs=0;
        for(HypothesesGraph::ArcIt a(*(graph[m])); a != lemon::INVALID; ++a)
        {
            LabelType arc_label = arc_labels[a];

            sltDataset.setGTS(
                m,
                (size_t) indexArcs,
                (size_t) arc_label);
            //LOG(logINFO) << " arc node: " << traxel_map_gm[graph[m]->source(a)].Id << "---> " << traxel_map_gm[graph[m]->target(a)].Id << "     arc num" << indexArcs;
            ++indexArcs;
        }
        assert ( indexArcs == number_of_transition_nodes);

        LOG(logINFO) << "[StructuredLearningTracking] Control print: appearances";
        size_t number_of_appearance_nodes = boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model[m])->get_number_of_appearance_nodes();
        size_t indexAppNodes=0;
        for (HypothesesGraph::NodeIt n(*(graph[m])); n != lemon::INVALID; ++n){

            sltDataset.setGTS(
                m,
                (size_t) indexArcs + indexAppNodes,
                (size_t)appearance_labels[n]);
            //LOG(logINFO) << "app. node: " << traxel_map_gm[n].Id << "    t " << timestep_map_gm[n] << "    node num " << indexArcs + indexAppNodes;
            ++indexAppNodes;
        }
        assert ( indexAppNodes == number_of_appearance_nodes );

        LOG(logINFO) << "[StructuredLearningTracking] Control print: disappearances";
        size_t number_of_disappearance_nodes = boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model[m])->get_number_of_disappearance_nodes();
        size_t indexDisAppNodes=0;
        for (HypothesesGraph::NodeIt n(*(graph[m])); n != lemon::INVALID; ++n){

            sltDataset.setGTS(
                m,
                (size_t) indexArcs + indexAppNodes + indexDisAppNodes,
                (size_t)disappearance_labels[n]);
            //LOG(logINFO) << "disapp. node: " << traxel_map_gm[n].Id << "    t " << timestep_map_gm[n] <<"    node num " << indexArcs + indexAppNodes + indexDisAppNodes;
            ++indexDisAppNodes;
        }
        assert ( indexDisAppNodes == number_of_disappearance_nodes);

        LOG(logINFO) << "[StructuredLearningTracking] Control print: divisions";
        if(boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model[m])->param_.with_divisions){
            size_t number_of_division_nodes = boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model[m])->get_number_of_division_nodes();
            std::vector<size_t> division_var_to_node_fun (number_of_division_nodes);
            size_t indexDivNodes=0;
            for (HypothesesGraph::NodeIt n(*(graph[m])); n != lemon::INVALID; ++n){
                size_t number_of_outarcs = 0;
                for (HypothesesGraph::OutArcIt a(*(graph[m]), n); a != lemon::INVALID; ++a)
                    ++number_of_outarcs;
                if (number_of_outarcs > 1){
                    division_var_to_node_fun[indexDivNodes] = traxel_map[n].Id;

                    sltDataset.setGTS(
                        m,
                        (size_t) indexArcs + indexAppNodes + indexDisAppNodes + indexDivNodes,
                        (size_t)division_labels[n]);
                    //LOG(logINFO) << "div. node: " << traxel_map_gm[n].Id << "    t " << timestep_map_gm[n] <<"    node num " << indexArcs + indexAppNodes + indexDisAppNodes + indexDivNodes;
                    ++indexDivNodes;
                }
            }
            assert ( indexDivNodes == number_of_division_nodes);
        }

        LOG(logINFO) << "[StructuredLearningTracking] Control print: build model with loss";
        sltDataset.build_model_with_loss(m);

        LOG(logINFO) << "[StructuredLearningTracking] Control print: unlock model";
        // unlock the model
        #ifdef WITH_OPENMP
        omp_set_lock(&modelLock);
        sltDataset.unlockModel(m);
        omp_unset_lock(&modelLock);
        #else
        sltDataset.unlockModel(m);
        #endif

        LOG(logINFO) << "[StructuredLearningTracking] Control print: end model";
    } // for model m

    opengm::learning::StructMaxMargin<DSS>::Parameter para;
    para.optimizerParameter_.lambda = 1.00;

    opengm::learning::StructMaxMargin<DSS> learner(sltDataset,para);

    typedef StructuredLearningTrackingInferenceModel::cplex2_optimizer INFCPLEX;
    INFCPLEX::Parameter infPara;

    //lpcplex2
    infPara.integerConstraintNodeVar_ = true;

    infPara.relaxation_ = infPara.TightPolytope;

//    infPara.relaxation_ = infPara.LocalPolytope;
//    infPara.maxNumIterations_ = 100;
//    infPara.maxNumConstraintsPerIter_ = 10;

    infPara.epGap_ = ep_gap;
    infPara.verbose_ = verbose;
    infPara.numberOfThreads_ = num_threads;
    infPara.challengeHeuristic_ = infPara.Weighted;
    infPara.useSoftConstraints_ = false;

    learner.learn<INFCPLEX>(infPara);
    const DSS::Weights& learnedWeights = learner.getWeights();
    
    for (size_t i=0; i<learnedWeights.numberOfWeights(); ++i){
        trackingWeights_.setWeight(i, learnedWeights.getWeight(i));
    }

//////////////////////////////////////////////
// running inference on a LEARNABLE model --- NOT to be done in practice, just a test
///////////////////////////////////////////////////////////
//#include <opengm/inference/lpcplex.hxx>

//    std::cout << "LPCplex START" << std::endl;
//    typedef opengm::LPCplex<StructuredLearningTrackingInferenceModel::GraphicalModelType, opengm::Minimizer> cplex_optimizer;
//    cplex_optimizer::Parameter cplex_param_;

//    cplex_param_.verbose_ = true;
//    cplex_param_.integerConstraint_ = true;
//    cplex_param_.epGap_ = 0.01;
//    cplex_param_.timeLimit_ = 1e+75;

//    cplex_optimizer  newinf  = cplex_optimizer(sltDataset.getModel(0), cplex_param_);
//    newinf.infer();
//    std::cout << "LPCplex END" << std::endl;

}

double StructuredLearningTracking::weight(int index){
    return trackingWeights_.getWeight(index);
}

void StructuredLearningTracking::setWeight(int index, double val){
    trackingWeights_.setWeight(index,val);
}

} // namespace tracking
