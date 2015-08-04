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

// include the LPDef symbols only once!
#undef OPENGM_LPDEF_NO_SYMBOLS
#include <opengm/inference/auxiliary/lpdef.hxx>

#include "pgmlink/randomforest.h"
#include "pgmlink/features/feature.h"
#include "pgmlink/pgm.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_pgm.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/reasoner_constracking_explicit.h"
#include "pgmlink/merger_resolving.h"
#include "pgmlink/structured_learning_tracking_dataset.h"
#include "pgmlink/tracking.h"
#include <boost/python.hpp>

#include <stdio.h>

using namespace std;
using boost::shared_ptr;
using boost::shared_array;

//Quick and dity utilities

// from http://stackoverflow.com/questions/478898/how-to-execute-a-command-and-get-output-of-command-within-c
std::string exec(const char* cmd)
{
    FILE* pipe = popen(cmd, "r");
    if (!pipe)
    {
        return "ERROR";
    }
    char buffer[128];
    std::string result = "";
    while(!feof(pipe))
    {
        if(fgets(buffer, 128, pipe) != NULL)
        {
            result += buffer;
        }
    }
    pclose(pipe);
    return result;
}

void transpose_matrix_in_file(std::string filename)
{
    //from http://stackoverflow.com/questions/1729824/transpose-a-file-in-bash
    std::string awk_program = "gawk '{\n for (i=1; i<=NF; i++)  {\n a[NR,i] = $i \n} \n}\n NF>p { p = NF } \nEND {\n for(j=1; j<=p; j++) {\n str=a[1,j]\n for(i=2; i<=NR; i++){\n str=str\" \"a[i,j];\n }\n print str\n }\n }' ";
    system( (awk_program + filename + "> tmp.txt").c_str() ) ;
    system( (std::string("rm ") + filename).c_str() ) ;
    system( (std::string("cp tmp.txt ") + filename).c_str() ) ;
    system( "rm tmp.txt") ;
}

//from  http://stackoverflow.com/questions/392981/how-can-i-convert-string-to-double-in-c
double string_to_double( const std::string& s )
{
    std::istringstream i(s);
    double x;
    if (!(i >> x))
    {
        return 0;
    }
    return x;
}

namespace pgmlink
{
////
//// class ChaingraphTracking
////

void ChaingraphTracking::set_with_divisions(bool state)
{
    with_divisions_ = state;
}

void ChaingraphTracking::set_cplex_timeout(double seconds)
{
    cplex_timeout_ = seconds;
}

std::vector<std::vector<Event> > ChaingraphTracking::operator()(TraxelStore& ts)
{
    LOG(logINFO) << "Calling chaingraph tracking with the following parameters:\n"
                 << "\trandom forest filename: " << rf_fn_ << "\n"
                 << "\tappearance: " << app_ << "\n"
                 << "\tdisappearance: " << dis_ << "\n"
                 << "\tdetection: " << det_ << "\n"
                 << "\tmisdetection: " << mis_  << "\n"
                 << "\tcellness_by_random_forest: " << use_rf_  << "\n"
                 << "\topportunity cost: " << opportunity_cost_ << "\n"
                 << "\tforbidden cost: " << forbidden_cost_ << "\n"
                 << "\twith constraints: " << with_constraints_ << "\n"
                 << "\tfixed detections: " << fixed_detections_ << "\n"
                 << "\tmean division distance: " << mean_div_dist_ << "\n"
                 << "\tminimal division angle: " << min_angle_  << "\n"
                 << "\tcplex ep gap: " << ep_gap_ << "\n"
                 << "\tn neighbors: " <<  n_neighbors_ << "\n"
                 << "\twith divisions: " << with_divisions_  << "\n"
                 << "\tcplex timeout: " << cplex_timeout_ << "\n"
                 << "\talternative builder: " << alternative_builder_;



    cout << "-> building feature functions " << endl;
    SquaredDistance move;
    BorderAwareConstant appearance(app_, earliest_timestep(ts), true, 0);
    BorderAwareConstant disappearance(dis_, latest_timestep(ts), false, 0);
    GeometryDivision2 division(mean_div_dist_, min_angle_);

    Traxels empty;
    // random forest?
    boost::function<double(const Traxel&)> detection, misdetection;
    if (use_rf_)
    {
        LOG(logINFO) << "Loading Random Forest";
        vigra::RandomForest<RF::RF_LABEL_TYPE> rf = RF::getRandomForest(rf_fn_);
        std::vector<std::string> rf_features;
        rf_features.push_back("volume");
        rf_features.push_back("bbox");
        rf_features.push_back("position");
        rf_features.push_back("com");
        rf_features.push_back("pc");
        rf_features.push_back("intensity");
        rf_features.push_back("intminmax");
        rf_features.push_back("pair");
        rf_features.push_back("sgf");
        rf_features.push_back("lcom");
        rf_features.push_back("lpc");
        rf_features.push_back("lintensity");
        rf_features.push_back("lintminmax");
        rf_features.push_back("lpair");
        rf_features.push_back("lsgf");

        LOG(logINFO) << "Predicting cellness";
        RF::predict_traxels(ts, rf, rf_features, 1, "cellness");

        detection = NegLnCellness(det_);
        misdetection = NegLnOneMinusCellness(mis_);
    }
    else if (ts.begin()->features.find("detProb") != ts.begin()->features.end())
    {
        for (TraxelStore::iterator it = ts.begin(); it != ts.end(); ++it)
        {
            Traxel trax = *it;
            trax.features["cellness"] = trax.features["detProb"];
            assert(trax.features["detProb"].size() == 2);
            ts.replace(it, trax);
        }
        detection = NegLnCellness(det_);
        misdetection = NegLnOneMinusCellness(mis_);
    }
    else
    {
        detection = ConstantFeature(det_);
        misdetection = ConstantFeature(mis_);
    }

    cout << "-> building hypotheses" << endl;
    SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(n_neighbors_, 50);
    SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
    boost::shared_ptr<HypothesesGraph> graph = boost::shared_ptr<HypothesesGraph>(hyp_builder.build());

    cout << "-> init MRF reasoner" << endl;
    std::auto_ptr<Chaingraph> mrf;

    if(alternative_builder_)
    {
        pgm::chaingraph::TrainableModelBuilder b(appearance,
                disappearance,
                move,
                opportunity_cost_,
                forbidden_cost_);

        if (with_divisions_)
        {
            b.with_divisions(division);
        }

        b.with_detection_vars(detection, misdetection);
        mrf = std::auto_ptr<Chaingraph>(new Chaingraph(b,
                                                       with_constraints_,
                                                       ep_gap_,
                                                       fixed_detections_,
                                                       cplex_timeout_));
    }
    else
    {
        pgm::chaingraph::ECCV12ModelBuilder b(appearance,
                                              disappearance,
                                              move,
                                              opportunity_cost_,
                                              forbidden_cost_);

        if (with_divisions_)
        {
            b.with_divisions(division);
        }

        b.with_detection_vars(detection, misdetection);
        mrf = std::auto_ptr<Chaingraph>(new Chaingraph(b,
                                                       with_constraints_,
                                                       ep_gap_,
                                                       fixed_detections_,
                                                       cplex_timeout_));
    }

    cout << "-> formulate MRF model" << endl;
    mrf->formulate(*graph);

    cout << "-> infer" << endl;
    mrf->infer();

    cout << "-> conclude" << endl;
    mrf->conclude(*graph);

    cout << "-> storing state of detection vars" << endl;
    last_detections_ = state_of_nodes(*graph);

    cout << "-> pruning inactive hypotheses" << endl;
    prune_inactive(*graph);

    cout << "-> constructing events" << endl;

    return *events(*graph);
}

std::vector<std::map<unsigned int, bool> > ChaingraphTracking::detections()
{
    std::vector<std::map<unsigned int, bool> > res;
    if (last_detections_)
    {
        return *last_detections_;
    }
    else
    {
        throw std::runtime_error(
            "ChaingraphTracking::detections(): previous tracking result required");
    }
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


////
//// class ConsTracking
////
EventVectorVectorVector ConsTracking::operator()(TraxelStore& ts,
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

    build_hypo_graph(ts);

    // TODO need solution without copying the event vector
    EventVectorVectorVector events = track(
                                         forbidden_cost,
                                         ep_gap,
                                         with_tracklets,
                                         10./*detection*/,
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
}

boost::shared_ptr<HypothesesGraph> ConsTracking::build_hypo_graph(TraxelStore& ts)
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

boost::shared_ptr<HypothesesGraph> ConsTracking::get_hypo_graph()
{
    return hypotheses_graph_;
}

boost::shared_ptr<HypothesesGraph> ConsTracking::get_resolved_hypotheses_graph()
{
    if(!resolved_graph_)
    {
        throw std::runtime_error("Merger Resolving was not run, cannot get resolved graph.");
    }
    return resolved_graph_;
}

EventVectorVectorVector ConsTracking::track(double forbidden_cost,
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

    return ConsTracking::track_from_param(param);
}

void ConsTracking::prepareTracking(ConservationTracking& pgm, ConservationTracking::Parameter& param)
{}

EventVectorVectorVector ConsTracking::track_from_param(ConservationTracking::Parameter& param,
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

    prepareTracking(pgm, param);

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

ConservationTracking::Parameter ConsTracking::get_conservation_tracking_parameters(
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

void ConsTracking::setParameterWeights(ConservationTracking::Parameter& param,std::vector<double> weights)
{

    param.detection_weight  =weights[0];
    param.division_weight   =weights[1];
    param.transition_weight =weights[2];
    param.appearance_weight = weights[3];
    param.disappearance_weight = weights[4];

    size_t tmin = hypotheses_graph_->earliest_timestep();
    size_t tmax = hypotheses_graph_->latest_timestep();

    if (use_classifier_prior_)
    {
        LOG(logINFO) << "Using classifier prior";
        param.detection = NegLnDetection(weights[0]);
        param.detectionNoWeight = NegLnDetectionNoWeight(weights[0]);
    }
    else if (use_size_dependent_detection_)
    {
        LOG(logINFO) << "Using size dependent prior";
        param.detection = NegLnDetection(weights[0]); // weight
        param.detectionNoWeight = NegLnDetectionNoWeight(weights[0]);
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
        param.detectionNoWeight = boost::bind<double>(NegLnConstantNoWeight(weights[0], prob_vector), _2);
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

EventVectorVector ConsTracking::resolve_mergers(
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

std::vector<std::map<unsigned int, bool> > ConsTracking::detections()
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

void ConsTracking::save_ilp_solutions(const std::string& filename)
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

void ConsTracking::setTrackLabelingExportFile(std::string file_name)
{
    tracking_labels_export_file_name_ = file_name;
}

void ConsTracking::createStructuredLearningFiles(std::string feature_file_name,
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

void ConsTracking::writeStructuredLearningFiles(std::string feature_file_name,
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

std::vector<double> ConsTracking::learnTrackingWeights(std::string feature_file_name,
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

double ConsTracking::hammingloss_of_files(std::string f1, std::string f2)
{
    ifstream in(f1);
    ifstream in2(f2);
    double loss = 0.;

    while ((!in.eof()) && (!in2.eof()))
    {
        string line, line2;
        getline(in, line);
        getline(in2, line2);
        LOG(logDEBUG4) << "ConsTracking::hammingloss_of_files: comparing  " << line[0] << " and " << line2[0] ;
        
        if(line[0] != line2[0])
            loss += 1;

    }
    return loss;
}

/////////////////////////////////////////////////
//ConsExplicitTracking

////
//// class ConsExplicitTracking
////
EventVectorVectorVector ConsExplicitTracking::operator()(TraxelStore& ts,
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
    build_hypo_graph(ts);

    // TODO need solution without copying the event vector
    EventVectorVectorVector events = track(
                                         forbidden_cost,
                                         ep_gap,
                                         with_tracklets,
                                         10./*detection*/,
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
}

boost::shared_ptr<HypothesesGraph> ConsExplicitTracking::build_hypo_graph(TraxelStore& ts)
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

boost::shared_ptr<HypothesesGraph> ConsExplicitTracking::get_hypo_graph()
{
    return hypotheses_graph_;
}

boost::shared_ptr<HypothesesGraph> ConsExplicitTracking::get_resolved_hypotheses_graph()
{
    if(!resolved_graph_)
    {
        throw std::runtime_error("Merger Resolving was not run, cannot get resolved graph.");
    }
    return resolved_graph_;
}

EventVectorVectorVector ConsExplicitTracking::track(double forbidden_cost,
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
    ConservationExplicitTracking::Parameter param = get_conservation_tracking_parameters(
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

    return ConsExplicitTracking::track_from_param(param);
}

void ConsExplicitTracking::prepareTracking(ConservationExplicitTracking& pgm, ConservationExplicitTracking::Parameter& param)
{
}

EventVectorVectorVector ConsExplicitTracking::track_from_param(ConservationExplicitTracking::Parameter& param,
                                                       bool fixLabeledNodes)
{

    original_hypotheses_graph_ = boost::make_shared<HypothesesGraph>();
    HypothesesGraph::copy(*hypotheses_graph_, *original_hypotheses_graph_);

//	PyEval_InitThreads();
//	PyGILState_STATE gilstate = PyGILState_Ensure();

    ConservationExplicitTracking pgm(param);

    pgm.labels_export_file_name_ = tracking_labels_export_file_name_;
    if(fixLabeledNodes)
    {
        pgm.enableFixingLabeledAppearanceNodes();
    }

    prepareTracking(pgm, param);

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

ConservationExplicitTracking::Parameter ConsExplicitTracking::get_conservation_tracking_parameters(
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
        boost::python::api::object transition_classifier,
        ConservationExplicitTracking::SolverType solver)
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


    ConservationExplicitTracking::Parameter param(
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

void ConsExplicitTracking::setParameterWeights(ConservationExplicitTracking::Parameter& param,std::vector<double> weights)
{
    param.detection_weight  =weights[0];
    param.division_weight   =weights[1];
    param.transition_weight =weights[2];
    param.appearance_weight = weights[3];
    param.disappearance_weight = weights[4];

    size_t tmin = hypotheses_graph_->earliest_timestep();
    size_t tmax = hypotheses_graph_->latest_timestep();

    if (use_classifier_prior_)
    {
        LOG(logINFO) << "Using classifier prior";
        param.detection = NegLnDetection(weights[0]);
        param.detectionNoWeight = NegLnDetectionNoWeight(weights[0]);
    }
    else if (use_size_dependent_detection_)
    {
        LOG(logINFO) << "Using size dependent prior";
        param.detection = NegLnDetection(weights[0]); // weight
        param.detectionNoWeight = NegLnDetectionNoWeight(weights[0]);
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
        param.detectionNoWeight = boost::bind<double>(NegLnConstantNoWeight(weights[0], prob_vector), _2);
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

EventVectorVector ConsExplicitTracking::resolve_mergers(
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
        resolve_graph_explicit(*resolved_graph_,
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

std::vector<std::map<unsigned int, bool> > ConsExplicitTracking::detections()
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

void ConsExplicitTracking::save_ilp_solutions(const std::string& filename)
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

void ConsExplicitTracking::setTrackLabelingExportFile(std::string file_name)
{
    tracking_labels_export_file_name_ = file_name;
}

void ConsExplicitTracking::createStructuredLearningFiles(std::string feature_file_name,
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

void ConsExplicitTracking::writeStructuredLearningFiles(std::string feature_file_name,
                                      std::string constraints_file_name,
                                      std::string ground_truth_file_name,
                                      ConservationExplicitTracking::Parameter param)
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

        ConservationExplicitTracking pgm(param);

        pgm.features_file_      = feature_file_name;
        pgm.constraints_file_   = constraints_file_name;

        //during the model creation the feature and constraint files are filled
        HypothesesGraph *graph = pgm.get_prepared_graph(*hypotheses_graph_);
        boost::shared_ptr<InferenceModel> inference_model = pgm.create_inference_model();
        inference_model->build_from_graph(*graph);

        boost::static_pointer_cast<StructuredLearningTrackingInferenceModel>(inference_model)->set_inference_params(
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

std::vector<double> ConsExplicitTracking::learnTrackingWeights(std::string feature_file_name,
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

double ConsExplicitTracking::hammingloss_of_files(std::string f1, std::string f2)
{
    ifstream in(f1);
    ifstream in2(f2);
    double loss = 0.;

    while ((!in.eof()) && (!in2.eof()))
    {
        string line, line2;
        getline(in, line);
        getline(in2, line2);
        LOG(logDEBUG4) << "ConsExplicitTracking::hammingloss_of_files: comparing  " << line[0] << " and " << line2[0] ;

        if(line[0] != line2[0])
            loss += 1;

    }
    return loss;
}

} // namespace tracking
