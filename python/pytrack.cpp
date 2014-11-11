#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY
#define BOOST_PYTHON_MAX_ARITY 30

#include <vector>

#include "../include/pgmlink/field_of_view.h"
#include "../include/pgmlink/tracking_feature_extractor.h"
#include "../include/pgmlink/feature_extraction.h"
#include "../include/pgmlink/reasoner_constracking.h"
#include "../include/pgmlink/tracking.h"
#include <boost/utility.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python.hpp>

#include "pytemplated_pickle_suite.h"

using namespace std;
using namespace pgmlink;
using namespace boost::python;

vector<vector<Event> > pythonChaingraphTracking(ChaingraphTracking& tr, TraxelStore& ts) {
    vector<vector<Event> > result;
	// release the GIL
	Py_BEGIN_ALLOW_THREADS
	try {
	  result = tr(ts);
	} catch (std::exception& e) {
		Py_BLOCK_THREADS
		throw;
	}
	Py_END_ALLOW_THREADS
	return result;
}

vector<vector<vector<Event> > > pythonConsTracking(ConsTracking& tr, TraxelStore& ts, TimestepIdCoordinateMapPtr& coordinates,
					  double forbidden_cost,
					  double ep_gap,
					  bool   with_tracklets,
					  double division_weight,
					  double transition_weight,
					  double disappearance_cost,
					  double appearance_cost,
					  bool   with_merger_resolution,
					  int    n_dim,
					  double transition_parameter,
					  double border_width,
                      UncertaintyParameter uncertaintyParam,
					  bool   with_constraints,
					  double cplex_timeout,
                      object transition_classifier) {
    vector<vector<vector<Event> > > result;
	// release the GIL
	Py_BEGIN_ALLOW_THREADS
	try {
		result = tr(ts,
			    forbidden_cost,
			    ep_gap,
			    with_tracklets,
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
			    coordinates,
                transition_classifier);
	} catch (std::exception& e) {
		Py_BLOCK_THREADS
		throw;
	}
	Py_END_ALLOW_THREADS
	return result;
	}

void pywrite_all_funkey_features(ConsTracking& tracking,
                                 TraxelStore& ts,
                                 UncertaintyParameter uncertaintyParam,
                                 double forbidden_cost,
                                 int ndim,
                                 bool with_tracklets,
                                 double transition_parameter,
                                 double border_width)
{
    vector<vector<double>> parameterlist;

    for(size_t i = 0; i < 5; i++)
    {
        vector<double> params(5);
        params[i] = 1.0;
        parameterlist.push_back(params);
    }

    tracking.write_funkey_features(ts,
                                   parameterlist,
                                   uncertaintyParam,
                                   forbidden_cost,
                                   ndim,
                                   with_tracklets,
                                   transition_parameter,
                                   border_width);
}

feature_array    (pgmlink::feature_extraction::FeatureExtractor::*extract1)(const Traxel& t1) const = &pgmlink::feature_extraction::FeatureExtractor::extract;
feature_array    (pgmlink::feature_extraction::FeatureExtractor::*extract2)(const Traxel& t1, const Traxel& t2) const = &pgmlink::feature_extraction::FeatureExtractor::extract;
feature_array    (pgmlink::feature_extraction::FeatureExtractor::*extract3)(const Traxel& t1, const Traxel& t2, const Traxel& t3) const = &pgmlink::feature_extraction::FeatureExtractor::extract;

std::vector<double> pyextractor_get_feature_vector(pgmlink::features::TrackingFeatureExtractor& fe)
{
    std::vector<double> feature_vector;
    fe.get_feature_vector(feature_vector);
    return feature_vector;
}

void export_track() {
    class_<vector<Event> >("EventVector")
	.def(vector_indexing_suite<vector<Event> >())
    .def_pickle(TemplatedPickleSuite<EventVector>() )
    ;

    class_<vector<vector<Event> > >("NestedEventVector")
	.def(vector_indexing_suite<vector<vector<Event> > >())
    .def_pickle(TemplatedPickleSuite<EventVectorVector>() )
    ;

    class_<vector<vector<vector<Event> > > >("NestedNestedEventVector")
	.def(vector_indexing_suite<vector<vector<vector<Event> > > >())
    .def_pickle(TemplatedPickleSuite<EventVectorVectorVector>() )
    ;

    class_<map<unsigned int, bool> >("DetectionMap")
      .def(map_indexing_suite<map<unsigned int, bool> >())
      .def_pickle(TemplatedPickleSuite< map<unsigned int, bool> >() )
    ;

    class_<vector<map<unsigned int, bool> > >("DetectionMapsVector")
      .def(vector_indexing_suite<vector<map<unsigned int, bool> > >())
      .def_pickle(TemplatedPickleSuite< vector< map<unsigned int, bool> > >() )
    ;

    class_<vector<double>>("WeightVector")
      .def(vector_indexing_suite<vector<double>>())
      .def_pickle(TemplatedPickleSuite< vector<double> >() )
    ;

    class_<ChaingraphTracking>("ChaingraphTracking", 
			       init<string,double,double,double,double,
			       	   bool,double,double,bool,
			       	   bool,double,double,double,double>(
							  args("random_forest_filename", "appearance", "disappearance", "detection", "misdetection",
									  "use_random_forest", "opportunity_cost", "forbidden_cost", "with_constraints",
									  "fixed_detections", "mean_div_dist", "min_angle", "ep_gap", "n_neighbors"
									  )))
      .def("__call__", &pythonChaingraphTracking)
      .def("detections", &ChaingraphTracking::detections)
      .def("set_with_divisions", &ChaingraphTracking::set_with_divisions)
      .def("set_cplex_timeout", &ChaingraphTracking::set_cplex_timeout)
    ;

    class_<ConsTracking>("ConsTracking",
                         init<int,bool,double,double,bool,double,string,FieldOfView, string>(
                                            args("max_number_objects",
                                                 "size_dependent_detection_prob",
                                                 "avg_obj_size",
                                                 "max_neighbor_distance",
                                                 "with_division",
                                                 "division_threshold",
                                                 "detection_rf_filename",
                                                 "fov",
                                                 "event_vector_dump_filename")))
      .def("__call__", &pythonConsTracking)
      .def("buildGraph", &ConsTracking::build_hypo_graph)
      .def("track", &ConsTracking::track)
      .def("resolve_mergers", &ConsTracking::resolve_mergers)
      .def("detections", &ConsTracking::detections)
      .def("get_hypotheses_graph", &ConsTracking::get_hypo_graph)
      .def("get_resolved_hypotheses_graph", &ConsTracking::get_resolved_hypotheses_graph)
      .def("SetFunkeyOutputFiles",&ConsTracking::write_funkey_set_output_files)
      .def("writeFunkeyFeatures",&ConsTracking::write_funkey_features)
      .def("writeAllFunkeyFeatures",&pywrite_all_funkey_features)
      .def("writeFunkeyFiles",&ConsTracking::write_funkey_files)
      .def("HamminglossOfFiles",&ConsTracking::hammingloss_of_files)
      .def("LearnWithFunkey",&ConsTracking::learn_from_funkey_files)
      .def("SetFunkeyExportLabeledGraph",&ConsTracking::set_export_labeled_graph)
      .def("save_ilp_solutions", &ConsTracking::save_ilp_solutions)
	;

    enum_<Event::EventType>("EventType")
	.value("Move", Event::Move)
	.value("Division", Event::Division)
	.value("Appearance", Event::Appearance)
	.value("Disappearance", Event::Disappearance)
	.value("Merger", Event::Merger)
	.value("ResolvedTo", Event::ResolvedTo)
	.value("MultiFrameMove", Event::MultiFrameMove)
	.value("Void", Event::Void)
    ;

	class_<pgmlink::feature_extraction::FeatureExtractor>("FeatureExtractor",
                                                        init<string, string>(args("Calculator, Feature")))
            .def("extract", extract1)
            .def("extract", extract2)
            .def("extract", extract3)
    ;

    class_<pgmlink::features::TrackingFeatureExtractor>("TrackingFeatureExtractor",
                                                        init<boost::shared_ptr<HypothesesGraph>, FieldOfView>(args("HypothesesGraph, FieldOfView")))
            .def("compute_features", &pgmlink::features::TrackingFeatureExtractor::compute_features)
            .def("append_feature_vector_to_file", &pgmlink::features::TrackingFeatureExtractor::append_feature_vector_to_file)
            .def("train_track_svm", &pgmlink::features::TrackingFeatureExtractor::train_track_svm)
            .def("get_track_svm", &pgmlink::features::TrackingFeatureExtractor::get_track_svm)
            .def("set_track_svm", &pgmlink::features::TrackingFeatureExtractor::set_track_svm)
            .def("get_feature_vector", &pyextractor_get_feature_vector)
            .def("set_track_feature_output_file", &pgmlink::features::TrackingFeatureExtractor::set_track_feature_output_file);

    class_<pgmlink::features::SVMOutlierCalculator, boost::shared_ptr<pgmlink::features::SVMOutlierCalculator> >("SVMOutlierCalculator")
            .def_pickle(TemplatedPickleSuite<pgmlink::features::SVMOutlierCalculator>());


    class_<vector<vigra::UInt64> >("IdVector")
    .def(vector_indexing_suite<vector<vigra::UInt64> >())
    ;

    class_<Event>("Event")
	.def_readonly("type", &Event::type)
	.def_readonly("traxel_ids", &Event::traxel_ids)
	.add_property("energy", &Event::energy)
    .def_pickle(TemplatedPickleSuite<Event>() )
    ;
}
