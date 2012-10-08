#include <cassert>
#include <memory>
#include <set>
#include <string>
#include <iostream>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>

#include "pgmlink/feature.h"
#include "pgmlink/graphical_model.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/log.h"
#include "pgmlink/reasoner_opengm.h"
#include "pgmlink/tracking.h"
#include "pgmlink/reasoner_nearestneighbor.h"

using namespace std;
using boost::shared_ptr;
using boost::shared_array;

namespace pgmlink {
////
//// class ChaingraphTracking
////
vector<vector<Event> > ChaingraphTracking::operator()(TraxelStore& ts) {
	cout << "-> building feature functions " << endl;
	SquaredDistance move;
	BorderAwareConstant appearance(app_, earliest_timestep(ts), true, 0);
	BorderAwareConstant disappearance(dis_, latest_timestep(ts), false, 0);
	GeometryDivision2 division(mean_div_dist_, min_angle_);

	Traxels empty;
	// random forest?
	boost::function<double(const Traxel&)> detection, misdetection;
	if (use_rf_) {
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
	} else {
	  detection = ConstantFeature(det_);
	  misdetection = ConstantFeature(mis_);
	}

	cout << "-> building hypotheses" << endl;
	SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(6, 50);
	SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
	shared_ptr<HypothesesGraph> graph = shared_ptr<HypothesesGraph>(hyp_builder.build());

	cout << "-> init MRF reasoner" << endl;
	std::auto_ptr<Chaingraph> mrf;

	if(alternative_builder_) {
	  pgm::TrainableChaingraphModelBuilder b(appearance,
						 disappearance,
						 move,
						 opportunity_cost_,
						 forbidden_cost_);
	  
	  b.with_divisions(division)
	   .with_detection_vars(detection, misdetection);
	  mrf = std::auto_ptr<Chaingraph>(new Chaingraph(b, with_constraints_, ep_gap_, fixed_detections_));
	} else {
	  pgm::ChaingraphModelBuilderECCV12 b(appearance,
					      disappearance,
					      move,
					      opportunity_cost_,
					      forbidden_cost_);
	  
	  b.with_divisions(division)
	   .with_detection_vars(detection, misdetection);
	  mrf = std::auto_ptr<Chaingraph>(new Chaingraph(b, with_constraints_, ep_gap_, fixed_detections_));
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

vector<map<unsigned int, bool> > ChaingraphTracking::detections() {
	vector<map<unsigned int, bool> > res;
	if (last_detections_) {
		return *last_detections_;
	} else {
		throw std::runtime_error(
				"ChaingraphTracking::detections(): previous tracking result required");
	}
}


////
//// class NNTracking
////
vector<vector<Event> > NNTracking::operator()(TraxelStore& ts) {
	cout << "-> building hypotheses" << endl;
	SingleTimestepTraxel_HypothesesBuilder::Options builder_opts(2, max(divDist_,movDist_));
	SingleTimestepTraxel_HypothesesBuilder hyp_builder(&ts, builder_opts);
	HypothesesGraph* graph = hyp_builder.build();
	HypothesesGraph& g = *graph;

	LOG(logDEBUG1) << "NNTracking: adding offered property to nodes";
	// adding 'offered' property and set it true for each node
	g.add(node_offered());
	property_map<node_offered, HypothesesGraph::base_graph>::type& offered_nodes = g.get(node_offered());
	for(HypothesesGraph::NodeIt n(g); n!=lemon::INVALID; ++n) {
		offered_nodes.set(n,true);
	}
	LOG(logDEBUG1) << "NNTracking: adding distance property to edges";
	g.add(arc_distance());
	property_map<arc_distance, HypothesesGraph::base_graph>::type& arc_distances = g.get(arc_distance());
	property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
	for(HypothesesGraph::ArcIt a(g); a!=lemon::INVALID; ++a) {
		HypothesesGraph::Node from = g.source(a);
		HypothesesGraph::Node to = g.target(a);
		Traxel from_tr = traxel_map[from];
		Traxel to_tr = traxel_map[to];
		double dist = from_tr.distance_to(to_tr);
		LOG(logDEBUG2) << "NNTracking:: distance from " << from_tr.Id << " to " << to_tr.Id << " = " << dist;
		arc_distances.set(a, dist);
	}


	cout << "-> init NN reasoner" << endl;
	SingleTimestepTraxelNN nn_reasoner(divDist_,movDist_);

	cout << "-> formulate NN model" << endl;
	nn_reasoner.formulate(*graph);

	cout << "-> infer" << endl;
	nn_reasoner.infer();

	cout << "-> conclude" << endl;
	nn_reasoner.conclude(*graph);

	cout << "-> storing state of detection vars" << endl;
	last_detections_ = state_of_nodes(*graph);

	cout << "-> pruning inactive hypotheses" << endl;
	prune_inactive(*graph);

	cout << "-> constructing events" << endl;

	return *events(*graph);
}

vector<map<unsigned int, bool> > NNTracking::detections() {
	vector<map<unsigned int, bool> > res;
	if (last_detections_) {
		return *last_detections_;
	} else {
		throw std::runtime_error(
				"NNTracking::detections(): previous tracking result required");
	}
}



} // namespace tracking
