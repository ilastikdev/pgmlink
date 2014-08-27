#define BOOST_TEST_MODULE reasoner_constracking_test

#include <vector>
#include <iostream>
#include <set>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/bind.hpp>

#include <lemon/color.h>
#include <lemon/graph_to_eps.h>

#include <lemon/core.h>
#include <lemon/concepts/digraph.h>
#include <lemon/list_graph.h>
#include <lemon/maps.h>

#include "pgmlink/graph.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/feature.h"
#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/traxels.h"
#include "pgmlink/tracking.h"
#include "pgmlink/field_of_view.h"


using namespace pgmlink;
using namespace std;
using namespace boost;


BOOST_AUTO_TEST_CASE( Tracking_ConservationTracking_Tracklets ) {

	std::cout << "Constructing HypothesesGraph" << std::endl;
	std::cout << std::endl;

	using lemon::INVALID;

	std::cout << "Adding Traxels to TraxelStore" << std::endl;
	std::cout << std::endl;

	/*
	//  t=1      2       3
	//  1
	//    |
	//    |
	//      ---- 1 ----  1
	//    |
	//  2
	//  1 ----   2 ----- 1
	//  0 ------ 1 ----- 0
	TraxelStore ts;
	Traxel n11, n12, n21, n31, n13, n22, n32;
	Traxel n14, n24, n34;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	feature_array count(feature_array::difference_type(1));
	n11.Id = 1; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1; count[0] = 1;
	n11.features["com"] = com; n11.features["divProb"] = divProb; n11.features["count"] = count;
	add(ts,n11);
	n12.Id = 3; n12.Timestep = 1; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1; count[0] = 3;
	n12.features["com"] = com; n12.features["divProb"] = divProb; n12.features["count"] = count;
	add(ts,n12);
	n21.Id = 10; n21.Timestep = 2; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1; count[0] = 1;
	n21.features["com"] = com; n21.features["divProb"] = divProb; n21.features["count"] = count;
	add(ts,n21);
	n31.Id = 11; n31.Timestep = 3; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1; count[0] = 1;
	n31.features["com"] = com; n31.features["divProb"] = divProb; n31.features["count"] = count;
	add(ts,n31);

	n13.Id = 12; n13.Timestep = 1; com[0] = 100; com[1] = 100; com[2] = 100; divProb[0] = 0.1; count[0] = 1;
	n13.features["com"] = com; n13.features["divProb"] = divProb; n13.features["count"] = count;
	add(ts,n13);
	n22.Id = 13; n22.Timestep = 2; com[0] = 100; com[1] = 100; com[2] = 100; divProb[0] = 0.1; count[0] = 3;
	n22.features["com"] = com; n22.features["divProb"] = divProb; n22.features["count"] = count;
	add(ts,n22);
	n32.Id = 14; n32.Timestep = 3; com[0] = 100; com[1] = 100; com[2] = 100; divProb[0] = 0.1; count[0] = 1;
	n32.features["com"] = com; n32.features["divProb"] = divProb; n32.features["count"] = count;
	add(ts,n32);

	n14.Id = 15; n14.Timestep = 1; com[0] = 200; com[1] = 100; com[2] = 100; divProb[0] = 0.1; count[0] = 0.1;
	n14.features["com"] = com; n14.features["divProb"] = divProb; n14.features["count"] = count;
	add(ts,n14);
	n24.Id = 16; n24.Timestep = 2; com[0] = 200; com[1] = 100; com[2] = 100; divProb[0] = 0.1; count[0] = 1;
	n24.features["com"] = com; n24.features["divProb"] = divProb; n24.features["count"] = count;
	add(ts,n24);
	n34.Id = 17; n34.Timestep = 3; com[0] = 200; com[1] = 100; com[2] = 100; divProb[0] = 0.1; count[0] = 0.1;
	n34.features["com"] = com; n34.features["divProb"] = divProb; n34.features["count"] = count;
	add(ts,n34);
	*/
	
		//  t=1      2      3       4
	//  o                       o
	//    |                    |
	//      ---- o ---- o ----
	//    |                    |
	//  o                       o
	TraxelStore ts;
	Traxel n11, n12, n21, n31, n41, n42;
	feature_array com(feature_array::difference_type(3));
	feature_array divProb(feature_array::difference_type(1));
	n11.Id = 1; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
	n11.features["com"] = com; n11.features["divProb"] = divProb;
	add(ts,n11);
	n12.Id = 3; n12.Timestep = 1; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1;
	n12.features["com"] = com; n12.features["divProb"] = divProb;
	add(ts,n12);
	n21.Id = 10; n21.Timestep = 2; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.1;
	n21.features["com"] = com; n21.features["divProb"] = divProb;
	add(ts,n21);
	n31.Id = 11; n31.Timestep = 3; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.1;
	n31.features["com"] = com; n31.features["divProb"] = divProb;
	add(ts,n31);
	n41.Id = 12; n41.Timestep = 4; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
	n41.features["com"] = com; n41.features["divProb"] = divProb;
	add(ts,n41);
	n42.Id = 13; n42.Timestep = 4; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
	n42.features["com"] = com; n42.features["divProb"] = divProb;
	add(ts,n42);


	std::cout << "Initialize Conservation tracking" << std::endl;
	std::cout << std::endl;

	std::string sbrmr_binary = "/home/wolf/Machine_Learning/sbmrm/build/binaries/sbmrm";

	FieldOfView fov(0, 0, 0, 0, 4, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
	ConsTracking tracking = ConsTracking(
					     3, // max_number_objects
					     false, // detection_by_volume
					     double(1.1), // avg_obj_size
					     20, // max_neighbor_distance
					     true, //with_divisions
					     0.3, // division_threshold
					     "none" // random_forest_filename
					     );

	std::cout << "Write Funkey Files" << std::endl;
	std::cout << std::endl;

	tracking.write_funkey_files(ts,"features.txt","constraints.txt");
	
	//create lables from tracking with all weights = 10.
	tracking.write_funkey_files(ts,"","","labels_1.txt",vector<double>(5,10.));

	vector<double> weights = tracking.learn_from_funkey_files(sbrmr_binary,"features.txt","constraints.txt","labels.txt");
	
	for (int i=0; i<weights.size();i++){
	  cout << weights[i] << endl;
	}

	tracking.write_funkey_files(ts,"","","labels_2.txt",weights);


	// check if labels created by tracking match for 
	// 1)given   weights
	// 2)learned weights 
	// weights might differ but tracking result must be indentical !
	ifstream in("labels_1.txt");
	ifstream in2("labels_2.txt");
	while ((!in.eof()) && (!in2.eof())) {
		string line,line2;
		getline(in,line); 
		getline(in2,line2);
		BOOST_CHECK_EQUAL(line,line2);
        }



}


