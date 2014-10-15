#define BOOST_TEST_MODULE tracking_feature_extractor_test

#include <vector>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "pgmlink/reasoner_constracking.h"
#include "pgmlink/tracking.h"
#include "pgmlink/tracking_feature_extractor.h"

using namespace pgmlink;
using namespace pgmlink::features;

BOOST_AUTO_TEST_CASE( TrackingFeatureExtractor_SimpleMove ) {

    std::cout << "Constructing HypothesesGraph" << std::endl;
    std::cout << std::endl;

    using lemon::INVALID;

    std::cout << "Adding Traxels to TraxelStore" << std::endl;
    std::cout << std::endl;

    //  t=1      2
    //  o ------ o
    TraxelStore ts;
    Traxel n11, n21;
    feature_array com(feature_array::difference_type(3));
    feature_array divProb(feature_array::difference_type(1));
    n11.Id = 1; n11.Timestep = 1; com[0] = 1.0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
    n11.features["com"] = com; n11.features["divProb"] = divProb;
    add(ts,n11);
    n21.Id = 10; n21.Timestep = 2; com[0] = 0; com[1] = 1.0; com[2] = 0; divProb[0] = 0.1;
    n21.features["com"] = com; n21.features["divProb"] = divProb;
    add(ts,n21);

    std::cout << "Initialize Conservation tracking" << std::endl;
    std::cout << std::endl;

    FieldOfView fov(0, 0, 0, 0, 4, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
    ConsTracking tracking = ConsTracking(
                         2, // max_number_objects
                         false, // detection_by_volume
                         double(1.1), // avg_obj_size
                         20, // max_neighbor_distance
                         true, //with_divisions
                         0.3, // division_threshold
                         "none", // random_forest_filename
                         fov
                  );

    std::cout << "Run Conservation tracking" << std::endl;
    std::cout << std::endl;

    boost::shared_ptr<HypothesesGraph> hypotheses_graph = tracking.build_hypo_graph(ts);

    EventVectorVector events = tracking.track(
                                0, // forbidden_cost
                                0.0, // ep_gap
                                false, // with_tracklets
                                10.0, //division_weight
                                10.0, //transition_weight
                                1500., // disappearance_cost,
                                1500., // appearance_cost
                                false, //with_merger_resolution
                                3, //n_dim
                                5, //transition_parameter
                                0 //border_width for app/disapp costs
                                )[0];


    std::cout << "Computing higher order features" << std::endl;

    TrackingFeatureExtractor extractor(ts, events);
    extractor.compute_features();
    TrackingFeatureExtractor::JointFeatureVector joint_feature_vector;
    extractor.get_feature_vector(joint_feature_vector);

    for(size_t i = 0; i < joint_feature_vector.size(); i++)
    {
        std::cout << "Feature:\t" << extractor.get_feature_description(i) << "\t\tValue:\t" << joint_feature_vector[i] << std::endl;
        BOOST_CHECK_EQUAL(joint_feature_vector[i], 2.0);
    }
    std::cout << "done" << std::endl;
}

// EOF
