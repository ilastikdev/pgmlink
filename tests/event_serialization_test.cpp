#define BOOST_TEST_MODULE event_serialization_test


#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/test/unit_test.hpp>

#include "pgmlink/graph.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/feature.h"
#include "pgmlink/traxels.h"
#include "pgmlink/tracking.h"
#include "pgmlink/field_of_view.h"
#include "pgmlink/reasoner_constracking.h"

using namespace pgmlink;

BOOST_AUTO_TEST_CASE( Event_Serialization )
{

    // build example and run conservation tracking

    //  t=1      2      3       4
    //  o                       o
    //    |                    |
    //      ---- o ---- o ----
    //    |                    |
    //  o                       o
    TraxelStore ts;
    boost::shared_ptr<FeatureStore> fs = boost::make_shared<FeatureStore>();
    Traxel n11, n12, n21, n31, n41, n42;
    feature_array com(feature_array::difference_type(3));
    feature_array divProb(feature_array::difference_type(1));
    n11.Id = 1; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
    n11.features["com"] = com; n11.features["divProb"] = divProb;
    add(ts,fs, n11);
    n12.Id = 3; n12.Timestep = 1; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1;
    n12.features["com"] = com; n12.features["divProb"] = divProb;
    add(ts,fs, n12);
    n21.Id = 10; n21.Timestep = 2; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.1;
    n21.features["com"] = com; n21.features["divProb"] = divProb;
    add(ts,fs, n21);
    n31.Id = 11; n31.Timestep = 3; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.1;
    n31.features["com"] = com; n31.features["divProb"] = divProb;
    add(ts,fs, n31);
    n41.Id = 12; n41.Timestep = 4; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
    n41.features["com"] = com; n41.features["divProb"] = divProb;
    add(ts,fs, n41);
    n42.Id = 13; n42.Timestep = 4; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
    n42.features["com"] = com; n42.features["divProb"] = divProb;
    add(ts,fs, n42);

    std::cout << "Initialize Conservation tracking" << std::endl;    

    FieldOfView fov(0, 0, 0, 0, 4, 5, 5, 5); // tlow, xlow, ylow, zlow, tup, xup, yup, zup
    ConsTracking tracking = ConsTracking(
          2, // max_number_objects
          false, // detection_by_volume
          double(1.1), // avg_obj_size
          20, // max_neighbor_distance
          true, //with_divisions
          0.3, // division_threshold
          "none", // random_forest_filename
          fov,
          "none" // event_vector_dump_filename
             );


    std::cout << "run conservation tracking" << std::endl;

    EventVectorVectorVector events = tracking(ts,
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
                            0, //border_width for app/disapp costs
                            true); // with_constraints

    std::cout << "serialize events" << std::endl;

    EventVectorVectorVector events_loaded;

    {
        // now store the results
        std::ofstream ofs("temp_filename.evt");
        boost::archive::text_oarchive out_archive(ofs);
        out_archive << events;
    }

    std::cout << "read events" << std::endl;

    {
        // read again
        std::ifstream ifs("temp_filename.evt");
        BOOST_CHECK(ifs.good());
        boost::archive::text_iarchive in_archive(ifs);
        in_archive >> events_loaded;
    }

    BOOST_CHECK_EQUAL(events.size(), events_loaded.size());

    for(size_t i = 0; i < events.size(); i++)
    {
        BOOST_CHECK_EQUAL(events[i].size(), events_loaded[i].size());

        for(size_t j = 0; j < events[i].size(); j++)
        {
            for(size_t k = 0; k < events[i][j].size(); k++)
            {
                BOOST_CHECK(events[i][j][k] == events_loaded[i][j][k]);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( Traxelstore_Serialization_Test )
{
    std::cout << "Traxelstore_Serialization_Test" << std::endl;

    TraxelStore ts;
    boost::shared_ptr<FeatureStore> fs = boost::make_shared<FeatureStore>();
    Traxel n11, n12, n21, n31, n41, n42;
    feature_array com(feature_array::difference_type(3));
    feature_array divProb(feature_array::difference_type(1));
    n11.Id = 1; n11.Timestep = 1; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
    n11.features["com"] = com; n11.features["divProb"] = divProb;
    add(ts,fs, n11);
    n12.Id = 3; n12.Timestep = 1; com[0] = 2; com[1] = 2; com[2] = 2; divProb[0] = 0.1;
    n12.features["com"] = com; n12.features["divProb"] = divProb;
    add(ts, fs, n12);
    n21.Id = 10; n21.Timestep = 2; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.1;
    n21.features["com"] = com; n21.features["divProb"] = divProb;
    add(ts, fs, n21);
    n31.Id = 11; n31.Timestep = 3; com[0] = 1; com[1] = 1; com[2] = 1; divProb[0] = 0.1;
    n31.features["com"] = com; n31.features["divProb"] = divProb;
    add(ts, fs, n31);
    n41.Id = 12; n41.Timestep = 4; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
    n41.features["com"] = com; n41.features["divProb"] = divProb;
    add(ts, fs, n41);
    n42.Id = 13; n42.Timestep = 4; com[0] = 0; com[1] = 0; com[2] = 0; divProb[0] = 0.1;
    n42.features["com"] = com; n42.features["divProb"] = divProb;
    add(ts, fs, n42);

    TraxelStore ts_loaded;
    boost::shared_ptr<FeatureStore> fs_loaded;

    std::cout << "store traxelstore and featurestore in file" << std::endl;
    {
        // now store the results
        std::ofstream ofs("temp_filename.evt");
        boost::archive::text_oarchive out_archive(ofs);
        out_archive << ts;
        out_archive << fs;
    }

    std::cout << "read traxelstore and featurestore from file" << std::endl;
    {
        // read again
        std::ifstream ifs("temp_filename.evt");
        BOOST_CHECK(ifs.good());
        boost::archive::text_iarchive in_archive(ifs);
        std::cout << "reading traxelstore" << std::endl;
        in_archive >> ts_loaded;
        std::cout << "reading featurestore" << std::endl;
        in_archive >> fs_loaded;
        std::cout << "setting feature store" << std::endl;
        set_feature_store(ts_loaded, fs_loaded);
    }

    std::cout << "compare if the traxelstores are equal" << std::endl;

    BOOST_CHECK_EQUAL(ts.size(), ts_loaded.size());
    TraxelStore::iterator it_loaded = ts_loaded.begin();

    for(TraxelStore::iterator it = ts.begin();
        it != ts.end();
        it++)
    {
        BOOST_CHECK(it_loaded != ts_loaded.end());

        const Traxel& a = *it;
        const Traxel& b = *it_loaded;

        BOOST_CHECK(b.get_feature_store());

        BOOST_CHECK(a == b);
        BOOST_CHECK(fabsf(a.X() - b.X()) < 0.001f);
        BOOST_CHECK(fabsf(a.Y() - b.Y()) < 0.001f);
        BOOST_CHECK(fabsf(a.Z() - b.Z()) < 0.001f);

        it_loaded++;
    }

    // should not contain more elements than before
    BOOST_CHECK(it_loaded == ts_loaded.end());
}
