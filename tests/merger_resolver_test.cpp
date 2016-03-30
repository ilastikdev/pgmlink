#define BOOST_TEST_MODULE merger_resolver_test

#include <stdexcept>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <set>
#include <vector>
#include <iterator>

#include <boost/test/unit_test.hpp>

#include <vigra/multi_array.hxx>
#include <vigra/tinyvector.hxx>

#include "pgmlink/hypotheses.h"
#include "pgmlink/traxels.h"
// enable tests of private class members
#define private public
#include "pgmlink/merger_resolving.h"
#undef private

#include <armadillo>
#include <mlpack/core.hpp>
#include <mlpack/methods/kmeans/kmeans.hpp>
#include <lemon/maps.h>
#include <lemon/concepts/digraph.h>

using namespace pgmlink;
using namespace std;
using namespace boost;


BOOST_AUTO_TEST_CASE( MergerResolver_no_mergers )
{
    HypothesesGraph src;
    HypothesesGraph dest;
    TraxelStore ts;

    src.add(node_active())
    .add(arc_active())
    .add(node_active2())
    .add(node_traxel())
    .add(node_originated_from())
    .add(node_resolution_candidate())
    .add(arc_resolution_candidate())
    .add(arc_distance());

    MergerResolver m(&src,2);
    FeatureExtractorMCOMsFromGMM extractor(2);
    DistanceFromCOMs distance;
    FeatureHandlerFromTraxels handler(extractor, distance, &ts);
    m.resolve_mergers(handler);

    HypothesesGraph::Node n1 = src.add_node(1);
    HypothesesGraph::Node n2 = src.add_node(2);

    src.addArc(n1, n2);

    Parameter consTrackingParams = Parameter();

    resolve_graph(src,
                  dest,
                  consTrackingParams,
                  0.05, // ep gap
                  true, // with_tracklets
                  100, // transition_parameter
                  true, // with_constraint
                  boost::python::object(),
                  SolverType::CplexSolver,
                  2
                 );

    int n_arcs = lemon::countArcs(dest);
    int n_node = lemon::countNodes(dest);

    BOOST_CHECK_EQUAL(n_arcs, 0);
    BOOST_CHECK_EQUAL(n_node, 0);

}


BOOST_AUTO_TEST_CASE( MergerResolver_subgraph )
{
    HypothesesGraph g1, g2;
    g1.add(node_active()).add(arc_active()).add(node_active2()).add(node_traxel()).add(node_originated_from());
    HypothesesGraph::Node n1 = g1.add_node(1);
    HypothesesGraph::Node n2 = g1.add_node(2);
    HypothesesGraph::Node n3 = g1.add_node(2);
    property_map<node_active, HypothesesGraph::base_graph>::type& na_map = g1.get(node_active());
    property_map<node_active2, HypothesesGraph::base_graph>::type& na2_map = g1.get(node_active2());
    na_map.set(n1, true);
    na_map.set(n2, false);
    na_map.set(n3, true);
    na2_map.set(n1, 2);
    na2_map.set(n1, 3);
    na2_map.set(n1, 1);
    std::map<HypothesesGraph::Node, HypothesesGraph::Node> nr;
    std::map<HypothesesGraph::Arc, HypothesesGraph::Arc> ar;
    std::map<HypothesesGraph::Node, HypothesesGraph::Node> ncr;
    std::map<HypothesesGraph::Arc, HypothesesGraph::Arc> acr;
    copy_hypotheses_graph_subset<node_active, arc_active>(g1, g2, nr, ar, ncr, acr);
    g2.add(node_active()).add(node_active2());
    translate_property_bool_map<node_active, HypothesesGraph::Node>(g1, g2, nr);
    translate_property_value_map<node_active2, HypothesesGraph::Node>(g1, g2, nr);
    property_map<node_active, HypothesesGraph::base_graph>::type& n2a_map = g2.get(node_active());
    property_map<node_active2, HypothesesGraph::base_graph>::type& n2a2_map = g2.get(node_active2());

    std::map<HypothesesGraph::Node, HypothesesGraph::Node>::iterator it;
    for (it = nr.begin(); it != nr.end(); ++it)
    {
        std::cout << "Id (old): " << g1.id(it->first) <<  "(," <<  na_map[it->first] << "," << na2_map[it->first]
                  << "), Id (new): " << g2.id(it->second) << "," << n2a_map[it->second] << "," << n2a2_map[it->second] << "\n";
    }

    for (it = ncr.begin(); it != ncr.end(); ++it)
    {
        std::cout << "Id (old): " << g1.id(it->second) << "(," << na_map[it->second] << "," << na2_map[it->second]
                  << "), Id (new): " << g2.id(it->first) << ",(" << n2a_map[it->first] << "," << n2a2_map[it->first] << ")\n";
    }



}


BOOST_AUTO_TEST_CASE( MergerResolver_constructor )
{
    HypothesesGraph g;
    size_t ndim=2;
    BOOST_CHECK_THROW(MergerResolver m(&g,ndim), std::runtime_error);
    g.add(node_active2());
    BOOST_CHECK_THROW(MergerResolver m(&g,ndim), std::runtime_error);
    g.add(arc_active());
    BOOST_CHECK_THROW(MergerResolver m(&g,ndim), std::runtime_error);
    g.add(arc_distance());
    MergerResolver m(&g,ndim);
    BOOST_CHECK_EQUAL(m.g_, &g);
    // check that merger_resolved_to property has been added
    BOOST_CHECK(m.g_->has_property(merger_resolved_to()));
    BOOST_CHECK(m.g_->has_property(node_originated_from()));

    // check exception on intialization with null pointer
    HypothesesGraph* G = 0; // = hyp_builder.build();
    BOOST_CHECK_THROW(MergerResolver M(G,ndim), std::runtime_error);
}


BOOST_AUTO_TEST_CASE( MergerResolver_resolve_mergers_3 )
{
    LOG(logINFO) << "Starting test MergerResolver_resolve_mergers_3";
    //  t=1       2
    //       --- (2)
    //     |
    //   (3) --- (1)

    // -> each of the nodes in timesteps t has a possible arc to all nodes in t+1
    //    o ----- o
    //     |     |
    //      -----
    //     | | | |
    //    o --x-- o
    //     | | | |
    //      -----
    //     |     |
    //    o ----- o

    HypothesesGraph g;
    g.add(node_traxel()).add(arc_distance()).add(arc_active()).add(node_active2()).add(arc_active_count()).add(node_active_count());
    feature_array com(3, 0);
    feature_array pCOM(6 * 3, 0);
    TraxelStore ts;
    boost::shared_ptr<FeatureStore> fs = boost::make_shared<FeatureStore>();

    pCOM[0]  = 3;
    pCOM[3]  = 1;
    pCOM[6]  = 6;
    pCOM[9]  = 1;
    pCOM[12] = 3;
    pCOM[15] = 6;

    Traxel t11;
    t11.Timestep = 1;
    t11.Id = 11;
    com[0] = 3;
    t11.features["com"] = com;
    t11.features["possibleCOMs"] = pCOM;
    add(ts, fs, t11);

    Traxel t21;
    t21.Timestep = 2;
    t21.Id = 21;
    com[0] = 1.5;
    t21.features["com"] = com;
    t21.features["possibleCOMs"] = pCOM;
    add(ts, fs, t21);

    Traxel t22;
    t22.Timestep = 2;
    t22.Id = 22;
    com[0] = 3;
    t22.features["com"] = com;
    add(ts,fs,t22);

    HypothesesGraph::Node n11 = g.add_node(1);
    HypothesesGraph::Node n21 = g.add_node(2);
    HypothesesGraph::Node n22 = g.add_node(2);

    HypothesesGraph::Arc a11_21 = g.addArc(n11, n21);
    HypothesesGraph::Arc a11_22 = g.addArc(n11, n22);

    //initialize vectors
    property_map<arc_active_count, HypothesesGraph::base_graph>::type& active_arcs_count =
        g.get(arc_active_count());
    property_map<node_active_count, HypothesesGraph::base_graph>::type& active_nodes_count =
        g.get(node_active_count());
    
    for (HypothesesGraph::ArcIt a(g); a != lemon::INVALID; ++a)
    {
        std::cout << "Setting arc " << g.id(a) << " to {false}" << std::endl;
        active_arcs_count.set(a, {false});
    }
    for (HypothesesGraph::NodeIt n(g); n != lemon::INVALID; ++n)
    {
        active_nodes_count.set(n, {0});
    }

    property_map<node_traxel, HypothesesGraph::base_graph>::type& traxel_map = g.get(node_traxel());
    traxel_map.set(n11, t11);
    traxel_map.set(n21, t21);
    traxel_map.set(n22, t22);

    property_map<arc_active, HypothesesGraph::base_graph>::type& arc_map = g.get(arc_active());
    g.set_arc_active(a11_21, true);
    g.set_arc_active(a11_22, true);

    property_map<node_active2, HypothesesGraph::base_graph>::type& active_map = g.get(node_active2());
    g.set_node_active(n11, 3);
    g.set_node_active(n21, 2);
    g.set_node_active(n22, 1);

    property_map<arc_distance, HypothesesGraph::base_graph>::type& dist_map = g.get(arc_distance());

    size_t ndim=2;
    MergerResolver m(&g,ndim);
    FeatureExtractorMCOMsFromPCOMs extractor;
    DistanceFromCOMs distance;
    FeatureHandlerFromTraxels handler(extractor, distance, &ts);
    m.resolve_mergers(handler);
    prune_inactive(g);


    // check that arcs and nodes have been deactivated

    BOOST_CHECK(!g.valid(a11_22) || (g.source(a11_22) != n11 || g.target(a11_22) != n22));
    BOOST_CHECK(!g.valid(a11_21) || (g.source(a11_21) != n11 || g.target(a11_21) != n21));

    BOOST_CHECK(!g.valid(n11));
    BOOST_CHECK(!g.valid(n21));


    // check that traxel ids have been set correctly
    set<int> trx_ids_1;
    set<int> trx_ids_2;
    trx_ids_1.insert(12);
    trx_ids_1.insert(13);
    trx_ids_1.insert(14);
    trx_ids_1.insert(22);
    trx_ids_1.insert(23);
    trx_ids_1.insert(24);
    int trx_count = 0;
    property_map<node_active2, HypothesesGraph::base_graph>::type::ItemIt active_it(active_map, 1);
    for (; active_it != lemon::INVALID; ++active_it, ++trx_count)
    {
        trx_ids_2.insert(traxel_map[active_it].Id);

        HypothesesGraph::InArcIt IAIT(g, active_it);
        int count_n = 0;
        for(; IAIT != lemon::INVALID; ++IAIT, ++count_n)
        {
            BOOST_CHECK(arc_map[IAIT]);
        }
        HypothesesGraph::OutArcIt OAIT(g, active_it);
        for(; OAIT != lemon::INVALID; ++OAIT, ++count_n)
        {
            BOOST_CHECK(arc_map[OAIT]);
        }
        BOOST_CHECK_EQUAL(count_n, 3);
    }
    BOOST_CHECK_EQUAL_COLLECTIONS(trx_ids_1.begin(), trx_ids_1.end(), trx_ids_2.begin(), trx_ids_2.end());
    BOOST_CHECK_EQUAL(trx_count, 6);


    // check that distances are calculated correctly and active arcs are valid
    set<double> arc_dist_1;
    set<double> arc_dist_2;
    arc_dist_1.insert(0);
    arc_dist_1.insert(2);
    arc_dist_1.insert(3);
    arc_dist_1.insert(5);
    int arc_count = 0;

    property_map<arc_active, HypothesesGraph::base_graph>::type::ItemIt arc_it(arc_map, true);
    for(; arc_it != lemon::INVALID; ++arc_it, ++arc_count)
    {
        arc_dist_2.insert(dist_map[arc_it]);
        BOOST_CHECK(g.valid(arc_it));
    }

    BOOST_CHECK_EQUAL_COLLECTIONS(arc_dist_1.begin(), arc_dist_1.end(), arc_dist_2.begin(), arc_dist_2.end());
    BOOST_CHECK_EQUAL(arc_count, 9);

    // check that deactivated nodes are pruned, i.e. ItemIt(active_map, 0) should be equal  to lemon::INVALID
    property_map<node_active2, HypothesesGraph::base_graph>::type::ItemIt deactive_it(active_map, 0);
    BOOST_CHECK(!(deactive_it != lemon::INVALID));

    // check that deactivated arcs are pruned, i.e. FalseIt should be equal to lemon::INVALID
    property_map<arc_active, HypothesesGraph::base_graph>::type::FalseIt f_it(arc_map);
    BOOST_CHECK(!(f_it != lemon::INVALID));

    g.add(division_active());

    HypothesesGraph g_res;
    Parameter param=Parameter();
    resolve_graph(g, g_res, param, 0.05, false,5,true,boost::python::object(),SolverType::CplexSolver,ndim);
    prune_inactive(g);

    vector<vector<Event> > ev = *(events(g));
    unsigned resolve_count = 0;
    for (vector<vector<Event> >::iterator t_it = ev.begin(); t_it != ev.end(); ++t_it)
    {
        for (vector<Event>::iterator e_it = t_it->begin(); e_it != t_it->end(); ++ e_it)
        {
            cout << *e_it << "\n";
            if (e_it->type == Event::ResolvedTo)
            {
                ++resolve_count;
            }
        }
    }
    BOOST_CHECK_EQUAL(resolve_count, 1);

}


BOOST_AUTO_TEST_CASE( arma_mat_serialization ) {
    arma::mat orig(3,3);
    orig.randn();

    // save to string
    string s;
    {
      stringstream ss;
      boost::archive::text_oarchive oa(ss);
      oa & orig;
      s = ss.str();
    }

    // load from string and compare
    arma::mat loaded;
    {
      stringstream ss(s);
      boost::archive::text_iarchive ia(ss);
      ia & loaded;
    }

    BOOST_CHECK(arma::sum(arma::sum(loaded - orig)) < 0.0001);
}

// EOF
