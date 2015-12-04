#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

#include <string>
#include <sstream>

#include "../include/pgmlink/hypotheses.h"
#include "../include/pgmlink/features/higher_order_features.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/python.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/utility.hpp>

#include <lemon/core.h>

using namespace pgmlink;
using namespace boost::python;

typedef property_map<node_traxel, HypothesesGraph::base_graph>::type node_traxel_m;
typedef property_map<node_tracklet, HypothesesGraph::base_graph>::type NodeTrackletMap;
typedef property_map<arc_active, HypothesesGraph::base_graph>::type ArcActiveMap;
typedef property_map<node_active2, HypothesesGraph::base_graph>::type NodeActiveMap;
typedef property_map<division_active, HypothesesGraph::base_graph>::type DivisionActiveMap;
typedef property_map<node_timestep, HypothesesGraph::base_graph>::type NodeTimestepMap;
typedef property_map<node_origin_reference, HypothesesGraph::base_graph>::type NodeOriginReferenceMap;
typedef property_map<arc_origin_reference, HypothesesGraph::base_graph>::type ArcOriginReferenceMap;

std::vector< std::vector<Event> > get_events_of_graph(const HypothesesGraph& g)
{
    return *events(g, 0);
}

node_traxel_m& addNodeTraxelMap(HypothesesGraph* g)
{
    g->add(node_traxel());
    return g->get(node_traxel());
}
node_traxel_m& getNodeTraxelMap(HypothesesGraph* g)
{
    return g->get(node_traxel());
}

NodeTrackletMap& getNodeTrackletMap(HypothesesGraph* g)
{
    return g->get(node_tracklet());
}

std::vector<Traxel> get_item_NodeTrackletMap(NodeTrackletMap& map, const NodeTrackletMap::Key& k)
{
    return map[k];
}

NodeActiveMap& getNodeActiveMap(HypothesesGraph* g)
{
    return g->get(node_active2());
}

size_t get_item_NodeActiveMap(NodeActiveMap& map, const NodeActiveMap::Key& k)
{
    return map[k];
}

ArcActiveMap& getArcActiveMap(HypothesesGraph* g)
{
    return g->get(arc_active());
}

bool get_item_ArcActiveMap(ArcActiveMap& map, const ArcActiveMap::Key& k)
{
    return map[k];
}

DivisionActiveMap& getDivisionActiveMap(HypothesesGraph* g)
{
    return g->get(division_active());
}

bool get_item_DivisionActiveMap(DivisionActiveMap& map, const DivisionActiveMap::Key& k)
{
    return map[k];
}

NodeOriginReferenceMap& getNodeOriginReferenceMap(HypothesesGraph* g)
{
    return g->get(node_origin_reference());
}

ArcOriginReferenceMap& getArcOriginReferenceMap(HypothesesGraph* g)
{
    return g->get(arc_origin_reference());
}

NodeTimestepMap& getNodeTimestepMap(HypothesesGraph* g)
{
    return g->get(node_timestep());
}

int get_item_NodeTimestepMap(NodeTimestepMap& map, const NodeTimestepMap::Key& k)
{
    return map[k];
}

struct HypothesesGraph_pickle_suite : pickle_suite
{
    static std::string getstate( const HypothesesGraph& g )
    {
        std::stringstream ss;
        boost::archive::text_oarchive oa(ss);
        oa & g;
        return ss.str();
    }

    static void setstate( HypothesesGraph& g, const std::string& state )
    {
        std::stringstream ss(state);
        boost::archive::text_iarchive ia(ss);
        ia & g;
    }
};

inline object pass_through(object const& o)
{
    return o;
}
template < typename ITERABLE_VALUE_MAP >
struct IterableValueMap_ValueIterator
{
    typedef ITERABLE_VALUE_MAP map_type;
    IterableValueMap_ValueIterator( const map_type& m )
    {
        it_ = m.beginValue();
        end_ = m.endValue();
    }
    typename map_type::Value next()
    {
        if( it_ == end_)
        {
            PyErr_SetString(PyExc_StopIteration, "No more data.");
            boost::python::throw_error_already_set();
        }
        const typename map_type::Value& result = *it_;
        ++it_;
        return result;
    }

    static IterableValueMap_ValueIterator<map_type>
    values ( const map_type& m )
    {
        return IterableValueMap_ValueIterator<map_type>( m );
    }

    static void
    wrap( const char* python_name)
    {
        class_<IterableValueMap_ValueIterator<map_type> >( python_name, init<const map_type&>(args("iterable_value_map")) )
        .def("next", &IterableValueMap_ValueIterator::next)
        .def("__iter__", &pass_through)
        ;
    }

    typename map_type::ValueIt it_;
    typename map_type::ValueIt end_;
};

template <typename GRAPH_ELEMENT_ITERATOR , typename GRAPH_ELEMENT >
struct GraphIterator
{
    typedef GRAPH_ELEMENT_ITERATOR iter_type;
    typedef GRAPH_ELEMENT elem_type;
    iter_type it_;

    GraphIterator( const HypothesesGraph& g )
    {
        it_ = iter_type(g);
    }
    elem_type next()
    {
        if( it_ == lemon::INVALID)
        {
            PyErr_SetString(PyExc_StopIteration, "No more data.");
            boost::python::throw_error_already_set();
        }
        const elem_type result = it_;
        ++it_;
        return result;
    }

    // static GraphIterator<iter_type,elem_type>values ( const elem_type& n ) {
    // return GraphIterator<iter_type,map_type>( n );
    // }

    static void
    wrap( const char* python_name)
    {
        class_<GraphIterator<iter_type, elem_type> >( python_name, init<const HypothesesGraph&>(args("hypotheses_graph")))
        .def("next", &GraphIterator::next)
        .def("__iter__", &pass_through)
        ;
    }
};

int id1(const HypothesesGraph& g, const HypothesesGraph::Node& n)
{
    return g.id(n);
}
int id2(const HypothesesGraph& g, const HypothesesGraph::Arc& a)
{
    return g.id(a);
}

size_t num_active_incoming_arcs(const HypothesesGraph& g, const HypothesesGraph::Node& n)
{
    ArcActiveMap& arc_active_map = g.get(arc_active());
    size_t num_active = 0;

    for(HypothesesGraph::InArcIt a(g, n); a != lemon::INVALID; ++a)
    {
        if(arc_active_map[a])
        {
            num_active++;
        }
    }

    return num_active;
}

/**
 * @brief create a tracklet graph and return the new graph
 */
boost::shared_ptr<HypothesesGraph> pyGenerateTrackletGraph(const HypothesesGraph& traxel_graph)
{
    boost::shared_ptr<HypothesesGraph> tracklet_graph(new HypothesesGraph());
    generateTrackletGraph2(traxel_graph, *tracklet_graph);
    return tracklet_graph;
}

void export_hypotheses()
{
    class_<HypothesesGraph::Arc>("Arc");
    class_<HypothesesGraph::Node>("Node");
    class_<HypothesesGraph::InArcIt>("InArcIt");
    class_<HypothesesGraph::OutArcIt>("OutArcIt");

    GraphIterator<HypothesesGraph::NodeIt , HypothesesGraph::Node>::wrap("NodeIt" );
    // GraphIterator<HypothesesGraph::OutArcIt,HypothesesGraph::Arc >::wrap("OutArcIt");
    // GraphIterator<HypothesesGraph::InArcIt ,HypothesesGraph::Arc >::wrap("InArcIt" );
    GraphIterator<HypothesesGraph::ArcIt , HypothesesGraph::Arc >::wrap("ArcIt" );

    IterableValueMap_ValueIterator<node_traxel_m>::wrap("NodeTraxelMap_ValueIt");

    class_<node_traxel_m, boost::noncopyable>("NodeTraxelMap", init<const HypothesesGraph&>(args("hypotheses_graph")))
    .def("__getitem__", &node_traxel_m::operator[],
         return_internal_reference<>())
    .def("__setitem__", &node_traxel_m::set)
    .def("values", &IterableValueMap_ValueIterator<node_traxel_m>::values)
    ;

    // node/arc active maps
    class_< ArcActiveMap, boost::noncopyable >("ArcActiveMap", init<const HypothesesGraph&>(args("hypotheses_graph")))
    .def("__getitem__", &get_item_ArcActiveMap)
    .def("__setitem__", &ArcActiveMap::set);

    class_< NodeActiveMap, boost::noncopyable >("NodeActiveMap", init<const HypothesesGraph&>(args("hypotheses_graph")))
    .def("__getitem__", &get_item_NodeActiveMap)
    .def("__setitem__", &NodeActiveMap::set);

    class_< DivisionActiveMap, boost::noncopyable >("DivisionActiveMap", init<const HypothesesGraph&>(args("hypotheses_graph")))
    .def("__getitem__", &get_item_DivisionActiveMap)
    .def("__setitem__", &DivisionActiveMap::set);

    class_< NodeTimestepMap, boost::noncopyable >("NodeTimestepMap", init<const HypothesesGraph&>(args("hypotheses_graph")))
    .def("__getitem__", &get_item_NodeTimestepMap)
    .def("__setitem__", &NodeTimestepMap::set);

    class_< NodeTrackletMap, boost::noncopyable >("NodeTrackletMap", init<const HypothesesGraph&>(args("hypotheses_graph")))
    .def("__getitem__", &get_item_NodeTrackletMap)
    .def("__setitem__", &NodeTrackletMap::set);

    // node origin reference map
    class_< NodeOriginReferenceMap, boost::noncopyable >("NodeOriginReferenceMap", init<const HypothesesGraph&>(args("hypotheses_graph")))
    .def("__getitem__", &NodeOriginReferenceMap::operator[], return_internal_reference<>())
    .def("__setitem__", &NodeOriginReferenceMap::set);
    class_< ArcOriginReferenceMap, boost::noncopyable >("ArcOriginReferenceMap", init<const HypothesesGraph&>(args("hypotheses_graph")))
    .def("__getitem__", &ArcOriginReferenceMap::operator[], return_internal_reference<>())
    .def("__setitem__", &ArcOriginReferenceMap::set);

    // handle function overloading
    HypothesesGraph::Node (HypothesesGraph::*addnode1)(const int)
        = &HypothesesGraph::add_node;
    HypothesesGraph::Node (HypothesesGraph::*addnode2)(const std::vector<int>)
        = &HypothesesGraph::add_node;
    void (HypothesesGraph::*erase1)(const HypothesesGraph::Node) = &HypothesesGraph::erase;
    void (HypothesesGraph::*erase2)(const HypothesesGraph::Arc) = &HypothesesGraph::erase;
    bool (HypothesesGraph::*valid1)(const HypothesesGraph::Node) const = &HypothesesGraph::valid;
    bool (HypothesesGraph::*valid2)(const HypothesesGraph::Arc) const = &HypothesesGraph::valid;
//  int (*id1)(HypothesesGraph::Node) = &HypothesesGraph::id;
//  int (*id2)(HypothesesGraph::Arc) = &HypothesesGraph::id;
    HypothesesGraph::Node (HypothesesGraph::*baseNode1)(const HypothesesGraph::InArcIt&) const =
        &HypothesesGraph::baseNode;
    HypothesesGraph::Node (HypothesesGraph::*baseNode2)(const HypothesesGraph::OutArcIt&) const =
        &HypothesesGraph::baseNode;
    HypothesesGraph::Node (HypothesesGraph::*runningNode1)(const HypothesesGraph::InArcIt&) const =
        &HypothesesGraph::runningNode;
    HypothesesGraph::Node (HypothesesGraph::*runningNode2)(const HypothesesGraph::OutArcIt&) const =
        &HypothesesGraph::runningNode;

    class_<HypothesesGraph, shared_ptr<HypothesesGraph>, boost::noncopyable>("HypothesesGraph")
    .def("addNode", addnode1)
    .def("addNode", addnode2)
    //.def("timesteps", &HypothesesGraph::timesteps,
    //   return_internal_reference<>())
    .def("earliest_timestep", &HypothesesGraph::earliest_timestep)
    .def("latest_timestep", &HypothesesGraph::latest_timestep)

    // lemon graph interface
    .def("addArc", &HypothesesGraph::addArc)
    .def("erase", erase1)
    .def("erase", erase2)
    .def("valid", valid1)
    .def("valid", valid2)
    .def("changeTarget", &HypothesesGraph::changeTarget)
    .def("changeSource", &HypothesesGraph::changeSource)
    .def("target", &HypothesesGraph::target)
    .def("source", &HypothesesGraph::source)
    .def("id", &id1)
    .def("id", &id2)
    .def("nodeFromId", &HypothesesGraph::nodeFromId)
    .def("arcFromId", &HypothesesGraph::arcFromId)
    .def("maxNodeId", &HypothesesGraph::maxNodeId)
    .def("maxArcId", &HypothesesGraph::maxArcId)
    .def("baseNode", baseNode1)
    .def("runningNode", runningNode1)
    .def("baseNode", baseNode2)
    .def("runningNode", runningNode2)
    .def("oppositeNode", &HypothesesGraph::oppositeNode)

    // hypotheses graph specifics
    .def("addTraxel", &HypothesesGraph::add_traxel)
    .def("initLabelingMaps", &HypothesesGraph::init_labeling_maps)
    .def("addArcLabel" , &HypothesesGraph::add_arc_label )
    .def("addAppearanceLabel", &HypothesesGraph::add_appearance_label)
    .def("addDisappearanceLabel", &HypothesesGraph::add_disappearance_label)
    .def("addDivisionLabel", &HypothesesGraph::add_division_label)
    .def("set_solution", &features::set_solution)
    .def("set_injected_solution", &features::set_injected_solution)
    .def("write_hypotheses_graph_state", &HypothesesGraph::write_hypotheses_graph_state)
    .def("num_active_incoming_arcs", &num_active_incoming_arcs)
    .def("generate_tracklet_graph", &pyGenerateTrackletGraph)

    // extensions
    .def("addNodeTraxelMap", &addNodeTraxelMap, return_internal_reference<>())
    .def("getNodeTraxelMap", &getNodeTraxelMap, return_internal_reference<>())
    .def("getNodeTrackletMap", &getNodeTrackletMap, return_internal_reference<>())
    .def("getNodeActiveMap", &getNodeActiveMap, return_internal_reference<>())
    .def("getArcActiveMap", &getArcActiveMap, return_internal_reference<>())
    .def("getDivisionActiveMap", &getDivisionActiveMap, return_internal_reference<>())
    .def("getNodeOriginReferenceMap", &getNodeOriginReferenceMap, return_internal_reference<>())
    .def("getArcOriginReferenceMap", &getArcOriginReferenceMap, return_internal_reference<>())
    .def_pickle(HypothesesGraph_pickle_suite())
    ;

    def("getEventsOfGraph", &get_events_of_graph);

    // selector maps and copy functions for hypotheses graph, need to specify type of operator[] to be the const version or it wont compile
    typedef HypothesesGraph::base_graph::NodeMap<bool> NodeMask;
    class_< NodeMask, boost::noncopyable >("NodeMask", init<const HypothesesGraph&>(args("hypotheses_graph")))
    .def("__getitem__", (NodeMask::ConstReference (NodeMask::*)(const NodeMask::Key&) const)&NodeMask::operator[])
    .def("__setitem__", &NodeMask::set);

    typedef HypothesesGraph::base_graph::ArcMap<bool> ArcMask;
    class_< ArcMask, boost::noncopyable >("ArcMask", init<const HypothesesGraph&>(args("hypotheses_graph")))
    .def("__getitem__", (ArcMask::ConstReference (ArcMask::*)(const ArcMask::Key&) const)&ArcMask::operator[])
    .def("__setitem__", &ArcMask::set);

    def("copy_hypotheses_graph", &HypothesesGraph::copy, args("source_graph", "dest_graph"));
    def("copy_hypotheses_subgraph", &HypothesesGraph::copy_subgraph, args("source_graph", "dest_graph", "node_mask", "arc_mask"));

    //
    // lemon
    //
    int (*countNodes)(const HypothesesGraph&) = lemon::countNodes;
    def("countNodes", countNodes);

    int (*countArcs)(const HypothesesGraph&) = lemon::countArcs;
    def("countArcs", countArcs);
}
