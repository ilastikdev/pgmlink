/**
   @file
   @ingroup matching
   @brief Hypotheses Graph
*/

#ifndef HYPOTHESES_H
#define HYPOTHESES_H

#include <vector>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <map>
#include <boost/serialization/set.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>
#include <lemon/list_graph.h>
#include <lemon/maps.h>

#include "pgmlink/event.h"
#include "pgmlink/graph.h"
#include "pgmlink/log.h"
#include "pgmlink/pgmlink_export.h"
#include "pgmlink/traxels.h"

namespace pgmlink
{

////
//// IterableEditableValueMap
////
template <typename Graph, typename Key, typename Value>
class IterableEditableValueMap : public lemon::IterableValueMap<Graph, Key, Value>
{
private:
public:
    Value& get_value(const Key& key);
    explicit IterableEditableValueMap(const Graph& graph,
                                      const Value& value = Value());
};

template <typename Graph, typename Key, typename Value>
IterableEditableValueMap<Graph, Key, Value>::IterableEditableValueMap(const Graph& graph,
        const Value& value) :
    lemon::IterableValueMap<Graph, Key, Value>(graph, value)
{

}

template <typename Graph, typename Key, typename Value>
Value& IterableEditableValueMap<Graph, Key, Value>::get_value(const Key& key)
{
    return lemon::IterableValueMap<Graph, Key, Value>::Parent::operator[](key).value;
}


////
//// HypothesesGraph
////

// Properties of a HypothesesGraph

// node_timestep
struct node_timestep {};
template <typename Graph>
struct property_map<node_timestep, Graph>
{
    typedef lemon::IterableValueMap< Graph, typename Graph::Node, int > type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<node_timestep, Graph>::name = "node_timestep";

// node_traxel
struct node_traxel {};
class Traxel;
template <typename Graph>
struct property_map<node_traxel, Graph>
{
    typedef IterableEditableValueMap< Graph, typename Graph::Node, Traxel > type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<node_traxel, Graph>::name = "node_traxel";

// node_traxel
struct node_tracklet {};
template <typename Graph>
struct property_map<node_tracklet, Graph>
{
    typedef lemon::IterableValueMap< Graph, typename Graph::Node, std::vector<Traxel> > type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<node_tracklet, Graph>::name = "node_tracklet";


// tracklet_arcs
struct tracklet_intern_dist {};
template <typename Graph>
struct property_map<tracklet_intern_dist, Graph>
{
    typedef lemon::IterableValueMap< Graph, typename Graph::Node, std::vector<double> > type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<tracklet_intern_dist, Graph>::name = "tracklet_intern_dist";

// tracklet_arcs
struct tracklet_intern_arc_ids {};
template <typename Graph>
struct property_map<tracklet_intern_arc_ids, Graph>
{
    typedef lemon::IterableValueMap< Graph, typename Graph::Node, std::vector<int> > type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<tracklet_intern_arc_ids, Graph>::name = "tracklet_intern_arc_ids";

// node_active_count
struct node_active_count {};
template <typename Graph>
struct property_map<node_active_count, Graph>
{
    typedef IterableEditableValueMap< Graph, typename Graph::Node, std::vector<size_t> > type;
    //typedef lemon::IterableIntMap< Graph, typename Graph::Node> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<node_active_count, Graph>::name = "node_active_count";

// arc_active_count
struct arc_active_count {};
template <typename Graph>
struct property_map<arc_active_count, Graph>
{
    typedef IterableEditableValueMap< Graph, typename Graph::Arc, std::vector<bool> > type;
    //typedef lemon::IterableIntMap< Graph, typename Graph::Arc> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<arc_active_count, Graph>::name = "arc_active_count";

// arc_value_count
struct arc_value_count {};
template <typename Graph>
struct property_map<arc_value_count, Graph>
{
    typedef IterableEditableValueMap< Graph, typename Graph::Arc, std::vector<size_t> > type;
    //typedef lemon::IterableIntMap< Graph, typename Graph::Arc> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<arc_value_count, Graph>::name = "arc_value_count";

// division_active_count
struct division_active_count {};
template <typename Graph>
struct property_map<division_active_count, Graph>
{
    typedef IterableEditableValueMap< Graph, typename Graph::Node, std::vector<bool> > type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<division_active_count, Graph>::name = "division_active_count";

// relative_uncertainty
struct relative_uncertainty {};
template <typename Graph>
struct property_map<relative_uncertainty, Graph>
{
    typedef IterableEditableValueMap< Graph, typename Graph::Node, double > type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<relative_uncertainty, Graph>::name = "relative_uncertainty";


// node_active
struct node_active {};
template <typename Graph>
struct property_map<node_active, Graph>
{
    typedef lemon::IterableBoolMap< Graph, typename Graph::Node> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<node_active, Graph>::name = "node_active";

// node_active2
struct node_active2 {};
template <typename Graph>
struct property_map<node_active2, Graph>
{
    typedef lemon::IterableValueMap< Graph, typename Graph::Node, std::size_t> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<node_active2, Graph>::name = "node_active2";

// node_offered
struct node_offered {};
template <typename Graph>
struct property_map<node_offered, Graph>
{
    typedef lemon::IterableBoolMap< Graph, typename Graph::Node> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<node_offered, Graph>::name = "node_offered";

// arc_distance
struct arc_distance {};
template <typename Graph>
struct property_map<arc_distance, Graph>
{
    typedef lemon::IterableValueMap< Graph, typename Graph::Arc, double> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<arc_distance, Graph>::name = "arc_distance";

// traxel_arc_id
struct traxel_arc_id {};
template <typename Graph>
struct property_map<traxel_arc_id, Graph>
{
    typedef lemon::IterableValueMap< Graph, typename Graph::Arc, int> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<traxel_arc_id, Graph>::name = "traxel_arc_id";

struct arc_vol_ratio {};
template <typename Graph>
struct property_map<arc_vol_ratio, Graph>
{
    typedef lemon::IterableValueMap< Graph, typename Graph::Arc, double> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<arc_vol_ratio, Graph>::name = "arc_vol_ratio";

// split_into
struct split_from {};
template <typename Graph>
struct property_map<split_from, Graph>
{
    typedef lemon::IterableValueMap< Graph, typename Graph::Node, int> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<split_from, Graph>::name = "split_from";

// arc_from_timestep
struct arc_from_timestep {};
template <typename Graph>
struct property_map<arc_from_timestep, Graph>
{
    typedef typename Graph::template ArcMap<int> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<arc_from_timestep, Graph>::name = "arc_from_timestep";

// arc_to_timestep
struct arc_to_timestep {};
template <typename Graph>
struct property_map<arc_to_timestep, Graph>
{
    typedef typename Graph::template ArcMap<int> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<arc_to_timestep, Graph>::name = "arc_to_timestep";

// arc_active
struct arc_active {};
template <typename Graph>
struct property_map<arc_active, Graph>
{
    typedef lemon::IterableBoolMap< Graph, typename Graph::Arc> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<arc_active, Graph>::name = "arc_active";

// division_active
struct division_active {};
template <typename Graph>
struct property_map<division_active, Graph>
{
    typedef lemon::IterableBoolMap< Graph, typename Graph::Node> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<division_active, Graph>::name = "division_active";

// merger_resolved_to
struct merger_resolved_to {};
template <typename Graph>
struct property_map<merger_resolved_to, Graph>
{
    // typedef std::map<typename Graph::Node, std::vector<unsigned int> > type;
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, std::vector<unsigned int> > type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<merger_resolved_to, Graph>::name = "merger_resolved_to";

// node_originated_from
struct node_originated_from {};
template <typename Graph>
struct property_map<node_originated_from, Graph>
{
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, std::vector<unsigned int> > type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<node_originated_from, Graph>::name = "node_originated_from";

// node_resolution_candidate
struct node_resolution_candidate {};
template <typename Graph>
struct property_map<node_resolution_candidate, Graph>
{
    typedef lemon::IterableBoolMap<Graph, typename Graph::Node> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<node_resolution_candidate, Graph>::name = "node_resolution_candidate";

// arc_resolution_candidate
struct arc_resolution_candidate {};
template <typename Graph>
struct property_map<arc_resolution_candidate, Graph>
{
    typedef lemon::IterableBoolMap<Graph, typename Graph::Arc> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<arc_resolution_candidate, Graph>::name = "arc_resolution_candidate";

typedef size_t label_type;

// appearance ground truth label
struct appearance_label {};
template <typename Graph>
struct property_map< appearance_label, Graph>
{
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, label_type> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<appearance_label, Graph>::name = "appearance_label";

// disappearance ground truth label
struct disappearance_label {};
template <typename Graph>
struct property_map< disappearance_label, Graph>
{
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, label_type> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<disappearance_label, Graph>::name = "disappearance_label";

// division ground truth label
struct division_label {};
template <typename Graph>
struct property_map< division_label, Graph>
{
    typedef lemon::IterableValueMap<Graph, typename Graph::Node, label_type> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<division_label, Graph>::name = "division_label";

// arc ground truth label
struct arc_label {};
template <typename Graph>
struct property_map< arc_label, Graph>
{
    typedef lemon::IterableValueMap<Graph, typename Graph::Arc, label_type> type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<arc_label, Graph>::name = "arc_label";


struct node_origin_reference {};
template <typename Graph>
struct property_map<node_origin_reference, Graph>
{
    typedef lemon::IterableValueMap< Graph, typename Graph::Node, typename Graph::Node > type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<node_origin_reference, Graph>::name = "node_origin_reference";

struct arc_origin_reference {};
template <typename Graph>
struct property_map<arc_origin_reference, Graph>
{
    typedef lemon::IterableValueMap< Graph, typename Graph::Arc, typename Graph::Arc > type;
    static const std::string name;
};
template <typename Graph>
const std::string property_map<arc_origin_reference, Graph>::name = "arc_origin_reference";

class HypothesesGraph
    : public PropertyGraph<lemon::ListDigraph>
{
public:
    typedef property_map<node_timestep, HypothesesGraph::base_graph>::type node_timestep_map;
    typedef PropertyGraph<lemon::ListDigraph>::base_graph base_graph;

    PGMLINK_EXPORT HypothesesGraph()
    {
        // Properties attached to every HypothesesGraph
        add(node_timestep());
        add(arc_from_timestep());
        add(arc_to_timestep());
    }

    // use this instead of calling the parent graph directly
    PGMLINK_EXPORT HypothesesGraph::Node add_node(node_timestep_map::Value timestep);
    // call this function to add a multi-temporal node (e.g. for tracklets)
    PGMLINK_EXPORT HypothesesGraph::Node add_node(std::vector<node_timestep_map::Value> timesteps);

    // override
    PGMLINK_EXPORT HypothesesGraph::Arc addArc(HypothesesGraph::Node s, HypothesesGraph::Node t);

    // add node with associated traxel
    PGMLINK_EXPORT HypothesesGraph::Node add_traxel(Traxel ts);

    // assign ground truth to a node
    PGMLINK_EXPORT void add_appearance_label(HypothesesGraph::Node, label_type label);

    // assign ground truth to a node
    PGMLINK_EXPORT void add_disappearance_label(HypothesesGraph::Node, label_type label);

    // assign division ground truth to a node
    PGMLINK_EXPORT void add_division_label(HypothesesGraph::Node, label_type label);
    PGMLINK_EXPORT label_type get_division_label(HypothesesGraph::Node);

    // assign ground truth to a node
    PGMLINK_EXPORT void add_arc_label( HypothesesGraph::Arc  , label_type label);

    // getters and setters for active states
    PGMLINK_EXPORT size_t get_node_active(HypothesesGraph::Node n, int iteration = -1) const;
    PGMLINK_EXPORT bool get_division_active(HypothesesGraph::Node n, int iteration = -1) const;
    PGMLINK_EXPORT bool get_arc_active(HypothesesGraph::Arc a, int iteration = -1) const;
    PGMLINK_EXPORT void set_node_active(HypothesesGraph::Node n, size_t newState, int iteration = -1);
    PGMLINK_EXPORT void set_division_active(HypothesesGraph::Node n, bool newState, int iteration = -1);
    PGMLINK_EXPORT void set_arc_active(HypothesesGraph::Arc a, bool newState, int iteration = -1);

    PGMLINK_EXPORT const std::set<HypothesesGraph::node_timestep_map::Value>& timesteps() const;
    PGMLINK_EXPORT node_timestep_map::Value earliest_timestep() const;
    PGMLINK_EXPORT node_timestep_map::Value latest_timestep() const;

    PGMLINK_EXPORT void init_labeling_maps();
    PGMLINK_EXPORT void write_hypotheses_graph_state(const std::string out_fn);

    static void copy(HypothesesGraph& src, HypothesesGraph& dest);
    static void copy_subgraph(HypothesesGraph &src, HypothesesGraph &dest,
                              HypothesesGraph::base_graph::NodeMap<bool> &selected_nodes,
                              HypothesesGraph::base_graph::ArcMap<bool> &selected_arcs);

    void save_to_graphviz_dot_file(const std::string& filename,
                                   bool with_tracklets,
                                   bool with_divisions,
                                   boost::function<double (const Traxel&, const size_t)> detection,
                                   boost::function<double (const Traxel&, const size_t)> division,
                                   boost::function<double (const double)> transition,
                                   //boost::function<double (const Traxel&, const Traxel&, const size_t)> transition,
                                   boost::function<double (const Traxel&)> disappearance_cost,
                                   boost::function<double (const Traxel&)> appearance_cost,
                                   size_t max_number_objects,
                                   double transition_parameter) const;
private:
    void initialize_node(HypothesesGraph::Node n);

    // boost serialize
    friend class boost::serialization::access;
    template< typename Archive >
    void save( Archive&, const unsigned int /*version*/ ) const;
    template< typename Archive >
    void load( Archive&, const unsigned int /*version*/ );
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    std::set<node_timestep_map::Value> timesteps_;
};

PGMLINK_EXPORT void generateTrackletGraph(const HypothesesGraph& traxel_graph, HypothesesGraph& tracklet_graph);
PGMLINK_EXPORT std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > generateTrackletGraph2(
    const HypothesesGraph& traxel_graph, HypothesesGraph& tracklet_graph);
PGMLINK_EXPORT HypothesesGraph& prune_inactive(HypothesesGraph&);
PGMLINK_EXPORT boost::shared_ptr<std::vector< std::vector<Event> > > events(const HypothesesGraph& g, int iterationStep = 0);
PGMLINK_EXPORT boost::shared_ptr<std::vector< std::vector<Event> > > multi_frame_move_events(const HypothesesGraph& g);
PGMLINK_EXPORT boost::shared_ptr<std::vector< std::vector<Event> > > resolved_to_events(const HypothesesGraph& g);
PGMLINK_EXPORT EventVectorVector merge_event_vectors(const EventVectorVector& ev1, const EventVectorVector& ev2);
PGMLINK_EXPORT boost::shared_ptr<std::vector< std::map<unsigned int, bool> > > state_of_nodes(const HypothesesGraph&);
// prune to specific start nodes
PGMLINK_EXPORT HypothesesGraph& prune_to_node_descendants(HypothesesGraph& graph, const std::vector<HypothesesGraph::Node>& start_nodes);
PGMLINK_EXPORT HypothesesGraph& set_descendants_active(HypothesesGraph& graph, const HypothesesGraph::Node& start_node);

// lemon graph format (lgf) serialization
PGMLINK_EXPORT void write_lgf(const HypothesesGraph&, std::ostream& os,
                              std::map<std::string, bool> &config);
PGMLINK_EXPORT void read_lgf(HypothesesGraph&, std::istream& is,
                             std::map<std::string, bool> &config);

////
//// HypothesesBuilder
////
class HypothesesBuilder
{
public:
    PGMLINK_EXPORT virtual HypothesesGraph* build() const;

protected:
    // template methods
    PGMLINK_EXPORT virtual HypothesesGraph* construct() const = 0;
    PGMLINK_EXPORT virtual HypothesesGraph* add_nodes(HypothesesGraph*) const = 0;
    PGMLINK_EXPORT virtual HypothesesGraph* add_edges(HypothesesGraph*) const = 0;
};



////
//// SingleTimestepTraxel_HypothesesBuilder
////
class SingleTimestepTraxel_HypothesesBuilder
    : public HypothesesBuilder
{
public:
    struct Options
    {
        PGMLINK_EXPORT Options(unsigned int mnn = 6, double dt = 50,
                               bool forward_backward = false, bool consider_divisions = false,
                               double division_threshold = 0.5)
            : max_nearest_neighbors(mnn), distance_threshold(dt), forward_backward(forward_backward),
              consider_divisions(consider_divisions),
              division_threshold(division_threshold)
        {}

        unsigned int max_nearest_neighbors;
        double distance_threshold;
        bool forward_backward, consider_divisions;
        double division_threshold;
    };

    PGMLINK_EXPORT SingleTimestepTraxel_HypothesesBuilder(const TraxelStore* ts, const Options& o = Options())
        : ts_(ts), options_(o)
    {}

protected:
    // builder method implementations
    PGMLINK_EXPORT virtual HypothesesGraph* construct() const;
    PGMLINK_EXPORT virtual HypothesesGraph* add_nodes(HypothesesGraph*) const;
    PGMLINK_EXPORT virtual HypothesesGraph* add_edges(HypothesesGraph*) const;

    const TraxelStore* ts_;
    Options options_;
private:
    HypothesesGraph* add_edges_at(HypothesesGraph*, int timestep, bool reverse = false) const;
};




/**/
/* implementation */
/**/
template< typename Archive >
void HypothesesGraph::save( Archive& ar, const unsigned int /*version*/ ) const
{
    ar & timesteps_;

    // store config
    std::map<std::string, bool> config;
    config["node_timestep"] = this->has_property(node_timestep());
    config["node_active"] = this->has_property(node_active());
    config["node_active2"] = this->has_property(node_active2());
    config["node_active_count"] = this->has_property(node_active_count());
    config["node_offered"] = this->has_property(node_offered());
    config["split_from"] = this->has_property(split_from());
    config["division_active"] = this->has_property(division_active());
    config["merger_resolved_to"] = this->has_property(merger_resolved_to());
    config["node_originated_from"] = this->has_property(node_originated_from());
    config["node_resolution_candidate"] = this->has_property(node_resolution_candidate());
    config["arc_distance"] = this->has_property(arc_distance());
    config["traxel_arc_id"] = this->has_property(traxel_arc_id());
    config["arc_vol_ratio"] = this->has_property(arc_vol_ratio());
    config["arc_from_timestep"] = this->has_property(arc_from_timestep());
    config["arc_to_timestep"] = this->has_property(arc_to_timestep());
    config["arc_active"] = this->has_property(arc_active());
    config["arc_resolution_candidate"] = this->has_property(arc_resolution_candidate());
    config["tracklet_intern_dist"] = this->has_property(tracklet_intern_dist());
    config["tracklet_intern_arc_ids"] = this->has_property(tracklet_intern_arc_ids());
    config["arc_active_count"] = this->has_property(arc_active_count());
    config["node_traxel"] = this->has_property(node_traxel());
    config["node_tracklet"] = this->has_property(node_tracklet());

    ar & config;

    std::string lgf;
    {
        std::stringstream ss;
        write_lgf(*this, ss, config);
        lgf = ss.str();
    }
    ar & lgf;
}

template< typename Archive >
void HypothesesGraph::load( Archive& ar, const unsigned int /*version*/ )
{
    ar & timesteps_;

    // read config first
    std::map<std::string, bool> config;
    ar & config;

    std::string lgf;
    ar & lgf;
    {
        std::stringstream ss(lgf);
        read_lgf(*this, ss, config);
    }
}

}
#endif /* HYPOTHESES_H */
