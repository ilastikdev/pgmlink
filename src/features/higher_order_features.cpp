#include "pgmlink/features/higher_order_features.h"
#include <cmath> /* for sqrt */

#include <vigra/multi_math.hxx> /* for operator+ */
#include <vigra/matrix.hxx> /* for covariance calculation */
#include <vigra/linear_algebra.hxx> /* for matrix inverse calculation */

#include <numeric> /* for accumulate */

#include <sstream> /* for printing a matrix on the debug level */

#include <algorithm> /* for std::copy, std::min */

namespace pgmlink
{
namespace features
{
////
//// Some useful typedefinitions
////
typedef typename
property_map<node_active, HypothesesGraph::base_graph>::type
node_active_map_type;
typedef typename
property_map<node_active2, HypothesesGraph::base_graph>::type
node_active2_map_type;
typedef typename
property_map<arc_active, HypothesesGraph::base_graph>::type
arc_active_map_type;
typedef typename
property_map<node_timestep, HypothesesGraph::base_graph>::type
node_timestep_map_type;
typedef typename
property_map<node_active_count, HypothesesGraph::base_graph>::type
nodes_active_map_type;
typedef typename
property_map<arc_active_count, HypothesesGraph::base_graph>::type
arcs_active_map_type;
typedef typename
property_map<division_active_count, HypothesesGraph::base_graph>::type
divs_active_count_map_type;
typedef typename
property_map<division_active, HypothesesGraph::base_graph>::type
division_active_map_type;
typedef typename
property_map<node_traxel, HypothesesGraph::base_graph>::type
node_traxel_map_type;
typedef typename
property_map<node_tracklet, HypothesesGraph::base_graph>::type
node_tracklet_map_type;
typedef typename HypothesesGraph::NodeIt NodeIt;
typedef typename HypothesesGraph::ArcIt ArcIt;
typedef typename HypothesesGraph::InArcIt InArcIt;
typedef typename HypothesesGraph::OutArcIt OutArcIt;

typedef typename node_active_map_type::TrueIt NodeActiveIt;
typedef typename arc_active_map_type::TrueIt ArcActiveIt;


/*=============================================================================
 functions
=============================================================================*/
////
//// function set_solution
////
void set_solution(HypothesesGraph& graph, const size_t solution_index)
{
    // check if the graph has the necessary property maps
    if (not graph.has_property(node_active_count()))
    {
        throw std::runtime_error(
            "Graph doesn't have a \"node_active_count\" property map"
        );
    }
    if (not graph.has_property(arc_active_count()))
    {
        throw std::runtime_error(
            "Graph doesn't have an \"arc_active_count\" property map"
        );
    }

    // Get the property maps
    nodes_active_map_type& nodes_active_map = graph.get(node_active_count());
    arcs_active_map_type& arcs_active_map = graph.get(arc_active_count());
    divs_active_count_map_type& divisions_active_map = graph.get(division_active_count());

    // create the the node_active, node_active2 and arc_active map
    if (not graph.has_property(node_active()))
    {
        graph.add(node_active());
    }
    if (not graph.has_property(node_active2()))
    {
        graph.add(node_active2());
    }
    if (not graph.has_property(arc_active()))
    {
        graph.add(arc_active());
    }
    if(not graph.has_property(division_active()))
    {
        graph.add(division_active());
    }

    // Get the property maps (caution: "node_active_map" not "nodes_active_map")
    node_active_map_type& node_active_map = graph.get(node_active());
    node_active2_map_type& node_active2_map = graph.get(node_active2());
    arc_active_map_type& arc_active_map = graph.get(arc_active());
    division_active_map_type& div_active_map = graph.get(division_active());

    // Now we can start the with writing the solution with the index
    // solution_index into the node_active_map
    for (NodeIt n_it(graph); n_it != lemon::INVALID; ++n_it)
    {
        if (nodes_active_map[n_it].size() <= solution_index)
        {
            throw std::runtime_error(
                "In set_solution(): Solution index out of range"
            );
        }
        node_active_map[n_it] = nodes_active_map[n_it][solution_index];
        node_active2_map.set(n_it, nodes_active_map[n_it][solution_index]);
        div_active_map[n_it] = divisions_active_map[n_it][solution_index];
    }
    for (ArcIt a_it(graph); a_it != lemon::INVALID; ++a_it)
    {
        if (arcs_active_map[a_it].size() <= solution_index)
        {
            throw std::runtime_error(
                "In set_solution(): Solution index out of range"
            );
        }
        arc_active_map[a_it] = arcs_active_map[a_it][solution_index];
    }
}

void set_injected_solution(HypothesesGraph& graph)
{
    // check if the graph has the necessary property maps
    if (!graph.has_property(appearance_label()))
    {
        throw std::runtime_error("Graph doesn't have a \"appearance_label\" property map");
    }

    if (not graph.has_property(disappearance_label()))
    {
        throw std::runtime_error("Graph doesn't have an \"disappearance_label\" property map");
    }

    if (not graph.has_property(division_label()))
    {
        throw std::runtime_error("Graph doesn't have an \"division_label\" property map");
    }

    if (not graph.has_property(arc_label()))
    {
        throw std::runtime_error("Graph doesn't have an \"arc_label\" property map");
    }

    // get ground truth label maps
    property_map<appearance_label, HypothesesGraph::base_graph>::type& gt_appearance = graph.get(appearance_label());
    property_map<disappearance_label, HypothesesGraph::base_graph>::type& gt_disappearance = graph.get(disappearance_label());
    property_map<division_label, HypothesesGraph::base_graph>::type& gt_division = graph.get(division_label());
    property_map<arc_label, HypothesesGraph::base_graph>::type& gt_arc = graph.get(arc_label());

    // if not present, add active maps
    if (not graph.has_property(node_active()))
    {
        graph.add(node_active());
    }
    if (not graph.has_property(node_active2()))
    {
        graph.add(node_active2());
    }
    if (not graph.has_property(arc_active()))
    {
        graph.add(arc_active());
    }
    if(not graph.has_property(division_active()))
    {
        graph.add(division_active());
    }

    if(not graph.has_property(division_active_count()))
    {
        graph.add(division_active_count());
    }
    if (not graph.has_property(node_active_count()))
    {
        graph.add(node_active_count());
    }
    if (not graph.has_property(arc_active_count()))
    {
        graph.add(arc_active_count());
    }

    // Get the property maps (caution: "node_active_map" not "nodes_active_map")
    node_active_map_type& node_active_map = graph.get(node_active());
    node_active2_map_type& node_active2_map = graph.get(node_active2());
    arc_active_map_type& arc_active_map = graph.get(arc_active());
    division_active_map_type& div_active_map = graph.get(division_active());

    property_map<node_active_count, HypothesesGraph::base_graph>::type& node_active_count_map = graph.get(node_active_count());
    property_map<arc_active_count, HypothesesGraph::base_graph>::type& arc_active_count_map = graph.get(arc_active_count());
    property_map<division_active_count, HypothesesGraph::base_graph>::type& div_active_count_map = graph.get(division_active_count());

    // insert solution
    bool is_active;
    for (NodeIt n_it(graph); n_it != lemon::INVALID; ++n_it)
    {
        is_active = gt_appearance[n_it] || gt_disappearance[n_it];
        node_active_map[n_it] = is_active;
        node_active2_map.set(n_it, std::max(gt_appearance[n_it], gt_disappearance[n_it]));
        div_active_map[n_it] = gt_division[n_it];
        div_active_count_map.set(n_it, { gt_division[n_it] });
        node_active_count_map.set(n_it, { std::max(gt_appearance[n_it], gt_disappearance[n_it]) });
    }

    for (ArcIt a_it(graph); a_it != lemon::INVALID; ++a_it)
    {
        arc_active_map[a_it] = gt_arc[a_it];
        arc_active_count_map.set(a_it, { gt_arc[a_it] });
    }
}

////
//// function get_out_nodes
////
void get_out_nodes(
    const HypothesesGraph::Node& node,
    const HypothesesGraph& graph,
    Nodevector& out_nodes
)
{
    out_nodes.clear();
    for (OutArcIt oa_it(graph, node); oa_it != lemon::INVALID; ++oa_it)
    {
        if (graph.get(arc_active())[oa_it])
        {
            out_nodes.push_back(graph.target(oa_it));
        }
    }
}

////
//// function get_in_nodes
////
void get_in_nodes(
    const HypothesesGraph::Node& node,
    const HypothesesGraph& graph,
    Nodevector& in_nodes
)
{
    in_nodes.clear();
    for (InArcIt ia_it(graph, node); ia_it != lemon::INVALID; ++ia_it)
    {
        if (graph.get(arc_active())[ia_it])
        {
            in_nodes.push_back(graph.source(ia_it));
        }
    }
}

/*=============================================================================
  virtual classes
=============================================================================*/
FeatureMatrix TraxelsFeatureExtractor::extract(
    const ConstTraxelRefVector& traxelrefs
) const
{
    FeatureMatrix ret_;
    extract(traxelrefs, ret_);
    return ret_;
}

FeatureMatrix TraxelsFeatureCalculator::calculate(
    const FeatureMatrix& feature_matrix
) const
{
    FeatureMatrix ret_;
    calculate(feature_matrix, ret_);
    return ret_;
}

void TraxelsFeatureCalculator::calculate_for_all(
    const ConstTraxelRefVectors& traxelrefs,
    FeatureMatrix& return_matrix,
    boost::shared_ptr<TraxelsFeatureExtractor> feature_extractor_ref
) const
{
    size_t col_count = traxelrefs.size();
    // check if the vector of traxels is empty
    if (col_count == 0)
    {
        LOG(logDEBUG) << "In " << name() << ": vector of traxels is empty";
        LOG(logDEBUG) << "Returning zero";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0.0);
    }
    else
    {
        // allocate memory for result of feature extractor and feature calculator
        FeatureMatrix extr_features;
        FeatureMatrix calc_features;
        // row count will contain the size of the row_count of the return matrix
        size_t row_count = 0;
        size_t current_col = 0;
        // iterate over all subsets of traxels and write the flattened result of
        // the feature calculator for the n-th subset into the n-th column of the
        // return matrix
        for (
            ConstTraxelRefVectors::const_iterator traxelrefs_it = traxelrefs.begin();
            traxelrefs_it != traxelrefs.end();
            traxelrefs_it++, current_col++
        )
        {
            // do the feature extraction and calculation
            feature_extractor_ref->extract(*traxelrefs_it, extr_features);
            calculate(extr_features, calc_features);
            if (row_count == 0)
            {
                // if we never had an feature matrix with a size larger than 0 get the
                // current feature size
                row_count = calc_features.size();
                if (row_count != 0)
                {
                    // if this feature size is not zero set the row count of the return
                    // matrix to this size
                    return_matrix.reshape(vigra::Shape2(col_count, row_count));
                    return_matrix.init(0.0);
                }
            }
            if (row_count != 0)
            {
                // write the result if the row count of the return matrix is set
                // only write the first <row_max> values from the flattened feature
                // matrix into the current row of the return matrix
                size_t row_max = std::min(row_count, size_t(calc_features.size()));
                FeatureMatrix::iterator calc_features_it = calc_features.begin();
                for (
                    size_t current_row = 0;
                    current_row < row_max;
                    current_row++, calc_features_it++
                )
                {
                    return_matrix(current_col, current_row) = *calc_features_it;
                }
            }
        }
        if (row_count == 0)
        {
            LOG(logDEBUG) << "In " << name()
                          << ": all feature calculators returned empty matrices";
            LOG(logDEBUG) << "Returning zero";
            return_matrix.reshape(vigra::Shape2(1, 1));
            return_matrix.init(0.0);
        }
    }
}

/*=============================================================================
  specific classes
=============================================================================*/
////
//// class GraphFeatureCalculator
////
void GraphFeatureCalculator::calculate(
    const HypothesesGraph& graph,
    FeatureMatrix& return_matrix
) const
{
    // Get all interesing subsets into the vector trx_vecs
    const std::vector<ConstTraxelRefVector> trx_vecs =
        subsets_extractor_ptr_->operator()(graph);

    return_matrix.reshape(vigra::Shape2(0, 0));
    size_t subset_count = trx_vecs.size();
    if (subset_count >= 1 )
    {
        std::vector<ConstTraxelRefVector>::const_iterator trx_it = trx_vecs.begin();
        std::vector<FeatureScalar> return_vector;
        for (; trx_it != trx_vecs.end(); trx_it++)
        {
            // Calculate the features for this subset and store it in the variable
            // feature_matrix
            FeatureMatrix feature_matrix;
            feature_calculator_ptr_->calculate(
                feature_extractor_ptr_->extract(*trx_it),
                feature_matrix
            );
            // Write the feature_matrix into an vector
            return_vector.insert(
                return_vector.end(),
                feature_matrix.begin(),
                feature_matrix.end()
            );
        }
        return_matrix.reshape(vigra::Shape2(return_vector.size(), 1));
        std::copy(return_vector.begin(), return_vector.end(), return_matrix.begin());
    }
}

////
//// class TraxelsFeaturesIdentity
////
const std::string TraxelsFeaturesIdentity::name_ = "TraxelsFeaturesIdentity";

TraxelsFeaturesIdentity::TraxelsFeaturesIdentity(
    const std::vector<std::string>& feature_names
) : feature_names_(feature_names)
{
}

TraxelsFeaturesIdentity::TraxelsFeaturesIdentity(
    const std::string& feature_name
)
{
    feature_names_.resize(1);
    feature_names_[0] = feature_name;
}

const std::string& TraxelsFeaturesIdentity::name() const
{
    return name_;
}

void TraxelsFeaturesIdentity::extract(
    const ConstTraxelRefVector& traxelrefs,
    FeatureMatrix& feature_matrix
) const
{
    // get the size of the return matrix
    size_t x_size = traxelrefs.size();

    // get y-size
    size_t y_size = 0;

    // get the feature map of the first traxel and iterate over the feature names
    const FeatureMap& feature_map = traxelrefs.front()->features.get();
    for (
        std::vector<std::string>::const_iterator fname_it = feature_names_.begin();
        fname_it != feature_names_.end();
        ++fname_it
    )
    {
        FeatureMap::const_iterator feature_map_it = feature_map.find(*fname_it);
        if (feature_map_it != feature_map.end())
        {
            y_size += feature_map_it->second.size();
        }
        else
        {
            LOG(logDEBUG) << "In TraxelsFeaturesIdentity: Feature \"" << *fname_it << "\" not found";
        }
    }

    // initialize the return matrix
    feature_matrix.reshape(vigra::Shape2(x_size, y_size));

    // iterate over all traxel
    size_t column_index = 0;
    for(
        ConstTraxelRefVector::const_iterator tref_it = traxelrefs.begin();
        tref_it != traxelrefs.end();
        tref_it++, column_index++
    )
    {
        // fetch the features which are stored in a map
        assert(*tref_it != nullptr);
        const FeatureMap& feature_map = (*tref_it)->features.get();

        // get the current column as a view
        FeatureVectorView column = feature_matrix.bind<0>(column_index);

        // reset the row_index
        size_t row_index = 0;

        // iterate over all feature names
        for(
            std::vector<std::string>::const_iterator fname_it = feature_names_.begin();
            fname_it != feature_names_.end();
            fname_it++
        )
        {
            // check if the feature name exists
            FeatureMap::const_iterator fmap_it = feature_map.find(*fname_it);
            if ( fmap_it != feature_map.end())
            {
                // iterate over the elements of the feature vector
                for (
                    feature_array::const_iterator f_it = fmap_it->second.begin();
                    f_it != fmap_it->second.end();
                    f_it++, row_index++
                )
                {
                    if (row_index >= y_size)
                    {
                        throw std::runtime_error(
                            "In TraxelsFeaturesIdentity: Invalid row index"
                        );
                    }
                    column(row_index) = *f_it;
                }
            } // if (fmap_it != feature_map.end())"
        } // for (fname_it = ... )
    } // for(traxels_it = .. )
}

////
//// class TrackTraxels
////
const std::string TrackTraxels::name_ = "TrackTraxels";

const std::string& TrackTraxels::name() const
{
    return name_;
}

TrackTraxels::TrackTraxels(
    bool require_div_start,
    bool require_div_end
) :
    require_div_start_(require_div_start),
    require_div_end_(require_div_end)
{}

const std::vector<ConstTraxelRefVector>& TrackTraxels::operator()(
    const HypothesesGraph& graph
)
{
    ret_.clear();

    // Check if the graph has the necessary attributes
    if (not graph.has_property(node_active2()))
    {
        throw std::runtime_error(
            "Graph doesn't have a \"node_active2\" property map"
        );
    }
    if (not graph.has_property(arc_active()))
    {
        throw std::runtime_error(
            "Graph doesn't have an \"arc_active\" property map"
        );
    }

    // Get the property maps
    node_active2_map_type& node_active_map = graph.get(node_active2());
    arc_active_map_type& arc_active_map = graph.get(arc_active());

    // check if we have a tracklet graph
    bool with_tracklets = graph.has_property(node_tracklet());

    if(with_tracklets)
    {
        // check if the traxel vector for any tracklet is empty
        property_map<node_active2, HypothesesGraph::base_graph>::type::ValueIt active_valueIt = node_active_map.beginValue();
        ++active_valueIt; // ignore nodes at state 0
        for (; active_valueIt != node_active_map.endValue(); ++active_valueIt)
        {
            property_map<node_active2, HypothesesGraph::base_graph>::type::ItemIt active_itemIt(node_active_map, *active_valueIt);
            for (; active_itemIt != lemon::INVALID; ++active_itemIt)
            {
                const std::vector<Traxel>& t_vec = graph.get(node_tracklet())[active_itemIt];
                if (t_vec.size() == 0)
                {
                    with_tracklets = false;
                    LOG(logDEBUG3) << "In TrackTraxels::operator(): "
                                  << "Empty traxel vector in tracklet map for node " << graph.id(active_itemIt);
                    LOG(logDEBUG3) << "Use therefore traxel in traxel map";
                }
            }
        }
    }

    // check if the graph is legal
    if (not (graph.has_property(node_traxel()) or with_tracklets))
    {
        throw std::runtime_error(
            "HypothesesGraph has neither traxel nor complete tracklet property map"
        );
    }

    // Make maps from child to parent and parent to child
    typedef std::map<HypothesesGraph::Node, HypothesesGraph::Node> NodeNodeMap;
    NodeNodeMap parent;
    NodeNodeMap child;
    // Store for each node if it has a natural parent (simply checks if there
    // is an active arc that points to this node -> only correct if our solution
    // is valid)
    std::map<HypothesesGraph::Node, bool> has_natural_parent;
    // Store for each node if is has a natural child as well
    std::map<HypothesesGraph::Node, bool> has_natural_child;

    // Initialize
    property_map<node_active2, HypothesesGraph::base_graph>::type::ValueIt active_valueIt = node_active_map.beginValue();
    ++active_valueIt; // ignore nodes at state 0
    for (; active_valueIt != node_active_map.endValue(); ++active_valueIt)
    {
        property_map<node_active2, HypothesesGraph::base_graph>::type::ItemIt n_it(node_active_map, *active_valueIt);
        for (; n_it != lemon::INVALID; ++n_it)
        {
            assert(with_tracklets || graph.get(node_traxel())[n_it].Id > 0);
            parent[n_it] = n_it;
            child[n_it] = n_it;
            has_natural_parent[n_it] = false;
            has_natural_child[n_it] = false;
        }
    }

    // Set connections
    for (ArcActiveIt a_it(arc_active_map); a_it != lemon::INVALID; ++a_it)
    {
        assert(with_tracklets || graph.get(node_traxel())[graph.target(a_it)].Id > 0);
        assert(with_tracklets || graph.get(node_traxel())[graph.source(a_it)].Id > 0);
        assert(node_active_map[graph.source(a_it)] > 0);
        assert(node_active_map[graph.target(a_it)] > 0);
        has_natural_parent[graph.target(a_it)] = true;
        has_natural_child[graph.source(a_it)] = true;
        // count the active arcs with the same source
        size_t active_out_arcs = 0;
        for (
            OutArcIt o_it(graph, graph.source(a_it));
            o_it != lemon::INVALID;
            ++o_it
        )
        {
            active_out_arcs += (arc_active_map[o_it] ? 1 : 0);
        }
        // link those nodes if there are no other active arcs with the same source
        if (active_out_arcs == 1)
        {
            parent[graph.target(a_it)] = graph.source(a_it);
            child[graph.source(a_it)] = graph.target(a_it);
        }
    }

    // Compose return vector of traxel reference vectors
    for (
        NodeNodeMap::const_iterator nmap_it = parent.begin();
        nmap_it != parent.end();
        ++nmap_it
    )
    {
        HypothesesGraph::Node current_node = nmap_it->first;
        LOG(logDEBUG4) << "Is parent node invalid?";
        LOG(logDEBUG4) << (parent[current_node] == lemon::INVALID);
        assert(with_tracklets || graph.get(node_traxel())[current_node].Id > 0);
        // check for node if we have to start a new track
        bool start_track = (parent[current_node] == current_node);
        if (start_track and require_div_start_)
        {
            start_track = has_natural_parent[current_node];
        }
        if (start_track)
        {
            // resize the return vector
            ret_.resize(ret_.size() + 1);
            bool loop = true;
            // loop as long as the track isn't finished
            while (loop)
            {
                assert(with_tracklets || graph.get(node_traxel())[current_node].Id > 0);
                if (with_tracklets)
                {
                    // get the traxel vector of this node
                    const std::vector<Traxel>& t_vec = graph.get(node_tracklet())[current_node];
                    for (
                        std::vector<Traxel>::const_iterator t_it = t_vec.begin();
                        t_it != t_vec.end();
                        t_it++
                    )
                    {
                        ret_.back().push_back( &(*t_it) );
                    }
                }
                else
                {
                    assert(with_tracklets || graph.get(node_traxel())[current_node].Id > 0);
                    ret_.back().push_back( &(graph.get(node_traxel())[current_node]) );
                }
                if(child.find(current_node) == child.end())
                    loop = false;
                else
                    loop = current_node != child[current_node];
                // remove the last vector if a division is required at the end of a
                // track and this condition isn't fulfilled
                if ((not loop) and (require_div_end_))
                {
                    if (not has_natural_child[current_node])
                    {
                        ret_.erase(ret_.end() - 1);
                    }
                }
                if(child.find(current_node) != child.end())
                    current_node = child[current_node];
            }
        }
    }
    return ret_;
}

////
//// class DivisionTraxels
////
const std::string DivisionTraxels::name_ = "DivisionTraxels";

const std::string& DivisionTraxels::name() const
{
    return name_;
}

const std::vector<ConstTraxelRefVector>& DivisionTraxels::operator()(
    const HypothesesGraph& graph
)
{
    return operator()(graph, depth_);
}
const std::vector<ConstTraxelRefVector>& DivisionTraxels::operator()(
    const HypothesesGraph& graph,
    size_t depth
)
{
    // Check if the graph has the necessary attributes
    if (not graph.has_property(node_active()))
    {
        throw std::runtime_error(
            "Graph doesn't have a \"node_active\" property map"
        );
    }
    if (not graph.has_property(arc_active()))
    {
        throw std::runtime_error(
            "Graph doesn't have an \"arc_active\" property map"
        );
    }

    // Get the node_active property map
    node_active2_map_type& node_active_map = graph.get(node_active2());

    // check if we have a tracklet graph
    bool with_tracklets = graph.has_property(node_tracklet());

    if(with_tracklets)
    {
        // check if the traxel vector for any tracklet is empty
        property_map<node_active2, HypothesesGraph::base_graph>::type::ValueIt active_valueIt = node_active_map.beginValue();
        ++active_valueIt; // ignore nodes at state 0
        for (; active_valueIt != node_active_map.endValue(); ++active_valueIt)
        {
            property_map<node_active2, HypothesesGraph::base_graph>::type::ItemIt active_itemIt(node_active_map, *active_valueIt);
            for (; active_itemIt != lemon::INVALID; ++active_itemIt)
            {
                const std::vector<Traxel>& t_vec = graph.get(node_tracklet())[active_itemIt];
                if (t_vec.size() == 0)
                {
                    with_tracklets = false;
                    LOG(logDEBUG3) << "In TrackTraxels::operator(): "
                                  << "Empty traxel vector in tracklet map for node " << graph.id(active_itemIt);
                    LOG(logDEBUG3) << "Use therefore traxel in traxel map";
                }
            }
        }
    }

    // check if the graph is legal
    if (not (graph.has_property(node_traxel()) or with_tracklets))
    {
        throw std::runtime_error(
            "HypothesesGraph has neither traxel nor complete tracklet property map"
        );
    }

    // make shure the depth is a valid value
    if (depth <= 0 )
    {
        depth = 1;
    }

    // call calculation functions regarding what kind of graph we have
    if (with_tracklets)
    {
        return from_tracklet_graph(graph, depth);
    }
    else
    {
        return from_traxel_graph(graph, depth);
    }
}

const std::vector<ConstTraxelRefVector>& DivisionTraxels::from_traxel_graph(
    const HypothesesGraph& graph,
    size_t depth
)
{
    ret_.clear();
    // check if we have a graph with mergers
    bool with_mergers = graph.has_property(node_active2());
    // Get the property maps
    node_active_map_type& node_active_map = graph.get(node_active());
    node_traxel_map_type& node_traxel_map = graph.get(node_traxel());

    // Find the divisions
    for (NodeActiveIt n_it(node_active_map); n_it != lemon::INVALID; ++n_it)
    {
        bool is_single_object_node = true;
        if (with_mergers)
        {
            node_active2_map_type& node_active2_map = graph.get(node_active2());
            is_single_object_node = (node_active2_map[n_it] == 1);
        }
        if (is_single_object_node)
        {
            // Count the active outgoing arcs
            std::vector<HypothesesGraph::Node> out_nodes;
            get_out_nodes(n_it, graph, out_nodes);
            ConstTraxelRefVector parents;
            ConstTraxelRefVector l_children;
            ConstTraxelRefVector r_children;
            // Two outgoing arcs: division
            if (out_nodes.size() == 2)
            {
                // Initialize the variables that change during the while loop
                // Now follows some ugly code
                bool valid_division = true;
                size_t curr_depth = depth;
                // These variables are the current parent node and the child nodes
                HypothesesGraph::Node parent = n_it;
                HypothesesGraph::Node l_node = out_nodes[0];
                HypothesesGraph::Node r_node = out_nodes[1];
                // Reserve space for the vectors in which the incoming and outgoing
                // nodes are written
                Nodevector in;
                Nodevector l_out;
                Nodevector r_out;
                while ((curr_depth != 0) and valid_division)
                {
                    // Save the reference to the traxel corresponding to the nodes
                    parents.push_back( &(node_traxel_map[parent]) );
                    l_children.push_back( &(node_traxel_map[l_node]) );
                    r_children.push_back( &(node_traxel_map[r_node]) );

                    // Get the following nodes
                    get_in_nodes(parent, graph, in);
                    get_out_nodes(l_node, graph, l_out);
                    get_out_nodes(r_node, graph, r_out);

                    // Check if the track is long enough to return the division to the
                    // full depth
                    valid_division  = (in.size() == 1);
                    valid_division &= (l_out.size() == 1);
                    valid_division &= (r_out.size() == 1);

                    // If all sizes fit get the following nodes
                    if (valid_division)
                    {
                        parent = in.front();
                        l_node = l_out.front();
                        r_node = r_out.front();
                    }
                    curr_depth--;
                }

                // store the results if the while loop ran over the whole depth
                if (curr_depth == 0)
                {
                    ret_.resize(ret_.size() + 1);
                    ret_.back().insert(ret_.back().end(), parents.begin(), parents.end());
                    ret_.back().insert(
                        ret_.back().end(),
                        l_children.begin(),
                        l_children.end()
                    );
                    ret_.back().insert(
                        ret_.back().end(),
                        r_children.begin(),
                        r_children.end()
                    );
                } // if (curr_depth == 0)
            } // if (out_nodes.size() == 2)
        } // if (is_single_object_node)
    } // for (NodeActiveIt n_it(...) ... )
    return ret_;
}

const std::vector<ConstTraxelRefVector>& DivisionTraxels::from_tracklet_graph(
    const HypothesesGraph& graph,
    size_t depth
)
{
    ret_.clear();
    // Get the property maps
    node_active_map_type& node_active_map = graph.get(node_active());
    // check if we have a graph with mergers
    bool with_mergers = graph.has_property(node_active2());

    // Find the divisions
    for (NodeActiveIt n_it(node_active_map); n_it != lemon::INVALID; ++n_it)
    {
        bool is_single_object_node = true;
        if (with_mergers)
        {
            node_active2_map_type& node_active2_map = graph.get(node_active2());
            is_single_object_node = (node_active2_map[n_it] == 1);
        }
        if (is_single_object_node)
        {
            // Count the active outgoing arcs
            std::vector<HypothesesGraph::Node> out_nodes;
            get_out_nodes(n_it, graph, out_nodes);
            // Two outgoing arcs: division
            if (out_nodes.size() == 2)
            {
                ConstTraxelRefVector parents;
                ConstTraxelRefVector l_children;
                ConstTraxelRefVector r_children;
                bool valid_division = true;

                valid_division &= get_parents_to_depth(
                                      n_it,
                                      graph,
                                      depth,
                                      parents
                                  );
                valid_division &= get_children_to_depth(
                                      out_nodes[0],
                                      graph,
                                      depth,
                                      l_children
                                  );
                valid_division &= get_children_to_depth(
                                      out_nodes[1],
                                      graph,
                                      depth,
                                      r_children
                                  );

                if (valid_division)
                {
                    ret_.resize(ret_.size() + 1);
                    ret_.back().insert(
                        ret_.back().end(),
                        parents.begin(),
                        parents.end()
                    );
                    ret_.back().insert(
                        ret_.back().end(),
                        l_children.begin(),
                        l_children.end()
                    );
                    ret_.back().insert(
                        ret_.back().end(),
                        r_children.begin(),
                        r_children.end()
                    );
                }
            } // if (out_nodes.size() == 2)
        } // if (is_single_object_node)
    }
    return ret_;
}

bool DivisionTraxels::get_children_to_depth(
    const HypothesesGraph::Node& node,
    const HypothesesGraph& graph,
    size_t depth,
    ConstTraxelRefVector& traxelrefs
)
{
    HypothesesGraph::Node curr_node = node;
    size_t curr_tracklet_index = 0;
    bool valid_division = true;
    node_tracklet_map_type& tracklet_map = graph.get(node_tracklet());

    while (valid_division and (depth != 0))
    {
        traxelrefs.push_back( &(tracklet_map[curr_node][curr_tracklet_index]) );

        curr_tracklet_index++;
        depth--;
        if (curr_tracklet_index == tracklet_map[curr_node].size())
        {
            Nodevector out_nodes;
            get_out_nodes(curr_node, graph, out_nodes);
            valid_division = (out_nodes.size() == 1);
            if (valid_division)
            {
                curr_node = out_nodes[0];
                curr_tracklet_index = 0;
            }
        }
    }
    return (depth == 0);
}

bool DivisionTraxels::get_parents_to_depth(
    const HypothesesGraph::Node& node,
    const HypothesesGraph& graph,
    size_t depth,
    ConstTraxelRefVector& traxelrefs
)
{
    node_tracklet_map_type& tracklet_map = graph.get(node_tracklet());

    HypothesesGraph::Node curr_node = node;
    size_t curr_tracklet_index = tracklet_map[node].size() - 1;
    bool valid_division = true;

    while (valid_division and (depth != 0))
    {
        traxelrefs.push_back( &(tracklet_map[curr_node][curr_tracklet_index]) );

        if (curr_tracklet_index == 0)
        {
            Nodevector in_nodes;
            get_in_nodes(curr_node, graph, in_nodes);
            valid_division = (in_nodes.size() == 1);
            if (valid_division)
            {
                curr_node = in_nodes[0];
                curr_tracklet_index = tracklet_map[in_nodes[0]].size();
            }
        }
        curr_tracklet_index--;
        depth--;
    }
    return (depth == 0);
}

////
//// class AppearanceTraxels
////
const std::string AppearanceTraxels::name_ = "AppearanceTraxels";

const std::string& AppearanceTraxels::name() const
{
    return name_;
}

AppearanceTraxels::AppearanceTraxels(
    const AppearanceType appearance,
    FilterFunctionType traxel_filter_function
) :
    appearance_(appearance),
    filter_function_(traxel_filter_function)
{}

const std::vector<ConstTraxelRefVector>& AppearanceTraxels::operator()(
    const HypothesesGraph& graph
)
{
    ret_.clear();

    // Check if the graph has the necessary attributes
    if (not graph.has_property(node_active2()))
    {
        throw std::runtime_error(
            "Graph doesn't have a \"node_active2\" property map"
        );
    }
    if (not graph.has_property(arc_active()))
    {
        throw std::runtime_error(
            "Graph doesn't have an \"arc_active\" property map"
        );
    }

    // Get the property maps
    node_active2_map_type& node_active_map = graph.get(node_active2());
    arc_active_map_type& arc_active_map = graph.get(arc_active());

    // check if we have a tracklet graph
    bool with_tracklets = graph.has_property(node_tracklet());

    if(with_tracklets)
    {
        // check if the traxel vector for any tracklet is empty
        property_map<node_active2, HypothesesGraph::base_graph>::type::ValueIt active_valueIt = node_active_map.beginValue();
        ++active_valueIt; // ignore nodes at state 0
        for (; active_valueIt != node_active_map.endValue(); ++active_valueIt)
        {
            property_map<node_active2, HypothesesGraph::base_graph>::type::ItemIt active_itemIt(node_active_map, *active_valueIt);
            for (; active_itemIt != lemon::INVALID; ++active_itemIt)
            {
                const std::vector<Traxel>& t_vec = graph.get(node_tracklet())[active_itemIt];
                if (t_vec.size() == 0)
                {
                    with_tracklets = false;
                    LOG(logDEBUG3) << "In TrackTraxels::operator(): "
                                  << "Empty traxel vector in tracklet map for node " << graph.id(active_itemIt);
                    LOG(logDEBUG3) << "Use therefore traxel in traxel map";
                }
            }
        }
    }

    // check if the graph is legal
    if (not (graph.has_property(node_traxel()) or with_tracklets))
    {
        throw std::runtime_error(
            "HypothesesGraph has neither traxel nor complete tracklet property map"
        );
    }

    property_map<node_active2, HypothesesGraph::base_graph>::type::ValueIt active_valueIt = node_active_map.beginValue();
    ++active_valueIt; // ignore nodes at state 0
    for (; active_valueIt != node_active_map.endValue(); ++active_valueIt)
    {
        property_map<node_active2, HypothesesGraph::base_graph>::type::ItemIt n_it(node_active_map, *active_valueIt);
        for (; n_it != lemon::INVALID; ++n_it)
        {
            if (appearance_ == AppearanceType::Appearance)
            {
                size_t in_arcs = 0;
                for (InArcIt a_it(graph, n_it); a_it != lemon::INVALID; ++a_it)
                {
                    in_arcs += (arc_active_map[a_it] ? 1 : 0);
                }
                if (in_arcs == 0)
                {
                    const Traxel* traxel_ref;
                    if (with_tracklets)
                    {
                        traxel_ref = &(graph.get(node_tracklet())[n_it].front());
                    }
                    else
                    {
                        traxel_ref = &(graph.get(node_traxel())[n_it]);
                    }
                    if (filter_function_)
                    {
                        if (filter_function_(*traxel_ref))
                        {
                            ret_.resize(ret_.size() + 1);
                            ret_.back().push_back(traxel_ref);
                        }
                    }
                    else
                    {
                        ret_.resize(ret_.size() + 1);
                        ret_.back().push_back(traxel_ref);
                    }
                }
            }
            else if (appearance_ == AppearanceType::Disappearance)
            {
                size_t out_arcs = 0;
                for (OutArcIt a_it(graph, n_it); a_it != lemon::INVALID; ++a_it)
                {
                    out_arcs += (arc_active_map[a_it] ? 1 : 0);
                }
                if (out_arcs == 0)
                {
                    const Traxel* traxel_ref;
                    if (with_tracklets)
                    {
                        traxel_ref = &(graph.get(node_tracklet())[n_it].front());
                    }
                    else
                    {
                        traxel_ref = &(graph.get(node_traxel())[n_it]);
                    }
                    if (filter_function_)
                    {
                        if (filter_function_(*traxel_ref))
                        {
                            ret_.resize(ret_.size() + 1);
                            ret_.back().push_back(traxel_ref);
                        }
                    }
                    else
                    {
                        ret_.resize(ret_.size() + 1);
                        ret_.back().push_back(traxel_ref);
                    }
                }
            }
        }
    }
    return ret_;
}

////
//// class CompositionCalculator
////
const std::string CompositionCalculator::name_ = "CompositionCalculator";

const std::string& CompositionCalculator::name() const
{
    return name_;
}

void CompositionCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    FeatureMatrix temp;
    first_calculator_ptr_->calculate(feature_matrix, temp);
    second_calculator_ptr_->calculate(temp, return_matrix);
}

////
//// class TraxelsFCFromFC
////
const std::string TraxelsFCFromFC::name_ = "TraxelsFCFromFC";

const std::string& TraxelsFCFromFC::name() const
{
    return name_;
}

void TraxelsFCFromFC::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    if (col_count < order_)
    {
        LOG(logDEBUG) << "In TraxelsFCFromFC: not enough data to calculate feature";
        LOG(logDEBUG) << "Returning zero";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0.0);
    }
    else
    {
        // convert feature_matrix to vec<vec<FeatureScalar> >
        std::vector<std::vector<FeatureScalar> > feature_vectors;
        for (size_t i = 0; i < col_count; i++)
        {
            feature_vectors.push_back(std::vector<FeatureScalar>());
            feature_vectors.back().insert(
                feature_vectors.back().end(),
                feature_matrix.bind<0>(i).begin(),
                feature_matrix.bind<0>(i).end()
            );
        }
        if (order_ == 1)
        {
            std::vector<FeatureScalar> ret_ = feature_calculator_ptr_->calculate(
                                                  feature_vectors[0]
                                              );
            size_t ret_row_count = ret_.size();
            return_matrix.reshape(vigra::Shape2(col_count, ret_row_count));
            for (size_t j = 0; j < ret_row_count; j++)
            {
                return_matrix(0, j) = ret_[j];
            }
            for (size_t i = 1; i < col_count; i++)
            {
                ret_ = feature_calculator_ptr_->calculate(feature_vectors[i]);
                for (size_t j = 0; j < ret_row_count; j++)
                {
                    return_matrix(i, j) = ret_[j];
                }
            }
        }
        else if (order_ == 2)
        {
            std::vector<FeatureScalar> ret_ = feature_calculator_ptr_->calculate(
                                                  feature_vectors[0],
                                                  feature_vectors[1]
                                              );
            size_t ret_row_count = ret_.size();
            return_matrix.reshape(vigra::Shape2(col_count - 1, ret_row_count));
            for (size_t j = 0; j < ret_row_count; j++)
            {
                return_matrix(0, j) = ret_[j];
            }
            for (size_t i = 1; i + 1 < col_count; i++)
            {
                ret_ = feature_calculator_ptr_->calculate(
                           feature_vectors[i],
                           feature_vectors[i + 1]
                       );
                for (size_t j = 0; j < ret_row_count; j++)
                {
                    return_matrix(i, j) = ret_[j];
                }
            }
        }
        else if (order_ == 3)
        {
            std::vector<FeatureScalar> ret_ = feature_calculator_ptr_->calculate(
                                                  feature_vectors[0],
                                                  feature_vectors[1],
                                                  feature_vectors[2]
                                              );
            size_t ret_row_count = ret_.size();
            return_matrix.reshape(vigra::Shape2(col_count - 2, ret_row_count));
            for (size_t j = 0; j < ret_row_count; j++)
            {
                return_matrix(0, j) = ret_[j];
            }
            for (size_t i = 1; i + 2 < col_count; i++)
            {
                ret_ = feature_calculator_ptr_->calculate(
                           feature_vectors[i],
                           feature_vectors[i + 1],
                           feature_vectors[i + 2]
                       );
                for (size_t j = 0; j < ret_row_count; j++)
                {
                    return_matrix(i, j) = ret_[j];
                }
            }
        }
    }
}

////
//// class SumCalculator
////
template<int N>
const std::string SumCalculator<N>::name_ = "SumCalculator";

template<int N>
const std::string& SumCalculator<N>::name() const
{
    return name_;
}

template<>
void SumCalculator<0>::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In SumCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning zero";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0.0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(1, row_count));
        return_matrix.init(0.0);
        for (size_t col = 0; col < col_count; col++)
        {
            return_matrix.bind<0>(0) += feature_matrix.bind<0>(col);
        }
    }
}

template<>
void SumCalculator<1>::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In SumCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning zero";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0.0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(col_count, 1));
        return_matrix.init(0.0);
        for (size_t row = 0; row < row_count; row++)
        {
            return_matrix.bind<1>(0) += feature_matrix.bind<1>(row);
        }
    }
}

template<>
void SumCalculator < -1 >::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    return_matrix.reshape(vigra::Shape2(1, 1));
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In SumCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning zero";
        return_matrix.init(0.0);
    }
    else
    {
        return_matrix(0, 0) = std::accumulate(
                                  feature_matrix.begin(),
                                  feature_matrix.end(),
                                  0.0
                              );
    }
}

// Explicit instantiation
template class SumCalculator<0>;
template class SumCalculator<1>;
template class SumCalculator < -1 >;

////
//// class DiffCalculator
////
const std::string DiffCalculator::name_ = "DiffCalculator";

const std::string& DiffCalculator::name() const
{
    return name_;
}

void DiffCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{

    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if (col_count <= 1)
    {
        LOG(logDEBUG) << "In DiffCalculator: matrix in argument has less than one column";
        LOG(logDEBUG) << "Returning a 0-vector";
        return_matrix.reshape(vigra::Shape2(1, row_count));
        return_matrix.init(0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(col_count - 1, row_count));
        for (size_t col = 0; col < col_count - 1; col++)
        {
            FeatureVectorView a = feature_matrix.bind<0>(col);
            FeatureVectorView b = feature_matrix.bind<0>(col + 1);
            using namespace vigra::multi_math;
            return_matrix.bind<0>(col) = b - a;
        }
    }
}

////
//// class CurveCalculator
////
const std::string CurveCalculator::name_ = "CurveCalculator";

const std::string& CurveCalculator::name() const
{
    return name_;
}

void CurveCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if (col_count <= 2)
    {
        LOG(logDEBUG) << "In CurveCalculator: matrix in argument has less than three columns";
        LOG(logDEBUG) << "Returning a 0-vector";
        return_matrix.reshape(vigra::Shape2(1, row_count));
        return_matrix.init(0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(col_count - 2, row_count));
        for (size_t col = 0; col < col_count - 2; col++)
        {
            using namespace vigra::multi_math;
            FeatureVectorView a = feature_matrix.bind<0>(col);
            FeatureVectorView b = feature_matrix.bind<0>(col + 1);
            FeatureVectorView c = feature_matrix.bind<0>(col + 2);
            return_matrix.bind<0>(col) = a - 2 * b + c;
        }
    }
}

////
//// class MinCalculator
////
template<int N>
const std::string MinCalculator<N>::name_ = "MinCalculator";

template<int N>
const std::string& MinCalculator<N>::name() const
{
    return name_;
}

template<>
void MinCalculator<0>::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In MinCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning a 0-vector";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(1, row_count));
        return_matrix.init(0);
        for (size_t i = 0; i < row_count; i++)
        {
            const FeatureVectorView current_row = feature_matrix.bind<1>(i);
            return_matrix(0, i) = *(std::min_element(current_row.begin(), current_row.end()));
        }
    }
}

template<>
void MinCalculator<1>::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In MinCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning a 0-vector";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(col_count, 1));
        return_matrix.init(0);
        for (size_t i = 0; i < col_count; i++)
        {
            const FeatureVectorView current_col = feature_matrix.bind<0>(i);
            return_matrix(i, 0) = *(
                                      std::min_element(current_col.begin(), current_col.end())
                                  );
        }
    }
}

template<>
void MinCalculator < -1 >::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    return_matrix.reshape(vigra::Shape2(1, 1));
    return_matrix.init(0);
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In MinCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning a 0-vector";
    }
    else
    {
        return_matrix(0, 0) = *(
                                  std::min_element(feature_matrix.begin(), feature_matrix.end())
                              );
    }
}

// Explicit instantiation
template class MinCalculator<0>;
template class MinCalculator<1>;
template class MinCalculator < -1 >;

////
//// class MaxCalculator
////
template<int N>
const std::string MaxCalculator<N>::name_ = "MaxCalculator";

template<int N>
const std::string& MaxCalculator<N>::name() const
{
    return name_;
}
template <>
void MaxCalculator<0>::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In MaxCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning a 0-vector";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(1, row_count));
        return_matrix.init(0);
        for (size_t i = 0; i < row_count; i++)
        {
            const FeatureVectorView current_row = feature_matrix.bind<1>(i);
            return_matrix(0, i) = *(
                                      std::max_element(current_row.begin(), current_row.end())
                                  );
        }
    }
}

template <>
void MaxCalculator<1>::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In MaxCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning a 0-vector";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(col_count, 1));
        return_matrix.init(0);
        for (size_t i = 0; i < col_count; i++)
        {
            const FeatureVectorView current_col = feature_matrix.bind<0>(i);
            return_matrix(i, 0) = *(
                                      std::max_element(current_col.begin(), current_col.end())
                                  );
        }
    }
}

template <>
void MaxCalculator < -1 >::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    return_matrix.reshape(vigra::Shape2(1, 1));
    return_matrix.init(0.0);
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In MaxCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning a 0-vector";
    }
    else
    {
        return_matrix(0, 0) = *(
                                  std::max_element(feature_matrix.begin(), feature_matrix.end())
                              );
    }
}

// Explicit instantiation
template class MaxCalculator<0>;
template class MaxCalculator<1>;
template class MaxCalculator < -1 >;

////
//// class MeanCalculator
////
template<int N>
const std::string MeanCalculator<N>::name_ = "MeanCalculator";

template<int N>
const std::string& MeanCalculator<N>::name() const
{
    return name_;
}

template<>
void MeanCalculator<0>::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In MeanCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning a 0-vector";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(1, row_count));
        for (size_t i = 0; i < row_count; i++)
        {
            return_matrix(0, i) = vigra::multi_math::sum<FeatureScalar>(
                                      feature_matrix.bind<1>(i)
                                  ) / static_cast<FeatureScalar>(col_count);
        }
    }
}

template<>
void MeanCalculator<1>::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In MeanCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning a 0-vector";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(col_count, 1));
        for (size_t i = 0; i < col_count; i++)
        {
            return_matrix(i, 0) = vigra::multi_math::sum<FeatureScalar>(
                                      feature_matrix.bind<0>(i)
                                  ) / static_cast<FeatureScalar>(row_count);
        }
    }
}

template<>
void MeanCalculator < -1 >::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    return_matrix.reshape(vigra::Shape2(1, 1));
    return_matrix.init(0);
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In MeanCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning a 0-vector";
    }
    else
    {
        return_matrix(0, 0) = std::accumulate(
                                  feature_matrix.begin(),
                                  feature_matrix.end(),
                                  0.0
                              ) / static_cast<FeatureScalar>(feature_matrix.size());
    }
}

// Explicit instantiation
template class MeanCalculator<0>;
template class MeanCalculator<1>;
template class MeanCalculator < -1 >;

////
//// class SquareCalculator
////
const std::string SquareCalculator::name_ = "SquareCalculator";

const std::string& SquareCalculator::name() const
{
    return name_;
}

void SquareCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    return_matrix.reshape(feature_matrix.shape());
    FeatureMatrix::const_iterator f_it = feature_matrix.begin();
    FeatureMatrix::iterator r_it = return_matrix.begin();
    for (; f_it != feature_matrix.end(); f_it++, r_it++ )
    {
        *r_it = *f_it **f_it;
    }
}

////
//// class SquareRootCalculator
////
const std::string SquareRootCalculator::name_ = "SquareRootCalculator";

const std::string& SquareRootCalculator::name() const
{
    return name_;
}

void SquareRootCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    return_matrix.reshape(feature_matrix.shape());
    FeatureMatrix::const_iterator f_it = feature_matrix.begin();
    FeatureMatrix::iterator r_it = return_matrix.begin();
    for (; f_it != feature_matrix.end(); f_it++, r_it++ )
    {
        *r_it = sqrt(*f_it);
    }
}

////
//// class SquaredNormCalculator
////
template<int N>
const std::string SquaredNormCalculator<N>::name_ = "SquaredNormCalculator";

template<int N>
const std::string& SquaredNormCalculator<N>::name() const
{
    return name_;
}

template<>
void SquaredNormCalculator<0>::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    return_matrix.reshape(vigra::Shape2(col_count, 1));
    for (size_t i = 0; i < col_count; i++)
    {
        return_matrix(i, 0) = vigra::squaredNorm(feature_matrix.bind<0>(i));
    }
}

template<>
void SquaredNormCalculator<1>::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t row_count = feature_matrix.shape(1);
    return_matrix.reshape(vigra::Shape2(1, row_count));
    for (size_t i = 0; i < row_count; i++)
    {
        return_matrix(0, i) = vigra::squaredNorm(feature_matrix.bind<1>(i));
    }
}

// Explicit instantiation
template class SquaredNormCalculator<0>;
template class SquaredNormCalculator<1>;

////
//// class DotProductCalculator
////
const std::string DotProductCalculator::name_ = "DotProductCalculator";

const std::string& DotProductCalculator::name() const
{
    return name_;
}

void DotProductCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    // Calculate the dot product depending on the size of the input matrix
    if (col_count == 0)
    {
        // empty matrix
        LOG(logDEBUG) << "In DotProductCalculator: matrix is empty";
        LOG(logDEBUG) << "Return zero";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix(0, 0) = 0.0;
    }
    else if (col_count == 1)
    {
        // matrix with one column
        LOG(logDEBUG) << "In DotProductCalculator: matrix has only one column";
        LOG(logDEBUG) << "Calculate the norm of this vector";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix(0, 0) =  vigra::linalg::dot(feature_matrix, feature_matrix);
    }
    else if (col_count >= 2)
    {
        // matrix with two or more columns
        return_matrix.reshape(vigra::Shape2(col_count - 1, 1));
        for (size_t i = 0; i < (col_count - 1); i++)
        {
            return_matrix(i, 0) = vigra::linalg::dot(
                                      feature_matrix.bind<0>(i),
                                      feature_matrix.bind<0>(i + 1)
                                  );
        }
    }
}

////
//// class AngleCosineCalculator
////
const std::string AngleCosineCalculator::name_ = "AngleCosineCalculator";

const std::string& AngleCosineCalculator::name() const
{
    return name_;
}

void AngleCosineCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if (col_count <= 2 or row_count == 0)
    {
        LOG(logDEBUG) << "In AngleCosineCalculator: Not enough column vectors";
        LOG(logDEBUG) << "Returning zero";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(col_count - 2, 1));
        FeatureMatrix diff_matrix;
        FeatureMatrix dot_matrix;
        FeatureMatrix norm_matrix;
        diff_calculator_.calculate(feature_matrix, diff_matrix);
        dot_product_calculator_.calculate(diff_matrix, dot_matrix);
        norm_calculator_.calculate(diff_matrix, norm_matrix);
        for (size_t i = 0; i < col_count - 2; i++)
        {
            if (norm_matrix(i, 0) != 0 and norm_matrix(i + 1, 0) != 0)
            {
                return_matrix(i, 0) = dot_matrix(i, 0);
                return_matrix(i, 0) /= (norm_matrix(i, 0) * norm_matrix(i + 1, 0));
            }
            else
            {
                return_matrix(i, 0) = 0;
            }
        }
    }
}

////
//// class DivAngleCosineCalculator
////
const std::string DivAngleCosineCalculator::name_ = "DivAngleCosineCalculator";

const std::string& DivAngleCosineCalculator::name() const
{
    return name_;
}

void DivAngleCosineCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if (col_count != 3 or row_count == 0)
    {
        LOG(logDEBUG) << "In DivAngleCosineCalculator: Invalid division depth";
        LOG(logDEBUG) << "Returning zero";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(1, 1));
        FeatureMatrix diff_matrix;
        FeatureMatrix dot_matrix;
        FeatureMatrix norm_matrix;
        diff_calculator_.calculate(feature_matrix, diff_matrix);
        dot_product_calculator_.calculate(diff_matrix, dot_matrix);
        norm_calculator_.calculate(diff_matrix, norm_matrix);
        if (norm_matrix(0, 0) != 0 and norm_matrix(1, 0) != 0)
        {
            return_matrix(0, 0) = dot_matrix(0, 0);
            return_matrix(0, 0) /= (norm_matrix(0, 0) * norm_matrix(1, 0));
        }
        else
        {
            return_matrix(0, 0) = 0;
        }
    }
}


////
//// class ChildParentDiffCalculator
////
const std::string ChildParentDiffCalculator::name_ = "ChildParentDiffCalculator";

const std::string& ChildParentDiffCalculator::name() const
{
    return name_;
}

void ChildParentDiffCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    if(col_count % 3 != 0)
    {
        calculate(feature_matrix, return_matrix, 0);
    }
    else
    {
        calculate(feature_matrix, return_matrix, col_count / 3);
    }
}

void ChildParentDiffCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix,
    size_t depth
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    return_matrix.reshape(vigra::Shape2(2, row_count));
    if ((col_count < depth * 3) or (depth == 0))
    {
        LOG(logDEBUG) << "In ChildParentDiffCalculator: Invalid division depth";
        LOG(logDEBUG) << "Return two 0-vectors";
        return_matrix.init(0);
    }
    else
    {
        return_matrix.bind<0>(0) = vigra::multi_math::operator-(
                                       feature_matrix.bind<0>(depth),
                                       feature_matrix.bind<0>(0)
                                   );
        return_matrix.bind<0>(1) = vigra::multi_math::operator-(
                                       feature_matrix.bind<0>(2 * depth),
                                       feature_matrix.bind<0>(0)
                                   );
    }
}

////
//// class ChildDeceleration
////
const std::string ChildDeceleration::name_ = "ChildDeceleration";

const std::string& ChildDeceleration::name() const
{
    return name_;
}

void ChildDeceleration::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    if(col_count % 3 != 0)
    {
        calculate(feature_matrix, return_matrix, 0);
    }
    else
    {
        calculate(feature_matrix, return_matrix, col_count / 3);
    }
}

void ChildDeceleration::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix,
    size_t depth
) const
{
    return_matrix.reshape(vigra::Shape2(2, 1));
    return_matrix.init(0.0);
    if (depth < 2)
    {
        LOG(logDEBUG) << "In ChildDeceleration: not enough depth of data";
        LOG(logDEBUG) << "Return 0-vector";
    }
    else
    {
        FeatureVector temp_a1;
        FeatureVector temp_a2;
        FeatureVector temp_b1;
        FeatureVector temp_b2;
        temp_a1 = vigra::multi_math::operator-(
                      feature_matrix.bind<0>(depth),
                      feature_matrix.bind<0>(0)
                  );
        temp_a2 = vigra::multi_math::operator-(
                      feature_matrix.bind<0>(2 * depth),
                      feature_matrix.bind<0>(0)
                  );
        temp_b1 = vigra::multi_math::operator-(
                      feature_matrix.bind<0>(depth + 1),
                      feature_matrix.bind<0>(depth)
                  );
        temp_b2 = vigra::multi_math::operator-(
                      feature_matrix.bind<0>(2 * depth + 1),
                      feature_matrix.bind<0>(2 * depth)
                  );
        using namespace vigra::linalg;
        //return_matrix(0, 0) = dot(temp_b1, temp_b1) / dot(temp_a1, temp_a1);
        //return_matrix(1, 0) = dot(temp_b2, temp_b2) / dot(temp_a2, temp_a2);
	FeatureScalar denominator = dot(temp_a1, temp_a1);
	if(denominator != 0.0)
		return_matrix(0, 0) = dot(temp_b1, temp_b1) / denominator;
	else
		return_matrix(0, 0) = 0.0;

	denominator = dot(temp_a2, temp_a2);
	if(denominator != 0.0)
		return_matrix(1, 0) = dot(temp_b2, temp_b2) / denominator;
	else
		return_matrix(1, 0) = 0.0;
    }
}

////
//// class DeviationCalculator
////
const std::string DeviationCalculator::name_ = "DeviationCalculator";

const std::string& DeviationCalculator::name() const
{
    return name_;
}

void DeviationCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if ((col_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In DeviationCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning a 1x1 matrix filled with zeros";
        return_matrix.reshape(vigra::Shape2(0, 0));
        return_matrix.init(0.0);
    }
    else
    {
        FeatureMatrix mean_vector;
        mean_calculator_.calculate(feature_matrix, mean_vector);
        return_matrix.reshape(vigra::Shape2(col_count, row_count));
        for (size_t i = 0; i < col_count; i++)
        {
            return_matrix.bind<0>(i) = vigra::multi_math::operator-(
                                           feature_matrix.bind<0>(i),
                                           mean_vector.bind<0>(0)
                                       );
        }
    }
}

////
//// class CovarianceCalculator
////
template<bool INV>
const std::string CovarianceCalculator<INV>::name_ = "CovarianceCalculator";

template<bool INV>
const std::string& CovarianceCalculator<INV>::name() const
{
    return name_;
}

template<>
void CovarianceCalculator<true>::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    return_matrix.reshape(vigra::Shape2(row_count, row_count));
    return_matrix.init(0.0);
    if (col_count <= row_count)
    {
        LOG(logDEBUG) << "In CovarianceCalculator: too few data to calculate covariance matrix";
        LOG(logDEBUG) << "Returning a Matrix filled with zeros";
    }
    FeatureMatrix covariance_matrix(vigra::Shape2(row_count, row_count), 0.0);
    vigra::linalg::covarianceMatrixOfColumns(feature_matrix, covariance_matrix);
    bool invertible = vigra::linalg::inverse(covariance_matrix, return_matrix);
    if (not invertible)
    {
        return_matrix.init(0.0);
    }
}

template<>
void CovarianceCalculator<false>::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    return_matrix.reshape(vigra::Shape2(row_count, row_count));
    return_matrix.init(0.0);
    if (col_count <= row_count)
    {
        LOG(logDEBUG) << "In CovarianceCalculator: too few data to calculate covariance matrix";
        LOG(logDEBUG) << "Returning a Matrix filled with zeros";
    }
    vigra::linalg::covarianceMatrixOfColumns(feature_matrix, return_matrix);
}

template class CovarianceCalculator<true>;
template class CovarianceCalculator<false>;

////
//// class SquaredMahalanobisCalculator
////
const std::string SquaredMahalanobisCalculator::name_ = "SquaredMahalanobisCalculator";

const std::string& SquaredMahalanobisCalculator::name() const
{
    return name_;
}

void SquaredMahalanobisCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    FeatureMatrix deviation_matrix;
    deviation_calculator_.calculate(feature_matrix, deviation_matrix);
    FeatureMatrix inv_cov_matrix;
    inv_covariance_calculator_.calculate(feature_matrix, inv_cov_matrix);
    calculate(deviation_matrix, return_matrix, inv_cov_matrix);
}

void SquaredMahalanobisCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix,
    const FeatureMatrix& inv_covariance_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    if ((row_count == 0) or (row_count == 0))
    {
        LOG(logDEBUG) << "In SquaredMahalanobisCalculator: empty input matrix";
        LOG(logDEBUG) << "Returning a 1x1 matrix filled with zeros";
        return_matrix.reshape(vigra::Shape2(1, 1));
        return_matrix.init(0.0);
    }
    else
    {
        return_matrix.reshape(vigra::Shape2(col_count, 1));
        FeatureMatrixView temp1;
        FeatureMatrix temp2(vigra::Shape2(row_count, 1));
        FeatureMatrix temp3(vigra::Shape2(1, 1));
        for (size_t i = 0; i < col_count; i++)
        {
            temp1 = feature_matrix.subarray(
                        vigra::Shape2(i, 0),
                        vigra::Shape2(i + 1, row_count)
                    );
            vigra::linalg::mmul(inv_covariance_matrix, temp1.transpose(), temp2);
            vigra::linalg::mmul(temp1, temp2, temp3);
            return_matrix(i, 0) = temp3(0, 0);
        }
    }
}

////
//// class MVNOutlierCalculator
////
const std::string MVNOutlierCalculator::name_ = "MVNOutlierCalculator";

const std::string& MVNOutlierCalculator::name() const
{
    return name_;
}

void MVNOutlierCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix,
    const FeatureScalar& sigma_threshold
) const
{
    size_t col_count = feature_matrix.shape(0);
    FeatureMatrix sigma_matrix;
    mahalanobis_calculator_.calculate(feature_matrix, sigma_matrix);
    return_matrix.reshape(vigra::Shape2(1, 1));
    return_matrix(0, 0) = 0.0;
    for (size_t col = 0; col < col_count; col++)
    {
        if(sigma_matrix(col, 0) > sigma_threshold)
        {
            return_matrix(0, 0) += 1.0 / static_cast<FeatureScalar>(col_count);
        }
    }
}

void MVNOutlierCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    calculate(feature_matrix, return_matrix, 3.0);
}

#ifdef WITH_DLIB
////
//// class SVMOutlierCalculator
////
const std::string SVMOutlierCalculator::name_ = "SVMOutlierCalculator";

const std::string& SVMOutlierCalculator::name() const
{
    return name_;
}

bool SVMOutlierCalculator::is_trained() const
{
    return is_trained_;
}

void SVMOutlierCalculator::compute_feature_mean_var(const FeatureMatrix& feature_matrix)
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);

    if(col_count < 2 || row_count == 0)
    {
        throw std::runtime_error("Cannot compute mean and var for empty matrix or with less than two samples!");
    }

    feature_means_.clear();
    feature_vars_.clear();

    for(size_t row_idx = 0; row_idx < row_count; row_idx++)
    {
        // get the mean and variance of each feature and store them
        FeatureVectorView row_view = feature_matrix.bind<1>(row_idx);
        FeatureScalar row_mean, row_var;
        row_view.meanVariance(&row_mean, &row_var);
        feature_means_.push_back(row_mean);
        feature_vars_.push_back(row_var);
    }
}

void SVMOutlierCalculator::normalize_features(FeatureMatrix& feature_matrix) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);

    if(col_count < 2 || row_count == 0)
    {
        throw std::runtime_error("Cannot normalize empty matrix or with less than two samples!");
    }

    if(row_count != feature_means_.size() || row_count != feature_vars_.size())
    {
        throw std::runtime_error("Stored normalization info has different number of features than new data!");
    }

    for(size_t row_idx = 0; row_idx < row_count; row_idx++)
    {
        // normalize each row
        FeatureVectorView row_view = feature_matrix.bind<1>(row_idx);
        FeatureScalar row_mean = feature_means_[row_idx];
        FeatureScalar row_var  = feature_vars_[row_idx];
        row_view -= row_mean;
        row_view /= row_var;
    }
}

void SVMOutlierCalculator::train(
    const FeatureMatrix& feature_matrix,
    const FeatureScalar& kernel_width
)
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    std::vector<SampleType> samples;
    // TODO assert num_samples >> dim ?
    if(col_count == 0 || row_count == 0)
    {
        throw std::runtime_error("Cannot train outlier SVM without any samples!");
    }

    // normalize features using mean and variance of feature_matrix
    compute_feature_mean_var(feature_matrix);
    FeatureMatrix normalized_feature_matrix(feature_matrix);
    normalize_features(normalized_feature_matrix);

//    LOG(logINFO) << "Training outlier SVM from " << col_count << " samples:\n" << normalized_feature_matrix;

    for (size_t col = 0; col < col_count; col++)
    {
        FeatureVectorView column = normalized_feature_matrix.bind<0>(col);
        SampleType sample(row_count);
        std::copy(column.begin(), column.end(), sample.begin());
        samples.push_back(sample);
    }

    KernelType kernel(kernel_width);
    OneClassSVMTrainerType trainer;
    trainer.set_kernel(kernel);
    decision_function_ = trainer.train(samples);
    is_trained_ = true;
}

std::ostream& operator<<(std::ostream& lhs, std::vector<dlib::matrix<FeatureScalar, 0, 1>>& rhs)
{
    for(dlib::matrix<FeatureScalar, 0, 1>& sample : rhs)
    {
        lhs << dlib::trans(sample);
    }
    return lhs;
}

void SVMOutlierCalculator::calculate(
    const FeatureMatrix& feature_matrix,
    FeatureMatrix& return_matrix
) const
{
    size_t col_count = feature_matrix.shape(0);
    size_t row_count = feature_matrix.shape(1);
    std::vector<SampleType> samples;

    // normalize features
    FeatureMatrix normalized_feature_matrix(feature_matrix);
    normalize_features(normalized_feature_matrix);

    // TODO assert num_samples >> dim ?
    for (size_t col = 0; col < col_count; col++)
    {
        FeatureVectorView column = normalized_feature_matrix.bind<0>(col);
        SampleType sample(row_count);
        std::copy(column.begin(), column.end(), sample.begin());
        samples.push_back(sample);
    }
    if (is_trained_)
    {
//        LOG(logINFO) << "Predicting from SVM with " << decision_function_.alpha.size() << " weights for samples of size " << row_count;
//        LOG(logINFO) << "Using learned SVM weights: \n" << decision_function_.alpha;
        return_matrix.reshape(vigra::Shape2(col_count, 1));
//        LOG(logINFO) << "and support vectors: \n" << decision_function_.basis_vectors;

        for(size_t col = 0; col < col_count; col++)
        {
            return_matrix(col, 0) = decision_function_(samples[col]);
        }

//        LOG(logINFO) << "SVM predictions from samples\n" << samples << "\n yielded: " << return_matrix.transpose();
    }
    else
    {
        throw std::runtime_error("SVMOutlierCalculator not trained");
    }
}
#endif

} // end namespace features
} // end namespace pgmlink
