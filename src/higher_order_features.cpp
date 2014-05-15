#include "pgmlink/higher_order_features.h"
#include <vigra/multi_math.hxx> /* for operator+ */
#include <vigra/matrix.hxx> /* for covariance calculation */
#include <vigra/linear_algebra.hxx> /* for matrix inverse calculation */

namespace pgmlink {

////
//// Some useful typedefinitions
////
typedef typename
  property_map<node_active, HypothesesGraph::base_graph>::type
  node_active_map_type;
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
void set_solution(HypothesesGraph& graph, const size_t solution_index) {
  // check if the graph has the necessary property maps
  if (not graph.has_property(node_active_count())) {
    throw std::runtime_error(
      "Graph doesn't have a \"node_active_count\" property map"
    );
  }
  if (not graph.has_property(arc_active_count())) {
    throw std::runtime_error(
      "Graph doesn't have an \"arc_active_count\" property map"
    );
  }

  // Get the property maps
  nodes_active_map_type& nodes_active_map = graph.get(node_active_count());
  arcs_active_map_type& arcs_active_map = graph.get(arc_active_count());

  // check if the solution_index is legal
  if (nodes_active_map.beginValue()->size() <= solution_index) {
    throw std::runtime_error("Index of solution out of range");
  }

  // create the a node_active and arc_active map
  if (not graph.has_property(node_active())) {
    graph.add(node_active());
  }
  if (not graph.has_property(arc_active())) {
    graph.add(arc_active());
  }

  // Get the property maps (caution: "node_active_map" not "nodes_active_map")
  node_active_map_type& node_active_map = graph.get(node_active());
  arc_active_map_type& arc_active_map = graph.get(arc_active());

  // Now we can start the with writing the solution with the index
  // solution_index into the node_active_map
  for (NodeIt n_it(graph); n_it != lemon::INVALID; ++n_it) {
    node_active_map[n_it] = nodes_active_map[n_it][solution_index];
  }
  for (ArcIt a_it(graph); a_it != lemon::INVALID; ++a_it) {
    arc_active_map[a_it] = arcs_active_map[a_it][solution_index];
  }
}

////
//// function get_out_nodes
////
void get_out_nodes(
  const HypothesesGraph::Node& node,
  const HypothesesGraph& graph,
  Nodevector& out_nodes
) {
  out_nodes.clear();
  for (OutArcIt oa_it(graph, node); oa_it != lemon::INVALID; ++oa_it) {
    if (graph.get(arc_active())[oa_it]) {
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
) {
  in_nodes.clear();
  for (InArcIt ia_it(graph, node); ia_it != lemon::INVALID; ++ia_it) {
    if (graph.get(arc_active())[ia_it]) {
      in_nodes.push_back(graph.source(ia_it));
    }
  }
}

/*=============================================================================
  pure virtual classes
=============================================================================*/
////
//// class SubsetFeatureExtractor
////
const FeatureMatrix& SubsetFeatureExtractor::extract_matrix(
  const ConstTraxelRefVector& /*traxelrefs*/
) {
  throw std::runtime_error(
    "Feature extractor " + name() + " doesn't provide a matrix valued result"
  );
  return *(new FeatureMatrix);
}
const FeatureVector& SubsetFeatureExtractor::extract_vector(
  const ConstTraxelRefVector& /*traxelrefs*/
) {
  throw std::runtime_error(
    "Feature extractor " + name() + " doesn't provide a vector valued result"
  );
  return *(new FeatureVector);
}
FeatureScalar SubsetFeatureExtractor::extract_scalar(
  const ConstTraxelRefVector& /*traxelrefs*/
) {
  throw std::runtime_error(
    "Feature extractor " + name() + " doesn't provide a scalar valued result"
  );
  FeatureScalar ret_ = 0.;
  return ret_;
}

////
//// class SubsetFeatureCalculator
////
const FeatureMatrix& SubsetFeatureCalculator::calculate_matrix(
  const FeatureMatrix& /*feature_matrix*/
) {
  throw std::runtime_error(
    "Feature calculator " + name() + " doesn't provide a matrix valued result"
  );
  return *(new FeatureMatrix);
}
const FeatureVector& SubsetFeatureCalculator::calculate_vector(
  const FeatureMatrix& /*feature_matrix*/
) {
  throw std::runtime_error(
    "Feature calculator " + name() + " doesn't provide a vector valued result"
  );
  return *(new FeatureVector);
}
FeatureScalar SubsetFeatureCalculator::calculate_scalar(
  const FeatureMatrix& /*feature_matrix*/
) {
  throw std::runtime_error(
    "Feature calculator " + name() + " doesn't provide a scalar valued result"
  );
  FeatureScalar ret_ = 0.;
  return ret_;
}

/*=============================================================================
  specific classes
=============================================================================*/
////
//// class SubsetFeaturesIdentity
////
const std::string SubsetFeaturesIdentity::name_ = "SubsetFeaturesIdentity";

SubsetFeaturesIdentity::SubsetFeaturesIdentity(
  const std::vector<std::string>& feature_names
) : feature_names_(feature_names) {
}

SubsetFeaturesIdentity::SubsetFeaturesIdentity(
  const std::string& feature_name
) {
  feature_names_.resize(1);
  feature_names_[0] = feature_name;
}

const std::string& SubsetFeaturesIdentity::name() const {
  return name_;
}

const FeatureMatrix& SubsetFeaturesIdentity::extract_matrix(
  const ConstTraxelRefVector& traxelrefs
) {
  // get the size of the return matrix
  size_t x_size = traxelrefs.size();

  // get y-size
  size_t y_size = 0;

  // get the feature map of the first traxel and iterate over the feature names
  const FeatureMap& feature_map = traxelrefs.front()->features;
  for (
    std::vector<std::string>::const_iterator fname_it = feature_names_.begin();
    fname_it != feature_names_.end();
    ++fname_it
  ) {
    FeatureMap::const_iterator feature_map_it = feature_map.find(*fname_it);
    if (feature_map_it != feature_map.end()) {
      y_size += feature_map_it->second.size();
    } else {
      LOG(logDEBUG) << "In SubsetFeaturesIdentity: Feature \"" << *fname_it << "\" not found";
    }
  }
  
  // initialize the return matrix
  ret_matrix_.reshape(vigra::Shape2(x_size, y_size));

  // iterate over all traxel
  size_t column_index = 0;
  for(
    ConstTraxelRefVector::const_iterator tref_it = traxelrefs.begin();
    tref_it != traxelrefs.end();
    tref_it++, column_index++
  ) {
    // fetch the features which are stored in a map
    const FeatureMap feature_map = (*tref_it)->features;

    // get the current column as a view
    FeatureVectorView column = ret_matrix_.bind<0>(column_index);

    // reset the row_index
    size_t row_index = 0;

    // iterate over all feature names
    for(
      std::vector<std::string>::const_iterator fname_it = feature_names_.begin();
      fname_it != feature_names_.end();
      fname_it++
    ) {
      // check if the feature name exists
      FeatureMap::const_iterator fmap_it = feature_map.find(*fname_it);
      if ( fmap_it != feature_map.end()) {
        // iterate over the elements of the feature vector
        for (
          feature_array::const_iterator f_it = fmap_it->second.begin();
          f_it != fmap_it->second.end();
          f_it++, row_index++
        ) {
          if (row_index >= y_size) {
            throw std::runtime_error(
              "In SubsetFeaturesIdentity: Invalid row index"
            );
          }
          column(row_index) = *f_it;
        }
      } // if (fmap_it != feature_map.end())"
    } // for (fname_it = ... )
  } // for(traxels_it = .. )
  return ret_matrix_;
}

////
//// class TrackSubsets
////
const std::string TrackSubsets::name_ = "TrackSubsets";

const std::string& TrackSubsets::name() const {
  return name_;
}

const std::vector<ConstTraxelRefVector>& TrackSubsets::operator()(
  const HypothesesGraph& graph
) {
  ret_.clear();

  // Check if the graph has the necessary attributes
  if (not graph.has_property(node_active())) {
    throw std::runtime_error(
      "Graph doesn't have a \"node_active\" property map"
    );
  }
  if (not graph.has_property(arc_active())) {
    throw std::runtime_error(
      "Graph doesn't have an \"arc_active\" property map"
    );
  }

  // check if we have a tracklet graph
  bool with_tracklets = graph.has_property(node_tracklet());

  // check if the graph is legal
  if (not (graph.has_property(node_traxel()) or with_tracklets)) {
    throw std::runtime_error(
      "HypothesesGraph has neither traxel nor tracklet property map"
    );
  }

  // Get the property maps
  node_active_map_type& node_active_map = graph.get(node_active());
  arc_active_map_type& arc_active_map = graph.get(arc_active());

  // Make maps from child to parent and parent to child
  typedef std::map<HypothesesGraph::Node, HypothesesGraph::Node> NodeNodeMap;
  NodeNodeMap parent;
  NodeNodeMap child;
  // Initialize
  for (NodeActiveIt n_it(node_active_map); n_it != lemon::INVALID; ++n_it) {
    parent[n_it] = n_it;
    child[n_it] = n_it;
  }

  // Set connections
  for (ArcActiveIt a_it(arc_active_map); a_it != lemon::INVALID; ++a_it) {
    // count the active arcs with the same source
    size_t out_arcs = 0;
    for (
      OutArcIt o_it(graph, graph.source(a_it));
      o_it != lemon::INVALID;
      ++o_it
    ) {
      out_arcs += (arc_active_map[o_it] ? 1 : 0);
    }
    // link those nodes if there are no other active arcs with the same source
    if (out_arcs == 1) {
      parent[graph.target(a_it)] = graph.source(a_it);
      child[graph.source(a_it)] = graph.target(a_it);
    }
  }
  
  // Compose return vector of traxel reference vectors
  for (
    NodeNodeMap::const_iterator nmap_it = parent.begin();
    nmap_it != parent.end();
    ++nmap_it
  ) {
    HypothesesGraph::Node current_node = nmap_it->first;
    LOG(logDEBUG4) << "Is parent node invalid?";
    LOG(logDEBUG4) << (parent[current_node] == lemon::INVALID);
    if (parent[current_node] == current_node) {
      // resize the return vector
      ret_.resize(ret_.size()+1);
      bool loop = true;
      // loop as long as the track isn't finished
      while (loop) {
        if (with_tracklets) {
          // get the traxel vector of this node
          const std::vector<Traxel>& t_vec = graph.get(node_tracklet())[current_node];
          for (
            std::vector<Traxel>::const_iterator t_it = t_vec.begin();
            t_it != t_vec.end();
            t_it++
          ) {
            ret_.back().push_back( &(*t_it) );
          }
        } else {
          ret_.back().push_back( &(graph.get(node_traxel())[current_node]) );
        }
        loop = current_node != child[current_node];
        current_node = child[current_node];
      }
    }
  }
  return ret_;
}

////
//// class DivisionSubsets
////
const std::string DivisionSubsets::name_ = "DivisionSubsets";

const std::string& DivisionSubsets::name() const {
  return name_;
}

const std::vector<ConstTraxelRefVector>& DivisionSubsets::operator()(
  const HypothesesGraph& graph
) {
  return operator()(graph, depth_);
}
const std::vector<ConstTraxelRefVector>& DivisionSubsets::operator()(
  const HypothesesGraph& graph,
  size_t depth
) {
  // Check if the graph has the necessary attributes
  if (not graph.has_property(node_active())) {
    throw std::runtime_error(
      "Graph doesn't have a \"node_active\" property map"
    );
  }
  if (not graph.has_property(arc_active())) {
    throw std::runtime_error(
      "Graph doesn't have an \"arc_active\" property map"
    );
  }

  // check if we have a tracklet graph
  bool with_tracklets = graph.has_property(node_tracklet());

  // check if the graph is legal
  if (not (graph.has_property(node_traxel()) or with_tracklets)) {
    throw std::runtime_error(
      "HypothesesGraph has neither traxel nor tracklet property map"
    );
  }

  // make shure the depth is a valid value
  if (depth <=0 ) {
    depth = 1;
  }

  // call calculation functions regarding what kind of graph we have
  if (with_tracklets) {
    return from_tracklet_graph(graph, depth);
  } else {
    return from_traxel_graph(graph, depth);
  }
}

const std::vector<ConstTraxelRefVector>& DivisionSubsets::from_traxel_graph(
  const HypothesesGraph& graph,
  size_t depth
) {
  ret_.clear();
  // Get the property maps
  node_active_map_type& node_active_map = graph.get(node_active());
  node_traxel_map_type& node_traxel_map = graph.get(node_traxel());

  // Find the divisions
  for (NodeActiveIt n_it(node_active_map); n_it != lemon::INVALID; ++n_it) {
    // Count the active outgoing arcs
    std::vector<HypothesesGraph::Node> out_nodes;
    get_out_nodes(n_it, graph, out_nodes);
    ConstTraxelRefVector parents;
    ConstTraxelRefVector l_children;
    ConstTraxelRefVector r_children;
    // Two outgoing arcs: division
    if (out_nodes.size() == 2) {
      // Initialize the variables that change during the while loop
      // Now follows some ugly code
      bool valid_division = true;
      size_t curr_depth = depth;
      // These variables are the current parent node and the child nodes
      HypothesesGraph::Node parent = n_it;
      HypothesesGraph::Node l_node = out_nodes[0];
      HypothesesGraph::Node r_node = out_nodes[1];
      // Reserve some space for the vectors in which the incoming and out going
      // nodes are written
      Nodevector in;
      Nodevector l_out;
      Nodevector r_out;
      while ((curr_depth != 0) and valid_division) {
        // Save the reference to the traxel corresponding to the nodes
        parents.push_back( &(node_traxel_map[parent]) );
        l_children.push_back( &(node_traxel_map[l_node]) );
        r_children.push_back( &(node_traxel_map[r_node]) );

        // Get the following nodes
        get_in_nodes(parent, graph, in);
        get_out_nodes(l_node, graph, l_out);
        get_out_nodes(r_node, graph, r_out);

        // Check if the track is long enough to return the division to the full
        // depth
        valid_division  = (in.size() == 1);
        valid_division &= (l_out.size() == 1);
        valid_division &= (r_out.size() == 1);

        // If all sizes fit get the following nodes
        if (valid_division) {
          parent = in.front();
          l_node = l_out.front();
          r_node = r_out.front();
        }
        curr_depth--;
      }

      // store the results if the while loop ran over the whole depth
      if (curr_depth == 0) {
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
      }
    }
  }
  return ret_;
}

const std::vector<ConstTraxelRefVector>& DivisionSubsets::from_tracklet_graph(
  const HypothesesGraph& graph,
  size_t depth
) {
  ret_.clear();
  // Get the property maps
  node_active_map_type& node_active_map = graph.get(node_active());

  // Find the divisions
  for (NodeActiveIt n_it(node_active_map); n_it != lemon::INVALID; ++n_it) {
    // Count the active outgoing arcs
    std::vector<HypothesesGraph::Node> out_nodes;
    get_out_nodes(n_it, graph, out_nodes);
    // Two outgoing arcs: division
    if (out_nodes.size() == 2) {
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

      if (valid_division) {
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
    }
  }
  return ret_;
}

bool DivisionSubsets::get_children_to_depth(
  const HypothesesGraph::Node& node,
  const HypothesesGraph& graph,
  size_t depth,
  ConstTraxelRefVector& traxelrefs
) {
  HypothesesGraph::Node curr_node = node;
  size_t curr_tracklet_index = 0;
  bool valid_division = true;
  node_tracklet_map_type& tracklet_map = graph.get(node_tracklet());

  while (valid_division and (depth != 0)) {
    traxelrefs.push_back( &(tracklet_map[curr_node][curr_tracklet_index]) );

    curr_tracklet_index++;
    depth--;
    if (curr_tracklet_index == tracklet_map[curr_node].size()) {
      Nodevector out_nodes;
      get_out_nodes(curr_node, graph, out_nodes);
      valid_division = (out_nodes.size() == 1);
      if (valid_division) {
        curr_node = out_nodes[0];
        curr_tracklet_index = 0;
      }
    }
  }
  return (depth == 0);
}

bool DivisionSubsets::get_parents_to_depth(
  const HypothesesGraph::Node& node,
  const HypothesesGraph& graph,
  size_t depth,
  ConstTraxelRefVector& traxelrefs
) {
  node_tracklet_map_type& tracklet_map = graph.get(node_tracklet());

  HypothesesGraph::Node curr_node = node;
  size_t curr_tracklet_index = tracklet_map[node].size()-1;
  bool valid_division = true;

  while (valid_division and (depth != 0)) {
    traxelrefs.push_back( &(tracklet_map[curr_node][curr_tracklet_index]) );

    if (curr_tracklet_index == 0) {
      Nodevector in_nodes;
      get_in_nodes(curr_node, graph, in_nodes);
      valid_division = (in_nodes.size() == 1);
      if (valid_division) {
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
//// class SumCalculator
////
const std::string SumCalculator::name_ = "SumCalculator";

const std::string& SumCalculator::name() const {
  return name_;
}

const FeatureVector& SumCalculator::calculate_vector(
  const FeatureMatrix& features
) {
  size_t col_count = features.shape(0);
  size_t row_count = features.shape(1);
  ret_.reshape(vigra::Shape1(row_count));
  ret_.init(0);
  for (size_t col = 0; col < col_count; col++) {
    ret_ += features.bind<0>(col);
  }
  return ret_;
}

FeatureScalar SumCalculator::calculate_scalar(
  const FeatureMatrix& features
) {
  FeatureScalar ret_ = 0;
  for (int i = 0; i < features.size(); i++) {
    ret_ += features[i];
  }
  return ret_;
}

////
//// class DiffCalculator
////
const std::string DiffCalculator::name_ = "DiffCalculator";

const std::string& DiffCalculator::name() const {
  return name_;
}

const FeatureMatrix& DiffCalculator::calculate_matrix(
  const FeatureMatrix& feature_matrix
) {
  size_t col_count = feature_matrix.shape(0);
  size_t row_count = feature_matrix.shape(1);
  if (col_count <= 1) {
    LOG(logDEBUG) << "In DiffCalculator: matrix in argument has less than one column";
    LOG(logDEBUG) << "Returning a 0-vector";
    ret_.reshape(vigra::Shape2(1, row_count));
    ret_.init(0);
  } else {
    ret_.reshape(vigra::Shape2(col_count-1, row_count));
    for (size_t col = 0; col < col_count-1; col++) {
      FeatureVectorView a = feature_matrix.bind<0>(col);
      FeatureVectorView b = feature_matrix.bind<0>(col+1);
      using namespace vigra::multi_math;
      ret_.bind<0>(col) = b - a;
    }
  }
  return ret_;
}

////
//// class CurveCalculator
////
const std::string CurveCalculator::name_ = "CurveCalculator";

const std::string& CurveCalculator::name() const {
  return name_;
}

const FeatureMatrix& CurveCalculator::calculate_matrix(
  const FeatureMatrix& feature_matrix
) {
  size_t col_count = feature_matrix.shape(0);
  size_t row_count = feature_matrix.shape(1);
  if (col_count <= 2) {
    LOG(logDEBUG) << "In CurveCalculator: matrix in argument has less than two columns";
    LOG(logDEBUG) << "Returning a 0-vector";
    ret_.reshape(vigra::Shape2(1, row_count));
    ret_.init(0);
  } else {
    ret_.reshape(vigra::Shape2(col_count-2, row_count));
    for (size_t col = 0; col < col_count-2; col++) {
      using namespace vigra::multi_math;
      FeatureVectorView a = feature_matrix.bind<0>(col);
      FeatureVectorView b = feature_matrix.bind<0>(col+1);
      FeatureVectorView c = feature_matrix.bind<0>(col+2);
      ret_.bind<0>(col) = a - 2*b + c;
    }
  }
  return ret_;
}

////
//// class InvCovarianceCalculator
////
const std::string InvCovarianceCalculator::name_ = "InvCovaricaneCalculator";

const std::string& InvCovarianceCalculator::name() const {
  return name_;
}

const FeatureMatrix& InvCovarianceCalculator::calculate_matrix(
  const FeatureMatrix& feature_matrix
) {
  size_t col_count = feature_matrix.shape(0);
  size_t row_count = feature_matrix.shape(1);
  ret_.reshape(vigra::Shape2(row_count, row_count));
  if (col_count <= row_count) {
    LOG(logDEBUG) << "In InvCovarianceCalculator: too few data to calculate covariance matrix";
    LOG(logDEBUG) << "Returning a Matrix filled with zeros";
    ret_.reshape(vigra::Shape2(row_count, row_count));
    ret_.init(0);
  }
  FeatureMatrixView covariance_matrix;
  vigra::linalg::covarianceMatrixOfRows(feature_matrix, covariance_matrix);
  bool invertible = vigra::linalg::inverse(feature_matrix, ret_);
  if (not invertible) {
    ret_.init(0);
  }
  return ret_;
}

/*
////
//// class OutlierCountAggregator
////
const std::string OutlierCountAggregator::name_ = "OutlierCountAggregator";

const std::string& OutlierCountAggregator::name() const {
  return name_;
}

const feature_type& OutlierCountAggregator::operator()(
  const feature_arrays& features
) {
  std::vector<size_t> outlier_ids = outlier_calculator_->calculate(features);
  ret_ = static_cast<feature_type>(outlier_ids.size());
  // prevent divide by 0 error
  if (features.size()) {
    ret_ = ret_ / static_cast<feature_type>(features.size());
  }
  return ret_;
}

////
//// class OutlierBadnessAggregator
////
const std::string OutlierBadnessAggregator::name_ = "OutlierBadnessAggregator";

const std::string& OutlierBadnessAggregator::name() const {
  return name_;
}

const feature_type& OutlierBadnessAggregator::operator()(
  const feature_arrays& features
) {
  outlier_calculator_->calculate(features);
  feature_array outlier_badness = outlier_calculator_->get_measures();
  ret_ = 0.0;
  for (
    feature_array::const_iterator f_it = outlier_badness.begin();
    f_it != outlier_badness.end();
    f_it++
  ) {
    ret_ = *f_it > ret_ ? *f_it : ret_;
  }
  return ret_;
}


////
//// function to_arma_matrix
////
arma::Mat<feature_type> to_arma_matrix(const feature_arrays& features) {
  assert(features.empty() == false);
  size_t cols = features.size();
  size_t rows = features[0].size();
  arma::Mat<feature_type> ret(rows, cols);

  typename std::vector<feature_array>::const_iterator feature_array_it;
  feature_array_it = features.begin();
  for (
    size_t j=0;
    feature_array_it != features.end();
    feature_array_it++, j++
  ) {
    assert(feature_array_it->size() == rows);
    arma::Col<feature_type> column(*feature_array_it);
    ret.col(j) = column;
  }
  return ret;
}

////
//// class OutlierCalculator
////
const std::string OutlierCalculator::name_ = "OutlierCalculator";

const std::string& OutlierCalculator::name() const {
  return name_;
}

const feature_array& OutlierCalculator::get_measures() const {
  throw std::runtime_error(
    "OutlierCalculator \"" + name() + "\"doesn't provide a measure"
  );
  return measures_;
}

////
//// class MVNOutlierCalculator
////
const std::string MVNOutlierCalculator::name_ = "MVNOutlierCalculator";

MVNOutlierCalculator::MVNOutlierCalculator(const feature_type sigma_threshold) {
  sigma_threshold_ = sigma_threshold;
}

const std::string& MVNOutlierCalculator::name() const {
  return name_;
}

const feature_array& MVNOutlierCalculator::get_measures() const {
  return measures_;
}

const arma::Mat<feature_type>& MVNOutlierCalculator::get_covariance() const {
  return covariance_;
}

const arma::Mat<feature_type>& MVNOutlierCalculator::get_inverse_covariance() const {
  return inv_covariance_;
}

const arma::Col<feature_type>& MVNOutlierCalculator::get_mean() const {
  return mean_;
}

const std::vector<size_t>& MVNOutlierCalculator::calculate(
  const feature_arrays& features
) {
  measures_.clear();
  outlier_ids_.clear();
  mean_.clear();
  covariance_.clear();
  inv_covariance_.clear();

  bool good_data = true;
  if (features.size() == 0) {
    good_data = false;
  } else if (features.size() <= features[0].size()+1) {
    good_data = false;
  } 
  if (good_data) {
    // Get covariance and inverse covariance matrix
    arma::Mat<feature_type> features_mat(to_arma_matrix(features));
    arma::Mat<feature_type> features_mat_t(trans(features_mat));
    covariance_ = arma::cov(features_mat_t);
    bool invertible = arma::inv_sympd(inv_covariance_, covariance_);

    if (invertible) {
      // Get mean values
      mean_ = arma::mean(features_mat, 1);
  
      // Calculate the outliers
      outlier_ids_.clear();
      measures_.clear();
      feature_arrays::const_iterator features_it = features.begin();
      for(size_t id=0; features_it != features.end(); features_it++, id++) {
        arma::Col<feature_type> diff_vector(*features_it);
        diff_vector -= mean_;
        feature_type norm_residual = arma::dot(diff_vector, inv_covariance_*diff_vector);
        measures_.push_back(norm_residual);
        if (norm_residual >= sigma_threshold_) {
          outlier_ids_.push_back(id);
        }
      }
    }
  } // else
  return outlier_ids_;
}
*/
} // namespace pgmlink
