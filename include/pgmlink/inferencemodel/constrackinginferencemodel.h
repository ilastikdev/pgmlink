#ifndef CONSTRACKINGINFERENCEMODEL_H
#define CONSTRACKINGINFERENCEMODEL_H

#include <boost/function.hpp>

#include <opengm/graphicalmodel/graphicalmodel.hxx>
#include <opengm/inference/inference.hxx>
#include <opengm/inference/lpcplex.hxx>

#include "pgmlink/hypotheses.h"
#include "pgmlink/pgm.h"
#include "pgmlink/inferencemodel/constraint_pool.hxx"
#include "pgmlink/inferencemodel/inferencemodel.h"

namespace pgmlink
{

/**
 * @brief The ConsTrackingInferenceModel class builds the OpenGM model needed to run basic conservation tracking.
 * Derived classes such as PerturbedInferenceModel can extend the functionality to support more advanced models.
 * The general usage is to set up this inference model from a hypotheses graph, retrieve the OpenGM model, create
 * an optimizer, add the constraints, and run inference.
 */
class ConsTrackingInferenceModel : public InferenceModel
{
public: // typedefs
    typedef double ValueType;
    typedef pgm::OpengmModelDeprecated::ogmGraphicalModel::OperatorType OperatorType;
    typedef pgm::OpengmModelDeprecated::ogmGraphicalModel::LabelType LabelType;
    typedef pgm::OpengmModelDeprecated::ogmGraphicalModel::IndexType IndexType;
    typedef std::vector<LabelType> IlpSolution;
    typedef PertGmType GraphicalModelType;
    typedef opengm::LPCplex<PertGmType, pgm::OpengmModelDeprecated::ogmAccumulator> cplex_optimizer;
    typedef std::map<HypothesesGraph::Node, size_t> HypothesesGraphNodeMap;
    typedef std::map<HypothesesGraph::Arc, size_t> HypothesesGraphArcMap;

public: // API
    // constructor
    ConsTrackingInferenceModel(const Parameter& param, double ep_gap, double cplex_timeout);

    // build the inference model from the given graph
    virtual void build_from_graph(const HypothesesGraph&);

    virtual void fixFirstDisappearanceNodesToLabels(
            const HypothesesGraph& g,
            const HypothesesGraph &tracklet_graph,
            std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > &traxel2tracklet_map);

    // extract the model
    GraphicalModelType& get_model();

    // run inference
    virtual IlpSolution infer();
    void set_inference_params(size_t numberOfSolutions,
                              const std::string& feature_filename,
                              const std::string& constraints_filename,
                              const std::string& ground_truth_filename);

    IlpSolution extractSolution(size_t k, const std::string& ground_trugh_filename);

    // write results to hypotheses graph
    virtual void conclude(HypothesesGraph &g,
                          HypothesesGraph &tracklet_graph,
                          std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > &tracklet2traxel_node_map,
                          IlpSolution &solution);

    // output
    void printResults(const HypothesesGraph &g);
    void write_labeledgraph_to_file(const HypothesesGraph & g,
                                    const std::string &ground_truth_filename);

    // structured learning tracking inference model
    opengm::learning::Weights<double> weights_;
    // weights in the same order as in:
    // enum EnergyType {Appearance = 0, Disappearance = 1, Detection = 2, Transition = 3, Division = 4 };
    void setWeight ( size_t, double);
    GraphicalModelType model();
    unsigned int get_number_of_division_nodes();

protected: // methods
    void add_appearance_nodes( const HypothesesGraph& );
    void add_disappearance_nodes( const HypothesesGraph& );
    void add_transition_nodes( const HypothesesGraph& );
    void add_division_nodes(const HypothesesGraph& );

    void add_finite_factors(const HypothesesGraph& );
    virtual size_t add_division_factors(const HypothesesGraph &g, size_t factorIndex);
    virtual size_t add_transition_factors(const HypothesesGraph &g, size_t factorIndex);
    virtual size_t add_detection_factors(const HypothesesGraph &g, size_t factorIndex);

    void add_constraints_to_pool(const HypothesesGraph& );

    // retrieve node and arc maps
    HypothesesGraphNodeMap& get_division_node_map();
    HypothesesGraphNodeMap& get_appearance_node_map();
    HypothesesGraphNodeMap& get_disappearance_node_map();
    HypothesesGraphArcMap& get_arc_map();
    std::map<HypothesesGraph::Node, size_t>& get_detection_factor_node_map();

    // add constraints to the model or the optimizer, depending on whether INF supports hard constraints
    template<class INF>
    void add_constraints(INF& optimizer);

protected: // members
    GraphicalModelType model_;

    HypothesesGraphNodeMap div_node_map_;
    HypothesesGraphNodeMap app_node_map_;
    HypothesesGraphNodeMap dis_node_map_;
    HypothesesGraphArcMap arc_map_;
    //factor id maps
    std::map<HypothesesGraph::Node, size_t> detection_f_node_map_;

    // optimizer
    cplex_optimizer::Parameter cplex_param_;
    boost::shared_ptr<cplex_optimizer> optimizer_;
    pgm::ConstraintPool constraint_pool_;

    // remove?
    unsigned int number_of_transition_nodes_, number_of_division_nodes_;
    unsigned int number_of_appearance_nodes_, number_of_disappearance_nodes_;
    std::map< size_t, std::vector<size_t> > nodes_per_timestep_;

    // funky export maps
    bool export_from_labeled_graph_;
    std::string ground_truth_filename_;
    std::map<std::pair<size_t, size_t>, size_t > cplex_variable_id_map_;
    std::map<std::pair<size_t, std::pair<size_t, size_t> >, size_t> cplex_factor_id_map_;

};

template<class INF>
void ConsTrackingInferenceModel::add_constraints(INF &optimizer)
{
    constraint_pool_.add_constraints_to_problem(model_, optimizer);
}

} // namespace pgmlink

#endif // CONSTRACKINGINFERENCEMODEL_H
