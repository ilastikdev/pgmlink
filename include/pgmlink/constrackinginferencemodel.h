#ifndef CONSTRACKINGINFERENCEMODEL_H
#define CONSTRACKINGINFERENCEMODEL_H

#include <boost/function.hpp>
#include <opengm/graphicalmodel/graphicalmodel.hxx>

#include "pgmlink/hypotheses.h"
#include "pgmlink/pgm.h"
#include "pgmlink/constraint_pool.hxx"
#include <boost/python.hpp>

namespace pgmlink
{

enum EnergyType {Appearance=0, Disappearance=1, Detection=2, Transition=3, Division=4 };

/**
 * @brief The ConsTrackingInferenceModel class builds the OpenGM model needed to run basic conservation tracking.
 * Derived classes such as PerturbedInferenceModel can extend the functionality to support more advanced models.
 * The general usage is to set up this inference model from a hypotheses graph, retrieve the OpenGM model, create
 * an optimizer, add the constraints, and run inference.
 */
class ConsTrackingInferenceModel
{
public: // typedefs
    typedef double ValueType;
    typedef pgm::OpengmModelDeprecated::ogmGraphicalModel::OperatorType OperatorType;
    typedef pgm::OpengmModelDeprecated::ogmGraphicalModel::LabelType LabelType;
    typedef pgm::OpengmModelDeprecated::ogmGraphicalModel::IndexType IndexType;
    typedef PertGmType GraphicalModelType;
    typedef std::map<HypothesesGraph::Node, size_t> HypothesesGraphNodeMap;
    typedef std::map<HypothesesGraph::Arc, size_t> HypothesesGraphArcMap;
    typedef std::map<std::pair<Traxel, Traxel >, std::pair<double, double > > TransitionPredictionsMap;

public: // Parameter object
    class Parameter{
    public:
        size_t max_number_objects;

        bool with_constraints;
        bool with_tracklets;
        bool with_divisions;
        bool with_optical_correction;
        bool with_misdetections_allowed;
        bool with_appearance;
        bool with_disappearance;

        boost::function<double (const Traxel&, const size_t)> detection;
        boost::function<double (const Traxel&, const size_t)> division;
        boost::function<double (const double)> transition;
        double transition_parameter;
        boost::python::object transition_classifier;

        boost::function<double (const Traxel&)> disappearance_cost;
        boost::function<double (const Traxel&)> appearance_cost;
        double forbidden_cost;
    };

public: // API
    ConsTrackingInferenceModel(const Parameter& param);

    // build the inference model from the given graph
    void build_from_graph(const HypothesesGraph&);

    // add constraints to the model or the optimizer, depending on whether INF supports hard constraints
    template<class INF>
    void add_constraints(INF& optimizer);

    // extract the model
    GraphicalModelType& get_model();

    // retrieve node and arc maps
    HypothesesGraphNodeMap& get_division_node_map();
    HypothesesGraphNodeMap& get_appearance_node_map();
    HypothesesGraphNodeMap& get_disappearance_node_map();
    HypothesesGraphArcMap& get_arc_map();

    // output
    void printResults(const HypothesesGraph &g);

protected: // methods
    void add_appearance_nodes( const HypothesesGraph& );
    void add_disappearance_nodes( const HypothesesGraph& );
    void add_transition_nodes( const HypothesesGraph& );
    void add_division_nodes(const HypothesesGraph& );

    void add_finite_factors(const HypothesesGraph& );
    size_t add_division_factors(const HypothesesGraph &g, size_t factorIndex);
    size_t add_transition_factors(const HypothesesGraph &g, size_t factorIndex);
    size_t add_detection_factors(const HypothesesGraph &g, size_t factorIndex);

    void add_constraints_to_pool(const HypothesesGraph& );
    bool callable(boost::python::object object);

    double get_transition_prob(double distance, size_t state, double alpha);
    virtual double get_transition_probability(Traxel &tr1, Traxel &tr2, size_t state);
    virtual double generateRandomOffset(EnergyType parameterIndex,  double energy = 0, Traxel tr = 0, Traxel tr2 = 0, size_t state = 0);
    virtual size_t add_div_m_best_perturbation(marray::Marray<double> &energies, EnergyType energy_type, size_t factorIndex);

protected: // members
    Parameter param_;
    GraphicalModelType model_;
    pgm::ConstraintPool constraint_pool_;

    HypothesesGraphNodeMap div_node_map_;
    HypothesesGraphNodeMap app_node_map_;
    HypothesesGraphNodeMap dis_node_map_;
    HypothesesGraphArcMap arc_map_;

    // remove?
    unsigned int number_of_transition_nodes_, number_of_division_nodes_;
    unsigned int number_of_appearance_nodes_, number_of_disappearance_nodes_;
    std::map< size_t, std::vector<size_t> > nodes_per_timestep_;
};

template<class INF>
void ConsTrackingInferenceModel::add_constraints(INF &optimizer)
{
    constraint_pool_.add_constraints_to_problem(model_, optimizer);
}

} // namespace pgmlink

#endif // CONSTRACKINGINFERENCEMODEL_H
