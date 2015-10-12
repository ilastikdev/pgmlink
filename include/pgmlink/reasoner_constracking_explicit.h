#ifndef CONSTRACKING_REASONER_EXPLICIT_H
#define CONSTRACKING_REASONER_EXPLICIT_H

#include <map>
#include <boost/function.hpp>

#include <boost/python.hpp>

#include "pgmlink/pgm.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/reasoner.h"
#include "pgmlink/features/feature.h"
#include "pgmlink/uncertaintyParameter.h"
#include "pgmlink/inferencemodel/constrackingexplicitinferencemodel.h"
#include "pgmlink/inferencemodel/perturbation/perturbation.h"
#include "pgmlink/inferencemodel/perturbedinferencemodel.h"

#include <opengm/opengm.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>

namespace pgmlink
{

class Traxel;

class ConservationExplicitTracking : public Reasoner
{
public:
    enum SolverType
    {
        CplexSolver,
        DynProgSolver
    };

    class Parameter
    {
    public:

        Parameter(
            unsigned int max_number_objects,
            boost::function<double (const Traxel&, const size_t)> detection,
            //boost::function<double (const Traxel&, const size_t)> detectionNoWeight,
            boost::function<double (const Traxel&, const size_t)> division,
            //boost::function<double (const double)> transition,
            boost::function<double (const Traxel&, const Traxel&, const size_t)> transition,
            double forbidden_cost = 0,
            double ep_gap = 0.01,
            bool with_tracklets = false,
            bool with_divisions = true,
            boost::function<double (const Traxel&)> disappearance_cost_fn = ConstantFeature(500.0),
            boost::function<double (const Traxel&)> appearance_cost_fn = ConstantFeature(500.0),
            bool with_misdetections_allowed = true,
            bool with_appearance = true,
            bool with_disappearance = true,
            double transition_parameter = 5,
            bool with_constraints = true,
            UncertaintyParameter uncertainty_param = UncertaintyParameter(),
            double cplex_timeout = 1e75,
            double division_weight = 10,
            double detection_weight = 10,
            double transition_weight = 10,
            double border_width = 0,
            boost::python::object transition_classifier = boost::python::object(),
            bool with_optical_correction = false,
            SolverType solver = CplexSolver):
            max_number_objects(max_number_objects),
            detection(detection),
            detectionNoWeight(0),
            division(division),
            transition(transition),
            border_width(border_width),
            forbidden_cost(forbidden_cost),
            ep_gap(ep_gap),
            with_tracklets(with_tracklets),
            with_divisions(with_divisions),
            disappearance_cost_fn(disappearance_cost_fn),
            appearance_cost_fn(appearance_cost_fn),
            with_misdetections_allowed(with_misdetections_allowed),
            with_appearance(with_appearance),
            with_disappearance(with_disappearance),
            transition_parameter(transition_parameter),
            with_constraints(with_constraints),
            uncertainty_param(uncertainty_param),
            cplex_timeout(cplex_timeout),
            division_weight(division_weight),
            detection_weight(detection_weight),
            transition_weight(transition_weight),
            transition_classifier(transition_classifier),
            with_optical_correction(with_optical_correction),
            solver_(solver)
        {}

        // empty parameter needed for python
        Parameter() {}

        // settings
        unsigned int max_number_objects;
        boost::function<double (const Traxel&, const size_t)> detection;
        boost::function<double (const Traxel&, const size_t)> detectionNoWeight;
        boost::function<double (const Traxel&, const size_t)> division;
        //boost::function<double (const double)> transition;
        boost::function<double (const Traxel&, const Traxel&, const size_t)> transition;
        double forbidden_cost;
        double ep_gap;
        bool with_tracklets;
        bool with_divisions;
        boost::function<double (const Traxel&)> disappearance_cost_fn;
        boost::function<double (const Traxel&)> appearance_cost_fn;
        bool with_misdetections_allowed;
        bool with_appearance;
        bool with_disappearance;
        double transition_parameter;
        bool with_constraints;
        UncertaintyParameter uncertainty_param;
        double cplex_timeout;
        double division_weight;
        double detection_weight;
        double transition_weight;
        double appearance_weight;
        double disappearance_weight;
        double border_width;
        boost::python::object transition_classifier;
        bool with_optical_correction;
        SolverType solver_;
    };

public:
    typedef std::vector<pgm::OpengmModelDeprecated::ogmInference::LabelType> IlpSolution;

    ConservationExplicitTracking(const Parameter& param);
    ~ConservationExplicitTracking();

    virtual void infer();
    virtual void conclude(HypothesesGraph&);
    virtual void formulate( const HypothesesGraph& );
    virtual void perturbedInference(HypothesesGraph&);
    void enableFixingLabeledAppearanceNodes();

    double forbidden_cost() const;

    /// Return reference to all CPLEX solution vectors
    const std::vector<IlpSolution>& get_ilp_solutions() const;
    void set_ilp_solutions(const std::vector<IlpSolution>&);

    //cplex export file names
    std::string features_file_;
    std::string constraints_file_;
    std::string labels_export_file_name_;

    static std::string get_export_filename(size_t iteration, const std::string &orig_file_name);

    HypothesesGraph *get_prepared_graph(HypothesesGraph &hypotheses);
    //boost::shared_ptr<InferenceModel> create_inference_model();
    virtual boost::shared_ptr<InferenceModel> create_inference_model(ConservationExplicitTracking::Parameter& param);
    boost::shared_ptr<InferenceModel> create_inference_model();
    void setInferenceModel(boost::shared_ptr<InferenceModel> inference_model);
    boost::shared_ptr<InferenceModel> getInferenceModel();

protected: // methods
    void reset();
    void compute_relative_uncertainty(HypothesesGraph *graph);

    boost::shared_ptr<Perturbation> create_perturbation();
    boost::shared_ptr<InferenceModel> create_perturbed_inference_model(boost::shared_ptr<Perturbation> perturb);

protected: // members
    unsigned int max_number_objects_;

    // energy functions
    boost::function<double (const Traxel&, const size_t)> detection_;
    boost::function<double (const Traxel&, const size_t)> detectionNoWeight_;
    boost::function<double (const Traxel&, const size_t)> division_;
    //boost::function<double (const double)> transition_;
    boost::function<double (const Traxel&, const Traxel&, const size_t)> transition_;

    double forbidden_cost_;
    std::vector<IlpSolution> solutions_;

    double ep_gap_;
    bool with_tracklets_;
    bool with_divisions_;
    boost::function<double (const Traxel&)> disappearance_cost_;
    boost::function<double (const Traxel&)> appearance_cost_;

    bool with_misdetections_allowed_;
    bool with_appearance_;
    bool with_disappearance_;
    bool with_constraints_;

    double transition_parameter_;

    UncertaintyParameter uncertainty_param_;
    ConsTrackingInferenceModel::Parameter inference_model_param_;
    Perturbation::Parameter perturbed_inference_model_param_;

    double cplex_timeout_;
    bool isMAP_;
    SolverType solver_;

    double division_weight_;
    double detection_weight_;
    double transition_weight_;

    boost::python::object transition_classifier_;

    bool with_optical_correction_;
    bool use_app_node_labels_to_fix_values_;

    HypothesesGraph tracklet_graph_;
    std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > tracklet2traxel_node_map_;

    boost::shared_ptr<InferenceModel> inference_model_;
    bool with_structured_learning_;
};



/******************/
/* Implementation */
/******************/

// template< typename table_t, typename const_iter >
//   void ConservationExplicitTracking::add_factor( const table_t& table, const_iter first_idx, const_iter last_idx ){
//   OpengmModelDeprecated::FunctionIdentifier id=pgm_->Model()->addFunction(table);
//   pgm_->Model()->addFactor(id, first_idx, last_idx);
// }

} /* namespace pgmlink */
#endif /* MRF_REASONER_H */

