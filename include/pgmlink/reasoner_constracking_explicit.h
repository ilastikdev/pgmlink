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
            bool with_merger_resolution = false,
            int n_dim = 2,
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
            with_merger_resolution(with_merger_resolution),
            n_dim(n_dim),
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
        boost::function<double (const Traxel&, const Traxel&, const Traxel&)> motion_model3;
        boost::function<double (const Traxel&, const Traxel&, const Traxel&, const Traxel&)> motion_model4;
        double motion_model3_default;
        double motion_model4_default;
        double forbidden_cost;
        double ep_gap;
        bool with_tracklets;
        bool with_divisions;
        boost::function<double (const Traxel&)> disappearance_cost_fn;
        boost::function<double (const Traxel&)> appearance_cost_fn;
        bool with_merger_resolution;
        int n_dim;
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

    private:
        // python extensions:
        double python_caller_explicit_det_div(boost::python::object func, const Traxel& t, const size_t state)
        {
            assert(1 == PyCallable_Check(func.ptr()));
            // PyGILState_STATE pygilstate = PyGILState_Ensure();
            boost::python::object py_result = func(t, state);
            double result = boost::python::extract<double>(py_result);
            // PyGILState_Release(pygilstate);
            return result;
        }

        double python_caller_explicit_dis_appear(boost::python::object func, const Traxel& t)
        {
            assert(1 == PyCallable_Check(func.ptr()));
            // PyGILState_STATE pygilstate = PyGILState_Ensure();
            boost::python::object py_result = func(t);
            double result = boost::python::extract<double>(py_result);
            // PyGILState_Release(pygilstate);
            return result;
        }

        //double python_caller_trans(boost::python::object func, double distance)
        double python_caller_explicit_trans(boost::python::object func, const Traxel& a, const Traxel& b, const size_t state)
        {
            assert(1 == PyCallable_Check(func.ptr()));
            // PyGILState_STATE pygilstate = PyGILState_Ensure();
            //boost::python::object py_result = func(distance);
            boost::python::object py_result = func(a, b, state);
            double result = boost::python::extract<double>(py_result);
            // PyGILState_Release(pygilstate);
            return result;
        }

        double python_caller_explicit_motion_model3(boost::python::object func, const Traxel& a, const Traxel& b, const Traxel& c)
        {
            assert(1 == PyCallable_Check(func.ptr()));
            boost::python::object py_result = func(a, b, c);
            double result = boost::python::extract<double>(py_result);
            return result;
        }

        double python_caller_explicit_motion_model4(boost::python::object func, const Traxel& a, const Traxel& b, const Traxel& c, const Traxel& d)
        {
            assert(1 == PyCallable_Check(func.ptr()));
            boost::python::object py_result = func(a, b, c, d);
            double result = boost::python::extract<double>(py_result);
            return result;
        }

    public:
        /// Expects a function with signature (Traxel traxel, size_t state) -> double energy
        void register_explicit_detection_func(boost::python::object func)
        {
            detection = boost::bind(&ConservationExplicitTracking::Parameter::python_caller_explicit_det_div, this, func, _1, _2);
        }

        /// Expects a function with signature (Traxel traxel, size_t state) -> double energy
        void register_explicit_division_func(boost::python::object func)
        {
            division = boost::bind(&ConservationExplicitTracking::Parameter::python_caller_explicit_det_div, this, func, _1, _2);
        }

        /// Expects a function with signature (double distance) -> double energy
        void register_explicit_transition_func(boost::python::object func)
        {
            //transition = boost::bind(&ConservationTracking::Parameter::python_caller_trans, this, func, _1);
            transition = boost::bind(&ConservationExplicitTracking::Parameter::python_caller_explicit_trans, this, func, _1, _2, _3);
        }

        /// Expects a function with signature (Traxel traxel) -> double energy
        void register_explicit_appearance_func(boost::python::object func)
        {
            appearance_cost_fn = boost::bind(&ConservationExplicitTracking::Parameter::python_caller_explicit_dis_appear, this, func, _1);
        }

        /// Expects a function with signature (Traxel traxel) -> double energy
        void register_explicit_disappearance_func(boost::python::object func)
        {
            disappearance_cost_fn = boost::bind(&ConservationExplicitTracking::Parameter::python_caller_explicit_dis_appear, this, func, _1);
        }

        /// Expects a function with signature (Traxel, Traxel, Traxel) -> double energy
        void register_explicit_motion_model3_func(boost::python::object func, double default_value)
        {
            motion_model3 = boost::bind(&ConservationExplicitTracking::Parameter::python_caller_explicit_motion_model3, this, func, _1, _2, _3);
            motion_model3_default = default_value;
        }

        /// Expects a function with signature (Traxel, Traxel, Traxel, Traxel) -> double energy
        void register_explicit_motion_model4_func(boost::python::object func, double default_value)
        {
            motion_model4 = boost::bind(&ConservationExplicitTracking::Parameter::python_caller_explicit_motion_model4, this, func, _1, _2, _3, _4);
            motion_model4_default = default_value;
        }
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

    bool with_merger_resolution_;
    int n_dim_;
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

