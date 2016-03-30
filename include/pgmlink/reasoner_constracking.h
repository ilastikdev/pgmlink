#ifndef CONSTRACKING_REASONER_H
#define CONSTRACKING_REASONER_H

#include <map>
#include <boost/function.hpp>
#include <boost/python.hpp>

#include "pgmlink/pgm.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/reasoner.h"
#include "pgmlink/features/feature.h"
#include "pgmlink/uncertaintyParameter.h"
#include "pgmlink/inferencemodel/constrackinginferencemodel.h"
#include "pgmlink/inferencemodel/perturbation/perturbation.h"
#include "pgmlink/inferencemodel/perturbedinferencemodel.h"
#include "pgmlink/conservationtracking_parameter.h"

#include <opengm/opengm.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>

namespace pgmlink
{

class Traxel;

class ConservationTracking : public Reasoner
{
public:
    typedef std::vector<pgm::OpengmModelDeprecated::ogmInference::LabelType> IlpSolution;

    ConservationTracking(const Parameter& param);
    ~ConservationTracking();

    virtual void infer();
    virtual void conclude(HypothesesGraph&);
    virtual void formulate( const HypothesesGraph& );
    virtual void perturbedInference(HypothesesGraph&);//, Parameter&);

    /**
    Run the dynamic programming solver first and then use it as initialization
    for the globally optimal solver.
    */
    void twoStageInference(HypothesesGraph & hypotheses);

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


    //    virtual boost::shared_ptr<InferenceModel> create_inference_model(ConservationTracking::Parameter& param);
    virtual boost::shared_ptr<InferenceModel> create_inference_model(Parameter& param);
    boost::shared_ptr<InferenceModel> create_inference_model();
    void setInferenceModel(boost::shared_ptr<InferenceModel> inference_model);
    boost::shared_ptr<InferenceModel> getInferenceModel();
    //void fixNodeToAppearanceLabel( HypothesesGraph& hypothesesGraph, int node, double label);
    bool training_to_hard_constraints_;

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
    boost::function<double (const Traxel&, const size_t)> divisionNoWeight_;
    boost::function<double (const double)> transition_;
    //boost::function<double (const Traxel&, const Traxel&, const size_t)> transition_;

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

    Parameter param_;
    UncertaintyParameter uncertainty_param_;
    Perturbation::Parameter perturbed_inference_model_param_;

    double cplex_timeout_;
    unsigned int num_threads_;
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
//   void ConservationTracking::add_factor( const table_t& table, const_iter first_idx, const_iter last_idx ){
//   OpengmModelDeprecated::FunctionIdentifier id=pgm_->Model()->addFunction(table);
//   pgm_->Model()->addFactor(id, first_idx, last_idx);
// }

} /* namespace pgmlink */
#endif /* MRF_REASONER_H */

