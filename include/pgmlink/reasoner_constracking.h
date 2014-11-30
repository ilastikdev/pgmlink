#ifndef CONSTRACKING_REASONER_H
#define CONSTRACKING_REASONER_H

#include <map>
#include <boost/function.hpp>
#include <opengm/inference/inference.hxx>
#include <opengm/inference/lpcplex.hxx>

#include <boost/python.hpp>

#include "pgmlink/pgm.h"
#include "pgmlink/hypotheses.h"
#include "pgmlink/reasoner.h"
#include "pgmlink/feature.h"
#include "pgmlink/uncertaintyParameter.h"
#include "pgmlink/constrackinginferencemodel.h"
#include "pgmlink/perturbedinferencemodel.h"

#include <opengm/opengm.hxx>
#include <opengm/graphicalmodel/graphicalmodel.hxx>

namespace pgmlink {	

typedef pgm::OpengmModelDeprecated::ogmGraphicalModel MAPGmType;
		

typedef opengm::LPCplex
	<	PertGmType,
			pgm::OpengmModelDeprecated::ogmAccumulator		>
    cplex_optimizer;



typedef pgm::OpengmModelDeprecated::ogmGraphicalModel::FactorType factorType;

typedef typename boost::variate_generator<boost::mt19937, boost::normal_distribution<> > normalRNGType; 
typedef typename boost::variate_generator<boost::mt19937, boost::uniform_real<> > uniformRNGType;    

class Traxel;

class ConservationTracking : public Reasoner {
    public:
    typedef std::vector<pgm::OpengmModelDeprecated::ogmInference::LabelType> IlpSolution;

    ConservationTracking(
                             unsigned int max_number_objects,
                             boost::function<double (const Traxel&, const size_t)> detection,
                             boost::function<double (const Traxel&, const size_t)> division,
                             boost::function<double (const double)> transition,
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
                             UncertaintyParameter param = UncertaintyParameter(),
                             double cplex_timeout = 1e75,
                             double division_weight = 10,
                             double detection_weight = 10,
                             double transition_weight = 10,
                             boost::python::object transition_classifier = boost::python::object(),
                             bool with_optical_correction = false
                             );
    ~ConservationTracking();

    virtual void infer();
    virtual void conclude(HypothesesGraph&);
    virtual void conclude(HypothesesGraph&, boost::shared_ptr<ConsTrackingInferenceModel> inference_model);
    virtual void formulate( const HypothesesGraph& );
    virtual void perturbedInference(HypothesesGraph&, bool with_inference = true);
    
    double forbidden_cost() const;
    
    /// Return reference to all CPLEX solution vectors
    const std::vector<IlpSolution>& get_ilp_solutions() const;
    void set_ilp_solutions(const std::vector<IlpSolution>&);

    //cplex export file names
    std::string features_file_;
    std::string constraints_file_;
    std::string ground_truth_file_;

    bool export_from_labeled_graph_;

    void write_labeledgraph_to_file(const HypothesesGraph&, boost::shared_ptr<ConsTrackingInferenceModel> inference_model);

    static std::string get_export_filename(size_t iteration, const std::string &orig_file_name);
private:
    // copy and assingment have to be implemented, yet

//    ConservationTracking(const ConservationTracking&):
//		random_normal_(rng_,boost::normal_distribution<>(0, 1)),
//		random_uniform_(rng_,boost::uniform_real<>(0, 1))
//        {}
//    ConservationTracking& operator=(const ConservationTracking&) { return *this;}

    void reset();
    void add_constraints(const HypothesesGraph& );
    void add_appearance_nodes( const HypothesesGraph& );
    void add_disappearance_nodes( const HypothesesGraph& );
    void add_transition_nodes( const HypothesesGraph& );
    void add_division_nodes(const HypothesesGraph& );
    template <typename ModelType> void add_finite_factors( const HypothesesGraph&, ModelType* model, bool perturb= false, vector<vector<vector<size_t> > >* detoff=NULL );
    double getEnergyByEvent(EnergyType event, HypothesesGraph::NodeIt n,bool perturb=false,size_t state=0);
    void printResults(const HypothesesGraph &);
    double sample_with_classifier_variance(double mean, double variance);
    double generateRandomOffset(EnergyType parameterIndex,  double energy=0, Traxel tr=0, Traxel tr2=0, size_t state=0);
    double get_transition_probability(Traxel& tr1, Traxel& tr2, size_t state);
    double get_transition_variance(Traxel& tr1, Traxel& tr2);
    const marray::Marray<ValueType>  perturbFactor(const factorType* factor,size_t factorId,std::vector<marray::Marray<ValueType> >* detoffset);
    void write_hypotheses_graph_state(const HypothesesGraph& g, const std::string out_fn);
    void compute_relative_uncertainty(HypothesesGraph *graph);

    // funky export maps
    std::map<std::pair<size_t,size_t>,size_t > clpex_variable_id_map_;
    std::map<std::pair<size_t,std::pair<size_t,size_t> >,size_t> clpex_factor_id_map_;
    
    unsigned int max_number_objects_;

    // energy functions
    boost::function<double (const Traxel&, const size_t)> detection_;
    boost::function<double (const Traxel&, const size_t)> division_;
    boost::function<double (const double)> transition_;

    double forbidden_cost_;
    
    //opengm::LPCplex<pgm::OpengmModelDeprecated::ogmGraphicalModel, pgm::OpengmModelDeprecated::ogmAccumulator>* optimizer_;
	boost::shared_ptr<cplex_optimizer> optimizer_;
    std::vector<IlpSolution> solutions_;

    //factor id maps
    std::map<HypothesesGraph::Node, size_t> detection_f_node_map_;

    double ep_gap_;
    bool with_tracklets_, with_divisions_;
    boost::function<double (const Traxel&)> disappearance_cost_;
    boost::function<double (const Traxel&)> appearance_cost_;

    bool with_misdetections_allowed_;
    bool with_appearance_;
    bool with_disappearance_;
    bool with_constraints_;

    double transition_parameter_;

    UncertaintyParameter uncertainty_param_;
    cplex_optimizer::Parameter cplex_param_;
    ConsTrackingInferenceModel::Parameter inference_model_param_;
    PerturbedInferenceModel::Parameter perturbed_inference_model_param_;

    double cplex_timeout_;
    bool isMAP_;
    
    double division_weight_; // these cannot be read from the division/detection variable since
    double detection_weight_;// those were converted to boost::function objects in tracking
    double transition_weight_;

    boost::python::object transition_classifier_;

    bool with_optical_correction_;

	HypothesesGraph tracklet_graph_;
    std::map<HypothesesGraph::Node, std::vector<HypothesesGraph::Node> > tracklet2traxel_node_map_;
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
  
