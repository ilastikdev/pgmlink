/**
@file
@ingroup tracking
@brief tracking API
 */

#ifndef FUNKEY_BINARY_FILE
#define FUNKEY_BINARY_FILE @FUNKEY_BINARY_FILE@
#endif // FUNKEY_BINARY_FILE

#ifndef STRUCTURED_LEARNING_TRACKING_H
#define STRUCTURED_LEARNING_TRACKING_H

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "pgmlink/event.h"
#include "pgmlink/pgmlink_export.h"
#include "pgmlink/traxels.h"
#include "pgmlink/tracking.h"
#include "pgmlink/field_of_view.h"
#include "pgmlink/merger_resolving.h"
#include "pgmlink/hypotheses.h"
#include <boost/python.hpp>

namespace pgmlink
{

  class StructuredLearningTracking : public ConsTracking
{
public:
    PGMLINK_EXPORT
    StructuredLearningTracking(
                boost::shared_ptr<HypothesesGraph> hypotheses_graph,
                int max_number_objects = 3,
                bool size_dependent_detection_prob = false,
                double avg_obj_size = 30.0,
                double max_neighbor_distance = 20,
                bool with_divisions = true,
                double division_threshold = 0.3,
                const std::string& random_forest_filename = "none",
                FieldOfView fov = FieldOfView(),
                const std::string& event_vector_dump_filename = "none",
                ConservationTracking::SolverType solver = ConservationTracking::CplexSolver
                ):
                hypotheses_graph_ (hypotheses_graph),
                max_number_objects_(max_number_objects),
                use_size_dependent_detection_(size_dependent_detection_prob),
                avg_obj_size_(avg_obj_size),
                max_dist_(max_neighbor_distance),
                with_divisions_(with_divisions),
                division_threshold_(division_threshold),
                detection_rf_fn_(random_forest_filename),
                means_(std::vector<double>()),
                sigmas_(std::vector<double>()),
                fov_(fov),
                event_vector_dump_filename_(event_vector_dump_filename),
                with_optical_correction_(false),
                solver_(solver),
                traxel_store_(nullptr)
    {}


    PGMLINK_EXPORT
    StructuredLearningTracking(boost::shared_ptr<HypothesesGraph> g,
                 TraxelStore& ts,
                 ConservationTracking::Parameter param,
                 UncertaintyParameter uncertainty_param,
                 FieldOfView fov = FieldOfView(),
                 bool size_dependent_detection_prob = false,
                 double avg_obj_size = 30.0,
                 double max_neighbor_distance = 20,
                 double division_threshold = 0.3):
        max_number_objects_(param.max_number_objects),
        max_dist_(max_neighbor_distance),
        with_divisions_(param.with_divisions),
        division_threshold_(division_threshold),
        detection_rf_fn_("none"),
        use_size_dependent_detection_(size_dependent_detection_prob),
        avg_obj_size_(avg_obj_size),
        means_(std::vector<double>()),
        sigmas_(std::vector<double>()),
        fov_(fov),
        event_vector_dump_filename_("none"),
        with_optical_correction_(false),
        solver_(param.solver_),
        traxel_store_(&ts),
        hypotheses_graph_(g),
        uncertainty_param_(uncertainty_param)
    {}




    PGMLINK_EXPORT EventVectorVectorVector operator()(
        TraxelStore& ts,
        double forbidden_cost = 0,
        double ep_gap = 0.01,
        bool with_tracklets = true,
        double division_weight = 10.0,
        double transition_weight = 10.0,
        double disappearance_cost = 0,
        double appearance_cost = 0,
        bool with_merger_resolution = true,
        int n_dim = 3,
        double transition_parameter = 5.,
        double border_width = 0,
        bool with_constraints = true,
        UncertaintyParameter uncertaintyParam = UncertaintyParameter(),
        double cplex_timeout = 1e+75,
        TimestepIdCoordinateMapPtr coordinates = TimestepIdCoordinateMapPtr(),
        boost::python::object transition_classifier = boost::python::object());

    PGMLINK_EXPORT bool exportCrop(FieldOfView);//, const std::string& );

    PGMLINK_EXPORT void hypothesesGraphTest( const HypothesesGraph& );
    PGMLINK_EXPORT void addLabels( HypothesesGraph& );
    PGMLINK_EXPORT void addFirstLabels( HypothesesGraph&, int, int, double );
    PGMLINK_EXPORT void addLastLabels( HypothesesGraph&, int, int, double );
    PGMLINK_EXPORT void addSingletonLabels( HypothesesGraph&, int, int, double );
    PGMLINK_EXPORT void addIntermediateLabels( HypothesesGraph&, int, int, double );




    /**
     * refactoring of operator().
     */





    //PGMLINK_EXPORT boost::shared_ptr<HypothesesGraph> build_hypo_graph(TraxelStore& ts);

    /*
    PGMLINK_EXPORT boost::shared_ptr<HypothesesGraph> get_hypo_graph();
    PGMLINK_EXPORT boost::shared_ptr<HypothesesGraph> get_resolved_hypotheses_graph();





    PGMLINK_EXPORT virtual EventVectorVectorVector track(double forbidden_cost = 0,
            double ep_gap = 0.01,
            bool with_tracklets = true,
            double detection_weight = 10.,
            double division_weight = 10.0,
            double transition_weight = 10.0,
            double disappearance_cost = 0,
            double appearance_cost = 0,
            bool with_merger_resolution = true,
            int n_dim = 3,
            double transition_parameter = 5.,
            double border_width = 0,
            bool with_constraints = true,
            UncertaintyParameter uncertaintyParam = UncertaintyParameter(),
            double cplex_timeout = 1e+75,
            boost::python::object TransitionClassifier = boost::python::object());

    PGMLINK_EXPORT EventVectorVectorVector track_from_param(ConservationTracking::Parameter& param,
                                                            bool fixLabeledNodes = false);

    PGMLINK_EXPORT ConservationTracking::Parameter get_conservation_tracking_parameters(double forbidden_cost = 0,
            double ep_gap = 0.01,
            bool with_tracklets = true,
            double detection_weight = 10.,
            double division_weight = 10.0,
            double transition_weight = 10.0,
            double disappearance_cost = 0,
            double appearance_cost = 0,
            bool with_merger_resolution = true,
            int n_dim = 3,
            double transition_parameter = 5.,
            double border_width = 0,
            bool with_constraints = true,
            UncertaintyParameter uncertaintyParam = UncertaintyParameter(),
            double cplex_timeout = 1e+75,
            boost::python::object transition_classifier = boost::python::object(),
            ConservationTracking::SolverType solver = ConservationTracking::CplexSolver);

    PGMLINK_EXPORT void setTrackLabelingExportFile(std::string file_name);

    PGMLINK_EXPORT void setParameterWeights(ConservationTracking::Parameter& param,std::vector<double> weights);

    PGMLINK_EXPORT EventVectorVector resolve_mergers(EventVectorVector &events,
            TimestepIdCoordinateMapPtr coordinates = TimestepIdCoordinateMapPtr(),
            double ep_gap = 0.01,
            double transition_weight = 10.0,
            bool with_tracklets = true,
            int n_dim = 3,
            double transition_parameter = 5.,
            bool with_constraints = true,
            boost::python::object transitionClassifier = boost::python::object());
*/

    /**
     * Get state of detection variables after call to operator().
     */


    /*
    PGMLINK_EXPORT std::vector< std::map<unsigned int, bool> > detections();
    PGMLINK_EXPORT void createStructuredLearningFiles(std::string feature_file_name,
                                      std::string constraints_file_name,
                                      std::string ground_truth_file_name);
    PGMLINK_EXPORT void writeStructuredLearningFiles(std::string feature_file_name,
                                      std::string constraints_file_name,
                                      std::string ground_truth_file_name,
                                      ConservationTracking::Parameter param);
    PGMLINK_EXPORT vector<double> learnTrackingWeights(std::string feature_file_name,
                                      std::string constraints_file_name,
                                      std::string ground_truth_file_name,
                                      std::string lossweights = "",
                                      std::string options = "");
    PGMLINK_EXPORT double hammingloss_of_files(std::string f1, std::string f2);

    /// Return reference to the ilp solutions
    PGMLINK_EXPORT void save_ilp_solutions(const std::string& filename);
*/
protected:
    int numCrops_;
    std::vector<FieldOfView> crops_;
    int max_number_objects_;
    double max_dist_;
    bool with_divisions_;
    double division_threshold_;
    const std::string detection_rf_fn_;
    bool use_size_dependent_detection_;
    bool use_classifier_prior_;
    double avg_obj_size_;
    std::vector<double> means_, sigmas_;
//    boost::shared_ptr<std::vector< std::map<unsigned int, bool> > > last_detections_;
    FieldOfView fov_;
    std::string event_vector_dump_filename_;
//    std::vector<ConservationTracking::IlpSolution> ilp_solutions_;

//    std::string tracking_labels_export_file_name_;

    TraxelStore* traxel_store_;

    boost::shared_ptr<HypothesesGraph> hypotheses_graph_;
//    boost::shared_ptr<HypothesesGraph> original_hypotheses_graph_;
//    boost::shared_ptr<HypothesesGraph> resolved_graph_;
//    boost::shared_ptr<ConservationTracking> pgm_;

    bool with_optical_correction_;
    UncertaintyParameter uncertainty_param_;
    ConservationTracking::SolverType solver_;
};

} // end namespace pgmlink

#endif /* STRUCTURED_LEARNING_TRACKING_H */
