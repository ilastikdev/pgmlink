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

#include <opengm/graphicalmodel/weights.hxx>

#include "pgmlink/event.h"
#include "pgmlink/pgmlink_export.h"
#include "pgmlink/traxels.h"
#include "pgmlink/tracking.h"
#include "pgmlink/field_of_view.h"
#include "pgmlink/merger_resolving.h"
#include "pgmlink/hypotheses.h"
#include <boost/python.hpp>
#include "pgmlink/reasoner_constracking_explicit.h"
#include "pgmlink/inferencemodel/structuredlearningtrackinginferencemodel.h"
#include "pgmlink/structured_learning_tracking_dataset.h"

namespace pgmlink
{

  class StructuredLearningTracking : public ConsExplicitTracking
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
        ConservationExplicitTracking::SolverType solver = ConservationExplicitTracking::CplexSolver,
        int ndim = 2
        ):

        ConsExplicitTracking(
            max_number_objects,
            size_dependent_detection_prob,
            avg_obj_size,
            max_neighbor_distance,
            with_divisions,
            division_threshold,
            random_forest_filename,
            fov,
            event_vector_dump_filename,
            solver,
            ndim),
        numLabels_(max_number_objects),
        numWeights_(5),
        trackingWeights_((size_t)5)

    {
        hypotheses_graph_ = hypotheses_graph;
    }

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

    PGMLINK_EXPORT bool exportCrop(FieldOfView);

    PGMLINK_EXPORT void hypothesesGraphTest( const HypothesesGraph& );
    PGMLINK_EXPORT void addLabels();
    PGMLINK_EXPORT void addAppearanceLabel(int, int, double );
    PGMLINK_EXPORT void addDisappearanceLabel(int, int, double );
    PGMLINK_EXPORT void addDivisionLabel(int, int, double );
    PGMLINK_EXPORT bool addArcLabel(int, int, int, double );
    PGMLINK_EXPORT void addFirstLabels(int, int, double );
    PGMLINK_EXPORT void addLastLabels(int, int, double );
    PGMLINK_EXPORT void addIntermediateLabels(int, int, double );
    PGMLINK_EXPORT virtual boost::shared_ptr<InferenceModel> create_inference_model(
            ConservationExplicitTracking::Parameter& param,
            opengm::learning::Weights<double>& trackingWeights,
            bool withNormalization);
    PGMLINK_EXPORT virtual void prepareTracking(
            ConservationExplicitTracking& pgm,
            ConservationExplicitTracking::Parameter& param,
            opengm::learning::Weights<double>& trackingWeights,
            bool withNormalization,
            bool withClassifierPrior);
    PGMLINK_EXPORT double weight(int);
    PGMLINK_EXPORT void setWeight(int,double);
    PGMLINK_EXPORT void structuredLearning(
            double forbidden_cost = 0,
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
            boost::python::object TransitionClassifier = boost::python::object(),
            bool withNormalization = true,
            bool withClassifierPrior = true);

    PGMLINK_EXPORT void structuredLearningFromParam(ConservationExplicitTracking::Parameter& param);

    PGMLINK_EXPORT ConservationExplicitTracking::Parameter get_structured_learning_tracking_parameters(
            double forbidden_cost = 0,
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
            ConservationExplicitTracking::SolverType solver = ConservationExplicitTracking::CplexSolver,
            bool withNormalization = true,
            bool withClassifierPrior = true);

    PGMLINK_EXPORT void setParameterWeights(ConservationExplicitTracking::Parameter& param,std::vector<double> weights);

public:
    int numCrops_;
    std::vector<FieldOfView> crops_;
    int numWeights_;
    int numLabels_;
    opengm::learning::Weights<double> trackingWeights_;

protected:
    StructuredLearningTrackingInferenceModel::Parameter inference_model_param_;
    double ep_gap_;
    double cplex_timeout_;
    double division_weight_;
    double detection_weight_;
    double transition_weight_;

    boost::function<double (const Traxel&, const Traxel&, const size_t)> transition_;
};

} // end namespace pgmlink

#endif /* STRUCTURED_LEARNING_TRACKING_H */
