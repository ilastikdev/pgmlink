#ifndef TRACKING_FEATURE_EXTRACTOR_H
#define TRACKING_FEATURE_EXTRACTOR_H

#include "pgmlink/hypotheses.h"
#include "pgmlink/higher_order_features.h"

namespace pgmlink {
namespace features {

/**
 * @brief Takes a set of events from a tracking solution and computes a bunch of features for struct learning.
 * @author Carsten Haubold
 */
class TrackingFeatureExtractor
{
public:
    typedef std::vector<double> JointFeatureVector;
    typedef std::vector<std::string> FeatureDescription;

public:
    /// Forbid usage of default constructor
    TrackingFeatureExtractor() = delete;

    /// Create the extractor given a hypotheses graph
    TrackingFeatureExtractor(HypothesesGraph& graph);

    /// Get the complete vector of features computed for the currently set solution
    void get_feature_vector(JointFeatureVector& feature_vector) const;

    /// Return a short description of the feature at the given index in the feature vector
    const std::string get_feature_description(size_t feature_index) const;

    /// Dispatch computation of features here
    void compute_features();

private:
    void push_back_feature(std::string feature_name, double feature_value);
    /**
     * methods that compute each feature
     */
    void compute_velocity_features(ConstTraxelRefVectors&);
    void compute_acceleration_features(ConstTraxelRefVectors&);
    void compute_track_length_features(ConstTraxelRefVectors&);
    void compute_size_difference_features();
    // TODO: add many more

private:
    boost::shared_ptr<TraxelsFeaturesIdentity> position_extractor_ptr_;

    boost::shared_ptr<SquaredDiffCalculator> sq_diff_calc_ptr_;
    boost::shared_ptr<SquaredCurveCalculator> sq_curve_calc_ptr_;

    boost::shared_ptr<MinCalculator<0> > row_min_calc_ptr_;
    boost::shared_ptr<MaxCalculator<0> > row_max_calc_ptr_;
    boost::shared_ptr<SumCalculator<0> > row_sum_calc_ptr_;
    boost::shared_ptr<VarianceCalculator> row_var_calc_ptr_;

    JointFeatureVector joint_feature_vector_;
    FeatureDescription feature_descriptions_;
    HypothesesGraph& graph_;
};

} // end namespace features
} // end namespace pgmlink

#endif // TRACKING_FEATURE_EXTRACTOR_H
