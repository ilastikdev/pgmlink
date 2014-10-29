#ifndef TRACKING_FEATURE_EXTRACTOR_H
#define TRACKING_FEATURE_EXTRACTOR_H

#include "pgmlink/hypotheses.h"
#include "pgmlink/higher_order_features.h"

namespace pgmlink {
namespace features {
class MinMaxMeanVarCalculator
{
public:
    MinMaxMeanVarCalculator();
    void reset();
    void add_value(const double& value);
    void add_values(const FeatureMatrix& values);
    void set_min(const double& value);
    size_t get_count() const;
    double get_mean() const;
    double get_var() const;
    double get_min() const;
    double get_max() const;
    double get_sum_squared_diff() const;

private:
    size_t n;
    double mean;
    double sum_squared_diff;
    double min;
    double max;
};

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
    TrackingFeatureExtractor(boost::shared_ptr<HypothesesGraph> graph,
        const FieldOfView &fov);

    /// Create the extractor given a hypotheses graph with filter function
    TrackingFeatureExtractor(boost::shared_ptr<HypothesesGraph> graph,
        const FieldOfView &fov,
        boost::function<bool (const Traxel&)> margin_filter_function);

    /// Get the complete vector of features computed for the currently set solution
    void get_feature_vector(JointFeatureVector& feature_vector) const;

    /// Return a short description of the feature at the given index in the feature vector
    const std::string get_feature_description(size_t feature_index) const;

    /// Dispatch computation of features here
    void compute_features();

    /// Append features for this solution to the given file.
    /// If file does not exist, create it.
    /// Comments are ignored, and will not be copied to the edited file
    void append_feature_vector_to_file(const std::string &filename);

private:
    void push_back_feature(std::string feature_name, double feature_value);
    /**
     * methods that compute each feature
     */
    void compute_velocity_features(ConstTraxelRefVectors&);
    void compute_acceleration_features(ConstTraxelRefVectors&);
    void compute_angle_features(ConstTraxelRefVectors&);
    void compute_track_outlier_features(ConstTraxelRefVectors&);
    void compute_track_length_features(ConstTraxelRefVectors&);
    void compute_division_move_distance(ConstTraxelRefVectors&);
    void compute_division_move_outlier(ConstTraxelRefVectors&);
    void compute_child_deceleration_features(ConstTraxelRefVectors&);
    void compute_child_deceleration_outlier(ConstTraxelRefVectors&);
    void compute_border_distances(ConstTraxelRefVectors&, std::string);
    void compute_size_difference_features();
    // TODO: add many more

private:
    boost::shared_ptr<TraxelsFeaturesIdentity> position_extractor_ptr_;

    boost::shared_ptr<SquaredDiffCalculator> sq_diff_calc_ptr_;
    boost::shared_ptr<DiffCalculator> diff_calc_ptr_;
    boost::shared_ptr<SquaredCurveCalculator> sq_curve_calc_ptr_;

    boost::shared_ptr<MinCalculator<0> > row_min_calc_ptr_;
    boost::shared_ptr<MaxCalculator<0> > row_max_calc_ptr_;
    boost::shared_ptr<MVNOutlierCalculator> mvn_outlier_calc_ptr_;
    boost::shared_ptr<AngleCosineCalculator> angle_cos_calc_ptr_;
    boost::shared_ptr<ChildParentDiffCalculator> child_parent_diff_calc_ptr_;
    boost::shared_ptr<SquaredNormCalculator<0> > sq_norm_calc_ptr_;
    boost::shared_ptr<ChildDeceleration> child_decel_calc_ptr_;

    FieldOfView fov_;
    boost::function<bool (const Traxel&)> margin_filter_function_;

    JointFeatureVector joint_feature_vector_;
    FeatureDescription feature_descriptions_;
    boost::shared_ptr<HypothesesGraph> graph_;
};

class BorderDistanceFilter {
public:
    BorderDistanceFilter(
        const FieldOfView& field_of_view,
        double t_margin,
        double spatial_margin);
    bool is_out_of_margin(const Traxel& traxel) const;
private:
    FieldOfView fov_;
};

} // end namespace features
} // end namespace pgmlink

#endif // TRACKING_FEATURE_EXTRACTOR_H
