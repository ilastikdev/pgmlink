#ifndef TRACKING_FEATURE_EXTRACTOR_H
#define TRACKING_FEATURE_EXTRACTOR_H

#include "pgmlink/hypotheses.h"
#include "pgmlink/features/higher_order_features.h"

namespace pgmlink
{
namespace features
{
class MinMaxMeanVarCalculator
{
public:
    MinMaxMeanVarCalculator();
    void reset();
    void add_value(const double& value);
    void add_values(const FeatureMatrix& v);
    size_t get_count() const;
    double get_mean() const;
    double get_var() const;
    double get_min() const;
    double get_max() const;

private:
    std::vector<double> values;
};

class TrackFeatureExtractor
/* TODO:
    - add track length as a feature as well?
*/
{
public:
    typedef std::vector<double> JointFeatureVector;
    typedef std::vector<std::string> FeatureDescription;
public:
    TrackFeatureExtractor();
    void compute_features(
        ConstTraxelRefVector& traxelref_vec,
        FeatureVectorView return_vector);
    void compute_features(
        ConstTraxelRefVectors& traxelref_vecs,
        FeatureMatrix& return_matrix);
    size_t get_feature_vector_length() const;
    void get_feature_descriptions(FeatureDescription& feature_descriptions) const;
private:
    void compute_sq_id_features(ConstTraxelRefVector&, std::string);
    void compute_sq_diff_features(ConstTraxelRefVector&, std::string);
    void compute_sq_curve_features(ConstTraxelRefVector&, std::string);
    void compute_angle_features(ConstTraxelRefVector&, std::string);

    void push_back_feature(std::string feature_name, double feature_value);
    void push_back_feature(
        std::string feature_name,
        const MinMaxMeanVarCalculator& mmmv_calculator);
private:
    FeatureVectorView::iterator feature_vector_offset_it_;
    FeatureDescription::iterator feature_descriptions_offset_it_;
    FeatureDescription feature_descriptions_;
};

class DivisionFeatureExtractor
{
public:
    typedef std::vector<double> JointFeatureVector;
    typedef std::vector<std::string> FeatureDescription;
public:
    DivisionFeatureExtractor();
    void compute_features(
        ConstTraxelRefVector& traxelref_vec,
        FeatureVectorView return_vector);
    void compute_features(
        ConstTraxelRefVectors& traxelref_vecs,
        FeatureMatrix& return_matrix);
    size_t get_feature_vector_length() const;
    void get_feature_descriptions(FeatureDescription& feature_descriptions) const;
private:
    void compute_id_features(ConstTraxelRefVector&, std::string);
    void compute_sq_diff_features(ConstTraxelRefVector&, std::string);
    void compute_angle_features(ConstTraxelRefVector&, std::string);

    void push_back_feature(std::string feature_name, double feature_value);
private:
    FeatureVectorView::iterator feature_vector_offset_it_;
    FeatureDescription::iterator feature_descriptions_offset_it_;
    FeatureDescription feature_descriptions_;
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

    /// Train the svms with the given hypotheses graph assuming that the arc and
    /// node active map contain a correct labeling
    void train_track_svm();
    void train_division_svm();

    /// set and get for the outlier svms
    boost::shared_ptr<SVMOutlierCalculator> get_track_svm() const;
    boost::shared_ptr<SVMOutlierCalculator> get_division_svm() const;
    void set_track_svm(boost::shared_ptr<SVMOutlierCalculator> track_svm);
    void set_division_svm(boost::shared_ptr<SVMOutlierCalculator> division_svm);

    /// Get the complete vector of features computed for the currently set solution
    void get_feature_vector(JointFeatureVector& feature_vector) const;

    /// Return a short description of the feature at the given index in the feature vector
    const std::string get_feature_description(size_t feature_index) const;

    /// Dispatch computation of features here
    void compute_features();

    /// Set HDF5 filename to which the features of all tracks will be written
    void set_track_feature_output_file(const std::string& filename);

    /// Append features for this solution to the given file.
    /// If file does not exist, create it.
    /// Comments are ignored, and will not be copied to the edited file
    void append_feature_vector_to_file(const std::string &filename);

private:
    void push_back_feature(std::string feature_name, double feature_value);
    void push_back_feature(
        std::string feature_name,
        const MinMaxMeanVarCalculator& mmmv_calculator);

    void save_features_to_h5(size_t track_id, const std::string& feature_name, FeatureMatrix &matrix, bool tracks = true);
    void save_traxel_ids_to_h5(ConstTraxelRefVectors& track_traxels);
    void save_division_traxels_to_h5(ConstTraxelRefVectors &division_traxels);

    /**
     * methods that compute each feature
     */
    void compute_sq_diff_features(ConstTraxelRefVectors&, std::string);
    void compute_sq_accel_features(ConstTraxelRefVectors&, std::string);
    void compute_angle_features(ConstTraxelRefVectors&, std::string);
    void compute_track_length_features(ConstTraxelRefVectors&);
    void compute_track_id_outlier(ConstTraxelRefVectors&, std::string);
    void compute_track_diff_outlier(ConstTraxelRefVectors&, std::string);
    void compute_track_feature_outlier(ConstTraxelRefVectors&);
    void compute_division_sq_diff_features(ConstTraxelRefVectors&, std::string);
    void compute_division_sq_diff_outlier(ConstTraxelRefVectors&, std::string);
    void compute_child_deceleration_features(ConstTraxelRefVectors&, std::string);
    void compute_child_deceleration_outlier(ConstTraxelRefVectors&, std::string);
    void compute_division_feature_outlier(ConstTraxelRefVectors&);
    void compute_border_distances(ConstTraxelRefVectors&, std::string);
    // TODO: add many more

private:
    boost::shared_ptr<SquaredDiffCalculator> sq_diff_calc_ptr_;
    boost::shared_ptr<DiffCalculator> diff_calc_ptr_;
    boost::shared_ptr<SquaredCurveCalculator> sq_curve_calc_ptr_;

    boost::shared_ptr<MinCalculator<0> > row_min_calc_ptr_;
    boost::shared_ptr<MaxCalculator<0> > row_max_calc_ptr_;
    boost::shared_ptr<MVNOutlierCalculator> mvn_outlier_calc_ptr_;
    boost::shared_ptr<SVMOutlierCalculator> svm_track_outlier_calc_ptr_;
    boost::shared_ptr<SVMOutlierCalculator> svm_div_outlier_calc_ptr_;
    boost::shared_ptr<SquaredMahalanobisCalculator> sq_mahal_calc_ptr_;
    boost::shared_ptr<AngleCosineCalculator> angle_cos_calc_ptr_;
    boost::shared_ptr<ChildParentDiffCalculator> child_parent_diff_calc_ptr_;
    boost::shared_ptr<SquaredNormCalculator<0> > sq_norm_calc_ptr_;
    boost::shared_ptr<ChildDeceleration> child_decel_calc_ptr_;

    FieldOfView fov_;
    boost::function<bool (const Traxel&)> margin_filter_function_;

    JointFeatureVector joint_feature_vector_;
    FeatureDescription feature_descriptions_;
    boost::shared_ptr<HypothesesGraph> graph_;
    std::string track_feature_output_file_;
};

class BorderDistanceFilter
{
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
