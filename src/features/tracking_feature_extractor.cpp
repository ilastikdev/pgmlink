#include "pgmlink/features/tracking_feature_extractor.h"
#include <boost/algorithm/string.hpp>

#include <cmath> /* for std::abs */

namespace pgmlink
{
namespace features
{

MinMaxMeanVarCalculator::MinMaxMeanVarCalculator()
{
    reset();
}

void MinMaxMeanVarCalculator::reset()
{
    values.clear();
}

void MinMaxMeanVarCalculator::add_value(const double& value)
{
    values.push_back(value);
}

void MinMaxMeanVarCalculator::add_values(const FeatureMatrix& v)
{
    values.insert(values.end(), v.begin(), v.end());
}

size_t MinMaxMeanVarCalculator::get_count() const
{
    return values.size();
}

double MinMaxMeanVarCalculator::get_mean() const
{
    if(values.size() == 0)
        return 0.0;

    double sum = 0.0;
    for(auto v : values)
        sum += v;
    return sum / values.size();
}

double MinMaxMeanVarCalculator::get_var() const
{
    if (values.size() > 0)
    {
        double mean = get_mean();
        double sum_squared_diff = 0.0;
        for(auto v : values)
            sum_squared_diff += pow(v - mean, 2.0);
        return sum_squared_diff / values.size();
    }
    else
    {
        return 0.0;
    }
}

double MinMaxMeanVarCalculator::get_min() const
{
    if(values.size() == 0)
//        return std::numeric_limits<double>::max();
        return 0.0;
    return *std::min_element(values.begin(), values.end());
}

double MinMaxMeanVarCalculator::get_max() const
{
    if(values.size() == 0)
//        return std::numeric_limits<double>::lowest();
        return 0.0;
    return *std::max_element(values.begin(), values.end());
}

void TrackingFeatureExtractorBase::init_compute(FeatureVectorView return_vector)
{
    assert(return_vector.shape(0) == get_feature_vector_length());
    if (feature_descriptions_.size() == 0)
    {
        feature_descriptions_.resize(get_feature_vector_length());
    }
    feature_descriptions_it_ = feature_descriptions_.begin();
    feature_vector_it_ = return_vector.begin();
}

void TrackingFeatureExtractorBase::get_feature_descriptions(
    FeatureDescriptionVector& feature_descriptions) const
{
    feature_descriptions.clear();
    feature_descriptions.insert(
        feature_descriptions.begin(),
        feature_descriptions_.begin(),
        feature_descriptions_.end());
}

void TrackingFeatureExtractorBase::push_back_feature(
    std::string feature_name,
    double feature_value)
{
    LOG(logDEBUG4) << "Push back feature " << feature_name << ": " << feature_value;
    *feature_vector_it_ = feature_value;
    feature_vector_it_++;
    *feature_descriptions_it_ = feature_name;
    feature_descriptions_it_++;
}

void TrackingFeatureExtractorBase::push_back_feature(
    std::string feature_name,
    const MinMaxMeanVarCalculator& mmmv_calculator)
{
//    push_back_feature("min of " + feature_name, mmmv_calculator.get_min());
//    push_back_feature("max of " + feature_name, mmmv_calculator.get_max());
    push_back_feature("mean of " + feature_name, mmmv_calculator.get_mean());
    push_back_feature("var of " + feature_name, mmmv_calculator.get_var());
}

TrackFeatureExtractor::TrackFeatureExtractor()
{}

size_t TrackFeatureExtractor::get_feature_vector_length() const
{
    // TODO: UUUUGLY!
    return 3 * 2 // compute_sq_id_features (Count, Mean, Variance)
           + 4 * 2 // compute_sq_diff_features (com, Count, Mean, Variance)
           + 4 * 2 // compute_sq_accel_features (com, Count, Mean, Variance)
           + 1 * 2 // compute_angle_features (com)
           + 1; // length
}

std::ostream& operator<<(std::ostream& lhs, TrackFeatureExtractor::FeatureDescriptionVector& rhs)
{
    for(auto s: rhs)
    {
        lhs << s << ", ";
    }
    return lhs;
}

void TrackFeatureExtractor::compute_features(
    ConstTraxelRefVector& traxelref_vec,
    FeatureVectorView return_vector)
{
    init_compute(return_vector);
    compute_sq_id_features(traxelref_vec, "Count");
    compute_sq_id_features(traxelref_vec, "Mean");
    compute_sq_id_features(traxelref_vec, "Variance");
    compute_sq_diff_features(traxelref_vec, "RegionCenter");
    compute_sq_diff_features(traxelref_vec, "Count");
    compute_sq_diff_features(traxelref_vec, "Mean");
    compute_sq_diff_features(traxelref_vec, "Variance");
    compute_sq_curve_features(traxelref_vec, "RegionCenter");
    compute_sq_curve_features(traxelref_vec, "Count");
    compute_sq_curve_features(traxelref_vec, "Mean");
    compute_sq_curve_features(traxelref_vec, "Variance");
    compute_angle_features(traxelref_vec, "RegionCenter");
    push_back_feature("Length", (double)traxelref_vec.size());
}

void TrackFeatureExtractor::compute_features(
    ConstTraxelRefVectors& traxelref_vecs,
    FeatureMatrix& return_matrix)
{
    size_t f_dim = get_feature_vector_length();
    return_matrix.reshape(vigra::Shape2(traxelref_vecs.size(), f_dim));
    for (size_t i = 0; i < traxelref_vecs.size(); i++)
    {
        compute_features(traxelref_vecs[i], return_matrix.bind<0>(i));
    }
}

void TrackFeatureExtractor::compute_sq_id_features(
    ConstTraxelRefVector& traxelref_vec,
    std::string feature_name)
{
    MinMaxMeanVarCalculator sq_id_mmmv;
//    sq_id_mmmv.set_max(0.0);

    if (traxelref_vec.size() >= 1)
    {
        // get the features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(traxelref_vec, feature_matrix);
        // calculate the squared norms
        SquaredNormCalculator<0> sq_id_calculator;
        FeatureMatrix sq_id_matrix;
        sq_id_calculator.calculate(feature_matrix, sq_id_matrix);

        sq_id_mmmv.add_values(sq_id_matrix);
    }
    push_back_feature("sq_id of " + feature_name, sq_id_mmmv);
}

void TrackFeatureExtractor::compute_sq_diff_features(
    ConstTraxelRefVector& traxelref_vec,
    std::string feature_name)
{
    MinMaxMeanVarCalculator sq_diff_mmmv;
//    sq_diff_mmmv.set_max(0.0);

    if (traxelref_vec.size() > 1)
    {
        // get the features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(traxelref_vec, feature_matrix);
        // calculate the squared differences
        SquaredDiffCalculator sq_diff_calculator;
        FeatureMatrix sq_diff_matrix;
        sq_diff_calculator.calculate(feature_matrix, sq_diff_matrix);

        sq_diff_mmmv.add_values(sq_diff_matrix);
    }
    push_back_feature("sq_diff of " + feature_name, sq_diff_mmmv);
}

void TrackFeatureExtractor::compute_sq_curve_features(
    ConstTraxelRefVector& traxelref_vec,
    std::string feature_name)
{
    MinMaxMeanVarCalculator sq_curve_mmmv;
//    sq_curve_mmmv.set_max(0.0);

    if (traxelref_vec.size() > 2)
    {
        // get the features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(traxelref_vec, feature_matrix);
        // calculate the squared curvatures
        SquaredCurveCalculator sq_curve_calculator;
        FeatureMatrix sq_curve_matrix;
        sq_curve_calculator.calculate(feature_matrix, sq_curve_matrix);

        sq_curve_mmmv.add_values(sq_curve_matrix);
    }
    push_back_feature("sq_curve of " + feature_name, sq_curve_mmmv);
}

void TrackFeatureExtractor::compute_angle_features(
    ConstTraxelRefVector& traxelref_vec,
    std::string feature_name)
{
    MinMaxMeanVarCalculator angle_mmmv;

    if (traxelref_vec.size() > 2)
    {
        // get the feature matrix
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(traxelref_vec, feature_matrix);
        // calculate the angle cosines
        AngleCosineCalculator angle_calculator;
        FeatureMatrix angle_matrix;
        angle_calculator.calculate(feature_matrix, angle_matrix);

        angle_mmmv.add_values(angle_matrix);
    }
    push_back_feature("mean of angle cosines of " + feature_name, angle_mmmv.get_mean());
    push_back_feature("var of angle cosines of " + feature_name, angle_mmmv.get_var());
}

//-----------------------------------------------------------------------------
// DivisionFeatureExtractor
//-----------------------------------------------------------------------------
DivisionFeatureExtractor::DivisionFeatureExtractor()
{}

size_t DivisionFeatureExtractor::get_feature_vector_length() const
{
    // TODO: still ugly
    return 3 * 3 // compute_sq_id_features (Count, Mean, Variance)
           + 2 * 4 // compute_sq_diff_features (com, Count, Mean, Variance)
           + 1 * 1; // compute_angle_features (com)
}

void DivisionFeatureExtractor::compute_features(
    ConstTraxelRefVector& traxelref_vec,
    FeatureVectorView return_vector)
{
    init_compute(return_vector);
    compute_id_features(traxelref_vec, "Count");
    compute_id_features(traxelref_vec, "Mean");
    compute_id_features(traxelref_vec, "Variance");
    compute_sq_diff_features(traxelref_vec, "RegionCenter");
    compute_sq_diff_features(traxelref_vec, "Count");
    compute_sq_diff_features(traxelref_vec, "Mean");
    compute_sq_diff_features(traxelref_vec, "Variance");
    compute_angle_features(traxelref_vec, "RegionCenter");
}

void DivisionFeatureExtractor::compute_features(
    ConstTraxelRefVectors& traxelref_vecs,
    FeatureMatrix& return_matrix)
{
    size_t f_dim = get_feature_vector_length();
    return_matrix.reshape(vigra::Shape2(traxelref_vecs.size(), f_dim));
    for (size_t i = 0; i < traxelref_vecs.size(); i++)
    {
        // TODO implement for divisions of higher order
        assert(traxelref_vecs[i].size() == 3);
        compute_features(traxelref_vecs[i], return_matrix.bind<0>(i));
    }
}

void DivisionFeatureExtractor::compute_id_features(
    ConstTraxelRefVector& traxelref_vec,
    std::string feature_name)
{
    // get the features
    FeatureMatrix feature_matrix;
    TraxelsFeaturesIdentity feature_extractor(feature_name);
    feature_extractor.extract(traxelref_vec, feature_matrix);
    assert(feature_matrix.shape(0) == 3);
    // the feature matrix does not have to have a feature dimension of one
    // but it's asserted anyway to indicate possible causes of bugs
    assert(feature_matrix.shape(1) == 1);
    double child_feature_sum = feature_matrix(1, 0) + feature_matrix(2, 0);
    double child_feature_diff = feature_matrix(1, 0) - feature_matrix(2, 0);
    child_feature_diff = std::abs(child_feature_diff);

    push_back_feature("parent " + feature_name, feature_matrix(0, 0));
    push_back_feature("child sum " + feature_name, child_feature_sum);
    push_back_feature("child diff " + feature_name, child_feature_diff);
}

void DivisionFeatureExtractor::compute_sq_diff_features(
    ConstTraxelRefVector& traxelref_vec,
    std::string feature_name)
{
    // get the features
    FeatureMatrix feature_matrix;
    TraxelsFeaturesIdentity feature_extractor(feature_name);
    feature_extractor.extract(traxelref_vec, feature_matrix);
    // calculate the squared differences
    TCompositionCalculator <
    ChildParentDiffCalculator,
    SquaredNormCalculator<0>
    > sq_diff_calculator;
    FeatureMatrix sq_diff_matrix;
    sq_diff_calculator.calculate(feature_matrix, sq_diff_matrix);
    assert(sq_diff_matrix.shape(0) == 2);
    assert(sq_diff_matrix.shape(1) == 1);

    double feature_sum = sq_diff_matrix(0, 0) + sq_diff_matrix(1, 0);
    double feature_diff = sq_diff_matrix(0, 0) - sq_diff_matrix(1, 0);
    feature_diff = std::abs(feature_diff);

    push_back_feature(
        "sum of squared child-parent diff of " + feature_name,
        feature_sum);
    push_back_feature(
        "diff of squared child-parent diff of " + feature_name,
        feature_diff);
}

void DivisionFeatureExtractor::compute_angle_features(
    ConstTraxelRefVector& traxelref_vec,
    std::string feature_name)
{
    // get the feature matrix
    FeatureMatrix feature_matrix;
    TraxelsFeaturesIdentity feature_extractor(feature_name);
    feature_extractor.extract(traxelref_vec, feature_matrix);
    // calculate the angle cosines (dot product)
    ChildParentDiffCalculator diff_calculator;
    FeatureMatrix diff_matrix;
    diff_calculator.calculate(feature_matrix, diff_matrix);
    assert(diff_matrix.shape(0) == 2);
    EuclideanNormCalculator norm_calculator;
    FeatureMatrix norm_matrix;
    norm_calculator.calculate(diff_matrix, norm_matrix);
    assert(norm_matrix.shape(0) == 2);
    assert(norm_matrix.shape(1) == 1);

    DotProductCalculator dot_calculator;
    FeatureMatrix dot_matrix;
    dot_calculator.calculate(diff_matrix, dot_matrix);
    assert(dot_matrix.shape(0) == 1);
    assert(dot_matrix.shape(1) == 1);

    push_back_feature(
        "angle cosine of " + feature_name,
        dot_matrix(0, 0) / (norm_matrix(0, 0) * norm_matrix(1, 0)));
}

//-----------------------------------------------------------------------------
// TrackingFeatureExtractor
//-----------------------------------------------------------------------------
TrackingFeatureExtractor::TrackingFeatureExtractor(boost::shared_ptr<HypothesesGraph> graph,
        const FieldOfView& fov,
        boost::function<bool (const Traxel&)> margin_filter_function):
    graph_(graph),
    fov_(fov),
    margin_filter_function_(margin_filter_function),
    sq_diff_calc_ptr_(new SquaredDiffCalculator),
    diff_calc_ptr_(new DiffCalculator),
    sq_curve_calc_ptr_(new SquaredCurveCalculator),
    row_min_calc_ptr_(new MinCalculator<0>),
    row_max_calc_ptr_(new MaxCalculator<0>),
    mvn_outlier_calc_ptr_(new MVNOutlierCalculator),
    svm_track_outlier_calc_ptr_(new SVMOutlierCalculator),
    svm_div_outlier_calc_ptr_(new SVMOutlierCalculator),
    sq_mahal_calc_ptr_(new SquaredMahalanobisCalculator),
    angle_cos_calc_ptr_(new AngleCosineCalculator),
    child_parent_diff_calc_ptr_(new ChildParentDiffCalculator),
    sq_norm_calc_ptr_(new SquaredNormCalculator<0>),
    child_decel_calc_ptr_(new ChildDeceleration)
{}

TrackingFeatureExtractor::TrackingFeatureExtractor(boost::shared_ptr<HypothesesGraph> graph,
        const FieldOfView& fov):
    graph_(graph),
    fov_(fov),
    margin_filter_function_(NULL),
    sq_diff_calc_ptr_(new SquaredDiffCalculator),
    diff_calc_ptr_(new DiffCalculator),
    sq_curve_calc_ptr_(new SquaredCurveCalculator),
    row_min_calc_ptr_(new MinCalculator<0>),
    row_max_calc_ptr_(new MaxCalculator<0>),
    mvn_outlier_calc_ptr_(new MVNOutlierCalculator),
    svm_track_outlier_calc_ptr_(new SVMOutlierCalculator),
    svm_div_outlier_calc_ptr_(new SVMOutlierCalculator),
    sq_mahal_calc_ptr_(new SquaredMahalanobisCalculator),
    angle_cos_calc_ptr_(new AngleCosineCalculator),
    child_parent_diff_calc_ptr_(new ChildParentDiffCalculator),
    sq_norm_calc_ptr_(new SquaredNormCalculator<0>),
    child_decel_calc_ptr_(new ChildDeceleration)
{}

void TrackingFeatureExtractor::get_feature_vector(TrackingFeatureExtractor::JointFeatureVector &feature_vector) const
{
    feature_vector.clear();
    feature_vector.insert(feature_vector.begin(), joint_feature_vector_.begin(), joint_feature_vector_.end());
}

void TrackingFeatureExtractor::train_track_svm()
{
    // get the track traxels
    // TODO train only with complete tracks? If yes it should also only be
    // applied to complete tracks, see comment in compute_features()
    TrackTraxels track_extractor(true, true);
    ConstTraxelRefVectors track_traxels = track_extractor(*graph_);
    // calculate the feature vectors for all tracks
    FeatureMatrix track_feature_matrix;
    TrackFeatureExtractor track_feature_extractor;
//    TrackFeatureExtractor::FeatureDescription fd;
//    track_feature_extractor.get_feature_descriptions(fd);
//    LOG(logINFO) << "Training SVM from features: \n" << fd;
    track_feature_extractor.compute_features(track_traxels, track_feature_matrix);
    // train the svm (feature normalization is performed internally)
    svm_track_outlier_calc_ptr_->train(track_feature_matrix, sqrt(track_feature_matrix.shape(1)));
}

void TrackingFeatureExtractor::train_division_svm()
{
    // get the division traxels
    DivisionTraxels div_extractor;
    ConstTraxelRefVectors div_traxels = div_extractor(*graph_);
    // calculate the feature vectors for all divisions
    FeatureMatrix div_feature_matrix;
    DivisionFeatureExtractor div_feature_extractor;
    div_feature_extractor.compute_features(div_traxels, div_feature_matrix);
    // train the svm (feature normalization is performed internally)
    svm_div_outlier_calc_ptr_->train(div_feature_matrix, sqrt(div_feature_matrix.shape(1)));
}

boost::shared_ptr<SVMOutlierCalculator> TrackingFeatureExtractor::get_track_svm() const
{
    return svm_track_outlier_calc_ptr_;
}

boost::shared_ptr<SVMOutlierCalculator> TrackingFeatureExtractor::get_division_svm() const
{
    return svm_div_outlier_calc_ptr_;
}

void TrackingFeatureExtractor::set_track_svm(
    boost::shared_ptr<SVMOutlierCalculator> track_svm)
{
    svm_track_outlier_calc_ptr_ = track_svm;
}

void TrackingFeatureExtractor::set_division_svm(
    boost::shared_ptr<SVMOutlierCalculator> division_svm)
{
    svm_div_outlier_calc_ptr_ = division_svm;
}

const std::string TrackingFeatureExtractor::get_feature_description(size_t feature_index) const
{
    return feature_descriptions_[feature_index];
}

void TrackingFeatureExtractor::save_features_to_h5(size_t track_id, const std::string& feature_name, FeatureMatrix &matrix, bool tracks)
{
    if(track_feature_output_file_.size() > 0)
    {
        std::stringstream dataset_name;
        if(tracks)
        {
            dataset_name << "tracks/";
        }
        else
        {
            dataset_name << "divisions/";
        }
        dataset_name << track_id << "/" << feature_name;
        vigra::writeHDF5(track_feature_output_file_.c_str(), dataset_name.str().c_str(), matrix);
    }
}

void TrackingFeatureExtractor::save_traxel_ids_to_h5(ConstTraxelRefVectors &track_traxels)
{
    size_t track_id = 0;
    for(auto track : track_traxels)
    {
        if(track.size() == 0)
        {
            LOG(logWARNING) << "Cannot output empty track!";
            continue;
        }
        FeatureMatrix traxels(vigra::Shape2(track.size(), 2));

        for(size_t i = 0; i < track.size(); ++i)
        {
            traxels(i, 0) = track[i]->Timestep;
            traxels(i, 1) = track[i]->Id;
        }
        save_features_to_h5(track_id++, "traxels", traxels);
    }
}

void TrackingFeatureExtractor::save_division_traxels_to_h5(ConstTraxelRefVectors &division_traxels)
{
    size_t division_id = 0;
    for(auto div : division_traxels)
    {
        if(div.size() < 3)
        {
            LOG(logWARNING) << "Cannot output empty/invalid division!";
            continue;
        }
        FeatureMatrix traxels(vigra::Shape2(div.size(), 2));

        for(size_t i = 0; i < div.size(); ++i)
        {
            traxels(i, 0) = div[i]->Timestep;
            traxels(i, 1) = div[i]->Id;
        }
        save_features_to_h5(division_id++, "traxels", traxels, false);
    }
}

void TrackingFeatureExtractor::compute_all_track_features()
{
    LOG(logDEBUG) << "Extract all tracks";
    TrackTraxels track_extractor;
    ConstTraxelRefVectors track_traxels = track_extractor(*graph_);
    compute_sq_diff_features(track_traxels, "RegionCenter");
    compute_sq_diff_features(track_traxels, "Count");
    compute_sq_diff_features(track_traxels, "Mean");
    compute_sq_diff_features(track_traxels, "Variance");
    compute_sq_accel_features(track_traxels, "RegionCenter");
    compute_sq_accel_features(track_traxels, "Count");
    compute_sq_accel_features(track_traxels, "Mean");
    compute_sq_accel_features(track_traxels, "Variance");
    compute_angle_features(track_traxels, "RegionCenter");
    compute_track_length_features(track_traxels);
    compute_track_id_outlier(track_traxels, "RegionCenter");
    compute_track_id_outlier(track_traxels, "Count");
    compute_track_id_outlier(track_traxels, "Mean");
    compute_track_id_outlier(track_traxels, "Variance");
    compute_track_diff_outlier(track_traxels, "RegionCenter");
    compute_track_diff_outlier(track_traxels, "Count");
    compute_track_diff_outlier(track_traxels, "Mean");
    compute_track_diff_outlier(track_traxels, "Variance");
    //TODO filter the tracks for the following? (division start / division end)
    compute_svm_track_feature_outlier(track_traxels);

    save_traxel_ids_to_h5(track_traxels);
}

void TrackingFeatureExtractor::compute_all_division_features()
{
    LOG(logDEBUG) << "Extract all divisions to depth 1";
    DivisionTraxels div_1_extractor(1);
    ConstTraxelRefVectors div_1_traxels = div_1_extractor(*graph_);
    LOG(logDEBUG) << "Extract all divisions to depth 2";
    DivisionTraxels div_2_extractor(2);
    ConstTraxelRefVectors div_2_traxels = div_2_extractor(*graph_);
    compute_division_sq_diff_features(div_1_traxels, "RegionCenter");
    compute_division_sq_diff_features(div_1_traxels, "Count");
    compute_division_sq_diff_features(div_1_traxels, "Mean");
    compute_division_sq_diff_features(div_1_traxels, "Variance");
    compute_division_sq_diff_outlier(div_1_traxels, "RegionCenter");
    compute_division_sq_diff_outlier(div_1_traxels, "Count");
    compute_division_sq_diff_outlier(div_1_traxels, "Mean");
    compute_division_sq_diff_outlier(div_1_traxels, "Variance");
    compute_child_deceleration_features(div_2_traxels, "RegionCenter");
    compute_child_deceleration_features(div_2_traxels, "Count");
    compute_child_deceleration_features(div_2_traxels, "Mean");
    compute_child_deceleration_features(div_2_traxels, "Variance");
    compute_child_deceleration_outlier(div_2_traxels, "RegionCenter");
    compute_child_deceleration_outlier(div_2_traxels, "Count");
    compute_child_deceleration_outlier(div_2_traxels, "Mean");
    compute_child_deceleration_outlier(div_2_traxels, "Variance");
    compute_svm_division_feature_outlier(div_1_traxels);

    save_division_traxels_to_h5(div_1_traxels);
}

void TrackingFeatureExtractor::compute_all_app_dis_features()
{
    typedef AppearanceTraxels::AppearanceType AppearanceType;

    LOG(logDEBUG) << "Extract all appearances";
    AppearanceTraxels appearance_extractor(AppearanceType::Appearance);
    ConstTraxelRefVectors all_app_traxels = appearance_extractor(*graph_);
    LOG(logDEBUG) << "Extract all disappearances";
    AppearanceTraxels disappearance_extractor(AppearanceType::Disappearance);
    ConstTraxelRefVectors all_disapp_traxels = disappearance_extractor(*graph_);
    LOG(logDEBUG) << "Extract filtered appearances";
    AppearanceTraxels appearance_extractor_f(
        AppearanceType::Appearance,
        margin_filter_function_);
    ConstTraxelRefVectors filtered_app_traxels = appearance_extractor_f(*graph_);
    LOG(logDEBUG) << "Extract filtered disappearances";
    AppearanceTraxels disappearance_extractor_f(
        AppearanceType::Disappearance,
        margin_filter_function_);
    ConstTraxelRefVectors filtered_disapp_traxels = disappearance_extractor_f(*graph_);
    push_back_feature(
        "Share of appearances within margin",
        static_cast<double>(filtered_app_traxels.size() / all_app_traxels.size()));
    push_back_feature(
        "Share of disappearances within margin",
        static_cast<double>(filtered_disapp_traxels.size() / all_disapp_traxels.size()));
    compute_border_distances(all_app_traxels, "appearance");
    compute_border_distances(all_disapp_traxels, "disappearance");
}

void TrackingFeatureExtractor::compute_features()
{
    compute_all_track_features();
    compute_all_division_features();
    compute_all_app_dis_features();
}

void TrackingFeatureExtractor::set_track_feature_output_file(const std::string &filename)
{
    track_feature_output_file_ = filename;
}

void TrackingFeatureExtractor::append_feature_vector_to_file(const std::string& filename)
{
    std::vector< std::vector <double> > other_proposal_features;

    {
        std::ifstream feature_vector_file(filename.c_str());
        if(feature_vector_file.good())
        {
            while(!feature_vector_file.eof())
            {
                // read line and remove comments and whitespace
                std::string line;
                std::getline(feature_vector_file, line);
                std::string::size_type comment_start = line.find('#');
                if(comment_start != std::string::npos)
                {
                    line = line.substr(comment_start);
                }
                boost::algorithm::trim(line);

                // skip lines without features
                if(line.size() == 0)
                {
                    continue;
                }

                // read features
                std::stringstream linestream(line);
                other_proposal_features.push_back(std::vector<double>());
                while(!linestream.eof())
                {
                    double f;
                    linestream >> f;
                    other_proposal_features.back().push_back(f);
                }
            }
        }

        if(other_proposal_features.size() > 0 && other_proposal_features.size() != joint_feature_vector_.size())
        {
            throw std::runtime_error("Feature vector already stored in file has different dimension!");
        }

        // stream is closed on scope exit
    }

    std::ofstream feature_vector_file(filename.c_str());

    for(size_t feature_idx = 0; feature_idx < joint_feature_vector_.size(); feature_idx++)
    {
        // add old features if there were any
        if(other_proposal_features.size() > 0)
        {
            for(size_t proposal_idx = 0; proposal_idx < other_proposal_features[feature_idx].size(); proposal_idx++)
            {
                feature_vector_file << other_proposal_features[feature_idx][proposal_idx] << " ";
            }
        }

        // add new feature and end line
        feature_vector_file << joint_feature_vector_[feature_idx] << "\n";
    }
}

void TrackingFeatureExtractor::compute_sq_diff_features(
    ConstTraxelRefVectors& track_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator sq_diff_mmmv;
//    sq_diff_mmmv.set_max(0.0);

    // for each track:
    size_t track_id = 0;
    for(auto track : track_traxels)
    {
        // only compute velocities if track is longer than 1 element
        if(track.size() < 2)
        {
            continue;
        }

        // extract features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(track, feature_matrix);

        // compute the squared norm of the differences of the column vectors
        FeatureMatrix sq_diff_matrix;
        sq_diff_calc_ptr_->calculate(feature_matrix, sq_diff_matrix);

        // add values to min/max/mean/var calculator
        sq_diff_mmmv.add_values(sq_diff_matrix);
        save_features_to_h5(track_id++, "sq_diff_" + feature_name, sq_diff_matrix);
    }
    push_back_feature("squared difference of " + feature_name, sq_diff_mmmv);
}

void TrackingFeatureExtractor::compute_sq_accel_features(
    ConstTraxelRefVectors& track_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator sq_accel_mmmv;
//    sq_accel_mmmv.set_max(0.0);

    // for each track:
    size_t track_id = 0;
    for(auto track : track_traxels)
    {
        // only compute accelerations if track is longer than 2 elements
        if(track.size() < 3)
        {
            continue;
        }

        // extract features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(track, feature_matrix);

        // compute the squared norm of the accelerations of the column vectors
        FeatureMatrix sq_accel_matrix;
        sq_curve_calc_ptr_->calculate(feature_matrix, sq_accel_matrix);

        // add values to min/max/mean/var calculator
        sq_accel_mmmv.add_values(sq_accel_matrix);
        save_features_to_h5(track_id++, "sq_accel_" + feature_name, sq_accel_matrix);
    }
    push_back_feature("squared acceleration of " + feature_name, sq_accel_mmmv);
}

void TrackingFeatureExtractor::compute_angle_features(
    ConstTraxelRefVectors& track_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator angle_mmmv;
//    angle_mmmv.set_max(0.0);

    // for each track:
    size_t track_id = 0;
    for(auto track : track_traxels)
    {
        // only compute angles in track if track is longer than 2 elements
        if(track.size() < 3)
        {
            continue;
        }

        // extract features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(track, feature_matrix);

        // compute for all triples of features the angle of change of direction
        FeatureMatrix angles;
        angle_cos_calc_ptr_->calculate(feature_matrix, angles);

        angle_mmmv.add_values(angles);
        save_features_to_h5(track_id++, "angles_" + feature_name, angles);
    }
    push_back_feature(
        "Mean of all angle cosines of feature " + feature_name,
        angle_mmmv.get_mean());
    push_back_feature(
        "Variance of all angle cosines of features " + feature_name,
        angle_mmmv.get_var());
}

void TrackingFeatureExtractor::compute_track_length_features(
    ConstTraxelRefVectors& track_traxels)
{
    MinMaxMeanVarCalculator track_length_mmmv;
//    track_length_mmmv.set_max(0.0);

    size_t track_id = 0;
    for (auto track : track_traxels)
    {
        track_length_mmmv.add_value(static_cast<double>(track.size()));
        FeatureMatrix m(vigra::Shape2(1, 1));
        m(0, 0) = track.size();
        save_features_to_h5(track_id++, "track_length", m);
    }
    push_back_feature("track length", track_length_mmmv);
}

void TrackingFeatureExtractor::compute_track_id_outlier(
    ConstTraxelRefVectors& track_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator id_out_mmmv;
//    id_out_mmmv.set_max(0.0);
    size_t track_id = 0;
    for (auto track : track_traxels)
    {
        // extract features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(track, feature_matrix);

        // calculate identity outlier if num_samples > dim
        FeatureMatrix temp;
        if (feature_matrix.size(0) > feature_matrix.size(1))
        {
            mvn_outlier_calc_ptr_->calculate(feature_matrix, temp);
            id_out_mmmv.add_value(temp(0, 0));
            save_features_to_h5(track_id++, "outlier_id_" + feature_name, temp);
        }
    }
    push_back_feature("track outlier share of " + feature_name, id_out_mmmv);
}

void TrackingFeatureExtractor::compute_track_diff_outlier(
    ConstTraxelRefVectors& track_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator diff_out_mmmv;
//    diff_out_mmmv.set_max(0.0);
    size_t track_id = 0;
    for (auto track : track_traxels)
    {
        // extract features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(track, feature_matrix);

        // calculate differences
        FeatureMatrix diff_matrix;
        diff_calc_ptr_->calculate(feature_matrix, diff_matrix);

        // calculate diff outlier if num_samples > dim
        FeatureMatrix temp;
        if (diff_matrix.size(0) > diff_matrix.size(1))
        {
            mvn_outlier_calc_ptr_->calculate(diff_matrix, temp);
            diff_out_mmmv.add_value(temp(0, 0));
            save_features_to_h5(track_id++, "diff_outlier_" + feature_name, temp);
        }
    }
    push_back_feature(
        "track outlier share of " + feature_name + " diff",
        diff_out_mmmv);
}

void TrackingFeatureExtractor::compute_svm_track_feature_outlier(
    ConstTraxelRefVectors& track)
{
    FeatureMatrix track_feature_matrix;
    TrackFeatureExtractor track_feature_extractor;
    track_feature_extractor.compute_features(track, track_feature_matrix);
    FeatureMatrix score_matrix;
    if (svm_track_outlier_calc_ptr_->is_trained())
    {
        svm_track_outlier_calc_ptr_->calculate(track_feature_matrix, score_matrix);
    }
    else
    {
        LOG(logWARNING) << "in TrackingFeatureExtractor::compute_track_feature_outlier()";
        LOG(logWARNING) << "SVMOutlierCalculator not trained";
        score_matrix.reshape(vigra::Shape2(1, 1));
        score_matrix.init(0.0);
    }
    MinMaxMeanVarCalculator mmmv_score;
//    mmmv_score.set_max(0.0);
    mmmv_score.add_values(score_matrix);
    push_back_feature("track feature outlier score", mmmv_score);

    if(track_feature_output_file_.size() > 0)
    {
        vigra::writeHDF5(track_feature_output_file_.c_str(), "track_outliers_svm", score_matrix);
    }
}

void TrackingFeatureExtractor::compute_division_sq_diff_features(
    ConstTraxelRefVectors& div_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator sq_diff_mmmv;
//    sq_diff_mmmv.set_max(0.0);

    size_t division_id = 0;
    for (auto div : div_traxels)
    {
        // extract features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(div, feature_matrix);

        // calculate the squared child-parent difference
        FeatureMatrix temp;
        FeatureMatrix sq_diff_matrix;
        child_parent_diff_calc_ptr_->calculate(feature_matrix, temp);
        sq_norm_calc_ptr_->calculate(temp, sq_diff_matrix);
        sq_diff_mmmv.add_values(sq_diff_matrix);
        save_features_to_h5(division_id++, "sq_diff_" + feature_name, sq_diff_matrix, false);
    }
    push_back_feature(
        "squared child-parent " + feature_name + " difference",
        sq_diff_mmmv);
}

void TrackingFeatureExtractor::compute_division_sq_diff_outlier(
    ConstTraxelRefVectors& div_traxels,
    std::string feature_name)
{
    // extract all squared child-parent differences
    boost::shared_ptr<TraxelsFeaturesIdentity> feature_extractor_ptr(
        new TraxelsFeaturesIdentity(feature_name));
    CompositionCalculator sq_diff_calc(
        child_parent_diff_calc_ptr_,
        sq_norm_calc_ptr_);
    FeatureMatrix sq_diff_matrix;
    sq_diff_calc.calculate_for_all(
        div_traxels,
        sq_diff_matrix,
        feature_extractor_ptr);
    double outlier = 0.0;
    if (sq_diff_matrix.size(0) > sq_diff_matrix.size(1))
    {
        FeatureMatrix score_matrix;
        sq_mahal_calc_ptr_->calculate(sq_diff_matrix, score_matrix);
        assert(score_matrix.size(0) == sq_diff_matrix.size(0));
        size_t division_count = score_matrix.size(0);

        // get the outlier count and save the outlier scores
        for (size_t col = 0; col < division_count; col++)
        {
            if (score_matrix(col, 0) > 3.0 * 3.0)
            {
                outlier += 1.0;
            }
            FeatureMatrix score(vigra::Shape2(1, 1), score_matrix(col, 0));
            save_features_to_h5(
                col,
                "sq_diff_" + feature_name + "_outlier_score",
                score,
                false);
        }
        outlier /= static_cast<double>(division_count);
    }
    push_back_feature(
        "outlier of squared child-parent " + feature_name + " difference",
        outlier);
}

void TrackingFeatureExtractor::compute_child_deceleration_features(
    ConstTraxelRefVectors& div_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator child_decel_mmmv;
//    child_decel_mmmv.set_max(0.0);

    size_t division_id = 0;
    for (auto div : div_traxels)
    {
        // extract features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(div, feature_matrix);

        // calculate the child decelerations
        FeatureMatrix child_decel_mat;
        child_decel_calc_ptr_->calculate(feature_matrix, child_decel_mat);

        child_decel_mmmv.add_values(child_decel_mat);
        save_features_to_h5(division_id++, "child_decel_" + feature_name, child_decel_mat, false);
    }
    push_back_feature("child " + feature_name + " deceleration", child_decel_mmmv);
}

void TrackingFeatureExtractor::compute_child_deceleration_outlier(
    ConstTraxelRefVectors& div_traxels,
    std::string feature_name)
{
    // extract all child decelerations
    boost::shared_ptr<TraxelsFeaturesIdentity> feature_extractor_ptr(
        new TraxelsFeaturesIdentity(feature_name));
    FeatureMatrix child_decel_mat;
    child_decel_calc_ptr_->calculate_for_all(
        div_traxels,
        child_decel_mat,
        feature_extractor_ptr);
    double outlier = 0.0;
    if (child_decel_mat.size(0) > child_decel_mat.size(1))
    {
        FeatureMatrix score_matrix;
        sq_mahal_calc_ptr_->calculate(child_decel_mat, score_matrix);
        assert(score_matrix.size(0) == child_decel_mat.size(0));
        size_t division_count = score_matrix.size(0);

        // get the outlier count and save the outlier scores
        for (size_t col = 0; col < division_count; col++)
        {
            if (score_matrix(col, 0) > 3.0 * 3.0)
            {
                outlier += 1.0;
            }
            FeatureMatrix score(vigra::Shape2(1, 1), score_matrix(col, 0));
            save_features_to_h5(
                col,
                "child_decel_" + feature_name + "_outlier_score",
                score,
                false);
        }
        outlier /= static_cast<double>(division_count);
    }
    push_back_feature("Outlier in child " + feature_name + " decelerations", outlier);
}

void TrackingFeatureExtractor::compute_svm_division_feature_outlier(
    ConstTraxelRefVectors& divisions)
{
    FeatureMatrix div_feature_matrix;
    DivisionFeatureExtractor div_feature_extractor;
    div_feature_extractor.compute_features(divisions, div_feature_matrix);
    FeatureMatrix score_matrix;
    if (svm_div_outlier_calc_ptr_->is_trained())
    {
        svm_div_outlier_calc_ptr_->calculate(div_feature_matrix, score_matrix);
    }
    else
    {
        LOG(logWARNING) << "in TrackingFeatureExtractor::compute_div_feature_outlier()";
        LOG(logWARNING) << "SVMOutlierCalculator not trained";
        score_matrix.reshape(vigra::Shape2(1, 1));
        score_matrix.init(0.0);
    }
    MinMaxMeanVarCalculator mmmv_score;
//    mmmv_score.set_max(0.0);
    mmmv_score.add_values(score_matrix);
    push_back_feature("division feature outlier score", mmmv_score);

    if(track_feature_output_file_.size() > 0)
    {
        vigra::writeHDF5(track_feature_output_file_.c_str(), "division_outliers_svm", score_matrix);
    }
}

void TrackingFeatureExtractor::compute_border_distances(
    ConstTraxelRefVectors& appearance_traxels,
    std::string description)
{
    MinMaxMeanVarCalculator border_dist_mmmv;
//    border_dist_mmmv.set_max(0.0);

    for (auto appearance : appearance_traxels)
    {
        FeatureMatrix position;
        TraxelsFeaturesIdentity position_extractor("RegionCenter");
        position_extractor.extract(appearance, position);
        double border_dist = 0;
        if (position.shape(1) == 3)
        {
            border_dist = fov_.spatial_distance_to_border(
                              0.0,
                              static_cast<double>(position(0, 0)),
                              static_cast<double>(position(0, 1)),
                              static_cast<double>(position(0, 2)),
                              false
                          );
        }
        else if (position.shape(1) == 2)
        {
            border_dist = fov_.spatial_distance_to_border(
                              0.0,
                              static_cast<double>(position(0, 0)),
                              static_cast<double>(position(0, 1)),
                              0.0,
                              false
                          );
        }
        else
        {
            LOG(logDEBUG) << "In TrackingFeatureExtractor::"
                          << "compute_appearance_border_distances:";
            LOG(logDEBUG) << "data is neither 2d nor 3d";
        }
        border_dist_mmmv.add_value(border_dist);
    }
    push_back_feature(description + " border distances", border_dist_mmmv);
}

void TrackingFeatureExtractor::push_back_feature(
    std::string feature_name,
    double feature_value)
{
    joint_feature_vector_.push_back(feature_value);
    feature_descriptions_.push_back(feature_name);
}

void TrackingFeatureExtractor::push_back_feature(
    std::string feature_name,
    const MinMaxMeanVarCalculator& mmmv_calculator)
{
    push_back_feature("Mean of " + feature_name, mmmv_calculator.get_mean());
    push_back_feature("Variance of " + feature_name, mmmv_calculator.get_var());
    push_back_feature("Min of " + feature_name, mmmv_calculator.get_min());
    push_back_feature("Max of " + feature_name, mmmv_calculator.get_max());
}

BorderDistanceFilter::BorderDistanceFilter(
    const FieldOfView& field_of_view,
    double t_margin,
    double spatial_margin)
{
    std::vector<double> lower_bound = field_of_view.lower_bound();
    std::vector<double> upper_bound = field_of_view.upper_bound();
    if ((upper_bound[0] - lower_bound[0]) < 2 * t_margin)
    {
        t_margin = (upper_bound[0] - lower_bound[0]) / 2.0;
    }
    lower_bound[0] += t_margin;
    upper_bound[0] -= t_margin;
    for (size_t i = 1; i < lower_bound.size(); i++)
    {
        double spatial_margin_temp = spatial_margin;
        if ((upper_bound[i] - lower_bound[i]) < 2 * spatial_margin)
        {
            spatial_margin_temp = (upper_bound[i] - lower_bound[i]) / 2.0;
        }
        lower_bound[i] += spatial_margin_temp;
        upper_bound[i] -= spatial_margin_temp;
    }
    fov_.set_boundingbox(
        lower_bound[0],
        lower_bound[1],
        lower_bound[2],
        lower_bound[3],
        upper_bound[0],
        upper_bound[1],
        upper_bound[2],
        upper_bound[3]
    );
}

bool BorderDistanceFilter::is_out_of_margin(const Traxel& traxel) const
{
    bool ret = false;
    const FeatureMap& feature_map = traxel.features.get();
    FeatureMap::const_iterator com_it = feature_map.find("RegionCenter");
    if (com_it == feature_map.end())
    {
        LOG(logDEBUG) << "in BorderDistanceFilter::is_out_of_margin(): "
                      << "feature map of traxel has no key \"com\"";
        LOG(logDEBUG) << "return \"false\"";
    }
    else
    {
        double t = static_cast<double>(traxel.Timestep);
        const std::vector<feature_type>& pos = com_it->second;
        if (pos.size() == 2)
        {
            ret = fov_.contains(t, pos[0], pos[1], 0.0);
        }
        else if (pos.size() == 3)
        {
            ret = fov_.contains(t, pos[0], pos[1], pos[2]);
        }
        else
        {
            LOG(logDEBUG) << "in BorderDistanceFilter::is_out_of_margin(): "
                          << "dimension is neither 2d+t nor 3d+t";
            LOG(logDEBUG) << "return \"false\"";
        }
    }
    return ret;
}

} // end namespace features
} // end namespace pgmlink
