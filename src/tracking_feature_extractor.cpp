#include "pgmlink/tracking_feature_extractor.h"

namespace pgmlink {
namespace features {

TrackingFeatureExtractor::TrackingFeatureExtractor(HypothesesGraph &graph):
    graph_(graph),
    position_extractor_ptr_(new TraxelsFeaturesIdentity("com")),
    sq_diff_calc_ptr_(new SquaredDiffCalculator),
    sq_curve_calc_ptr_(new SquaredCurveCalculator),
    row_min_calc_ptr_(new MinCalculator<0>),
    row_max_calc_ptr_(new MaxCalculator<0>),
    angle_cos_calc_ptr_(new AngleCosineCalculator),
    child_parent_diff_calc_ptr_(new ChildParentDiffCalculator),
    sq_norm_calc_ptr_(new SquaredNormCalculator<0>)
{
}

void TrackingFeatureExtractor::get_feature_vector(TrackingFeatureExtractor::JointFeatureVector &feature_vector) const
{
    feature_vector.clear();
    feature_vector.insert(feature_vector.begin(), joint_feature_vector_.begin(), joint_feature_vector_.end());
}

const std::string TrackingFeatureExtractor::get_feature_description(size_t feature_index) const
{
    return feature_descriptions_[feature_index];
}

void TrackingFeatureExtractor::compute_features()
{
    LOG(logDEBUG) << "Extract all tracks";
    TrackTraxels track_extractor;
    ConstTraxelRefVectors track_traxels = track_extractor(graph_);
    LOG(logDEBUG) << "Extract all divisions to depth 1";
    DivisionTraxels div_1_extractor(1);
    ConstTraxelRefVectors div_1_traxels = div_1_extractor(graph_);

    compute_velocity_features(track_traxels);
    compute_acceleration_features(track_traxels);
    compute_angle_features(track_traxels);
    compute_track_length_features(track_traxels);
    compute_division_move_distance(div_1_traxels);
    //compute_size_difference_features();
}

void TrackingFeatureExtractor::compute_velocity_features(
    ConstTraxelRefVectors& track_traxels)
{

    size_t n = 0;
    double mean = 0.0;
    double prev_mean = 0.0;
    double sum_squared_diff = 0.0;
    double min_squared_velocity = std::numeric_limits<double>::max();
    double max_squared_velocity = 0.0;

    // for each track:
    for(auto track : track_traxels)
    {
        // only compute velocities if track is longer than 1 element
        if(track.size() < 2)
            continue;

        // extract positions
        FeatureMatrix positions;
        position_extractor_ptr_->extract(track, positions);

        // compute velocity vector's squared magnitude for all pairs of positions
        FeatureMatrix velocities;
        sq_diff_calc_ptr_->calculate(positions, velocities);

        // compute per track min/max/sum of velocity
        FeatureMatrix min_velocity;
        row_min_calc_ptr_->calculate(velocities, min_velocity);
        min_squared_velocity = std::min(min_squared_velocity, double(min_velocity(0,0)));

        FeatureMatrix max_velocity;
        row_max_calc_ptr_->calculate(velocities, max_velocity);
        max_squared_velocity = std::max(max_squared_velocity, double(max_velocity(0,0)));

        for(FeatureMatrix::iterator v_it = velocities.begin();
            v_it != velocities.end();
            v_it++)
        {
            n++;
            prev_mean = mean;
            mean += (*v_it - prev_mean) / n;
            sum_squared_diff += (*v_it - prev_mean) * (*v_it - mean);
        }
    }

    if (n != 0)
        sum_squared_diff /= n;
    push_back_feature("Mean of all velocities (squared)", mean);
    push_back_feature("Variance of all velocities (squared)", sum_squared_diff);
    push_back_feature("Min of all velocities (squared)", min_squared_velocity);
    push_back_feature("Max of all velocities (squared)", max_squared_velocity);
}

void TrackingFeatureExtractor::compute_acceleration_features(
    ConstTraxelRefVectors& track_traxels)
{
    size_t n = 0;
    double mean = 0.0;
    double prev_mean = 0.0;
    double sum_squared_diff = 0.0;
    double sq_accel_min = std::numeric_limits<double>::max();
    double sq_accel_max = 0.0;

    // for each track:
    for(auto track : track_traxels)
    {
        // only compute accelerations if track is longer than 2 elements
        if(track.size() < 3)
            continue;

        // extract positions
        FeatureMatrix positions;
        position_extractor_ptr_->extract(track, positions);

        // compute acceleration vector's squared magnitude for all triples of positions
        FeatureMatrix sq_accel;
        sq_curve_calc_ptr_->calculate(positions, sq_accel);

        // compute per track min/max/sum of acceleration
        FeatureMatrix temp;

        row_min_calc_ptr_->calculate(sq_accel, temp);
        sq_accel_min = std::min(sq_accel_min, double(temp(0,0)));

        row_max_calc_ptr_->calculate(sq_accel, temp);
        sq_accel_max = std::max(sq_accel_max, double(temp(0,0)));

        for(FeatureMatrix::iterator a_it = sq_accel.begin();
            a_it != sq_accel.end();
            a_it++)
        {
            n++;
            prev_mean = mean;
            mean += (*a_it - prev_mean) / n;
            sum_squared_diff += (*a_it - prev_mean) * (*a_it - mean);
        }
    }

    if (n != 0)
        sum_squared_diff /= n;
    push_back_feature("Mean of all accelerations (squared)", mean);
    push_back_feature("Variance of all accelerations (squared)", sum_squared_diff);
    push_back_feature("Min of all accelerations (squared)", sq_accel_min);
    push_back_feature("Max of all accelerations (squared)", sq_accel_max);
}

void TrackingFeatureExtractor::compute_angle_features(
    ConstTraxelRefVectors& track_traxels)
{
    size_t n = 0;
    double mean = 0.0;
    double prev_mean = 0.0;
    double sum_squared_diff = 0.0;

    // for each track:
    for(auto track : track_traxels)
    {
        // only compute angles in track if track is longer than 2 elements
        if(track.size() < 3)
            continue;

        // extract positions
        FeatureMatrix positions;
        position_extractor_ptr_->extract(track, positions);

        // compute for all triples of positions the angle of change of direction
        FeatureMatrix angles;
        angle_cos_calc_ptr_->calculate(positions, angles);

        for(FeatureMatrix::iterator a_it = angles.begin();
            a_it != angles.end();
            a_it++)
        {
            n++;
            prev_mean = mean;
            mean += (*a_it - prev_mean) / n;
            sum_squared_diff += (*a_it - prev_mean) * (*a_it - mean);
        }
    }
    if (n != 0)
        sum_squared_diff /= n;
    push_back_feature("Mean of all angle cosines", mean);
    push_back_feature("Variance of all angle cosines", sum_squared_diff);
}

void TrackingFeatureExtractor::compute_track_length_features(
    ConstTraxelRefVectors& track_traxels)
{
    size_t n = 0;
    double mean = 0.0;
    double prev_mean = 0.0;
    double sum_squared_diff = 0.0;
    double min_track_length = std::numeric_limits<double>::max();
    double max_track_length = 0.0;
    for (auto track : track_traxels)
    {
        double track_size = static_cast<double>(track.size());
        min_track_length = std::min(min_track_length, track_size);
        max_track_length = std::max(max_track_length, track_size);
        n++;
        prev_mean = mean;
        mean += (track_size - prev_mean) / n;
        sum_squared_diff += (track_size - prev_mean) * (track_size - mean);
    }
    if (n != 0) {
      sum_squared_diff /= n;
    }
    push_back_feature("Mean of track length", mean);
    push_back_feature("Variance of track length", sum_squared_diff);
    push_back_feature("Min of track length", min_track_length);
    push_back_feature("Max of track length", max_track_length);
}

void TrackingFeatureExtractor::compute_division_move_distance(
    ConstTraxelRefVectors& div_traxels)
{
    size_t n = 0;
    double mean = 0.0;
    double prev_mean = 0.0;
    double sum_squared_diff = 0.0;
    double min = std::numeric_limits<double>::max();
    double max = 0.0;

    for (auto div : div_traxels)
    {
        // extract positions
        FeatureMatrix positions;
        position_extractor_ptr_->extract(div, positions);

        // calculate the squared move distance
        FeatureMatrix temp;
        FeatureMatrix sq_move_dist;
        child_parent_diff_calc_ptr_->calculate(positions, temp);
        sq_norm_calc_ptr_->calculate(temp, sq_move_dist);
        for(FeatureMatrix::iterator md_it = sq_move_dist.begin();
            md_it != sq_move_dist.end();
            md_it++)
        {
            min = std::min(min, double(*md_it));
            max = std::max(max, double(*md_it));
            n++;
            prev_mean = mean;
            mean += (*md_it - prev_mean) / n;
            sum_squared_diff += (*md_it - prev_mean) * (*md_it - mean);
        }
    }
    if (n != 0)
      sum_squared_diff /= n;
    push_back_feature("Mean of child-parent velocities (squared)", mean);
    push_back_feature("Variance of child-parent velocities (squared)", sum_squared_diff);
    push_back_feature("Min of child-parent velocities (squared)", min);
    push_back_feature("Max of child-parent velocities (squared)", max);
}

void TrackingFeatureExtractor::compute_size_difference_features()
{
    throw std::runtime_error("not yet implemented");
}

void TrackingFeatureExtractor::push_back_feature(
    std::string feature_name,
    double feature_value)
{
    joint_feature_vector_.push_back(feature_value);
    feature_descriptions_.push_back(feature_name);
}

} // end namespace features
} // end namespace pgmlink
