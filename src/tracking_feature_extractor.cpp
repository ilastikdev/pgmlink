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
    row_sum_calc_ptr_(new SumCalculator<0>),
    row_var_calc_ptr_(new VarianceCalculator)
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
    // extract all tracks
    TrackTraxels track_extractor;
    ConstTraxelRefVectors track_traxels = track_extractor(graph_);
    // extract all divisions to depth 1
    DivisionTraxels div_1_extractor(1);
    ConstTraxelRefVectors div_1_traxels = div_1_extractor(graph_);

    compute_velocity_features(track_traxels);
    compute_acceleration_features(track_traxels);
    compute_track_length_features(track_traxels);
    //compute_size_difference_features();
}

void TrackingFeatureExtractor::compute_velocity_features(
    ConstTraxelRefVectors& track_traxels)
{

    size_t num_velocity_entries = 0;
    double sum_of_squared_velocities = 0;
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

        FeatureMatrix sum_velocity;
        row_sum_calc_ptr_->calculate(velocities, sum_velocity);

        sum_of_squared_velocities += sum_velocity(0,0);
        num_velocity_entries += track.size() - 1;
    }

    double mean_squared_velocity = sum_of_squared_velocities / num_velocity_entries;

    push_back_feature("Mean of all velocities (squared)", mean_squared_velocity);
    push_back_feature("Min of all velocities (squared)", min_squared_velocity);
    push_back_feature("Max of all velocities (squared)", max_squared_velocity);
}

void TrackingFeatureExtractor::compute_acceleration_features(
    ConstTraxelRefVectors& track_traxels)
{
    size_t num_accel_entries = 0;
    double sq_accel_sum = 0;
    double sq_accel_min = std::numeric_limits<double>::max();
    double sq_accel_max = 0.0;

    // for each track:
    for(auto track : track_traxels)
    {
        // only compute velocities if track is longer than 2 elements
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

        row_sum_calc_ptr_->calculate(sq_accel, temp);
        sq_accel_sum+= temp(0,0);

        num_accel_entries += track.size() - 2;
    }

    if (num_accel_entries != 0)
        sq_accel_sum/= num_accel_entries;
    else
        sq_accel_sum = 0.0;

    push_back_feature("Mean of all accelerations (squared)", sq_accel_sum);
    push_back_feature("Min of all accelerations (squared)", sq_accel_min);
    push_back_feature("Max of all accelerations (squared)", sq_accel_max);
}

void TrackingFeatureExtractor::compute_track_length_features(
    ConstTraxelRefVectors& track_traxels)
{
    double sum_track_length = 0;
    double min_track_length = std::numeric_limits<double>::max();
    double max_track_length = 0.0;
    for (auto track : track_traxels)
    {
        double track_size = static_cast<double>(track.size());
        min_track_length = std::min(min_track_length, track_size);
        max_track_length = std::max(max_track_length, track_size);
        sum_track_length+= track_size;
    }
    if (track_traxels.size() != 0) {
      sum_track_length/= track_traxels.size();
    }
    push_back_feature("Mean of track length", sum_track_length);
    push_back_feature("Min of track length", min_track_length);
    push_back_feature("Max of track length", max_track_length);
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
