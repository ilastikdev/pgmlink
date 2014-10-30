#include "pgmlink/tracking_feature_extractor.h"
#include <boost/algorithm/string.hpp>

namespace pgmlink {
namespace features {

MinMaxMeanVarCalculator::MinMaxMeanVarCalculator()
{
    reset();
}

void MinMaxMeanVarCalculator::reset()
{
    n = 0;
    mean = 0;
    sum_squared_diff = 0;
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::lowest();
}

void MinMaxMeanVarCalculator::add_value(const double& value)
{
    n++;
    double delta = value - mean;
    mean += delta / n;
    sum_squared_diff += delta * (value - mean);

    min = std::min(min, value);
    max = std::max(max, value);
}

void MinMaxMeanVarCalculator::add_values(const FeatureMatrix& values)
{
    for(FeatureMatrix::const_iterator it = values.begin();
        it != values.end();
        it++)
    {
        add_value(*it);
    }
}

void MinMaxMeanVarCalculator::set_min(const double& value)
{
    min = value;
}

void MinMaxMeanVarCalculator::set_max(const double& value)
{
    max = value;
}

size_t MinMaxMeanVarCalculator::get_count() const { return n; }

double MinMaxMeanVarCalculator::get_mean() const { return mean; }

double MinMaxMeanVarCalculator::get_var() const
{
    if (n)
        return sum_squared_diff / n;
    else
        return 0.0;
}

double MinMaxMeanVarCalculator::get_min() const { return min; }

double MinMaxMeanVarCalculator::get_max() const { return max; }

TrackingFeatureExtractor::TrackingFeatureExtractor(boost::shared_ptr<HypothesesGraph> graph,
        const FieldOfView& fov,
        boost::function<bool (const Traxel&)> margin_filter_function):
    graph_(graph),
    fov_(fov),
    margin_filter_function_(margin_filter_function),
    position_extractor_ptr_(new TraxelsFeaturesIdentity("com")),
    sq_diff_calc_ptr_(new SquaredDiffCalculator),
    diff_calc_ptr_(new DiffCalculator),
    sq_curve_calc_ptr_(new SquaredCurveCalculator),
    row_min_calc_ptr_(new MinCalculator<0>),
    row_max_calc_ptr_(new MaxCalculator<0>),
    mvn_outlier_calc_ptr_(new MVNOutlierCalculator),
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
    position_extractor_ptr_(new TraxelsFeaturesIdentity("com")),
    sq_diff_calc_ptr_(new SquaredDiffCalculator),
    diff_calc_ptr_(new DiffCalculator),
    sq_curve_calc_ptr_(new SquaredCurveCalculator),
    row_min_calc_ptr_(new MinCalculator<0>),
    row_max_calc_ptr_(new MaxCalculator<0>),
    mvn_outlier_calc_ptr_(new MVNOutlierCalculator),
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

const std::string TrackingFeatureExtractor::get_feature_description(size_t feature_index) const
{
    return feature_descriptions_[feature_index];
}

void TrackingFeatureExtractor::compute_features()
{
    typedef AppearanceTraxels::AppearanceType AppearanceType;
    // extract traxels of interest
    LOG(logDEBUG) << "Extract all tracks";
    TrackTraxels track_extractor;
    ConstTraxelRefVectors track_traxels = track_extractor(*graph_);
    LOG(logDEBUG) << "Extract all divisions to depth 1";
    DivisionTraxels div_1_extractor(1);
    ConstTraxelRefVectors div_1_traxels = div_1_extractor(*graph_);
    LOG(logDEBUG) << "Extract all divisions to depth 2";
    DivisionTraxels div_2_extractor(2);
    ConstTraxelRefVectors div_2_traxels = div_2_extractor(*graph_);
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

    compute_velocity_features(track_traxels);
    compute_acceleration_features(track_traxels);
    compute_angle_features(track_traxels);
    compute_track_length_features(track_traxels);
    compute_track_outlier_features(track_traxels);
    compute_division_move_distance(div_1_traxels);
    compute_division_move_outlier(div_1_traxels);
    compute_child_deceleration_features(div_2_traxels);
    compute_child_deceleration_outlier(div_2_traxels);
    push_back_feature(
        "Count of all appearances",
        static_cast<double>(all_app_traxels.size()));
    push_back_feature(
        "Count of all disappearances",
        static_cast<double>(all_disapp_traxels.size()));
    push_back_feature(
        "Count of appearances within margin",
        static_cast<double>(filtered_app_traxels.size()));
    push_back_feature(
        "Count of disappearances within margin",
        static_cast<double>(filtered_disapp_traxels.size()));
    compute_border_distances(all_app_traxels, "appearance");
    compute_border_distances(all_disapp_traxels, "disappearance");
    //compute_size_difference_features();
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
                    line = line.substr(comment_start);
                boost::algorithm::trim(line);

                // skip lines without features
                if(line.size() == 0)
                    continue;

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

void TrackingFeatureExtractor::compute_velocity_features(ConstTraxelRefVectors& track_traxels)
{
    MinMaxMeanVarCalculator sq_velocity_mmmv;
    sq_velocity_mmmv.set_max(0.0);

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

        // add values to min/max/mean/var calculator
        sq_velocity_mmmv.add_values(velocities);
    }
    push_back_feature("all velocities (squared)", sq_velocity_mmmv);
}

void TrackingFeatureExtractor::compute_acceleration_features(
    ConstTraxelRefVectors& track_traxels)
{
    MinMaxMeanVarCalculator sq_accel_mmmv;
    sq_accel_mmmv.set_max(0.0);

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

        // add values to min/max/mean/var calculator
        sq_accel_mmmv.add_values(sq_accel);
    }
    push_back_feature("all accelerations (squared)", sq_accel_mmmv);
}

void TrackingFeatureExtractor::compute_angle_features(
    ConstTraxelRefVectors& track_traxels)
{
    MinMaxMeanVarCalculator angle_mmmv;
    angle_mmmv.set_max(0.0);

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

        angle_mmmv.add_values(angles);
    }
    push_back_feature("Mean of all angle cosines", angle_mmmv.get_mean());
    push_back_feature("Variance of all angle cosines", angle_mmmv.get_var());
}

void TrackingFeatureExtractor::compute_track_length_features(
    ConstTraxelRefVectors& track_traxels)
{
    MinMaxMeanVarCalculator track_length_mmmv;
    track_length_mmmv.set_max(0.0);

    for (auto track : track_traxels)
    {
        track_length_mmmv.add_value(static_cast<double>(track.size()));
    }
    push_back_feature("track length", track_length_mmmv);
}

void TrackingFeatureExtractor::compute_track_outlier_features(
    ConstTraxelRefVectors& track_traxels)
{
    MinMaxMeanVarCalculator pos_out_mmmv;
    MinMaxMeanVarCalculator vel_out_mmmv;
    for (auto track : track_traxels)
    {
        // extract positions
        FeatureMatrix positions;
        position_extractor_ptr_->extract(track, positions);

        // calculate velocities
        FeatureMatrix velocities;
        diff_calc_ptr_->calculate(positions, velocities);

        // calculate position outlier if num_samples > dim
        FeatureMatrix temp;
        if (positions.size(0) > positions.size(1))
        {
            mvn_outlier_calc_ptr_->calculate(positions, temp);
            pos_out_mmmv.add_value(temp(0, 0));
        }

        // calculate velocity outlier if num_samples > dim
        if (velocities.size(0) > velocities.size(1))
        {
            mvn_outlier_calc_ptr_->calculate(velocities, temp);
            vel_out_mmmv.add_value(temp(0, 0));
        }
    }
    push_back_feature("position outlier count", pos_out_mmmv);
    push_back_feature("velocity outlier count", vel_out_mmmv);
}

void TrackingFeatureExtractor::compute_division_move_distance(
    ConstTraxelRefVectors& div_traxels)
{
    MinMaxMeanVarCalculator move_dist_mmmv;
    move_dist_mmmv.set_max(0.0);

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
        move_dist_mmmv.add_values(sq_move_dist);
    }
    push_back_feature("child-parent velocity (squared)", move_dist_mmmv);
}

void TrackingFeatureExtractor::compute_division_move_outlier(
    ConstTraxelRefVectors& div_traxels)
{
    // extract all squared move distances
    CompositionCalculator sq_move_dist_calc(
        child_parent_diff_calc_ptr_,
        sq_norm_calc_ptr_);
    FeatureMatrix sq_move_distances;
    sq_move_dist_calc.calculate_for_all(
        div_traxels,
        sq_move_distances,
        position_extractor_ptr_);
    double outlier = 0.0;
    if (sq_move_distances.size(0) > sq_move_distances.size(1))
    {
        FeatureMatrix outlier_mat;
        mvn_outlier_calc_ptr_->calculate(sq_move_distances, outlier_mat);
        outlier = outlier_mat(0,0);
    }
    push_back_feature("Outlier in division move distance", outlier);
}

void TrackingFeatureExtractor::compute_child_deceleration_features(
    ConstTraxelRefVectors& div_traxels)
{
    MinMaxMeanVarCalculator child_decel_mmmv;
    child_decel_mmmv.set_max(0.0);

    for (auto div : div_traxels)
    {
        // extract positions
        FeatureMatrix positions;
        position_extractor_ptr_->extract(div, positions);

        // calculate the child decelerations
        FeatureMatrix child_decel_mat;
        child_decel_calc_ptr_->calculate(positions, child_decel_mat);

        child_decel_mmmv.add_values(child_decel_mat);
    }
    push_back_feature("child deceleration", child_decel_mmmv);
}

void TrackingFeatureExtractor::compute_child_deceleration_outlier(
    ConstTraxelRefVectors& div_traxels)
{
    // extract all child decelerations
    FeatureMatrix child_decelerations;
    child_decel_calc_ptr_->calculate_for_all(
        div_traxels,
        child_decelerations,
        position_extractor_ptr_);
    double outlier = 0.0;
    if (child_decelerations.size(0) > child_decelerations.size(1))
    {
        FeatureMatrix outlier_mat;
        mvn_outlier_calc_ptr_->calculate(child_decelerations, outlier_mat);
        outlier = outlier_mat(0,0);
    }
    push_back_feature("Outlier in child decelerations", outlier);
}

void TrackingFeatureExtractor::compute_border_distances(
    ConstTraxelRefVectors& appearance_traxels,
    std::string description)
{
    MinMaxMeanVarCalculator border_dist_mmmv;
    border_dist_mmmv.set_max(0.0);

    for (auto appearance : appearance_traxels)
    {
        FeatureMatrix position;
        position_extractor_ptr_->extract(appearance, position);
        double border_dist = 0;
        if (position.shape(1) == 3)
        {
            border_dist = fov_.spatial_distance_to_border(
                0.0,
                static_cast<double>(position(0,0)),
                static_cast<double>(position(0,1)),
                static_cast<double>(position(0,2)),
                false
            );
        }
        else if (position.shape(1) == 2)
        {
            border_dist = fov_.spatial_distance_to_border(
                0.0,
                static_cast<double>(position(0,0)),
                static_cast<double>(position(0,1)),
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
    FeatureMap::const_iterator com_it = feature_map.find("com");
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
