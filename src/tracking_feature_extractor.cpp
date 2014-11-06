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

    compute_sq_diff_features(track_traxels, "com");
    compute_sq_diff_features(track_traxels, "Count");
    compute_sq_diff_features(track_traxels, "Mean");
    compute_sq_diff_features(track_traxels, "Variance");
    compute_sq_accel_features(track_traxels, "com");
    compute_sq_accel_features(track_traxels, "Count");
    compute_sq_accel_features(track_traxels, "Mean");
    compute_sq_accel_features(track_traxels, "Variance");
    compute_angle_features(track_traxels, "com");
    compute_track_length_features(track_traxels);
    compute_track_id_outlier(track_traxels, "com");
    compute_track_id_outlier(track_traxels, "Count");
    compute_track_id_outlier(track_traxels, "Mean");
    compute_track_id_outlier(track_traxels, "Variance");
    compute_track_diff_outlier(track_traxels, "com");
    compute_track_diff_outlier(track_traxels, "Count");
    compute_track_diff_outlier(track_traxels, "Mean");
    compute_track_diff_outlier(track_traxels, "Variance");
    compute_division_sq_diff_features(div_1_traxels, "com");
    compute_division_sq_diff_features(div_1_traxels, "Count");
    compute_division_sq_diff_features(div_1_traxels, "Mean");
    compute_division_sq_diff_features(div_1_traxels, "Variance");
    compute_division_sq_diff_outlier(div_1_traxels, "com");
    compute_division_sq_diff_outlier(div_1_traxels, "Count");
    compute_division_sq_diff_outlier(div_1_traxels, "Mean");
    compute_division_sq_diff_outlier(div_1_traxels, "Variance");
    compute_child_deceleration_features(div_2_traxels, "com");
    compute_child_deceleration_features(div_2_traxels, "Count");
    compute_child_deceleration_features(div_2_traxels, "Mean");
    compute_child_deceleration_features(div_2_traxels, "Variance");
    compute_child_deceleration_outlier(div_2_traxels, "com");
    compute_child_deceleration_outlier(div_2_traxels, "Count");
    compute_child_deceleration_outlier(div_2_traxels, "Mean");
    compute_child_deceleration_outlier(div_2_traxels, "Variance");
    push_back_feature(
        "Share of appearances within margin",
        static_cast<double>(filtered_app_traxels.size() / all_app_traxels.size()));
    push_back_feature(
        "Share of disappearances within margin",
        static_cast<double>(filtered_disapp_traxels.size() / all_disapp_traxels.size()));
    compute_border_distances(all_app_traxels, "appearance");
    compute_border_distances(all_disapp_traxels, "disappearance");
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

void TrackingFeatureExtractor::compute_sq_diff_features(
    ConstTraxelRefVectors& track_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator sq_diff_mmmv;
    sq_diff_mmmv.set_max(0.0);

    // for each track:
    for(auto track : track_traxels)
    {
        // only compute velocities if track is longer than 1 element
        if(track.size() < 2)
            continue;

        // extract features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(track, feature_matrix);

        // compute the squared norm of the differences of the column vectors
        FeatureMatrix sq_diff_matrix;
        sq_diff_calc_ptr_->calculate(feature_matrix, sq_diff_matrix);

        // add values to min/max/mean/var calculator
        sq_diff_mmmv.add_values(sq_diff_matrix);
    }
    push_back_feature("squared difference of " + feature_name, sq_diff_mmmv);
}

void TrackingFeatureExtractor::compute_sq_accel_features(
    ConstTraxelRefVectors& track_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator sq_accel_mmmv;
    sq_accel_mmmv.set_max(0.0);

    // for each track:
    for(auto track : track_traxels)
    {
        // only compute accelerations if track is longer than 2 elements
        if(track.size() < 3)
            continue;

        // extract features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(track, feature_matrix);

        // compute the squared norm of the accelerations of the column vectors
        FeatureMatrix sq_accel_matrix;
        sq_curve_calc_ptr_->calculate(feature_matrix, sq_accel_matrix);

        // add values to min/max/mean/var calculator
        sq_accel_mmmv.add_values(sq_accel_matrix);
    }
    push_back_feature("squared acceleration of " + feature_name, sq_accel_mmmv);
}

void TrackingFeatureExtractor::compute_angle_features(
    ConstTraxelRefVectors& track_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator angle_mmmv;
    angle_mmmv.set_max(0.0);

    // for each track:
    for(auto track : track_traxels)
    {
        // only compute angles in track if track is longer than 2 elements
        if(track.size() < 3)
            continue;

        // extract features
        FeatureMatrix feature_matrix;
        TraxelsFeaturesIdentity feature_extractor(feature_name);
        feature_extractor.extract(track, feature_matrix);

        // compute for all triples of features the angle of change of direction
        FeatureMatrix angles;
        angle_cos_calc_ptr_->calculate(feature_matrix, angles);

        angle_mmmv.add_values(angles);
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
    track_length_mmmv.set_max(0.0);

    for (auto track : track_traxels)
    {
        track_length_mmmv.add_value(static_cast<double>(track.size()));
    }
    push_back_feature("track length", track_length_mmmv);
}

void TrackingFeatureExtractor::compute_track_id_outlier(
    ConstTraxelRefVectors& track_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator id_out_mmmv;
    id_out_mmmv.set_max(0.0);
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
        }
    }
    push_back_feature("track outlier share of " + feature_name, id_out_mmmv);
}

void TrackingFeatureExtractor::compute_track_diff_outlier(
    ConstTraxelRefVectors& track_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator diff_out_mmmv;
    diff_out_mmmv.set_max(0.0);
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
        }
    }
    push_back_feature(
        "track outlier share of " + feature_name + " diff",
        diff_out_mmmv);
}

void TrackingFeatureExtractor::compute_division_sq_diff_features(
    ConstTraxelRefVectors& div_traxels,
    std::string feature_name)
{
    MinMaxMeanVarCalculator sq_diff_mmmv;
    sq_diff_mmmv.set_max(0.0);

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
        FeatureMatrix outlier_mat;
        mvn_outlier_calc_ptr_->calculate(sq_diff_matrix, outlier_mat);
        outlier = outlier_mat(0,0);
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
    child_decel_mmmv.set_max(0.0);

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
        FeatureMatrix outlier_mat;
        mvn_outlier_calc_ptr_->calculate(child_decel_mat, outlier_mat);
        outlier = outlier_mat(0,0);
    }
    push_back_feature("Outlier in child " + feature_name + " decelerations", outlier);
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
        TraxelsFeaturesIdentity position_extractor("com");
        position_extractor.extract(appearance, position);
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
