#include "pgmlink/tracking_feature_extractor.h"
#include <boost/algorithm/string.hpp>

namespace pgmlink {
namespace features {

TrackingFeatureExtractor::TrackingFeatureExtractor(HypothesesGraph &graph):
    graph_(graph)
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
    compute_velocity_features();
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

void TrackingFeatureExtractor::compute_velocity_features()
{
    // extract all tracks
    TrackTraxels track_extractor;
    std::vector<ConstTraxelRefVector> track_traxels = track_extractor(graph_);

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
        TraxelsFeaturesIdentity position_extractor("com");
        FeatureMatrix positions;
        position_extractor.extract(track, positions);

        // compute velocity vector's squared magnitude for all pairs of positions
        SquaredDiffCalculator velocity_calculator;
        FeatureMatrix velocities;
        velocity_calculator.calculate(positions, velocities);

        // compute per track min/max/mean of velocity
        MinCalculator<0> min_calculator;
        FeatureMatrix min_velocity;
        min_calculator.calculate(velocities, min_velocity);
        min_squared_velocity = std::min(min_squared_velocity, double(min_velocity(0,0)));

        MaxCalculator<0> max_calculator;
        FeatureMatrix max_velocity;
        max_calculator.calculate(velocities, max_velocity);
        max_squared_velocity = std::max(max_squared_velocity, double(max_velocity(0,0)));

        MeanCalculator<0> mean_calculator;
        FeatureMatrix mean_velocity;
        mean_calculator.calculate(velocities, mean_velocity);

        // accumulate all velocities
        SumCalculator<0> sum_calculator;
        FeatureMatrix sum_velocity;
        sum_calculator.calculate(velocities, sum_velocity);

        sum_of_squared_velocities += sum_velocity(0,0);
        num_velocity_entries += track.size() - 1;
    }

    double mean_squared_velocity = sum_of_squared_velocities / num_velocity_entries;

    joint_feature_vector_.push_back(mean_squared_velocity);
    feature_descriptions_.push_back("Mean of all velocities (squared)");

    joint_feature_vector_.push_back(min_squared_velocity);
    feature_descriptions_.push_back("Min of all velocities (squared)");

    joint_feature_vector_.push_back(max_squared_velocity);
    feature_descriptions_.push_back("Max of all velocities (squared)");
}

void TrackingFeatureExtractor::compute_size_difference_features()
{
    throw std::runtime_error("not yet implemented");
}

} // end namespace features
} // end namespace pgmlink
