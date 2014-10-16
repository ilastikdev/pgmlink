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

    /// Append features for this solution to the given file.
    /// If file does not exist, create it.
    /// Comments are ignored, and will not be copied to the edited file
    void append_feature_vector_to_file(const std::string &filename);

private:
    /**
     * methods that compute each feature
     */
    void compute_velocity_features();
    void compute_size_difference_features();
    // TODO: add many more

private:
    JointFeatureVector joint_feature_vector_;
    FeatureDescription feature_descriptions_;
    HypothesesGraph& graph_;
};

} // end namespace features
} // end namespace pgmlink

#endif // TRACKING_FEATURE_EXTRACTOR_H
