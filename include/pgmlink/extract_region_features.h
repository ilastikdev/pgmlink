#ifndef EXTRACT_REGION_FEATURES_H
#define EXTRACT_REGION_FEATURES_H

#include <vigra/accumulator.hxx>
#include <boost/shared_ptr.hpp>
#include "pgmlink/featurestore.h"
#include "pgmlink/log.h"
#include <iostream>

namespace pgmlink {
namespace features{

template<typename T>
void set_feature(FeatureMap& feature_map, const std::string& name, T value)
{
    feature_map[name].clear();
    feature_map[name].push_back(feature_type(value));
}

template<>
void set_feature(FeatureMap &feature_map, const std::string &name, vigra::TinyVector<double,2> value)
{
    feature_map[name].clear();
    feature_map[name].push_back(feature_type(value[0]));
    feature_map[name].push_back(feature_type(value[1]));
}

template<>
void set_feature(FeatureMap &feature_map, const std::string &name, vigra::TinyVector<double,3> value)
{
    feature_map[name].clear();
    feature_map[name].push_back(feature_type(value[0]));
    feature_map[name].push_back(feature_type(value[1]));
    feature_map[name].push_back(feature_type(value[2]));
}

template<>
void set_feature(FeatureMap &feature_map, const std::string &name, vigra::linalg::Matrix<double> value)
{
    feature_map[name].clear();
    for(auto it = value.begin(); it != value.end(); ++it)
    {
        feature_map[name].push_back(*it);
    }
}

///
/// Extract features from all regions in the given MultiArrayView and insert them into the corresponding traxels at the current timestep.
/// \return the maximal label id
///
template<int N, typename DataType, typename LabelType>
int extract_region_features(const vigra::MultiArrayView<N, DataType> data,
                             const vigra::MultiArrayView<N, LabelType> labels,
                             boost::shared_ptr<pgmlink::FeatureStore> fs,
                             unsigned int timestep)
{
    // extract features using vigra
    using namespace vigra::acc;

    typedef AccumulatorChainArray<vigra::CoupledArrays<N, DataType, LabelType>,
            Select< // what statistics to compute (same as in joint seg & track, but without quantiles atm)
                RegionCenter,
                Count,
                Variance,
                Sum,
                Mean,
                RegionRadii,
                Central< PowerSum<2> >,
                Central< PowerSum<3> >,
                Central< PowerSum<4> >,
                Kurtosis,
                Maximum,
                Minimum,
                RegionAxes,
                Skewness,
                Weighted<PowerSum<0> >,
                Coord< Minimum >,
                Coord< Maximum >,
                DataArg<1>,
                LabelArg<2> // where to look for data and region labels
            > >
            FeatureAccumulator;
    FeatureAccumulator a;
    a.ignoreLabel(0); // do not compute features for the background
    LOG(pgmlink::logDEBUG1) << "Beginning feature extraction for frame " << timestep;
    extractFeatures(data, labels, a);
    LOG(pgmlink::logDEBUG1) << "Finished feature extraction for frame " << timestep;

    LabelType label_min, label_max;
    labels.minmax(&label_min, &label_max);

    // insert features into featurestore
    for(size_t label = label_min; label < label_max; label++)
    {
        if(label == 0) // ignore background
            continue;

        // get respective feature map from FeatureStore
        pgmlink::FeatureMap& feature_map = fs->get_traxel_features(timestep, label);

        // add features
        set_feature(feature_map, "Mean", get<Mean>(a, label));
        set_feature(feature_map, "Sum", get<Sum>(a, label));
        set_feature(feature_map, "Variance", get<Variance>(a, label));
        set_feature(feature_map, "Count", get<Count>(a, label));
        set_feature(feature_map, "RegionRadii", get<RegionRadii>(a, label));
        set_feature(feature_map, "RegionCenter", get<RegionCenter>(a, label));
        set_feature(feature_map, "Coord< Maximum >", get<Coord< Maximum > >(a, label));
        set_feature(feature_map, "Coord< Minimum >", get<Coord< Minimum > >(a, label));
        set_feature(feature_map, "RegionAxes", get<RegionAxes>(a, label));
        set_feature(feature_map, "Kurtosis", get<Kurtosis>(a, label));
        set_feature(feature_map, "Minimum", get<Minimum>(a, label));
        set_feature(feature_map, "Maximum", get<Maximum>(a, label));
        set_feature(feature_map, "Skewness", get<Skewness>(a, label));
        set_feature(feature_map, "Central< PowerSum<2> >", get<Central< PowerSum<2> > >(a, label));
        set_feature(feature_map, "Central< PowerSum<3> >", get<Central< PowerSum<3> > >(a, label));
        set_feature(feature_map, "Central< PowerSum<4> >", get<Central< PowerSum<4> > >(a, label));
        set_feature(feature_map, "Weighted<PowerSum<0> >", get<Weighted<PowerSum<0> > >(a, label));
    }

    return label_max;
}

} // namespace features
} // namespace pgmlink

#endif // EXTRACT_REGION_FEATURES_H
