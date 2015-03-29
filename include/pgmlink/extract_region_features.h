#ifndef EXTRACT_REGION_FEATURES_H
#define EXTRACT_REGION_FEATURES_H

#include <vigra/accumulator.hxx>
#include <boost/shared_ptr.hpp>
#include "pgmlink/featurestore.h"
#include "pgmlink/log.h"
#include <iostream>
#include <numeric>

namespace pgmlink
{
namespace features
{

template<typename T>
void set_feature(FeatureMap& feature_map, const std::string& name, T value)
{
    feature_map[name].clear();
    feature_map[name].push_back(feature_type(value));
}

template<>
void set_feature(FeatureMap &feature_map, const std::string &name, vigra::TinyVector<double, 2> value)
{
    feature_map[name].clear();
    feature_map[name].push_back(feature_type(value[0]));
    feature_map[name].push_back(feature_type(value[1]));
}

template<>
void set_feature(FeatureMap &feature_map, const std::string &name, vigra::TinyVector<double, 3> value)
{
    feature_map[name].clear();
    feature_map[name].push_back(feature_type(value[0]));
    feature_map[name].push_back(feature_type(value[1]));
    feature_map[name].push_back(feature_type(value[2]));
}

template<int N, typename T1, typename T2>
void set_feature_with_offset(FeatureMap &feature_map,
                             const std::string &name,
                             vigra::TinyVector<T1, N> value,
                             vigra::TinyVector<T2, N> offset)
{
    feature_map[name].clear();
    for(int i = 0; i < N; i++)
    {
        feature_map[name].push_back(feature_type(value[i] + offset[i]));
    }
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
/// Extract features from the selected regions in the given MultiArrayView
/// and insert them into the corresponding traxels (with offset) at the current timestep.
/// Coordinates are shifted by the given offset if this is only a ROI of the full image.
/// \return the maximal label id
///
template<int N, typename DataType, typename LabelType>
void extract_region_features_roi(const vigra::MultiArrayView<N, DataType>& data,
                                 const vigra::MultiArrayView<N, LabelType>& labels,
                                 const std::vector<size_t>& label_indices,
                                 unsigned int traxel_index_offset,
                                 const vigra::TinyVector<size_t, N>& coord_offsets,
                                 boost::shared_ptr<pgmlink::FeatureStore> fs,
                                 unsigned int timestep,
                                 bool verbose = true)
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

    // insert features into featurestore
    for(LabelType label : label_indices)
    {
        if(label == 0) // ignore background
        {
            continue;
        }

        // get respective feature map from FeatureStore
        pgmlink::FeatureMap& feature_map = fs->get_traxel_features(timestep, label + traxel_index_offset);

        // add features
        set_feature(feature_map, "Mean", get<Mean>(a, label));
        set_feature(feature_map, "Sum", get<Sum>(a, label));
        set_feature(feature_map, "Variance", get<Variance>(a, label));
        set_feature(feature_map, "Count", get<Count>(a, label));
        set_feature(feature_map, "RegionRadii", get<RegionRadii>(a, label));
        set_feature_with_offset(feature_map, "RegionCenter", get<RegionCenter>(a, label), coord_offsets);
        set_feature_with_offset(feature_map, "Coord< Maximum >", get<Coord< Maximum > >(a, label), coord_offsets);
        set_feature_with_offset(feature_map, "Coord< Minimum >", get<Coord< Minimum > >(a, label), coord_offsets);
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
}

///
/// Extract features from all regions in the given MultiArrayView
/// and insert them into the corresponding traxels at the current timestep.
/// \return the maximal label id
///
template<int N, typename DataType, typename LabelType>
int extract_region_features(const vigra::MultiArrayView<N, DataType>& data,
                            const vigra::MultiArrayView<N, LabelType>& labels,
                            boost::shared_ptr<pgmlink::FeatureStore> fs,
                            unsigned int timestep)
{
    // create list of all labels we want to extract features for
    LabelType label_min, label_max;
    labels.minmax(&label_min, &label_max);
    std::vector<size_t> label_indices(label_max + 1 - label_min);
    std::iota(label_indices.begin(), label_indices.end(), label_min);

    // coordinate offset is zero, use default constructor
    vigra::TinyVector<size_t, N> coord_offset;

    extract_region_features_roi<N, DataType, LabelType>(data, labels, label_indices, 0, coord_offset, fs, timestep, false);
    return label_max;
}

} // namespace features
} // namespace pgmlink

#endif // EXTRACT_REGION_FEATURES_H
