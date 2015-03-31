#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

// vigra
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

#include "pgmlink/features/extract_region_features.h"

#include <boost/python.hpp>

template<int N, typename DataType, typename LabelType>
int py_extract_region_features(
    const vigra::NumpyArray<N, DataType>& image,
    const vigra::NumpyArray<N, LabelType>& labels,
    boost::shared_ptr<pgmlink::FeatureStore> fs,
    const size_t timestep
)
{
    return pgmlink::features::extract_region_features<N, DataType, LabelType>(image, labels, fs, timestep);
}

template<int N, typename DataType, typename LabelType>
void py_extract_region_features_roi(
    const vigra::NumpyArray<N, DataType>& image,
    const vigra::NumpyArray<N, LabelType>& labels,
    const vigra::NumpyArray<1, vigra::Int64>& label_indices,
    int traxel_index_offset,
    const vigra::NumpyArray<1, vigra::Int64>& coord_offsets,
    boost::shared_ptr<pgmlink::FeatureStore> fs,
    const size_t timestep
)
{
    vigra::TinyVector<size_t, N> coord_offsets_tv;
    for (size_t idx = 0; idx < N; ++idx)
    {
        coord_offsets_tv[idx] = coord_offsets[idx];
    }

    std::vector<size_t> label_indices_std;
    for (auto it = label_indices.begin(); it != label_indices.end(); ++it)
    {
        label_indices_std.push_back(*it);
    }

    pgmlink::features::extract_region_features_roi<N, DataType, LabelType>(image,
            labels,
            label_indices_std,
            (unsigned int)traxel_index_offset,
            coord_offsets_tv,
            fs,
            (unsigned int)timestep);
}



void export_region_feature_extraction()
{
    using namespace boost::python;

    // uint32
    def("extract_region_features", vigra::registerConverters(&py_extract_region_features<3, float, vigra::UInt32>));
    def("extract_region_features", vigra::registerConverters(&py_extract_region_features<2, float, vigra::UInt32>));

    def("extract_region_features_roi", vigra::registerConverters(&py_extract_region_features_roi<3, float, vigra::UInt32>));
    def("extract_region_features_roi", vigra::registerConverters(&py_extract_region_features_roi<2, float, vigra::UInt32>));
}
