#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

// vigra
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

#include "pgmlink/extract_region_features.h"

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

void export_region_feature_extraction()
{
    using namespace boost::python;

    // uint32
    def("extract_region_features", vigra::registerConverters(&py_extract_region_features<3, float, vigra::UInt32>));
    def("extract_region_features", vigra::registerConverters(&py_extract_region_features<2, float, vigra::UInt32>));
}
