#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

// stl
// temp
#include <iostream>
#include <string>
#include <sstream>

// pgmlink
#include "pgmlink/features/feature.h"
#include "pgmlink/log.h"
#include "pgmlink/tracking_evaluation.h"

// vigra
#include <vigra/numpy_array.hxx>
#include <vigra/numpy_array_converters.hxx>

// boost
#include <boost/python.hpp>
#include <boost/python/return_internal_reference.hpp>


using namespace pgmlink;
using namespace boost::python;


// evaluation
template <int N, typename T, typename U>
typename IntersectCountMap<T, U>::type py_get_intersect_count(vigra::NumpyArray<N, T> image1, vigra::NumpyArray<N, U> image2)
{
    return get_intersect_count<N, T, U>(image1, image2);
}


template <typename T, typename U>
std::pair<IntersectCountType, IntersectCountType> py_calculate_intersect_union(const typename IntersectCountMap<T, U>::type& intersect_counts,
        T region1_label,
        U region2_label,
        IntersectCountType region1_count,
        IntersectCountType region2_count)
{
    return calculate_intersect_union(intersect_counts, region1_label, region2_label, region1_count, region2_count);
}


void export_evaluation()
{
    ////
    //// Evaluation
    ////
    class_<std::pair<IntersectCountType, IntersectCountType> >("CountPair")
    .def_readwrite("first", &std::pair<IntersectCountType, IntersectCountType>::first)
    .def_readwrite("second", &std::pair<IntersectCountType, IntersectCountType>::second)
    ;

    class_<IntersectCountMap<>::type >("IntersectCountMap");

    def("getIntersectCount", &py_get_intersect_count<2, LabelType, LabelType>);
    def("getIntersectCount", &py_get_intersect_count<3, LabelType, LabelType>);
    def("getIntersectCount", &py_get_intersect_count<4, LabelType, LabelType>);
    def("getIntersectCount", &py_get_intersect_count<2, int, int>);
    def("getIntersectCount", &py_get_intersect_count<3, int, int>);
    def("getIntersectCount", &py_get_intersect_count<4, int, int>);
    def("getIntersectCount", &py_get_intersect_count<2, long, long>);
    def("getIntersectCount", &py_get_intersect_count<3, long, long>);
    def("getIntersectCount", &py_get_intersect_count<4, long, long>);
    def("getIntersectCount", &py_get_intersect_count<2, short, short>);
    def("getIntersectCount", &py_get_intersect_count<3, short, short>);
    def("getIntersectCount", &py_get_intersect_count<4, short, short>);

    def("calculateIntersectUnion", &py_calculate_intersect_union<LabelType, LabelType>);

}
