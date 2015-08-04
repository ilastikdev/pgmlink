#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

#include <boost/python.hpp>
#include <boost/utility.hpp>

#include "../include/pgmlink/field_of_view.h"

using namespace pgmlink;
using namespace boost::python;

void set_time_bounds(FieldOfView& fov, int min_ts, int max_ts)
{
    assert(min_ts <= max_ts);
    auto lb = fov.lower_bound();
    auto ub = fov.upper_bound();
    fov.set_boundingbox(min_ts, lb[1], lb[2], lb[3], max_ts, ub[1], ub[2], ub[3]);
}

void export_field_of_view()
{
    class_< FieldOfView >("FieldOfView")
    .def(init<double, double, double, double, double, double, double, double>(
             args("lt", "lx", "ly", "lz", "ut", "ux", "uy", "uz")))
    .def("set_boundingbox", &FieldOfView::set_boundingbox, return_self<>())
    .def("set_time_bounds", &set_time_bounds)
    //.def("contains", &FieldOfView::contains)
    //.def("lower_bound", &FieldOfView::lower_bound, return_value_policy<copy_const_reference>())
    //.def("upper_bound", &FieldOfView::upper_bound, return_value_policy<copy_const_reference>())
    //.def("spatial_margin", &FieldOfView::spatial_margin)
    //.def("temporal_margin", &FieldOfView::temporal_margin)
    ;
}
