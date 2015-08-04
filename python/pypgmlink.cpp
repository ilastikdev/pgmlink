#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray

#include <boost/python.hpp>
#include <vigra/numpy_array.hxx>
#include <iostream>
#include <vector>
#include <complex>


//forward declarations
void export_field_of_view();
void export_uncertaintyParameter();
void export_hypotheses();
void export_track();
void export_traxels();
void export_gmm();
void export_region_feature_extraction();
void export_evaluation();

BOOST_PYTHON_MODULE( pgmlink )
{
    vigra::import_vigranumpy();
    export_field_of_view();
    export_uncertaintyParameter();
    export_hypotheses();
    export_track();
    export_traxels();
    export_gmm();
    export_region_feature_extraction();
    export_evaluation();
}
