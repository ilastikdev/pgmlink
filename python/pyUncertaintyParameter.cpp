#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

#include <vector>
#include <boost/python.hpp>
#include <boost/utility.hpp>

#include "../include/pgmlink/uncertaintyParameter.h"

using namespace pgmlink;
using namespace boost::python;

void export_uncertaintyParameter()
{
    enum_<DistrId>("DistrId")
    .value("GaussianPertubation", Gaussian)
    .value("PerturbAndMAP", PerturbAndMAP)
    .value("DiverseMbest", DiverseMbest)
    .value("MbestCPLEX", MbestCPLEX)
    .value("ClassifierUncertainty", ClassifierUncertainty)
    ;

    class_< UncertaintyParameter >("UncertaintyParameter")
    .def(init<int, DistrId, std::vector<double> >(
             args("number_of_iterations", "distribution_id", "distribution_parameters")))
    ;

}
