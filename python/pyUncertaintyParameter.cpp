#define PY_ARRAY_UNIQUE_SYMBOL pgmlink_pyarray
#define NO_IMPORT_ARRAY

#include <vector>
#include <boost/python.hpp>
#include <boost/utility.hpp>

#include "../include/pgmlink/uncertaintyParameter.h"

using namespace pgmlink;
using namespace boost::python;

void export_uncertaintyParameter() {
    class_< UncertaintyParameter >("UncertaintyParameter")
      .def(init<int,DistrId,std::vector<double> >(
	   args("nOI","dI","dP")))
    ;
    class_< std::vector<double> >("VectorOfDoubles")
              .def(init<>())
              .def("append",&std::vector<double>::push_back)
              ;

}
