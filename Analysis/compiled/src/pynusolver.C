#include <boost/python.hpp>
#include "boost/python/numpy.hpp"
#include "NuSolveLJ.h"
namespace py = boost::python;
namespace np = boost::python::numpy;

BOOST_PYTHON_MODULE(pynusolver) {
  Py_Initialize();
  np::initialize();
  def("run_nu_solver", &run_nu_solver);
}
