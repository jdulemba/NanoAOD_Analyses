#include <pybind11/pybind11.h>
#include "smoothing.h"

namespace py = pybind11;

////init.
PYBIND11_MODULE(pysmoother, m) {
  m.def("py_smoothhist", &py_smoothhist);
}
