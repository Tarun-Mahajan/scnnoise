#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>
// #include "../../lib/pybind11/pybind11/include/pybind11/pybind11.h"
// #include "../../lib/pybind11/pybind11/include/pybind11/stl.h"
#include "scnnoise.hpp"
#include "gillespieSSA.hpp"
#include "gillespieSDM.hpp"

namespace py = pybind11;

PYBIND11_MODULE(scnnoise, m) {
    py::bind_vector<std::vector<int>>(m, "VectorInt", py::buffer_protocol());
    py::class_<ScnnoiseInterface::scNNoiSE>(m, "scNNoiSE")
        .def("add_gene_state", &ScnnoiseInterface::scNNoiSE::add_gene_state)
        .def("add_GRN_edge", &ScnnoiseInterface::scNNoiSE::add_GRN_edge)
        .def("add_dependency_edge", &ScnnoiseInterface::scNNoiSE::add_dependency_edge);
    py::class_<ScnnoiseInterface::GillespieSSA,
              ScnnoiseInterface::scNNoiSE>(m, "GillespieSSA")
        .def("simulate", &ScnnoiseInterface::simulate);
    py::class_<ScnnoiseInterface::GillespieSDM,
              ScnnoiseInterface::GillespieSSA>(m, "GillespieSDM")
        .def(py::init<int, int, const std::vector<int>, const std::vector<int> \
             double, bool, int>());
}
