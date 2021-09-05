#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
// #include "../../lib/pybind11/pybind11/include/pybind11/pybind11.h"
// #include "../../lib/pybind11/pybind11/include/pybind11/stl.h"
#include <string>
#include "scnnoise.hpp"
#include "gillespieSSA.hpp"
#include "gillespieSDMnoCellCycle.hpp"
#include "gillespieSDMCellCycle.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

PYBIND11_MODULE(_scnnoise, m) {
    // py::bind_vector<std::vector<int>>(m, "VectorInt", py::buffer_protocol());
    py::class_<ScnnoiseInterface::scNNoiSE>(m, "scNNoiSE")
        .def("add_new_dependency_graph",
            &ScnnoiseInterface::scNNoiSE::add_new_dependency_graph)
        .def("set_simulation_params",
            &ScnnoiseInterface::scNNoiSE::set_simulation_params)
        .def("swap_rxn_rates",
            &ScnnoiseInterface::scNNoiSE::swap_rxn_rates);
    py::class_<ScnnoiseInterface::GillespieSSA,
              ScnnoiseInterface::scNNoiSE>(m, "GillespieSSA")
        .def("simulate", &ScnnoiseInterface::GillespieSSA::simulate);
    py::class_<ScnnoiseInterface::gillespieSDMnoCellCycle,
              ScnnoiseInterface::GillespieSSA>(m, "gillespieSDMnoCellCycle")
        .def(py::init<int, std::string,
            std::string, std::string, bool,
            std::string>());
    py::class_<ScnnoiseInterface::gillespieSDMCellCycle,
              ScnnoiseInterface::GillespieSSA>(m, "gillespieSDMCellCycle")
        .def(py::init<int, std::string,
            std::string, std::string, bool,
            std::string>())
        .def("set_cell_cycle_params",
            &ScnnoiseInterface::gillespieSDMCellCycle::set_cell_cycle_params)
        .def("set_dosage_compensation",
            &ScnnoiseInterface::gillespieSDMCellCycle::set_dosage_compensation);

#ifdef VERSION_INFO
 m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
 m.attr("__version__") = "dev";
#endif
}
