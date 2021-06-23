cmake_minimum_required(VERSION 3.10)

# set the project name and version
project(scnnoise VERSION 0.0.1)

# Set source directory
set(SOURCE_DIR "src/scnnoise")
# Tell CMake that headers are also in SOURCE_DIR
include_directories(${SOURCE_DIR})
set(SOURCES "${SOURCE_DIR}/graph.cpp" "${SOURCE_DIR}/GRN.cpp"
"${SOURCE_DIR}/graph_dependency.cpp" "${SOURCE_DIR}/scnnoise.cpp"
"${SOURCE_DIR}/gillespieSSA.cpp", "${SOURCE_DIR}/gillespieSDM.cpp")

# Generate Python module
add_subdirectory(lib/pybind11)
pybind11_add_module(snnoise ${SOURCES} "${SOURCE_DIR}/bindings.cpp")

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(cmake_example PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
