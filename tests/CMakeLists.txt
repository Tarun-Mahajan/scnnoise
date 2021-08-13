cmake_minimum_required(VERSION 3.7...3.21)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

# set the project name
project(scnnoise LANGUAGES CXX)

# Set source directory
set(SOURCE_DIR "../src/scnnoise")

# Tell CMake that headers are also in SOURCE_DIR
set(SOURCES "${SOURCE_DIR}/graph.cpp" "${SOURCE_DIR}/graph_derived.cpp"
"${SOURCE_DIR}/utils.cpp" )

# Set test directory
# SET(TEST_DIR "tests")
# SET(TESTS ${SOURCES} "test_graph.cpp")
SET(TESTS "test_main.cpp" "test_GRN.cpp" "test_utils.cpp" ${SOURCES})

# Generate a test executable for graph-related classes
# add_executable("scnnoise_graph_test" ${TESTS})
# target_include_directories("scnnoise_graph_test" PUBLIC ${SOURCE_DIR})
SET(INCLUDE_DIRS ${SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

# Generate a test executable
find_package(Catch2 REQUIRED)
add_executable("scnnoise_graph_test" ${TESTS})
target_include_directories("scnnoise_graph_test" PUBLIC ${INCLUDE_DIRS})
target_link_libraries("scnnoise_graph_test" PRIVATE Catch2::Catch2WithMain)
