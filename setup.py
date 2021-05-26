from glob import glob
from setuptools import setup

# Available at setup time due to pyproject.toml
from pybind11.setup_helpers import Pybind11Extension, build_ext, ParallelCompile
from pybind11 import get_cmake_dir

# Optional multithreaded build
# ParallelCompile("NPY_NUM_BUILD_JOBS").install()

import sys

__version__ = "0.0.1"

# The main interface is through Pybind11Extension.
# * You can add cxx_std=11/14/17, and then build_ext can be removed.
# * You can set include_pybind11=false to add the include directory yourself,
#   say from a submodule.
#
# Note:
#   Sort input source files if you glob sources to ensure bit-for-bit
#   reproducible builds (https://github.com/pybind/python_example/pull/53)

ext_modules = [
    Pybind11Extension("scnnoise",
        # sorted(glob("src/*.cpp")),
        sorted(glob(os.path.join('src', '*.cpp'))),
        # Example: passing in the version to the compiled code
        define_macros = [('VERSION_INFO', __version__)],
        ),
]

setup(
    name="scnnoise",
    version=__version__,
    author="Tarun Mahajan",
    author_email="tarunm3@illinois.edu",
    url="https://github.com/Tarun-Mahajan/scnnoise",
    description="Single cell Network-aware Noise Simulator for gene Expression",
    long_description="A simulator for generating single-cell RNA sequencing data \
    using the two-state model of gene expression. Further, the simulator also \
    takes account of gene expression noise propagation over gene regulatory \
    network while generating the simulated data.",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    ext_modules=ext_modules,
    extras_require={"test": "pytest"},
    # Currently, build_ext only provides an optional "highest supported C++
    # level" feature, but in the future it may provide more features.
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    options={"bdist_wheel": {"universal": True}}
)
