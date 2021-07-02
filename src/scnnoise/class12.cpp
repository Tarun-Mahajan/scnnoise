#include <pybind11/pybind11.h>
#include <iostream>
#include <vector>
#include <random>
#include "class12.hpp"

namespace py = pybind11;

namespace CLASS_ {
  class1::class1 (py::list vals) {
    std::cout << "Address for vals = " << &vals << std::endl;
    std::cout << "Address for val1 = " << &this->vals << std::endl;
    this->vals = vals;

    typedef std::mt19937 RNG;
    std::random_device rd;
    std::vector<std::uint_least32_t> rd_seeds = {rd(), rd(), rd(), rd()};
    std::seed_seq sd(rd_seeds.begin(), rd_seeds.end());
    RNG generator{sd};
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double rand_num = distribution(generator);
    std::cout << "Random number = " << rand_num << std::endl;
  }

  int class1::add () {
    int add_out = 0;
    for (auto &item : vals) {
      std::cout << "Address for element = " << &item << std::endl;
      add_out += item.attr("__str__")().cast<int>();
    }
    return add_out;
  }

  class2::class2 (py::list vals): class1(vals) {
  }

  int class2::subtract_ () {
    int tmp = add();
    return tmp - vals[0].attr("__str__")().cast<int>();
  }
}
