#ifndef CLS12
#define CLS12

#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace CLASS_ {
  class class1 {
  protected:
    /* data */
    py::list vals;


  public:
    class1 (py::list vals);
    int add ();
    virtual int subtract_ () = 0;
  };

  class class2: public class1 {
  private:
    /* data */

  public:
    class2 (py::list vals);
    int subtract_ ();
  };
}

#endif
