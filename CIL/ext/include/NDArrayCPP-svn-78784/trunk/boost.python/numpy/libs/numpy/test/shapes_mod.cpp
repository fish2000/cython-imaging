#include <boost/numpy.hpp>

namespace p = boost::python;
namespace np = boost::numpy;

np::ndarray reshape(np::ndarray old_array, p::tuple shape)
{
  np::ndarray local_shape =  old_array.reshape(shape);
  return local_shape;
}

BOOST_PYTHON_MODULE(shapes_mod)
{
  np::initialize();
  p::def("reshape", reshape);
}
