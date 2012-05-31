#ifndef BOOST_NUMPY_NDARRAY_HPP_INCLUDED
#define BOOST_NUMPY_NDARRAY_HPP_INCLUDED

/**
 *  @file boost/numpy/ndarray.hpp
 *  @brief Object manager and various utilities for numpy.ndarray.
 */

#include <boost/python.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_integral.hpp>
#include <boost/numpy/numpy_object_mgr_traits.hpp>
#include <boost/numpy/dtype.hpp>

#include <vector>

namespace boost 
{
namespace numpy 
{

/**
 *  @brief A boost.python "object manager" (subclass of object) for numpy.ndarray.
 *
 *  @todo This could have a lot more functionality (like boost::python::numeric::array).
 *        Right now all that exists is what was needed to move raw data between C++ and Python.
 */
class ndarray : public python::object 
{

  /**
   *  @brief An internal struct that's byte-compatible with PyArrayObject.
   *
   *  This is just a hack to allow inline access to this stuff while hiding numpy/arrayobject.h
   *  from the user.
   */
  struct array_struct 
  {
    PyObject_HEAD
    char * data;
    int nd;
    Py_intptr_t * shape;
    Py_intptr_t * strides;
    PyObject * base;
    PyObject * descr;
    int flags;
    PyObject * weakreflist;
  };
  
  /// @brief Return the held Python object as an array_struct.
  array_struct * get_struct() const { return reinterpret_cast<array_struct*>(this->ptr()); }

public:
  
  /**
   *  @brief Enum to represent (some) of Numpy's internal flags.
   *
   *  These don't match the actual Numpy flag values; we can't get those without including 
   *  numpy/arrayobject.h or copying them directly.  That's very unfortunate.
   *
   *  @todo I'm torn about whether this should be an enum.  It's very convenient to not
   *        make these simple integer values for overloading purposes, but the need to
   *        define every possible combination and custom bitwise operators is ugly.
   */
  enum bitflag 
  {
    NONE=0x0, C_CONTIGUOUS=0x1, F_CONTIGUOUS=0x2, V_CONTIGUOUS=0x1|0x2, 
    ALIGNED=0x4, WRITEABLE=0x8, BEHAVED=0x4|0x8,
    CARRAY_RO=0x1|0x4, CARRAY=0x1|0x4|0x8, CARRAY_MIS=0x1|0x8,
    FARRAY_RO=0x2|0x4, FARRAY=0x2|0x4|0x8, FARRAY_MIS=0x2|0x8,
    UPDATE_ALL=0x1|0x2|0x4, VARRAY=0x1|0x2|0x8, ALL=0x1|0x2|0x4|0x8
  };

  BOOST_PYTHON_FORWARD_OBJECT_CONSTRUCTORS(ndarray, object);

  /// @brief Return a view of the scalar with the given dtype.
  ndarray view(dtype const & dt) const;

  /// @brief Copy the scalar (deep for all non-object fields).
  ndarray copy() const;

  /// @brief Return the size of the nth dimension.
  int const shape(int n) const { return get_shape()[n]; }

  /// @brief Return the stride of the nth dimension.
  int const strides(int n) const { return get_strides()[n]; }
    
  /**
   *  @brief Return the array's raw data pointer.
   *
   *  This returns char so stride math works properly on it.  It's pretty much
   *  expected that the user will have to reinterpret_cast it.
   */
  char * get_data() const { return get_struct()->data; }

  /// @brief Return the array's data-type descriptor object.
  dtype get_dtype() const;
  
  /// @brief Return the object that owns the array's data, or None if the array owns its own data.
  python::object get_base() const;
  
  /// @brief Set the object that owns the array's data.  Use with care.
  void set_base(object const & base);
  
  /// @brief Return the shape of the array as an array of integers (length == get_nd()).
  Py_intptr_t const * get_shape() const { return get_struct()->shape; }
  
  /// @brief Return the stride of the array as an array of integers (length == get_nd()).
  Py_intptr_t const * get_strides() const { return get_struct()->strides; }
  
  /// @brief Return the number of array dimensions.
  int const get_nd() const { return get_struct()->nd; }
  
  /// @brief Return the array flags.
  bitflag const get_flags() const;
  
  /// @brief Reverse the dimensions of the array.
  ndarray transpose() const;
  
  /// @brief Eliminate any unit-sized dimensions.
  ndarray squeeze() const;
  
  /// @brief Equivalent to self.reshape(*shape) in Python.
  ndarray reshape(python::tuple const & shape) const;
  
  /**
   *  @brief If the array contains only a single element, return it as an array scalar; otherwise return
   *         the array.
   *
   *  @internal This is simply a call to PyArray_Return();
   */
  python::object scalarize() const;
};

/**
 *  @brief Construct a new array with the given shape and data type, with data initialized to zero.
 */
ndarray zeros(python::tuple const & shape, dtype const & dt);
ndarray zeros(int nd, Py_intptr_t const * shape, dtype const & dt);

/**
 *  @brief Construct a new array with the given shape and data type, with data left uninitialized.
 */
ndarray empty(python::tuple const & shape, dtype const & dt);
ndarray empty(int nd, Py_intptr_t const * shape, dtype const & dt);

/**
 *  @brief Construct a new array from an arbitrary Python sequence.
 *
 *  @todo This does't seem to handle ndarray subtypes the same way that "numpy.array" does in Python.
 */
ndarray array(python::object const & obj);
ndarray array(python::object const & obj, dtype const & dt);

namespace detail 
{

ndarray from_data_impl(void * data,
		       dtype const & dt,
		       std::vector<Py_intptr_t> const & shape,
		       std::vector<Py_intptr_t> const & strides,
		       python::object const & owner,
		       bool writeable);

template <typename Container>
ndarray from_data_impl(void * data,
		       dtype const & dt,
		       Container shape,
		       Container strides,
		       python::object const & owner,
		       bool writeable,
		       typename boost::enable_if< boost::is_integral<typename Container::value_type> >::type * enabled = NULL)
{
  std::vector<Py_intptr_t> shape_(shape.begin(),shape.end());
  std::vector<Py_intptr_t> strides_(strides.begin(), strides.end());
  return from_data_impl(data, dt, shape_, strides_, owner, writeable);    
}

ndarray from_data_impl(void * data,
		       dtype const & dt,
		       python::object const & shape,
		       python::object const & strides,
		       python::object const & owner,
		       bool writeable);

} // namespace boost::numpy::detail

/**
 *  @brief Construct a new ndarray object from a raw pointer.
 *
 *  @param[in] data    Raw pointer to the first element of the array.
 *  @param[in] dt      Data type descriptor.  Often retrieved with dtype::get_builtin().
 *  @param[in] shape   Shape of the array as STL container of integers; must have begin() and end().
 *  @param[in] strides Shape of the array as STL container of integers; must have begin() and end().
 *  @param[in] owner   An arbitray Python object that owns that data pointer.  The array object will
 *                     keep a reference to the object, and decrement it's reference count when the
 *                     array goes out of scope.  Pass None at your own peril.
 *
 *  @todo Should probably take ranges of iterators rather than actual container objects.
 */
template <typename Container>
inline ndarray from_data(void * data,
			 dtype const & dt,
			 Container shape,
			 Container strides,
			 python::object const & owner)
{
  return numpy::detail::from_data_impl(data, dt, shape, strides, owner, true);
}    

/**
 *  @brief Construct a new ndarray object from a raw pointer.
 *
 *  @param[in] data    Raw pointer to the first element of the array.
 *  @param[in] dt      Data type descriptor.  Often retrieved with dtype::get_builtin().
 *  @param[in] shape   Shape of the array as STL container of integers; must have begin() and end().
 *  @param[in] strides Shape of the array as STL container of integers; must have begin() and end().
 *  @param[in] owner   An arbitray Python object that owns that data pointer.  The array object will
 *                     keep a reference to the object, and decrement it's reference count when the
 *                     array goes out of scope.  Pass None at your own peril.
 *
 *  This overload takes a const void pointer and sets the "writeable" flag of the array to false.
 *
 *  @todo Should probably take ranges of iterators rather than actual container objects.
 */
template <typename Container>
inline ndarray from_data(void const * data,
			 dtype const & dt,
			 Container shape,
			 Container strides,
			 python::object const & owner)
{
  return numpy::detail::from_data_impl(const_cast<void*>(data), dt, shape, strides, owner, false);
}    

/**
 *  @brief Transform an arbitrary object into a numpy array with the given requirements.
 *
 *  @param[in] obj     An arbitrary python object to convert.  Arrays that meet the requirements
 *                     will be passed through directly.
 *  @param[in] dt      Data type descriptor.  Often retrieved with dtype::get_builtin().
 *  @param[in] nd_min  Minimum number of dimensions.
 *  @param[in] nd_max  Maximum number of dimensions.
 *  @param[in] flags   Bitwise OR of flags specifying additional requirements.
 */
ndarray from_object(python::object const & obj, dtype const & dt, 
                    int nd_min, int nd_max, ndarray::bitflag flags=ndarray::NONE);

inline ndarray from_object(python::object const & obj, dtype const & dt, 
                           int nd, ndarray::bitflag flags=ndarray::NONE)
{
  return from_object(obj, dt, nd, nd, flags);
}

inline ndarray from_object(python::object const & obj, dtype const & dt, ndarray::bitflag flags=ndarray::NONE)
{
  return from_object(obj, dt, 0, 0, flags);
}

ndarray from_object(python::object const & obj, int nd_min, int nd_max, 
                    ndarray::bitflag flags=ndarray::NONE);

inline ndarray from_object(python::object const & obj, int nd, ndarray::bitflag flags=ndarray::NONE)
{
  return from_object(obj, nd, nd, flags);
}

inline ndarray from_object(python::object const & obj, ndarray::bitflag flags=ndarray::NONE)
{
  return from_object(obj, 0, 0, flags);
}

inline ndarray::bitflag operator|(ndarray::bitflag a, ndarray::bitflag b)
{
  return ndarray::bitflag(int(a) | int(b));
}

inline ndarray::bitflag operator&(ndarray::bitflag a, ndarray::bitflag b)
{
  return ndarray::bitflag(int(a) & int(b));
}

} // namespace boost::numpy

namespace python
{
namespace converter 
{

NUMPY_OBJECT_MANAGER_TRAITS(numpy::ndarray);

} // namespace boost::python::converter
} // namespace boost::python
} // namespace boost

#endif // !BOOST_NUMPY_NDARRAY_HPP_INCLUDED
