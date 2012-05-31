#ifndef NDARRAY_initialization_hpp_INCLUDED
#define NDARRAY_initialization_hpp_INCLUDED

/** 
 *  \file ndarray/initialization.hpp @brief Construction functions for array.
 */

#include "ndarray/Array.hpp"
#include "ndarray/Manager.hpp"

namespace ndarray {
namespace detail {

struct NullOwner {};

template <int N, typename Derived>
class Initializer {
public:

    template <typename T, int C>
    operator Array<T,N,C> () const {
        return static_cast<Derived const *>(this)->template apply< Array<T,N,C> >();
    }

    template <typename T, int C>
    operator ArrayRef<T,N,C> () const {
        return static_cast<Derived const *>(this)->template apply< ArrayRef<T,N,C> >();
    }

};

template <int N>
class SimpleInitializer : public Initializer< N, SimpleInitializer<N> > {
public:

    template <typename Target>
    Target apply() const {
        typedef detail::ArrayAccess< Target > Access;
        typedef typename Access::Core Core;
        typedef typename Access::Element Element;
        int total = _shape.product();
        std::pair<Manager::Ptr,Element*> p = SimpleManager<Element>::allocate(total);
        return Access::construct(p.second, Core::create(_shape, p.first));
    }

    explicit SimpleInitializer(Vector<int,N> const & shape) : _shape(shape) {}

private:
    Vector<int,N> _shape;
};

template <typename T, int N, typename Owner>
class ExternalInitializer : public Initializer< N, ExternalInitializer<T,N,Owner> > {
public:

    template <typename Target>
    Target apply() const {
        typedef detail::ArrayAccess< Target > Access;
        typedef typename Access::Core Core;
        typedef typename Access::Element Element;
        Manager::Ptr manager;
        if (!boost::is_same<Owner,NullOwner>::value) {
            manager = ExternalManager<Owner>::make(_owner);
        }
        return Access::construct(_data, Core::create(_shape, _strides, manager));
    }

    ExternalInitializer(
        T * data, 
        Vector<int,N> const & shape,
        Vector<int,N> const & strides,
        Owner const & owner
    ) : _data(data), _owner(owner), _shape(shape), _strides(strides) {}

private:
    T * _data;
    Owner _owner;
    Vector<int,N> _shape;
    Vector<int,N> _strides;
};

} // namespace detail

/// @addtogroup MainGroup
/// @{

/** 
 *  @brief Create an expression that allocates uninitialized memory for an array.
 *
 *  @returns A temporary object convertible to an Array with fully contiguous row-major strides.
 */
template <int N>
inline detail::SimpleInitializer<N> allocate(Vector<int,N> const & shape) {
    return detail::SimpleInitializer<N>(shape); 
}

/** 
 *  @brief Create an expression that allocates uninitialized memory for a 1-d array.
 *
 *  @returns A temporary object convertible to an Array with fully contiguous row-major strides.
 */
inline detail::SimpleInitializer<1> allocate(int n) {
    return detail::SimpleInitializer<1>(ndarray::makeVector(n)); 
}

/** 
 *  @brief Create an expression that allocates uninitialized memory for a 2-d array.
 *
 *  @returns A temporary object convertible to an Array with fully contiguous row-major strides.
 */
inline detail::SimpleInitializer<2> allocate(int n1, int n2) {
    return detail::SimpleInitializer<2>(ndarray::makeVector(n1, n2)); 
}

/** 
 *  @brief Create an expression that allocates uninitialized memory for a 3-d array.
 *
 *  @returns A temporary object convertible to an Array with fully contiguous row-major strides.
 */
inline detail::SimpleInitializer<3> allocate(int n1, int n2, int n3) {
    return detail::SimpleInitializer<3>(ndarray::makeVector(n1, n2, n3)); 
}

/** 
 *  @brief Create a new Array by copying an Expression.
 */
template <typename Derived>
inline ArrayRef<typename boost::remove_const<typename Derived::Element>::type, 
                Derived::ND::value, Derived::ND::value>
copy(ExpressionBase<Derived> const & expr) {
    ArrayRef<typename boost::remove_const<typename Derived::Element>::type, 
        Derived::ND::value,Derived::ND::value> r(
            allocate(expr.getShape())
        );
    r = expr;
    return r;
}

/// @brief An enumeration for stride computation standards.
enum DataOrderEnum { ROW_MAJOR=1, COLUMN_MAJOR=2 };

/// @brief Compute row- or column-major strides for the given shape.
template <int N>
Vector<int,N> computeStrides(Vector<int,N> const & shape, DataOrderEnum order=ROW_MAJOR) {
    Vector<int,N> r(1);
    if (order == ROW_MAJOR) {
        for (int n=N-1; n > 0; --n) r[n-1] = r[n] * shape[n];
    } else {
        for (int n=1; n < N; ++n) r[n] = r[n-1] * shape[n-1];
    }
    return r;
}

/** 
 *  @brief Create an expression that initializes an Array with externally allocated memory.
 *
 *  No checking is done to ensure the shape, strides, and data pointers are sensible.
 *
 *  @param[in] data     A raw pointer to the first element of the Array.
 *  @param[in] shape    A Vector of dimensions for the new Array.
 *  @param[in] strides  A Vector of strides for the new Array.
 *  @param[in] owner    A copy-constructable object with an internal reference count
 *                      that owns the memory pointed at by 'data'.
 *
 *  @returns A temporary object convertible to an Array.
 */
template <typename T, int N, typename Owner>
inline detail::ExternalInitializer<T,N,Owner> external(
    T * data,
    Vector<int,N> const & shape,
    Vector<int,N> const & strides,
    Owner const & owner
) {
    return detail::ExternalInitializer<T,N,Owner>(data, shape, strides, owner);
}

/** 
 *  @brief Create an expression that initializes an Array with externally allocated memory.
 *
 *  No checking is done to ensure the shape, strides, and data pointers are sensible.  Memory will not
 *  be managed at all; the user must ensure the data pointer remains valid for the lifetime of the array.
 *
 *  @param[in] data     A raw pointer to the first element of the Array.
 *  @param[in] shape    A Vector of dimensions for the new Array.
 *  @param[in] strides  A Vector of strides for the new Array.
 *
 *  @returns A temporary object convertible to an Array.
 */
template <typename T, int N>
inline detail::ExternalInitializer<T,N,detail::NullOwner> external(
    T * data,
    Vector<int,N> const & shape,
    Vector<int,N> const & strides
) {
    return detail::ExternalInitializer<T,N,detail::NullOwner>(data, shape, strides, detail::NullOwner());
}

/** 
 *  @brief Create an expression that initializes an Array with externally allocated memory.
 *
 *  No checking is done to ensure the shape and data pointers are sensible.
 *
 *  @param[in] data     A raw pointer to the first element of the Array.
 *  @param[in] shape    A Vector of dimensions for the new Array.
 *  @param[in] order    Whether the strides are row- or column-major.
 *  @param[in] owner    A copy-constructable object with an internal reference count
 *                      that owns the memory pointed at by 'data'.
 *
 *  @returns A temporary object convertible to an Array.
 */
template <typename T, int N, typename Owner>
inline detail::ExternalInitializer<T,N,Owner> external(
    T * data,
    Vector<int,N> const & shape,
    DataOrderEnum order,
    Owner const & owner
) {
    return detail::ExternalInitializer<T,N,Owner>(data, shape, computeStrides(shape, order), owner);
}

/** 
 *  @brief Create an expression that initializes an Array with externally allocated memory.
 *
 *  No checking is done to ensure the shape and data pointers are sensible.  Memory will not
 *  be managed at all; the user must ensure the data pointer remains valid for the lifetime of the array.
 *
 *  @param[in] data     A raw pointer to the first element of the Array.
 *  @param[in] shape    A Vector of dimensions for the new Array.
 *  @param[in] order    Whether the strides are row- or column-major.
 *
 *  @returns A temporary object convertible to an Array.
 */
template <typename T, int N>
inline detail::ExternalInitializer<T,N,detail::NullOwner> external(
    T * data,
    Vector<int,N> const & shape,
    DataOrderEnum order = ROW_MAJOR
) {
    return detail::ExternalInitializer<T,N,detail::NullOwner>(
        data, shape, computeStrides(shape, order), detail::NullOwner()
    );
}

/// @}

} // namespace ndarray

#endif // !NDARRAY_initialization_hpp_INCLUDED
