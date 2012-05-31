#ifndef NDARRAY_tables_fwd_hpp_INCLUDED
#define NDARRAY_tables_fwd_hpp_INCLUDED

/**
 * @file ndarray/tables_fwd.hpp 
 *
 * @brief Forward declarations for ndarray Tables library.
 *
 *  \note This file is not included by the main "ndarray.hpp" header file.
 */

#include "ndarray_fwd.hpp"
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/remove_const.hpp>

namespace ndarray { namespace tables {

template <int N> struct Index {};

template <typename T, int N=0> struct Field;
template <typename T> class Layout;

template <typename T> class Row;
template <typename T> class Table;

namespace detail {

template <typename T, bool isConst = boost::is_const<T>::value> struct TraitsAccess;
template <typename Field_> struct FieldInfo;
template <typename T> class Columns;
template <typename T> class Iterator;

} // namespace detail

template <typename T>
struct Traits {
    typedef typename T::FieldSequence FieldSequence;
    typedef Row<T> Record;
    typedef Row<T> ConstRecord;
};

}} // namespace ndarray::tables

#endif // !NDARRAY_tables_fwd_hpp_INCLUDED
