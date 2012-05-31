#ifndef NDARRAY_TABLES_DETAIL_Iterator_hpp_INCLUDED
#define NDARRAY_TABLES_DETAIL_Iterator_hpp_INCLUDED

#include "ndarray_fwd.hpp"
#include "ndarray/tables/detail/TraitsAccess.hpp"

#include <boost/iterator/iterator_facade.hpp>

namespace ndarray { namespace tables { namespace detail {

template <typename T>
class Iterator : public boost::iterator_facade< 
    Iterator<T>, 
    typename TraitsAccess<T>::Row_,
    boost::random_access_traversal_tag, 
    typename TraitsAccess<T>::Row_
    > 
{
public:

    typedef typename TraitsAccess<T>::Row_ Value;
    typedef Value Reference;

    Reference operator[](int n) const {
        return Reference(getRow()._n + n, getRow()._columns);
    }

    Reference const & operator*() const { return _ref; }

    Reference const * operator->() { return &_ref; }

    Iterator() : _ref() {}

    explicit Iterator(Reference const & ref) : _ref(ref) {}

    Iterator(Iterator const & other) : _ref(other._ref) {}

    template <typename U>
    Iterator(Iterator<U> const & other) : _ref(other._ref) {}
    
    Iterator & operator=(Iterator const & other) {
        if (&other != this) _ref = other._ref;
        return *this;
    }

    template <typename U>
    Iterator & operator=(Iterator<U> const & other) {
        if (&other != this) _ref = other._ref;
        return *this;
    }

private:

    friend class boost::iterator_core_access;

    template <typename U> friend class Iterator;

    Row<T> & getRow() { return _ref; }
    Row<T> const & getRow() const { return _ref; }

    Reference const & dereference() const { return _ref; }

    void increment() { ++getRow()._n; }
    void decrement() { --getRow()._n; }
    void advance(int n) { getRow()._n += n; }

    template <typename U>
    int distance_to(Iterator<U> const & other) const {
        return other.getRow()._n - getRow()._n; 
    }

    template <typename U>
    bool equal(Iterator<U> const & other) const {
        return other.getRow()._n == getRow()._n && other.getRow()._columns == other.getRow()._columns;
    }

    Reference _ref;
};


}}} // namespace ndarray::tables::detail

#endif // !NDARRAY_TABLES_DETAIL_Iterator_hpp_INCLUDED
