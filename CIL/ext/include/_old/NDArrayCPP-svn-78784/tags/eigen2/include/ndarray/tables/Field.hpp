#ifndef NDARRAY_TABLES_Field_hpp_INCLUDED
#define NDARRAY_TABLES_Field_hpp_INCLUDED

#include "ndarray/tables_fwd.hpp"
#include "ndarray/Vector.hpp"

namespace ndarray { namespace tables {

template <typename T, int N>
struct Field {
    typedef T Element;
    typedef typename boost::mpl::int_<N>::type ND;
    typedef typename boost::mpl::bool_<(N==0)>::type IsScalar;

    std::string name;
    Vector<int,N> shape;
    int offset;

    template <typename T1, int N1>
    bool operator==(Field<T1,N1> const & other) const { return false; }

    template <typename T1, int N1>
    bool operator!=(Field<T1,N1> const & other) const { return true; }

    bool operator==(Field<T,N> const & other) const {
        return name == other.name && shape == other.shape 
            && (offset == other.offset || offset < 0 || other.offset < 0);
    }

    bool operator!=(Field<T,N> const & other) const {
        return !this->operator==(other);
    }

    Field() : name(), shape(), offset(-1) {}

    Field(Field const & other) : name(other.name), shape(other.shape), offset(other.offset) {}

    Field & operator=(Field const & other) {
        if (&other != this) {
            name = other.name;
            shape = other.shape;
            offset = other.offset;
        }
        return *this;
    }
};

}} // namespace ndarray::tables

#endif // !NDARRAY_TABLES_Field_hpp_INCLUDED
