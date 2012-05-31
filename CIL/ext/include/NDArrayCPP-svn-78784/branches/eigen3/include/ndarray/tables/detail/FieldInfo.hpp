#ifndef NDARRAY_TABLES_DETAIL_FieldInfo_hpp_INCLUDED
#define NDARRAY_TABLES_DETAIL_FieldInfo_hpp_INCLUDED

#include "ndarray/tables_fwd.hpp"

namespace ndarray { namespace tables { namespace detail {

template <typename Field>
struct FieldInfo {};

template <typename T, int N>
struct FieldInfo< Field<T,N> > {
    typedef typename boost::mpl::int_<N>::type ND;
    
    typedef Field<T,N> Field_;

    typedef typename boost::mpl::bool_<(N==0)>::type IsScalar;

    typedef ndarray::Array<T,N+1,N> ColumnValue;

    template <bool isConst>
    struct Qualified {
        typedef typename boost::mpl::if_c<isConst, T const, T>::type Element;
        typedef typename boost::mpl::if_<IsScalar, Element &, ndarray::ArrayRef<Element,N,N> >::type RowValue;
        typedef ndarray::ArrayRef<Element,N+1,N> TableValue;
    };

};

template <typename T, int N>
struct FieldInfo< Field<T,N> & > : public FieldInfo< Field<T,N> > {};

template <typename T, int N>
struct FieldInfo< Field<T,N> const & > : public FieldInfo< Field<T,N> > {};

}}} // namespace ndarray::tables::detail

#endif // !NDARRAY_TABLES_DETAIL_FieldInfo_hpp_INCLUDED
