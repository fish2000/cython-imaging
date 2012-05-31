#ifndef NDARRAY_TABLES_DETAIL_TraitsAccess_hpp_INCLUDED
#define NDARRAY_TABLES_DETAIL_TraitsAccess_hpp_INCLUDED

#include "ndarray/tables_fwd.hpp"
#include "ndarray/tables/detail/FieldInfo.hpp"

#include <boost/intrusive_ptr.hpp>

namespace ndarray { namespace tables { namespace detail {

template <typename T>
struct TraitsAccess<T,false> {
    typedef typename Traits<T>::FieldSequence FieldSequence;
    typedef Layout<T> Layout_;
    typedef Row<T> Row_;
    typedef Table<T> Table_;
    typedef unsigned char Raw;
    typedef detail::Columns<T> Columns_;
    typedef boost::intrusive_ptr< detail::Columns<T> > ColumnsPtr;

    template <int N>
    struct Fields {
        typedef typename boost::fusion::result_of::value_at_c<FieldSequence,N>::type Type;
        typedef FieldInfo<Type> Info;
        typedef typename Info::template Qualified<false> Qualified;
    };

};

template <typename T>
struct TraitsAccess<T,true> {
    typedef typename boost::remove_const<T>::type U;
    typedef typename Traits<T>::FieldSequence FieldSequence;
    typedef Layout<U> Layout_;
    typedef Row<T> Row_;
    typedef Table<T> Table_;
    typedef unsigned char const Raw;
    typedef detail::Columns<U> Columns_;
    typedef boost::intrusive_ptr< detail::Columns<U> > ColumnsPtr;

    template <int N>
    struct Fields {
        typedef typename boost::fusion::result_of::value_at_c<FieldSequence,N>::type Type;
        typedef FieldInfo<Type> Info;
        typedef typename Info::template Qualified<true> Qualified;
    };
};

}}} // namespace ndarray::tables::detail

#endif // !NDARRAY_TABLES_DETAIL_TraitsAccess_hpp_INCLUDED
