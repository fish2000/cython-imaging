#ifndef NDARRAY_TABLES_Layout_hpp_INCLUDED
#define NDARRAY_TABLES_Layout_hpp_INCLUDED

#include <boost/fusion/sequence/intrinsic/size.hpp>
#include <boost/fusion/sequence/intrinsic/value_at.hpp>
#include <boost/fusion/sequence/intrinsic/at.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>

#include "ndarray/tables_fwd.hpp"
#include "ndarray/tables/detail/TraitsAccess.hpp"
#include "ndarray/tables/Field.hpp"
#include "ndarray/tables/detail/functional.hpp"

namespace ndarray { namespace tables {

template <typename T>
class Layout {
public:
    typedef typename detail::TraitsAccess<T>::FieldSequence FieldSequence;
    typedef typename boost::fusion::result_of::size<FieldSequence>::type Size;

    template <int N>
    struct At {
        typedef typename boost::fusion::result_of::value_at_c<FieldSequence,N>::type Type;
    };

    template <int N>
    typename At<N>::Type const &
    operator[](Index<N> index) const {
        return boost::fusion::at_c<N>(_sequence);
    }

    template <int N>
    typename At<N>::Type &
    operator[](Index<N> index) {
        return boost::fusion::at_c<N>(_sequence);
    }

    void normalize(bool pack=false) {
        detail::SetOffsets function(pack, _bytes, _alignment);
        boost::fusion::for_each(_sequence, function);
    }

    FieldSequence const & getSequence() const { return _sequence; }

    FieldSequence & getSequence() { return _sequence; }

    int getMinStride() const { return _bytes + _bytes % _alignment; }

    int getBytes() const { return _bytes; }

    int getAlignment() const { return _alignment; }

    explicit Layout() : _sequence(), _bytes(0), _alignment(1) {}

    Layout(Layout const & other) :
        _sequence(other._sequence), _bytes(other._bytes), _alignment(other._alignment) {}

    bool operator==(Layout const & other) const {
        return _sequence == other._sequence;
    }

    bool operator!=(Layout const & other) const {
        return _sequence != other._sequence;
    }

private:
    FieldSequence _sequence;
    int _bytes;
    int _alignment;
};

}} // namespace ndarray::tables

#endif // !NDARRAY_TABLES_Layout_hpp_INCLUDED
