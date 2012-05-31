#ifndef NDARRAY_TABLES_DETAIL_Columns_hpp_INCLUDED
#define NDARRAY_TABLES_DETAIL_Columns_hpp_INCLUDED

#include <boost/intrusive_ptr.hpp>
#include <boost/fusion/algorithm/transformation/transform.hpp>

#include "ndarray/tables/detail/functional.hpp"
#include "ndarray/tables/Layout.hpp"

namespace ndarray { namespace tables { namespace detail {

template <typename T>
class Columns {
    BOOST_STATIC_ASSERT( !boost::is_const<T>::value );
public:

    typedef boost::intrusive_ptr<Columns> Ptr;

    typedef unsigned char Raw;

    typedef typename TraitsAccess<T>::FieldSequence FieldSequence;

    typedef typename boost::fusion::result_of::as_vector<
        typename boost::fusion::result_of::transform<FieldSequence,MakeColumns>::type
        >::type ColumnSequence;

    template <int N>
    struct At {
        typedef typename TraitsAccess<T>::template Fields<N>::Info::ColumnValue Type;
    };

    template <int N>
    typename At<N>::Type operator[](Index<N> index) const {
        return boost::fusion::at_c<N>(_sequence);
    }

    Ptr index(ndarray::index::Range const & dim) const {
        return Ptr(new Columns(*this, dim));
    }

    Ptr index(ndarray::index::Slice const & dim) const {
        return Ptr(new Columns(*this, dim));
    }

    int getSize() const { return _size; }

    int getStride() const { return _stride; }

    Raw * getRaw() const { return _raw; }

    Manager::Ptr getManager() const { return _manager; }

    Layout<T> const & getLayout() const { return _layout; }

    ColumnSequence const & getSequence() const { return _sequence; }

    static Ptr create(
        int size, int stride, Raw * raw, Manager::Ptr const & manager, Layout<T> const & layout
    ) {
        return Ptr(new Columns(size, stride, raw, manager, layout));
    }

    friend inline void intrusive_ptr_add_ref(Columns const * self) {
        ++self->_rc;
    }
 
    friend inline void intrusive_ptr_release(Columns const * self) {
        if ((--self->_rc)==0) delete self;
    }

private:

    static FieldSequence const & getNormalizedFieldSequence(Layout<T> & layout, int & stride) {
        layout.normalize();
        if (stride < 0) {
            stride = layout.getMinStride();
        } else {
            if (stride < layout.getBytes()) {
                throw std::logic_error("Table stride is smaller than layout size.");
            }
            if (stride % layout.getAlignment() != 0) {
                throw std::logic_error("Table stride is not evenly disible by its maximum element size.");
            }
        }
        return layout.getSequence();
    }

    Columns(Columns const & other, ndarray::index::Range const & dim) :
        _rc(1), _size(dim.stop - dim.start), _stride(other.getStride()),
        _raw(other.getRaw() + dim.start * other.getStride()),
        _manager(other.getManager()), _layout(other.getLayout()),
        _sequence(boost::fusion::transform(other.getSequence(), makeViewColumns(dim)))
    {}

    Columns(Columns const & other, ndarray::index::Slice const & dim) :
        _rc(1), _size(dim.stop - dim.start), _stride(other.getStride()),
        _raw(other.getRaw() + dim.start * other.getStride()),
        _manager(other.getManager()), _layout(other.getLayout()),
        _sequence(boost::fusion::transform(other.getSequence(), makeViewColumns(dim)))
    {}

    Columns(int size, int stride, Raw * raw, Manager::Ptr const & manager, Layout<T> const & layout) :
        _rc(1), _size(size), _stride(stride), _raw(raw), _manager(manager), _layout(layout),
        _sequence(
            boost::fusion::transform(
                getNormalizedFieldSequence(_layout, _stride),
                MakeColumns(_size, _stride, _raw, _manager)
            )
        )
    {}

    mutable int _rc;
    int _size;
    int _stride;
    Raw * _raw;
    Manager::Ptr _manager;
    Layout<T> _layout;
    ColumnSequence _sequence;
};

}}} // namespace ndarray::tables::detail

#endif // !NDARRAY_TABLES_DETAIL_Columns_hpp_INCLUDED
