#ifndef NDARRAY_TABLES_Row_hpp_INCLUDED
#define NDARRAY_TABLES_Row_hpp_INCLUDED

#include "ndarray.hpp"
#include "ndarray/tables/Layout.hpp"
#include "ndarray/tables/detail/Columns.hpp"

namespace ndarray { namespace tables {

template <typename T>
class Row {
    typedef typename detail::TraitsAccess<T>::Columns_ Columns_;
    typedef typename detail::TraitsAccess<T>::ColumnsPtr ColumnsPtr;
public:

    typedef typename detail::TraitsAccess<T>::Layout_ Layout_;
    typedef typename detail::TraitsAccess<T>::Raw Raw;


    /// @brief Metafunction for the return type of operator[].
    template <int N>
    struct At {
        typedef typename detail::TraitsAccess<T>::template Fields<N>::Qualified::RowValue Type;
    };

    /// @brief Return a single field value.
    template <int N>
    typename At<N>::Type operator[](Index<N> index) const { return (*_columns)[index][_n]; }

    /// @brief Get the Layout that describes the field types.
    Layout_ const & getLayout() const { return _columns->getLayout(); }

    /// @brief Shallow copy construction.
    Row(Row const & other) : _n(other._n), _columns(other._columns) {}

    /// @brief Shallow assignment.
    Row & operator=(Row const & other) {
        if (&other != this) {
            _n = other._n;
            _columns = other._columns;
        }
        return *this;
    }

    /// @brief Converting (non-const to const) shallow copy construction.
    template <typename U>
    Row(Row<U> const & other, typename boost::enable_if< boost::is_same<U const,T> >::type* e = 0) :
        _n(other._n), _columns(other._columns)
    {}

    /// @brief Converting (non-const to const) shallow assignment.
    template <typename U>
    typename boost::enable_if< boost::is_same<U const,T>, Row & >::type
    operator=(Row<U> const & other) {
        _n = other._n;
        _columns = other._columns;
        return *this;
    }

    /// @brief Construct a new row with the given layout.
    static Row allocate(Layout_ const & layout) {
        Layout_ normal_layout(layout);
        normal_layout.normalize();
        std::pair<Manager::Ptr,Raw*> p = SimpleManager<Raw>::allocate(normal_layout.getBytes());
        ColumnsPtr columns = Columns_::create(1, -1, p.second, p.first, normal_layout);
        return Row(0, columns);
    }

protected:

    Row() : _n(-1), _columns() {}

    Row(int n, ColumnsPtr const & columns) : _n(n), _columns(columns) {}

private:

    template <typename U> friend class Table;
    template <typename U> friend class Row;
    template <typename U> friend class detail::Iterator;

    int _n;
    ColumnsPtr _columns;
};

}} // namespace ndarray::tables

#endif // !NDARRAY_TABLES_Row_hpp_INCLUDED
