#ifndef NDARRAY_TABLES_Table_hpp_INCLUDED
#define NDARRAY_TABLES_Table_hpp_INCLUDED

#include "ndarray.hpp"
#include "ndarray/tables/Row.hpp"
#include "ndarray/tables/detail/Iterator.hpp"

namespace ndarray { namespace tables {

template <typename T>
class Table {
    typedef typename detail::TraitsAccess<T>::Columns_ Columns_;
    typedef typename detail::TraitsAccess<T>::ColumnsPtr ColumnsPtr;
public:

    typedef detail::Iterator<T> Iterator;
    typedef typename detail::TraitsAccess<T>::Row_ Record;
    typedef typename detail::TraitsAccess<T>::Layout_ Layout_;
    typedef typename detail::TraitsAccess<T>::Raw Raw;

    /// @brief Metafunction for the return type of operator[].
    template <int N>
    struct At {
        typedef typename detail::TraitsAccess<T>::template Fields<N>::Qualified::TableValue Type;
    };

    /// @brief Return a column array.
    template <int N>
    typename At<N>::Type operator[](Index<N> index) const {
        return typename At<N>::Type((*_columns)[index]);
    }

    Iterator begin() const { return Iterator(Record(0, _columns)); }
    Iterator end() const { return Iterator(Record(_columns->getSize(), _columns)); }

    /// @brief Return the nth record.
    Record const operator[](int n) const { return Record(n, _columns); }

    /// @brief Trivial indexing.
    Table const & operator[](View< boost::fusion::vector1<ndarray::index::Full> > const & def) const {
        return *this;
    }

    /// @brief Range indexing.
    Table const operator[](View< boost::fusion::vector1<ndarray::index::Range> > const & def) const {
        return Table(_columns->index(boost::fusion::front(def._seq)));
    }

    /// @brief Slice indexing.
    Table const operator[](View< boost::fusion::vector1<ndarray::index::Slice> > const & def) const {
        return Table(_columns->index(boost::fusion::front(def._seq)));
    }

    /// @brief Get the manager that determines the lifetime of the underlying memory block.
    Manager::Ptr getManager() const { return _columns->getManager(); }

    /// @brief Get a raw pointer to the beginning of the table.
    Raw * getRaw() const { return _columns->getRaw(); }

    /// @brief Return the number of rows in the table.
    int getSize() const { return _columns->getSize(); }

    /// @brief Return difference in bytes between two adjacent rows.
    int getStride() const { return _columns->getStride(); }

    /// @brief Get the Layout that describes the field types.
    Layout_ const & getLayout() const { return _columns->getLayout(); }

    /// @brief Shallow copy construction.
    Table(Table const & other) : _columns(other._columns) {}

    /// @brief Shallow assignment.
    Table & operator=(Table const & other)  {
        if (&other != this) _columns = other._columns;
        return *this;
    }

    /// @brief Converting (non-const to const) shallow copy construction.
    template <typename U>
    Table(Table<U> const & other, typename boost::enable_if< boost::is_same<U const,T> >::type* e = 0) :
        _columns(other._columns)
    {}

    /// @brief Converting (non-const to const) shallow assignment.
    template <typename U>
    typename boost::enable_if< boost::is_same<U const,T>, Table & >::type
    operator=(Table<U> const & other) {
        _columns = other._columns;
        return *this;
    }

    /// @brief Construct a new table with the given size and layout.
    static Table allocate(int size, Layout_ const & layout) {
        Layout_ normal_layout(layout);
        normal_layout.normalize();
        int stride = normal_layout.getMinStride();
        std::pair<Manager::Ptr,Raw*> p = SimpleManager<Raw>::allocate(size * stride);
        ColumnsPtr columns = Columns_::create(size, stride, p.second, p.first, normal_layout);
        return Table(columns);
    }

    /// @brief Construct a new table from existing data.
    template <typename Owner>
    static Table external(Raw * raw, int size, int stride, Layout_ const & layout, Owner const & owner) {
        Layout_ normal_layout(layout);
        normal_layout.normalize();
        Manager::Ptr manager = ExternalManager<Owner>::make(owner);
        ColumnsPtr columns = Columns_::create(size, stride, raw, manager, normal_layout);
        return Table(columns);
    }

private:
    
    template <typename U> friend class Table;

    explicit Table(ColumnsPtr const & columns) : _columns(columns) {}

    ColumnsPtr _columns;

};

}} // namespace ndarray::tables

#endif // !NDARRAY_TABLES_Table_hpp_INCLUDED
