#ifndef NDARRAY_Manager_hpp_INCLUDED
#define NDARRAY_Manager_hpp_INCLUDED

/** 
 *  @file ndarray/Manager.hpp
 *
 *  @brief Definition of Manager, which manages the ownership of array data.
 */

#include "ndarray_fwd.hpp"
#include <boost/noncopyable.hpp>
#include <boost/intrusive_ptr.hpp>
#include <boost/scoped_array.hpp>

namespace ndarray {

class Manager : private boost::noncopyable {
public:

    typedef boost::intrusive_ptr<Manager> Ptr;

    friend inline void intrusive_ptr_add_ref(Manager const * manager) {
        ++manager->_rc;
    }
 
    friend inline void intrusive_ptr_release(Manager const * manager) {
        if ((--manager->_rc)==0) delete manager;
    }

    int getRC() const { return _rc; }

    virtual bool isUnique() const { return false; }

protected:

    explicit Manager() : _rc(1) {}

    virtual ~Manager() {}

private:
    mutable int _rc;
};

template <typename T>
class SimpleManager : public Manager {
    typedef typename boost::remove_const<T>::type U;
public:
    
    static std::pair<Manager::Ptr,T*> allocate(int size) {
        boost::intrusive_ptr<SimpleManager> r(new SimpleManager(size), false);
        return std::pair<Manager::Ptr,T*>(r, r->_p.get());
    }

    virtual bool isUnique() const { return true; }

private:
    explicit SimpleManager(int size) : _p() {
        if (size > 0) _p.reset(new U[size]);
    }
    boost::scoped_array<U> _p;
};

template <typename U>
class ExternalManager : public Manager, private U {
public:
    typedef U Owner;

    static Manager::Ptr make(Owner const & owner) {
        return Manager::Ptr(new ExternalManager(owner), false);
    }

    Owner const & getOwner() const { return *static_cast<Owner const *>(this); }

private:
    explicit ExternalManager(Owner const & owner) : Owner(owner) {}
};

} // namespace ndarray

#endif // !NDARRAY_Manager_hpp_INCLUDED
