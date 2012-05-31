#ifndef NDARRAY_levmar_hpp_INCLUDED
#define NDARRAY_levmar_hpp_INCLUDED

#include <boost/scoped_array.hpp>
#include <vector>
#include "ndarray.hpp"

namespace ndarray {
namespace levmar {

template <typename T>
struct Function {
    virtual void operator()(Array<T const,1,1> const & p, Array<T,1,1> const & hx) = 0;
    virtual ~Function() {}
};

template <typename T>
struct Jacobian {
    virtual void operator()(Array<T const,1,1> const & p, Array<T,2,2> const & k) = 0;
    virtual ~Jacobian() {}
};

template <typename T, typename Wrapped>
struct WrapperFunction : public Function<T> {
    WrapperFunction(Wrapped wrapped) : _wrapped(wrapped) {}

    virtual void operator()(Array<T const,1,1> const & p, Array<T,1,1> const & hx) {
        _wrapped(p, hx);
    }

private:
    Wrapped _wrapped;
};

template <typename T, typename Wrapped>
struct WrapperJacobian : public Jacobian<T> {
    WrapperJacobian(Wrapped wrapped) : _wrapped(wrapped) {}

    virtual void operator()(Array<T const,1,1> const & p, Array<T,2,2> const & j) {
        _wrapped(p, j);
    }

private:
    Wrapped _wrapped;
};

template <typename T>
struct BoxConstraint {
    Array<T,1,1> const lower;
    Array<T,1,1> const upper;
    
    explicit BoxConstraint(int m, bool constrain_lower=true, bool constrain_upper=true) :
        lower(constrain_lower ? allocate(makeVector(m)) : Array<T,1,1>()),
        upper(constrain_upper ? allocate(makeVector(m)) : Array<T,1,1>())
    {}
};

template <typename T>
class Optimizer {
public:

    struct Options {
        enum Enum { TAU=0, EPSILON1=1, EPSILON2=2, EPSILON3=3, DELTA=4 };

        T & operator[](Enum n) { return _data[n]; }
        T const & operator[](Enum n) const { _data[n]; }

        Options();

        Options & set(Enum n, T value) { _data[n] = value; return *this; }

    private:
        template <typename T_> friend class Optimizer;
        Vector<T,5> _data;
    };

    explicit Optimizer(int itmax, bool compute_cov=true);

    template <typename Func, typename JacF>
    int operator()(Func func, JacF jacf, Array<T,1,1> const & p, Array<T,1,1> const & x) {
        WrapperFunction<T,Func> wfunc(func);
        WrapperJacobian<T,JacF> wjacf(jacf);
        return operator()(static_cast<Function<T>&>(wfunc), static_cast<Jacobian<T>&>(wjacf), p, x);
    }

    int operator()(Function<T> & func, Jacobian<T> & jacf, Array<T,1,1> const & p, Array<T,1,1> const & x);

    template <typename Func>
    int operator()(Func func, Array<T,1,1> const & p, Array<T,1,1> const & x) {
        WrapperFunction<T,Func> wfunc(func);
        return operator()(static_cast<Function<T>&>(wfunc), p, x);
    }

    int operator()(Function<T> & func, Array<T,1,1> const & p, Array<T,1,1> const & x);

    template <typename Func, typename JacF>
    int operator()(
        Func func, JacF jacf, Array<T,1,1> const & p, Array<T,1,1> const & x,
        BoxConstraint<T> const & bc
    ) {
        WrapperFunction<T,Func> wfunc(func);
        WrapperJacobian<T,JacF> wjacf(jacf);
        return operator()(static_cast<Function<T>&>(wfunc), static_cast<Jacobian<T>&>(wjacf), p, x, bc);
    }

    int operator()(
        Function<T> & func, Jacobian<T> & jacf, Array<T,1,1> const & p, Array<T,1,1> const & x,
        BoxConstraint<T> const & bc
    );

    template <typename Func>
    int operator()(
        Func func, Array<T,1,1> const & p, Array<T,1,1> const & x,
        BoxConstraint<T> const & bc
    ) {
        WrapperFunction<T,Func> wfunc(func);
        return operator()(static_cast<Function<T>&>(wfunc), p, x);
    }

    int operator()(
        Function<T> & func, Array<T,1,1> const & p, Array<T,1,1> const & x,
        BoxConstraint<T> const & bc
    );

    void setMaxIterations(int itmax) { _itmax = itmax; }
    void enableCovariance() { _compute_cov = true; }
    void disableCovariance() { _compute_cov = false; }

    T const operator[](int n) const { return _info[n]; }

    Array<T,2,2> getCovariance() const { return _covariance; }

    Options options;

private:

    static void _func(T * p, T * hx, int m, int n, void * adata);
    static void _jacf(T * p, T * j, int m, int n, void * adata);

    bool _compute_cov;
    int _itmax;
    boost::scoped_array<T> _info;
    Array<T,2,2> _covariance;
    std::vector<T> _workspace;
};

} // namespace ndarray::levmar
} // namespace ndarray

#ifndef NDARRAY_LEVMAR_MANUAL_INCLUDE
#include "ndarray/levmar.cc"
#endif

#endif // !NDARRAY_levmar_hpp_INCLUDED
