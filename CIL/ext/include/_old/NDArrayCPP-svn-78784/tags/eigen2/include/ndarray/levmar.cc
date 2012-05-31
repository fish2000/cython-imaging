#include "ndarray/levmar.hpp"
#include <levmar.h>

namespace ndarray {
namespace levmar {

template <typename T> struct Traits;

template <>
struct Traits<double> {

    typedef int (*der_t)(
        void (*func)(double *p, double *hx, int m, int n, void *adata),
        void (*jacf)(double *p, double *j, int m, int n, void *adata),
        double *p, double *x, int m, int n, int itmax, double *opts,
        double *info, double *work, double *covar, void *adata
    ); 
    static der_t der;
    
    typedef int (*dif_t)(
        void (*func)(double *p, double *hx, int m, int n, void *adata),
        double *p, double *x, int m, int n, int itmax, double *opts,
        double *info, double *work, double *covar, void *adata
    );
    static dif_t dif;
    
    typedef int (*bc_der_t)(
        void (*func)(double *p, double *hx, int m, int n, void *adata),
        void (*jacf)(double *p, double *j, int m, int n, void *adata),
        double *p, double *x, int m, int n, double *lb, double *ub, int itmax, double *opts,
        double *info, double *work, double *covar, void *adata
    );
    static bc_der_t bc_der;

    typedef int (*bc_dif_t)(
        void (*func)(double *p, double *hx, int m, int n, void *adata),
        double *p, double *x, int m, int n, double *lb, double *ub, int itmax, double *opts,
        double *info, double *work, double *covar, void *adata
    );
    static bc_dif_t bc_dif;

};

Traits<double>::der_t Traits<double>::der = &dlevmar_der;
Traits<double>::dif_t Traits<double>::dif = &dlevmar_dif;
Traits<double>::bc_der_t Traits<double>::bc_der = &dlevmar_bc_der;
Traits<double>::bc_dif_t Traits<double>::bc_dif = &dlevmar_bc_dif;

template <typename T>
Optimizer<T>::Options::Options() :
    _data(makeVector(LM_INIT_MU, LM_STOP_THRESH, LM_STOP_THRESH, LM_STOP_THRESH, LM_DIFF_DELTA))
{}

template <typename T>
Optimizer<T>::Optimizer(int itmax, bool compute_cov) :
    options(), _compute_cov(compute_cov), _itmax(itmax), _info(new T[LM_INFO_SZ])
{}

template <typename T>
int Optimizer<T>::operator()(
    Function<T> & func,
    Jacobian<T> & jacf,
    Array<T,1,1> const & p,
    Array<T,1,1> const & x
) {
    int m = p.template getSize<0>();
    int n = x.template getSize<0>();
    double * cov = 0;
    if (_compute_cov) {
        if (_covariance.template getSize<0>() != m || _covariance.template getSize<1>() != m)
            _covariance = allocate(makeVector(m, m));
        cov = _covariance.getData();
    }
    std::size_t work_size = LM_DER_WORKSZ(m, n);
    if (_workspace.size() < work_size) _workspace.resize(work_size);
    std::pair< Function<T> *, Jacobian<T> * > data(&func, &jacf);
    return Traits<T>::der(
        &_func, &_jacf, p.getData(), x.getData(), m, n, _itmax, options._data.begin(),
        _info.get(), &_workspace.front(), cov, &data
    );
}

template <typename T>
int Optimizer<T>::operator()(
    Function<T> & func,
    Array<T,1,1> const & p,
    Array<T,1,1> const & x
) {
    int m = p.template getSize<0>();
    int n = x.template getSize<0>();
    double * cov = 0;
    if (_compute_cov) {
        if (_covariance.template getSize<0>() != m || _covariance.template getSize<1>() != m)
            _covariance = allocate(makeVector(m, m));
        cov = _covariance.getData();
    }
    std::size_t work_size = LM_DIF_WORKSZ(m, n);
    if (_workspace.size() < work_size) _workspace.resize(work_size);
    std::pair< Function<T> *, Jacobian<T> * > data(&func, 0);
    return Traits<T>::dif(
        &_func, p.getData(), x.getData(), m, n, _itmax, options._data.begin(),
        _info.get(), &_workspace.front(), cov, &data
    );
}

template <typename T>
int Optimizer<T>::operator()(
    Function<T> & func,
    Jacobian<T> & jacf,
    Array<T,1,1> const & p,
    Array<T,1,1> const & x,
    BoxConstraint<T> const & bc
) {
    int m = p.template getSize<0>();
    int n = x.template getSize<0>();
    double * cov = 0;
    if (_compute_cov) {
        if (_covariance.template getSize<0>() != m || _covariance.template getSize<1>() != m)
            _covariance = allocate(makeVector(m, m));
        cov = _covariance.getData();
    }
    std::size_t work_size = LM_DER_WORKSZ(m, n);
    if (_workspace.size() < work_size) _workspace.resize(work_size);
    std::pair< Function<T> *, Jacobian<T> * > data(&func, &jacf);
    return Traits<T>::bc_der(
        &_func, &_jacf, p.getData(), x.getData(), m, n, bc.lower.getData(), bc.upper.getData(),
        _itmax, options._data.begin(),
        _info.get(), &_workspace.front(), cov, &data
    );
}

template <typename T>
int Optimizer<T>::operator()(
    Function<T> & func,
    Array<T,1,1> const & p,
    Array<T,1,1> const & x,
    BoxConstraint<T> const & bc
) {
    int m = p.template getSize<0>();
    int n = x.template getSize<0>();
    double * cov = 0;
    if (_compute_cov) {
        if (_covariance.template getSize<0>() != m || _covariance.template getSize<1>() != m)
            _covariance = allocate(makeVector(m, m));
        cov = _covariance.getData();
    }
    std::size_t work_size = LM_DIF_WORKSZ(m, n);
    if (_workspace.size() < work_size) _workspace.resize(work_size);
    std::pair< Function<T> *, Jacobian<T> * > data(&func, 0);
    return Traits<T>::dif(
        &_func, p.getData(), x.getData(), m, n, bc.lower.getData(), bc.upper.getData(),
        _itmax, options._data.begin(),
        _info.get(), &_workspace.front(), cov, &data
    );
}

template <typename T>
void Optimizer<T>::_func(T * p, T * hx, int m, int n, void * adata) {
    Array<T,1,1> p_array = external(p, makeVector(m));
    Array<T,1,1> hx_array = external(hx, makeVector(n));
    std::pair< Function<T> *, Jacobian<T> * > * data 
        = reinterpret_cast< std::pair< Function<T> *, Jacobian<T> * > *>(adata);
    (*data->first)(p_array, hx_array);
}

template <typename T>
void Optimizer<T>::_jacf(T * p, T * j, int m, int n, void * adata) {
    Array<T,1,1> p_array = external(p, makeVector(m));
    Array<T,2,2> j_array = external(j, makeVector(n, m));
    std::pair< Function<T> *, Jacobian<T> * > * data 
        = reinterpret_cast< std::pair< Function<T> *, Jacobian<T> * > *>(adata);
    (*data->second)(p_array, j_array);
}

} // namespace ndarray::levmar
} // namespace ndarray
