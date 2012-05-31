#ifndef NDARRAY_eigen_fwd_hpp_INCLUDED
#define NDARRAY_eigen_fwd_hpp_INCLUDED

/**
 *  @file ndarray/eigen_fwd.hpp
 *  @brief Forward declarations for ndarray/eigen interface.
 */

/** 
 * \defgroup EigenGroup Eigen
 * Interoperability with the Eigen 3 linear algebra library.
 */

namespace Eigen {

struct MatrixXpr;

} // namespace Eigen

namespace ndarray {
namespace detail {

template <int N, int C, int Rows, int Cols> struct EigenStrideTraits;

} // namespace detail

template <
    typename T, int N, int C, 
    typename XprKind_=Eigen::MatrixXpr,
    int Rows_=-1,
    int Cols_=((N == 1) ? 1 : -1)
    >
class EigenView;

} // namespace ndarray

#endif // !NDARRAY_eigen_fwd_hpp_INCLUDED
