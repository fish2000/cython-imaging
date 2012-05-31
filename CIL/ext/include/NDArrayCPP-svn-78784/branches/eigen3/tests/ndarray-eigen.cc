#include <ndarray/eigen.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ndarray-eigen
#include <boost/test/unit_test.hpp>

template <typename T, typename U>
void testElements2(T const & a, U const & b) {
    BOOST_CHECK( a.template getSize<0>() == b.rows() );
    BOOST_CHECK( a.template getSize<1>() == b.cols() );
    BOOST_CHECK( a.template getStride<0>() == b.rowStride() );
    BOOST_CHECK( a.template getStride<1>() == b.colStride() );
    for (int i = 0; i < b.rows(); ++i) {
        for (int j = 0; j < b.cols(); ++j) {
            BOOST_CHECK(&a[i][j] == &b(i,j));
        }
    }    
}

template <typename T, typename U>
void testElements1(T const & a, U const & b) {
    BOOST_CHECK( a.template getSize<0>() == b.size() );
    BOOST_CHECK( a.template getStride<0>() == b.innerStride() );
    for (int i = 0; i < b.size(); ++i) {
        BOOST_CHECK(&a[i] == &b[i]);
    }
}

template <int C, int Rows, int Cols>
void testEigenView(ndarray::EigenView<double,2,C,Eigen::ArrayXpr,Rows,Cols> b) {
    ndarray::Array<double,2,C> a(b.shallow());
    b.setRandom();
    testElements2(a, b);
    Eigen::Matrix<double,Rows,Eigen::Dynamic> m1(b.rows(), 6);
    m1.setRandom(b.rows(), 6);
    Eigen::Matrix<double,Eigen::Dynamic,Cols> m2(6, b.cols());
    m2.setRandom(6, b.cols());
    b.matrix() = m1 * m2;
    Eigen::Matrix<double,Rows,Cols> m3 = m1 * m2;
    for (int i = 0; i < b.rows(); ++i) {
        for (int j = 0; j < b.cols(); ++j) {
            BOOST_CHECK(a[i][j] == m3(i,j));
        }
    }
    Eigen::Array<double,Rows,Cols> m4(b.rows(), b.cols());
    m4.setRandom();
    Eigen::Array<double,Rows,Cols> m5 = m4 * b;
    Eigen::Array<double,Rows,Cols> m6 = m4 * m3.array();
    BOOST_CHECK( (m5 == m6).all() );
}

template <int C, int Rows, int Cols>
void testEigenView(ndarray::EigenView<double,1,C,Eigen::ArrayXpr,Rows,Cols> b) {
    ndarray::Array<double,1,C> a(b.shallow());
    b.setRandom();
    testElements1(a, b);
    Eigen::Matrix<double,Rows,Eigen::Dynamic> m1(b.rows(), 6);
    m1.setRandom(b.rows(), 6);
    Eigen::Matrix<double,Eigen::Dynamic,Cols> m2(6, b.cols());
    m2.setRandom(6, b.cols());
    b.matrix() = m1 * m2;
    Eigen::Matrix<double,Rows,Cols> m3 = m1 * m2;
    for (int i = 0; i < b.rows(); ++i) {
        BOOST_CHECK(a[i] == m3[i]);
    }
    Eigen::Array<double,Rows,Cols> m4(b.rows(), b.cols());
    m4.setRandom();
    Eigen::Array<double,Rows,Cols> m5 = m4 * b;
    Eigen::Array<double,Rows,Cols> m6 = m4 * m3.array();
    BOOST_CHECK( (m5 == m6).all() );
}

template <int C, int Rows, int Cols>
void testEigenView(ndarray::EigenView<double,2,C,Eigen::MatrixXpr,Rows,Cols> b) {
    ndarray::Array<double,2,C> a(b.shallow());
    b.setRandom();
    testElements2(a, b);
    Eigen::Matrix<double,Rows,Eigen::Dynamic> m1(b.rows(), 6);
    m1.setRandom(b.rows(), 6);
    Eigen::Matrix<double,Eigen::Dynamic,Cols> m2(6, b.cols());
    m2.setRandom(6, b.cols());
    b = m1 * m2;
    Eigen::Matrix<double,Rows,Cols> m3 = m1 * m2;
    for (int i = 0; i < b.rows(); ++i) {
        for (int j = 0; j < b.cols(); ++j) {
            BOOST_CHECK(a[i][j] == m3(i,j));
        }
    }
    Eigen::Array<double,Rows,Cols> m4(b.rows(), b.cols());
    m4.setRandom();
    Eigen::Array<double,Rows,Cols> m5 = m4 * b.array();
    Eigen::Array<double,Rows,Cols> m6 = m4 * m3.array();
    BOOST_CHECK( (m5 == m6).all() );
}

template <int C, int Rows, int Cols>
void testEigenView(ndarray::EigenView<double,1,C,Eigen::MatrixXpr,Rows,Cols> b) {
    ndarray::Array<double,1,C> a(b.shallow());
    b.setRandom();
    testElements1(a, b);
    Eigen::Matrix<double,Rows,Eigen::Dynamic> m1(b.rows(), 6);
    m1.setRandom(b.rows(), 6);
    Eigen::Matrix<double,Eigen::Dynamic,Cols> m2(6, b.cols());
    m2.setRandom(6, b.cols());
    b = m1 * m2;
    Eigen::Matrix<double,Rows,Cols> m3 = m1 * m2;
    for (int i = 0; i < b.rows(); ++i) {
        BOOST_CHECK(a[i] == m3[i]);
    }
    Eigen::Array<double,Rows,Cols> m4(b.rows(), b.cols());
    m4.setRandom();
    Eigen::Array<double,Rows,Cols> m5 = m4 * b.array();
    Eigen::Array<double,Rows,Cols> m6 = m4 * m3.array();
    BOOST_CHECK( (m5 == m6).all() );
}

template <typename XprKind>
void invokeEigenViewTests() {
    ndarray::Array<double,2,2> a22(ndarray::allocate(5,4));
    testEigenView(a22.asEigen<XprKind>());
    testEigenView(a22.asEigen<XprKind,5,4>());
    testEigenView(a22.transpose().asEigen<XprKind>());
    testEigenView(a22.transpose().asEigen<XprKind,4,5>());
    ndarray::Array<double,2,1> a21(a22[ndarray::view()(0,3)]);
    testEigenView(a21.asEigen<XprKind>());
    testEigenView(a21.asEigen<XprKind,5,3>());
    testEigenView(a21.transpose().asEigen<XprKind>());
    testEigenView(a21.transpose().asEigen<XprKind,3,5>());
    ndarray::Array<double,2,0> a20(a22[ndarray::view()(0,4,2)]);
    testEigenView(a20.asEigen<XprKind>());
    testEigenView(a20.asEigen<XprKind,5,2>());
    testEigenView(a20.transpose().asEigen<XprKind>());
    testEigenView(a20.transpose().asEigen<XprKind,2,5>());
    ndarray::Array<double,1,1> a11(ndarray::allocate(4));
    testEigenView(a11.asEigen<XprKind>());
    testEigenView(a11.asEigen<XprKind,4,1>());
    testEigenView(a11.asEigen<XprKind,1,4>());
    testEigenView(a11.transpose().asEigen<XprKind>());
    testEigenView(a11.transpose().asEigen<XprKind,4,1>());
    testEigenView(a11.transpose().asEigen<XprKind,1,4>());
    ndarray::Array<double,1,0> a10(a11[ndarray::view(0,4,2)]);
    testEigenView(a10.asEigen<XprKind>());
    testEigenView(a10.asEigen<XprKind,2,1>());
    testEigenView(a10.asEigen<XprKind,1,2>());
    testEigenView(a10.transpose().asEigen<XprKind>());
    testEigenView(a10.transpose().asEigen<XprKind,2,1>());
    testEigenView(a10.transpose().asEigen<XprKind,1,2>());
}

BOOST_AUTO_TEST_CASE(EigenView) {
    invokeEigenViewTests<Eigen::ArrayXpr>();
    invokeEigenViewTests<Eigen::MatrixXpr>();

    Eigen::MatrixXd m(Eigen::MatrixXd::Random(5,6));
    ndarray::SelectEigenView<Eigen::MatrixXd>::Type v(ndarray::copy(m));
    BOOST_CHECK( (v.array() == m.array()).all() );
}
