#include <boost/python/eigen.hpp>

namespace bp = boost::python;

template <typename M>
bool acceptMatrix(M m) {
    return (m(0,0) == 1) && (m(0,1) == 2) && (m(0,2) == 3) 
        && (m(1,0) == 4) && (m(1,1) == 5) && (m(1,2) == 6);
}

template <typename M>
bool acceptVector(M m) {
    return (m[0] == 1) && (m[1] == 2) && (m[2] == 3) && (m[3] == 4);
}

template <typename M>
void fillMatrix(M & m) {
    m(0,0) = 1;
    m(0,1) = 2;
    m(0,2) = 3;
    m(1,0) = 4;
    m(1,1) = 5;
    m(1,2) = 6;
}

template <typename M>
M returnMatrix() {
    static typename boost::remove_const<typename boost::remove_reference<M>::type>::type m(2,3);
    fillMatrix(m);
    return m;
}

template <typename M>
void fillVector(M & m) {
    m[0] = 1;
    m[1] = 2;
    m[2] = 3;
    m[3] = 4;
}

template <typename M>
M returnVector() {
    static typename boost::remove_const<typename boost::remove_reference<M>::type>::type m(4);
    fillVector(m);
    return m;
}

template <typename M>
bp::object returnObject() {
    static typename boost::remove_const<typename boost::remove_reference<M>::type>::type m(2,3);
    fillMatrix(m);
    bp::object o(m);
    return o;
}

template <typename M>
Eigen::Transpose<M> returnTranspose() {
    static typename boost::remove_const<typename boost::remove_reference<M>::type>::type m(2,3);
    fillMatrix(m);
    return m.transpose();
}

template <typename M>
class MatrixOwner {
    M _matrix;
public:
    MatrixOwner() : _matrix(2,3) { fillMatrix(_matrix); }

    M const & getMatrix_cref() const { return _matrix; }
    M & getMatrix_ref() { return _matrix; }

    bool compareData(boost::numpy::matrix const & mp) const {
        return _matrix.data() == reinterpret_cast<double const*>(mp.get_data());
    }

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    static void declare(char const * name) {
        bp::class_< MatrixOwner, boost::shared_ptr<MatrixOwner> >(name)
            .def("getMatrix_cref", &MatrixOwner::getMatrix_cref,
                 bp::return_internal_matrix<>())
            .def("getMatrix_ref", &MatrixOwner::getMatrix_ref,
                 bp::return_internal_matrix<>())
            .def("compareData", &MatrixOwner::compareData)
            ;
    }
};

static const int X = Eigen::Dynamic;

BOOST_PYTHON_MODULE(eigen_mod) {
    boost::numpy::initialize();
    bp::def("acceptMatrix_23d_cref", acceptMatrix< Eigen::Matrix<double,2,3> const &>);
    bp::def("acceptMatrix_X3d_cref", acceptMatrix< Eigen::Matrix<double,X,3> const &>);
    bp::def("acceptMatrix_2Xd_cref", acceptMatrix< Eigen::Matrix<double,2,X> const &>);
    bp::def("acceptMatrix_XXd_cref", acceptMatrix< Eigen::Matrix<double,X,X> const &>);
    bp::def("acceptVector_41d_cref", acceptVector< Eigen::Matrix<double,4,1> const &>);
    bp::def("acceptVector_X1d_cref", acceptVector< Eigen::Matrix<double,X,1> const &>);
    bp::def("acceptVector_14d_cref", acceptVector< Eigen::Matrix<double,1,4> const &>);
    bp::def("acceptVector_1Xd_cref", acceptVector< Eigen::Matrix<double,1,X> const &>);
    bp::def("returnVector_41d_sq", returnVector< Eigen::Matrix<double,4,1> >, bp::squeeze_matrix<>());
    bp::def("returnVector_14d_sq", returnVector< Eigen::Matrix<double,1,4> >, bp::squeeze_matrix<>());
    bp::def("returnVector_X1d_sq", returnVector< Eigen::Matrix<double,X,1> >, bp::squeeze_matrix<>());
    bp::def("returnVector_1Xd_sq", returnVector< Eigen::Matrix<double,1,X> >, bp::squeeze_matrix<>());
    bp::def("returnVector_41d", returnVector< Eigen::Matrix<double,4,1> >);
    bp::def("returnVector_14d", returnVector< Eigen::Matrix<double,1,4> >);
    bp::def("returnVector_X1d", returnVector< Eigen::Matrix<double,X,1> >);
    bp::def("returnVector_1Xd", returnVector< Eigen::Matrix<double,1,X> >);
    bp::def("returnMatrix_23d", returnMatrix< Eigen::Matrix<double,2,3> >);
    bp::def("returnMatrix_X3d", returnMatrix< Eigen::Matrix<double,X,3> >);
    bp::def("returnMatrix_2Xd", returnMatrix< Eigen::Matrix<double,2,X> >);
    bp::def("returnMatrix_XXd", returnMatrix< Eigen::Matrix<double,X,X> >);
    bp::def("returnMatrix_23d_c", returnMatrix< Eigen::Matrix<double,2,3> const>);
    bp::def("returnMatrix_X3d_c", returnMatrix< Eigen::Matrix<double,X,3> const>);
    bp::def("returnMatrix_2Xd_c", returnMatrix< Eigen::Matrix<double,2,X> const>);
    bp::def("returnMatrix_XXd_c", returnMatrix< Eigen::Matrix<double,X,X> const>);
    bp::def("returnObject_23d", returnObject< Eigen::Matrix<double,2,3> >);
    bp::def("returnObject_X3d", returnObject< Eigen::Matrix<double,X,3> >);
    bp::def("returnObject_2Xd", returnObject< Eigen::Matrix<double,2,X> >);
    bp::def("returnObject_XXd", returnObject< Eigen::Matrix<double,X,X> >);
    MatrixOwner< Eigen::Matrix<double,2,3> >::declare("MatrixOwner_23d");
    MatrixOwner< Eigen::Matrix<double,X,3> >::declare("MatrixOwner_X3d");
    MatrixOwner< Eigen::Matrix<double,2,X> >::declare("MatrixOwner_2Xd");
    MatrixOwner< Eigen::Matrix<double,X,X> >::declare("MatrixOwner_XXd");
}
