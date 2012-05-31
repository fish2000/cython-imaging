#include <ndarray/levmar.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ndarray-levmar
#include <boost/test/unit_test.hpp>

// largely copied from levmar demo code

struct Meyer : public ndarray::levmar::Function<double>, public ndarray::levmar::Jacobian<double> {

    virtual void operator()(ndarray::Array<double const,1,1> const & p, 
                            ndarray::Array<double,1,1> const & x) {
	for (register int i=0; i < x.getSize<0>(); ++i) {
            double ui = 0.45+0.05*i;
            x[i] = p[0]*std::exp(10.0*p[1]/(ui+p[2]) - 13.0);
	}
    }

    virtual void operator()(ndarray::Array<double const,1,1> const & p,
                            ndarray::Array<double,2,2> const & jac) {
        for (register int i=0; i < jac.getSize<0>(); ++i) {
            double ui = 0.45+0.05*i;
            double tmp = std::exp(10.0*p[1]/(ui+p[2]) - 13.0);
            jac[i][0] = tmp;
            jac[i][1] = 10.0*p[0]*tmp/(ui+p[2]);
            jac[i][2] = -10.0*p[0]*p[1]*tmp/((ui+p[2])*(ui+p[2]));
        }
    }

    static void initialize(ndarray::Array<double,1,1> & p, ndarray::Array<double,1,1> & x) {
        p = ndarray::allocate(ndarray::makeVector(3));
        x = ndarray::allocate(ndarray::makeVector(16));
        p[0]=8.85; p[1]=4.0; p[2]=2.5;
        x[0]=34.780;	x[1]=28.610; x[2]=23.650; x[3]=19.630;
        x[4]=16.370;	x[5]=13.720; x[6]=11.540; x[7]=9.744;
        x[8]=8.261;	x[9]=7.030; x[10]=6.005; x[11]=5.147;
        x[12]=4.427;	x[13]=3.820; x[14]=3.307; x[15]=2.872;
    }

};

BOOST_AUTO_TEST_CASE(meyer) {
    ndarray::Array<double,1,1> p;
    ndarray::Array<double,1,1> x;
    Meyer::initialize(p, x);
    ndarray::levmar::Optimizer<double> optimizer(1000, true);
    Meyer func;
    int result = optimizer(func, func, p, x);
    ndarray::Array<double,1,1> pf = ndarray::copy(p);
    BOOST_CHECK(result > 0);
    Meyer::initialize(p, x);
    result = optimizer(func, p, x);
    BOOST_CHECK(result > 0);
#ifndef GCC_45
    BOOST_CHECK(ndarray::allclose(p, pf, 1E-5));
#endif
}
