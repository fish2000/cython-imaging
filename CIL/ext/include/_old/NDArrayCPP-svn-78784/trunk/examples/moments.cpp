// Include this first (includes Python.h):
#include <boost/python.hpp>

// Include this AFTER Python.h; in projects with multiple source files, you may have to
// #define some things first (see the numpy C-API documentation).
#include <numpy/arrayobject.h>

// Include this AFTER numpy/arrayobject.h; provides automatic Boost.Python converters
// between ndarray::Array and numpy.ndarray.
#include <ndarray/python/boost/Array.hpp>

// Provides automatic Boost.Python converters between ndarray::Vector and Python's tuple.
#include <ndarray/python/boost/Vector.hpp>

#include <stdexcept>

typedef ndarray::Array<double const,2> Image;
typedef ndarray::Array<double const,2>::Reference Row;

ndarray::Vector<double,2> centroid(
    Image const & image,
    Image const & x,
    Image const & y
) {
    if (x.getShape() != image.getShape() || y.getShape() != image.getShape()) {
        throw std::invalid_argument("image and coordinate arrays do not have the same shape");
    }
    Image::Iterator const row_end = image.end();
    Image::Iterator row_iter = image.begin();
    Image::Iterator x_row_iter = x.begin();
    Image::Iterator y_row_iter = y.begin();
    double total = 0.0;
    ndarray::Vector<double,2> result;
    for (; row_iter != row_end; ++row_iter, ++x_row_iter, ++y_row_iter) {
        Row::Iterator const pix_end = row_iter->end();
        Row::Iterator pix_iter = row_iter->begin();
        Row::Iterator x_pix_iter = x_row_iter->begin();
        Row::Iterator y_pix_iter = y_row_iter->begin();
        for (; pix_iter != pix_end; ++pix_iter, ++x_pix_iter, ++y_pix_iter) {
            total += *pix_iter;
            result[0] += (*pix_iter) * (*y_pix_iter);
            result[1] += (*pix_iter) * (*x_pix_iter);
        }
    }
    result /= total;
    return result;
}

ndarray::Array<double,1,1> gaussian(int size, double sigma, double mu) {
    ndarray::Array<double,1,1> result = ndarray::allocate(ndarray::makeVector(size));
    ndarray::Array<double,1,1>::Iterator i = result.begin();
    for (int n = 0; n < size; ++n, ++i) {
        double x = (n - mu) / sigma;
        *i = std::exp(-0.5 * x * x);
    }
    return result;
}

BOOST_PYTHON_MODULE(moments) {
    import_array(); // need to call this at the beginning of the module
    boost::python::def(
        "centroid", &centroid, boost::python::args("image", "x", "y"),
        "xc,yc = centroid(image, x, y)\n\n"
        "Compute the centroid of an image, with x and y coordinates\n"
        "provided by a pair of 2-d coordinate grids.\n"
    );
    boost::python::def(
        "gaussian", &gaussian, boost::python::args("size", "sigma", "mu"),
        "array = gaussian(size, sigma, mu)\n\n"
        "Create a 1-d Gaussian function with scale sigma and offset mu.\n"
    );
}
