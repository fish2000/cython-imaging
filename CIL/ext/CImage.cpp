
extern "C" {
#include "Python.h"
#include "numpy/arrayobject.h"
#include "numpy/ndarraytypes.h"
}

#define cimg_display 0
#ifndef cimg_debug
    #define cimg_debug 1
#endif
#include "CImg.h"

#define cimg_use_magick
#include "include/CImage.h"
#include <cstdlib>

namespace CIL {

    CImage::CImage() {}
    CImage::CImage(int newDType) {}
    
    CImage::CImage(CImage *cimageFrom) {}
    CImage::CImage(CImage *cimageFrom, int newDType) {}
    CImage::CImage(CImage *cimageFrom, const bool cimageSharesData) {}
    
    CImage::CImage(cimg_library::CImg<uint> cimgFrom) {
        raw = cimgFrom;
        dtype = (char*)raw.pixel_type();
    }
    CImage::CImage(cimg_library::CImg<uint> cimgFrom, int newDType) {}
    CImage::CImage(cimg_library::CImg<uint> cimgFrom, const bool cimgIsShared) {}
    
    CImage::CImage(const char *const cimgFilename) {
        raw = cimg_library::CImg<uint>(cimgFilename);
        dtype = (char*)raw.pixel_type();
    }
    
    #if CIL_NUMPY
    CImage::CImage(PyObject *ndimg) {
        PyArrayObject *ndarrayimage;
        
        assert(npimg->nd==3);
        printf("YO DOGG: %p %d x %d x %d\\n", \
            PyArray_DATA(ndimg), \
            PyArray_DIM(ndimg, 1), PyArray_DIM(ndimg, 0), PyArray_DIM(ndimg, 2));
        
        raw = cimg_library::CImg<uint>();
        ndarrayimage = (PyArrayObject *)PyArray_ContiguousFromObject(ndimg, PyArray_UBYTE, 2, 4);
        raw.assign(cimg_library::CImg<uint>(ndarrayimage->data), \
            PyArray_DIM(ndimg, 1), PyArray_DIM(ndimg, 0), 1, PyArray_DIM(ndimg, 2));
        Py_DECREF(ndarrayimage);
    }
    #endif
    
    /*
    CImage::CImage(PyObject *pilFrom) {}
    */
    
    /* ---> DESTRUCTOR <--- */
    CImage::~CImage() {}
}


