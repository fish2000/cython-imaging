
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
    
    CImage::CImage(cimg_library::CImg<const char> cimgFrom) {
        raw = cimgFrom;
        dtype = const charFromType(atoi(raw.pixel_type()))[0];
    }
    CImage::CImage(cimg_library::CImg<const char> cimgFrom, int newDType) {}
    CImage::CImage(cimg_library::CImg<const char> cimgFrom, const bool cimgIsShared) {}
    
    CImage::CImage(const char *const cimgFilename) {
        raw = cimg_library::CImg<const char>();
        raw.assign(cimgFilename[0]);
        dtype = const charFromType(atoi(raw.pixel_type()))[0];
    }
    
    #if CIL_NUMPY
    CImage::CImage(PyArrayObject *ndimg) {
        assert(npimg->nd==3);
        printf("YO DOGG: %p %d x %d x %d\\n", \
            PyArray_DATA(ndimg), \
            PyArray_DIM(ndimg, 1), PyArray_DIM(ndimg, 0), PyArray_DIM(ndimg, 2));
        
        raw = cimg_library::CImg<const char>();
        ndraw = ndimg->PyArray_DATA(ndimg);
        raw.assign(cimg_library::CImg<const char>(ndraw->data), \
            PyArray_DIM(ndimg, 1), PyArray_DIM(ndimg, 0), 1, PyArray_DIM(ndimg, 2));
        
        
        /*image=image.blur(2.5);
        return 0;*/
        
    }
    #endif
    
    /*
    CImage::CImage(PyObject *pilFrom) {}
    */
    
    /* ---> DESTRUCTOR <--- */
    CImage::~CImage() {}
            
}


int main(void) {}

