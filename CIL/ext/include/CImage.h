
/* Python EXTMODULE_H flag */
#ifndef EXTMODULE_H
    #define EXTMODULE_H
    
    /* CIMG LIBRARY STUFF */
    #ifndef CIL_CIMG
        #define CIL_CIMG 1
    #endif
    
    #define cimg_display 0
    #ifndef cimg_debug
        #define cimg_debug 1
    #endif
    #include "CImg.h"
    
    #define cimg_use_magick
    
    #ifdef __cplusplus
    extern "C" {
    #endif
    #include "Python.h"
    
    /* NUMPY STUFF */
    
    #ifndef CIL_NUMPY
        #define CIL_NUMPY 0
    #else
        #define CIL_NUMPY 1
        #if !defined(EXTMODULE_IMPORT_ARRAY)
            #define NO_IMPORT_ARRAY
        #endif
        #include "numpy/arrayobject.h"
        #include "numpy/ndarraytypes.h"
    #endif
    
    #ifdef __cplusplus
    } /* end of 'extern "C"' */
    #endif
    
    /*#define cimg *cimg_library::CImg<PyArray_Descr>*/
    
    namespace CIL {
    
        class CImage {
        private:
            int ddtype;
        
        public:
            /* ---> THE WRAPPED INSTANCE <--- */
            cimg_library::CImg<PyArray_Descr> raw;
            PyArray_Descr dtype;
            
            /* ---> CONSTRUCTOR <--- */
            CImage();
            CImage(int newDType);
            
            CImage(CImage *cimageFrom);
            CImage(CImage *cimageFrom, int newDType);
            CImage(CImage *cimageFrom, const bool cimageSharesData);
            
            /* insert other constructors here:
            - with numpy instance
            - with PIL instance
            - with file path
            - with python data buffer (or whatever) etc */
            
            CImage(cimg_library::CImg<PyArray_Descr> cimgFrom);
            CImage(cimg_library::CImg<PyArray_Descr> cimgFrom, int newDType);
            CImage(cimg_library::CImg<PyArray_Descr> cimgFrom, const bool cimgIsShared);
            
            CImage(const char *const cimgFilename);
            
            #if CIL_NUMPY
            CImage(PyArrayObject *ndimg);
            #endif
            
            /* CImage(PyObject *pilFrom); */
            
            /* ---> DESTRUCTOR <--- */
            ~CImage();
            
            CImage *assign() {}/* in-place constructors */
            /* (et al...) */
            
            CImage *clear() {} /* same as assign w/ no args */
            CImage *move_to() {}
            
            /* ---> TYPE CONVERTERS <--- */
            #if CIL_NUMPY
            PyArrayObject* to_ndimg() {}
            PyArrayObject* to_ndimg(int newDType) {}
            #endif
            
        };
    
    }

#endif

