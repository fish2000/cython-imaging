
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

#include <CImg.h>

extern "C" {
    #include "Python.h"
    #include <boost/python.hpp>

    /* NUMPY STUFF */
    #ifndef CIL_NUMPY
        #define CIL_NUMPY 0
    #else
        #define CIL_NUMPY 1
        #if !defined(EXTMODULE_IMPORT_ARRAY)
            #define IMPORT_ARRAY 0
        #else
            #define IMPORT_ARRAY 1
            #include "numpy/arrayobject.h"
            #include "numpy/ndarraytypes.h"
        #endif
    #endif

} /* end of 'extern "C"' */


using namespace cimg_library;

template<typename T>
class PImg : public CImg<T> {
    public:
        Py_ssize_t view_count = 0;
    public:
        int PImg_getbuffer(PyObject *self, Py_buffer *view, int flags) {
            static Py_ssize_t suboffsets[2] = { 0, -1 };
            view->buf = this->data();
            view->len = this->size();
            view->readonly = 0;
            view->ndims = 1;
            view->shape = { this->height(), this->width() };
            view->strides = { sizeof(T*), cimg::type<T>::cut() };
            view->suboffsets = suboffsets;
            this->view_count++;
            return 0;
        }
        int PImg_releasebuffer(PyObject *self, Py_buffer *view) {
            this->view_count--;
            return 0;
        }

        static void boost_export(std::string name) {
            boost::python::class_<PImg<T>>(name.c_str());
        }

};

