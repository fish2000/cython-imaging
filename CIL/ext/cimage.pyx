
from libcpp cimport bool
from cimagestruct cimport CImg
cimport numpy
import numpy

from cython.operator cimport dereference as deref

numpy.import_array()

cdef class CWrappedImage:

    cdef __cythonbufferdefaults__ = {"ndim": 3, "mode": "c"}
    cdef CImg[int]* cimage

    def __cinit__(self):
        cdef CImg[int] ccimage = CImg[int]()
        self.cimage = &ccimage

    def fromarray(self,
        object[numpy.int32_t, ndim=3, mode="c"] ndimage_in):
        self.cimage.assign(<CImg[int]>ndimage_in)
        return self

    def blur(self,
        float sX=1.0, float sY=1.0, float sZ=1.0,
        bool boundary_conditions=True):
        if self.cimage:
            self.cimage.blur(sX, sY, sZ, boundary_conditions)
        return self

    '''
    def toarray(self not None,
        numpy.ndarray[numpy.uint8_t, ndim=3, mode="c"] ndimage_out not None):
    '''

    def toarray(self,
        object[numpy.int32_t, ndim=3, mode="c"] ndimage_out):
        if self.cimage:
            ndimage_out = <CImg[int]>self.cimage._data
        return self


