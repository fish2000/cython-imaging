
cdef extern from "CImg.h" namespace "cimg_library":
    
    # exceptions must be declared before throwing
    # like: cdef int CImgException except + 
    
    cdef cppclass CImgException:
        CImgException()
        CImgException(const char *const format, ...)
        const char *what() const throw()
    
    cdef cppclass CImgInstanceException(CImgException):
        CImgInstanceException(const char *const format, ...)
    
    cdef cppclass CImgArgumentException(CImgException):
        CImgArgumentException(const char *const format, ...)
    
    cdef cppclass CImgIOException(CImgException):
        CImgIOException(const char *const format, ...)
    
    cdef cppclass CImgDisplayException(CImgException):
        CImgDisplayException(const char *const format, ...)
    
    cdef cppclass CImgWarningException(CImgException):
        CImgWarningException(const char *const format, ...)
    
    cdef cppclass CImg:
        CImg()
        ~CImg()
        CImg(const unsigned int size_x, const unsigned int size_y=1, const unsigned int size_z=1, const unsigned int size_c=1)
        CImg(const unsigned int size_x, const unsigned int size_y, const unsigned int size_z, const unsigned int size_c, const T value)
        CImg(const unsigned int size_x, const unsigned int size_y, const unsigned int size_z, const unsigned int size_c, const int value0, const int value1, ...)
        CImg(const unsigned int size_x, const unsigned int size_y, const unsigned int size_z, const unsigned int size_c, const double value0, const double value1, ...)
        CImg(const unsigned int size_x, const unsigned int size_y, const unsigned int size_z, const unsigned int size_c, const char *const values, const bool repeat_values)
        