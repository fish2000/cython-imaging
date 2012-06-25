
from libcpp cimport bool

cdef extern from "CImg.h" namespace "cimg_library":

    cdef cppclass CImg[T]:

        T* _data
        T* iterator
        T* _iterator
        T value_type

        CImg()

        CImg( unsigned int size_x,  unsigned int size_y,  unsigned int size_z,  unsigned int size_c)
        CImg( unsigned int size_x,  unsigned int size_y,  unsigned int size_z,  unsigned int size_c,  T value)

        CImg( unsigned int size_x,  unsigned int size_y,  unsigned int size_z,  unsigned int size_c,  char * values,  bool repeat_values)

        # template
        CImg( T * values,  unsigned int size_x,  unsigned int size_y,  unsigned int size_z,  unsigned int size_c,  bool is_shared)
        #CImg( T * values,  unsigned int size_x,  unsigned int size_y,  unsigned int size_z,  unsigned int size_c,  bool is_shared)
        CImg( char * filename)

        CImg( CImg[T] &img)
        #CImg( CImg[T] &img)
        CImg( CImg[T] &img,  bool is_shared)
        #CImg( CImg[T] &img,  bool is_shared)

        CImg[T]& assign()
        CImg[T]& assign( unsigned int size_x,  unsigned int size_y,  unsigned int size_z,  unsigned int size_c)
        CImg[T]& assign( unsigned int size_x,  unsigned int size_y,  unsigned int size_z,  unsigned int size_c,  T value)
        CImg[T]& assign( unsigned int size_x,  unsigned int size_y,  unsigned int size_z,  unsigned int size_c,  char * values,  bool repeat_values)

        CImg[T]& assign ( CImg[T] &img)
        CImg[T]& assign ( CImg[T] &img,  bool is_shared)

        CImg[T]& blur( float sigma_x,  float sigma_y,  float sigma_z,  bool boundary_conditions)
