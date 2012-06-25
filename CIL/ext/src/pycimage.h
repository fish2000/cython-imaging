
extern "C" {
    #include <CImg.h>
    #include <Python.h>
}

typedef cimg_library::CImg<float>&  CImg_Float;
typedef cimg_library::CImg<int>&    CImg_Int32;

class PyRasterBuffer {
    public:
        PyObject * self;

    public:
        virtual int _get(PyObject *self, Py_buffer *view, int flags);
        virtual int _fill(PyObject *self, Py_buffer *view, int flags);
        virtual int _release(PyObject *self, Py_buffer *view);
};

class PyCImage : public virtual PyRasterBuffer {

    public:
        const CImg_Int32 cimg;

    public:
        virtual ~PyCImage();
        virtual PyCImage();

        virtual PyCImage(int w, int h);
        virtual PyCImage(int w, int h, float v);
        virtual PyCImage(int w, int h, float* v);
        virtual PyCImage(const PyCImage& im);

    public:
        virtual T* operator->();
        /*virtual Py_buffer* operator->*();*/

};


/* THE FRONT LINE IN THE THEATER OF FORWARD DECLARATION */


class PyRasterBuffer {
    public:

        virtual int _get(PyObject *self, Py_buffer *view, int flags) {
            static Py_ssize_t suboffsets[2] = { 0, -1 };

            view->buf = this->

        }

        virtual int _fill(PyObject *self, Py_buffer *view, int flags) {

        }

        virtual int _release(PyObject *self, Py_buffer *view) {

        }
};


class PyCImage : public virtual PyRasterBuffer {

    public:
        CImg_Int32 cimg;
        typedef

    public:
        virtual ~PyCImage();
        virtual PyCImage();

        virtual PyCImage(int w, int h);
        virtual PyCImage(int w, int h, float v);
        virtual PyCImage(int w, int h, float* v);
        virtual PyCImage(const PyCImage& im);

    public:
        virtual int* operator->() { return this.cimg.is_empty() ? 0 : this.cimg.data() };

};

int PyCImage::getbuffer(PyObject *self, Py_buffer *view, int flags) {
    PyCImage *that = (PyCImage *)self;
    static Py_ssize_t suboffsets[2] = { 0, -1 };
    view->buf = that
};



class ImageFactory {
    public:
        virtual ~ImageFactory() { }
        virtual std::auto_ptr<Image>
            create(int nbits, int d0, int d1, int d2, int d3=-1, int d4=-1) = 0;
};
