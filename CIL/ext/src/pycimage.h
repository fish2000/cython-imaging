extern "C" {
    #include <CImg.h>
    #include <Python.h>
}

typedef cimg_library::CImg<float>&  CImg_Float;
typedef cimg_library::CImg<int>&    CImg_Int32;

class PyCImage {
    public:
        virtual CImg_Int32 cimg;

    public:
        virtual ~PyCImage() { }

        virtual void* rowp(int r) = 0;

        virtual int nbits() const = 0;

        virtual int ndims() const = 0;
        virtual int dim(int) const = 0;

        virtual int dim_or(int dim, int def) const {
            if (dim >= this->ndims()) return def;
            return this->dim(dim);
        }

        template<typename T>
        T* rowp_as(const int r) {
            return static_cast<T*>(this->rowp(r));
        }
};

class ImageFactory {
    public:
        virtual ~ImageFactory() { }
        virtual std::auto_ptr<Image>
            create(int nbits, int d0, int d1, int d2, int d3=-1, int d4=-1) = 0;

    protected:
};
