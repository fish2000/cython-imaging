// Bradley Atcheson
// University of British Columbia
// 2011
//
// Note that NumPy arrays passed here are assumed to be in C order (row-major)
// A 3xN NumPy array has 3 rows, but since it's stored in C order, the second
// (N) dimension runs fastest, which then corresponds to the first (x) dimension
// of CImg arrays (or simple float* buffers) that are used in this C++ code.
// So beware of the transpose going on here, and follow the examples set by
// existing functions.
//
// XCode Instruments (time profiler module) is great for optimising.
// Compile in debug mode (with -g) to make it work.
//
// Note that according to Anderson & Kak "SART" [1984] it is common practice
// to have a matrix that is overdetermined by a factor of 4. So if you have
// N unknown voxels, be sure to trace at least 4N linearly independent rays.
//
// Be sure not to write to the output HDF5 file while this module still holds
// active handles. Since it's a singleton it will only release handles when
// the library is unloaded, or when you manually call the Close_File method.



/*****************************************************************************
 * Includes
 ****************************************************************************/
#include <iostream>
#include <vector>
#include <deque>
#include <algorithm>
#include <cmath>

#include <boost/utility.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

/*
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_minmax.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
*/

#include <flann/flann.hpp>

#define cimg_display 0
#define cimg_debug 1

/*#include <CImg.h>*/
#include <pHash.h>
#include <H5Cpp.h>

/*****************************************************************************
 * Typedefs
 ****************************************************************************/
typedef flann::Index< flann::L2<float> >                Flann_index_base;
typedef boost::scoped_ptr<Flann_index_base>             Flann_index;
typedef flann::Matrix<int>                              Matrixi;
typedef flann::Matrix<float>                            Matrixf;
typedef cimg_library::CImg<float>                       Numpy;
typedef cimg_library::CImg<int>                         Numpyi;


#define DCT_MATRIX_SIZE 32

cimg_library::CImg<float>* get_dct_matrix() {

    cimg_library::CImg<float> *dct_matrix = new cimg_library::CImg<float>(
        DCT_MATRIX_SIZE, DCT_MATRIX_SIZE,
        1, 1, 1/sqrt((float)DCT_MATRIX_SIZE));

    for (int x = 0; x < DCT_MATRIX_SIZE; x++) {
        for (int y = 1; y < DCT_MATRIX_SIZE; y++) {
            dct_matrix->at(x, y) = SQRT_TWO * cos((cimg::PI/2/DCT_MATRIX_SIZE) * y * (2 * (x + 1)));
        }
    }

    return dct_matrix;
}

struct Voxel
{
    int     label;          // unique. corresponds to column in matrix
    float   x, y, z;        // 3D world position
    float   rindex;         // refractive index
    float   gx, gy, gz;     // refractive index gradient

    Voxel(  int _label,
            float _x, float _y, float _z,
            float _rindex,
            float _gx, float _gy, float _gz):

            label(_label),
            x(_x), y(_y), z(_z),
            rindex(_rindex),
            gx(_gx), gy(_gy), gz(_gz)
            {}

    Voxel():
            label(-999), x(0),y(0),z(0), rindex(0), gx(0),gy(0),gz(0)
            {}

};

struct Ray
{
    float  x,  y,  z;
    float dx, dy, dz;

    Ray(
        const Numpy& data
    ):   x(data(0)),  y(data(1)),  z(data(2)),
        dx(data(3)), dy(data(4)), dz(data(5))
    {
        assert( fabs(dx*dx + dy*dy + dz*dz - 1) < 1e-4 );
    }
};

// Axis-aligned bounding box
struct BBox
{
    double xmin, ymin, zmin;
    double xmax, ymax, zmax;

    BBox(
        float x0, float y0, float z0,
        float x1, float y1, float z1
    ): xmin(x0), ymin(y0), zmin(z0),
       xmax(x1), ymax(y1), zmax(z1)
    {}

    bool
    Inside(
        const double x, const double y, const double z
    ) {
        return (x >= xmin) && (x <= xmax) &&
               (y >= ymin) && (y <= ymax) &&
               (z >= zmin) && (z <= zmax);
    }

    // if ray intersects, move its origin to the nearest intersection point
    // returns dist to nearest intersection point, else 0
    float
    Advance_to_Intersection(
        Ray&    ray,
        double* interval=NULL       // output dist nearest to farthest intersection
    ) {
        // http://people.csail.mit.edu/amy/papers/box-jgt.pdf
        float tmin, tmax, tymin, tymax, tzmin, tzmax;

        float idx=1/ray.dx, idy=1/ray.dy, idz=1/ray.dz;

        float xbounds[] = { (float)xmin, (float)xmax };

        int b = (idx >= 0) ? 0 : 1;
        tmin = ( xbounds[b  ] - ray.x ) * idx;
        tmax = ( xbounds[1-b] - ray.x ) * idx;

        float ybounds[] = { (float)ymin, (float)ymax };

        b = (idy >= 0) ? 0 : 1;
        tymin = ( ybounds[b  ] - ray.y ) * idy;
        tymax = ( ybounds[1-b] - ray.y ) * idy;
        if( (tmin > tymax) || (tymin > tmax) ) { return false; }
        tmin = std::max( tmin, tymin );
        tmax = std::min( tmax, tymax );

        float zbounds[] = { (float)zmin, (float)zmax };

        b = (idz >= 0) ? 0 : 1;
        tzmin = ( zbounds[b  ] - ray.z ) * idz;
        tzmax = ( zbounds[1-b] - ray.z ) * idz;
        if( (tmin > tzmax) || (tzmin > tmax) ) { return false; }
        tmin = std::max( tmin, tzmin );
        tmax = std::min( tmax, tzmax );
        // output distance to nearest intersection point
        if( tmin > 0 )
        {
            // go just a little bit further to ensure we don't lie too close
            // to the boundary for it to be numerically outside
            tmin  *= 1.00001;
            ray.x += tmin * ray.dx;
            ray.y += tmin * ray.dy;
            ray.z += tmin * ray.dz;
            // output one mandatory, plus one optional value
            if( interval )
            {
                *interval = tmax - tmin;
            }
            return tmin;
        }
        else
        {
            return 0;
        }
    }

};




/*****************************************************************************
 * Forward declarations
 ****************************************************************************/
int Y_Prime_Wrapper(double, const double*, double*, void*);

template<typename T>
void zero_vector( std::vector<T>& vec, size_t minSize );

void cross( float a[], const float b[], const float c[] );
float norm( const float x[] );
void normalise( float x[] );



/*****************************************************************************
 * Namespaces
 ****************************************************************************/
using namespace std;
using namespace H5;
using namespace boost;


/*****************************************************************************
 * Main class
 ****************************************************************************/
class RasterSystem : boost::noncopyable
{

private: // members

    vector<Voxel>               _voxels;
    Numpy                       _rbf_lut;
    Numpy                       _rbf_integral_lut;

    BBox                        _bbox;

    Flann_index                 _flann;
    /*
    vector<int>                 _flann_indices_buf;
    Matrixi                     _flann_indices;
    vector<float>               _flann_distances_sq_buf;
    Matrixf                     _flann_distances_sq;
    vector<float>               _flann_datapoints;
    */

    // the radius MUST be at least as large as the voxel spacing
    // otherwise interpolation will look very strange
    float                       _rbf_radius_sq;
    float                       _rbf_radius;

    // path of ray currently being traced stored here. This is a member of a
    // singleton object, so only one instance exists - this is NOT threadsafe!
    // Only trace one ray at a time (unless you run multiple processes, in
    // which case each one loads its own copy of the RasterSystem lib)
    // MAXSTEPS x 3, each row is x,y,z
    // In CImg parlance, the array has size 3 x MAXSTEPS
    // if pathlength is N<MAXSTEPS steps, then rows N..MAXSTEPS-1 are garbage
    // This will be resampled (cubic spline) to uniform spacing by Solve_ODE

    Numpy                       _path;
    // ODE step results, these have nonuniform spacing
    vector<double>              _path_t;
    vector<double>              _path_x;
    vector<double>              _path_y;
    vector<double>              _path_z;
    // XYZ ray exit direction vector after having been traced
    Numpy                       _exit_dir;
    // indices/weights for a matrix row (sparse vector in CSC format)
    Numpyi                      _row_indices;
    Numpy                       _row_weights;

    // input ray to be traced, 6 element vector: x,y,z,dx,dy,dz, unit-length d
    Numpy                       _ray;

    // absolute error tolerance for determining ODE adaptive step size
    float                       _tol;

    // ODE (spatial) step size limits
    float                       _hmin;
    float                       _hmax;

    // controls ODE termination criteria
    float                       _tpad;

    // HDF5 output file
    scoped_ptr<H5File>          _h5;
    const int                   _h5_chunk_size;
    deque<float>                _h5_data_buf;       // matrix value
    deque<int32_t>              _h5_cols_buf;       // matrix col coord
    deque<int32_t>              _h5_rows_buf;       // matrix row coord
    deque<float>                _h5_exit_buf;       // ray exit dir (xyz tuples)
    int                         _h5_curRow;         // current row index
    int                         _h5_nRows;          // num rows written to disk so far
    int                         _h5_nnz;            // num nonzeros written to disk so far



private: // ctor

    // real constructor
    RasterSystem():
        _flann(NULL),
        _tol(1e-4),
        _hmin(0),
        _hmax(0),
        _tpad(2.0),
        _bbox(0,0,0,0,0,0),
        _rbf_radius(0),
        _rbf_radius_sq(0),
        /*_ode_T(gsl_odeiv_step_rkf45),*/

        _h5_chunk_size(8192),                       // 32Kb @ 4bytes each
        _h5_curRow(0),
        _h5_nRows(0),
        _h5_nnz(0)
    {

        const int dim       = 6;

    }


public: // ctor, dtor

    // singleton accessor
    static RasterSystem& Get_Instance()
    {
        static RasterSystem instance;
        return instance;
    }


    ~RasterSystem()
    {
        // write final partial chunks to disk
        Close_HDF5();
    }



public: // accessors



    int                             /// Calculate a Perceptual Hash for an image
                                    /// Populate a Digest struct with the result
                                    /// _ph_image_digest()
                                    /// see pHash.h, ln# 284
    Calculate_PHash_Digest(
        const Numpyi&   img,        /// uint8_t
        const double    sigma,
        const double    gamma,
              Digest    &digest,    /// see pHash.h, ln# 193
              int       N = 180
    ) {
                int     outout      = _ph_image_digest(img, sigma, gamma, digest, N);
                return  outout;
    }


    int                             /// Calculate the DCT of a Perceptual Hash
    Calculate_PHash_DCT(
        const Numpyi&   srcimg,     /// uint8_t
              ulong64   &hash       /// see pHash.h, ln# 360
    ) {

        CImg<float> meanfilter(7, 7, 1, 1, 1);
        CImg<float> img;

        if (srcimg.spectrum() == 3) {
            img = srcimg.get_RGBtoYCbCr().channel(0).get_convolve(meanfilter);

        } else if (srcimg.spectrum() == 4) {
            int width = img.width();
            int height = img.height();
            int depth = img.depth();
            img = srcimg.get_crop(0, 0, 0, 0, width-1, height-1, depth-1, 2).\
                        RGBtoYCbCr().get_channel(0).get_convolve(meanfilter);

        } else {
            img = srcimg.get_channel(0).get_convolve(meanfilter);
        }

        img.resize(32, 32);
        /*CImg<float> *C = ph_dct_matrix(32);*/
        CImg<float> *C = get_dct_matrix();
        CImg<float> Ctransp = C->get_transpose();

        CImg<float> dctImage = (*C)*img*Ctransp;

        CImg<float> subsec = dctImage.crop(1, 1, 8, 8).unroll('x');

        float median = subsec.median();
        ulong64 one = 0x0000000000000001;
        hash = 0x0000000000000000;
        for (int i = 0; i < 64; i++) {
            float current = subsec(i);
            if (current > median) {
                hash |= one;
            }
            one = one << 1;
        }

        delete C;
        return 0;
    }


    /*
    int
    Range_Search(
        const Voxel& voxel      // input query point
    ) {
        float query[] = { voxel.x, voxel.y, voxel.z };
        flann::Matrix<float> queryMatrix( query, 1, 3 );

        return _flann->radiusSearch(
                queryMatrix,
                _flann_indices, _flann_distances_sq,
                _rbf_radius_sq,
                flann::SearchParams(FLANN_CHECKS_UNLIMITED, 0, false)
        );
    }


    bool
    Interpolate(
        Voxel& voxel            // input query point & output result
    ) {
        int neighbours = Range_Search( voxel );
        if( neighbours > 0 )
        {
            float normFactor = 0;
            float lutSize = static_cast<float>( _rbf_lut.size() );
            voxel.rindex = 0;
            voxel.gx = voxel.gy = voxel.gz = 0;
            for( int i=0; i<neighbours; i++ )
            {
                // sum contributions from all neighbours within range
                float d = _flann_distances_sq.data[i];
                // all points should have d<_rbf_rad_sq by construction
                if( d < _rbf_radius_sq )
                {
                    // map d from distance squared in [0.._rad_sq] to distance in [0..1]
                    float dNorm   = sqrtf( d ) / _rbf_radius;
                    // lookup weight in precomputed table
                    float weight  = _rbf_lut( static_cast<int>(dNorm*lutSize) );
                    normFactor   += weight;
                    const Voxel& neighbour = _voxels[ _flann_indices.data[i] ];
                    voxel.rindex += weight * neighbour.rindex;
                    voxel.gx     += weight * neighbour.gx;
                    voxel.gy     += weight * neighbour.gy;
                    voxel.gz     += weight * neighbour.gz;
                }
                else
                {
                    cerr << "Warning: d >= _rbf_radius_sq" << endl;
                    voxel.rindex = 1;
                }
            }
            // normalise
            if( normFactor > 0 )
            {
                normFactor    = 1/normFactor;
                voxel.rindex *=   normFactor;
                voxel.gx     *=   normFactor;
                voxel.gy     *=   normFactor;
            }
            return true;
        }
        else
        {
            voxel.rindex = 1;
            voxel.gx = voxel.gy = voxel.gz = 0;
            return false;
        }
    }


    // This writes the exit direction vectors to HDF5. The reason we don't combine
    // this function with the below Write_Buffers_to_HDF5 is that the buffers grow
    // at different rates. The main buffers fill up with an entry for each nonzero
    // in a matrix row, whereas the exitbuf only gains 3 entries per matrix row.
    // Since we only want to flush to disk rarely (when buffers are full) we do
    // these writes at potentially different times.
    */
    void
    Write_ExitBuf_to_HDF5(
        const int N                         // num float values to write, N=3k
    ) {
        if( N > 0 )
        {
            assert( N % 3 == 0 );
            const int k    = N/3;
            const int rank = 2;
            hsize_t totalSize[] = { _h5_nRows + k, 3 };
            hsize_t offset[]    = { _h5_nRows,     0 };
            hsize_t ext[]       = { k,             3 };

            // allocate space for new chunk
            DataSet exit = _h5->openDataSet( "exit" );
            exit.extend( totalSize );

            // define output slice of h5 array
            DataSpace mspace( rank, ext );
            DataSpace fspace = exit.getSpace();
            fspace.selectHyperslab( H5S_SELECT_SET, ext, offset );

            // copy deque contents to contiguous memory
            deque<float>::iterator itf = _h5_exit_buf.begin();
            vector<float> tmpExitBuf( itf, itf+N );
            _h5_exit_buf.erase( itf, itf+N );

            // write contiguous array to h5 file
            PredType f32 = PredType::NATIVE_FLOAT;
            exit.write( &tmpExitBuf[0], f32, mspace, fspace );

            _h5_nRows = _h5_curRow;

            // update matrix dimension attributes
            Group root = _h5->openGroup( "/" );
            Attribute size = root.openAttribute( "size" );
            int sizeVal[] = { _h5_nRows, (int)_voxels.size() };
            size.write( PredType::NATIVE_INT32, sizeVal );
        }
    }


    // This writes the main coefficient matrix to HDF5. Is is the matrix where
    // each row contains the coefficients of an individual ray for each voxel.
    void
    Write_Buffers_to_HDF5(
        const int N                     // num values to write from the buffer
    ) {
        if( N > 0 )
        {
            hsize_t totalSize[] = { _h5_nnz + N };
            hsize_t offset[]    = { _h5_nnz     };
            hsize_t ext[]       = { N           };

            // allocate space for new chunk
            DataSet data = _h5->openDataSet( "data" );
            DataSet cols = _h5->openDataSet( "cols" );
            DataSet rows = _h5->openDataSet( "rows" );
            data.extend( totalSize );
            cols.extend( totalSize );
            rows.extend( totalSize );

            // define output slice of h5 array
            int rank = 1;
            DataSpace mspace( rank, ext );
            DataSpace fspace = data.getSpace();
            fspace.selectHyperslab( H5S_SELECT_SET, ext, offset );
            _h5_nnz += N;

            // copy deque contents to contiguous memory
            deque<float>::iterator itf = _h5_data_buf.begin();
            vector<float> tmpDataBuf( itf, itf+N );
            _h5_data_buf.erase( itf, itf+N );

            deque<int32_t>::iterator iti = _h5_cols_buf.begin();
            vector<int32_t> tmpColsBuf( iti, iti+N );
            _h5_cols_buf.erase( iti, iti+N );

            iti = _h5_rows_buf.begin();
            vector<int32_t> tmpRowsBuf( iti, iti+N );
            _h5_rows_buf.erase( iti, iti+N );

            // write contiguous array to h5 file
            PredType f32 = PredType::NATIVE_FLOAT;
            PredType i32 = PredType::NATIVE_INT32;
            data.write( &tmpDataBuf[0], f32, mspace, fspace );
            cols.write( &tmpColsBuf[0], i32, mspace, fspace );
            rows.write( &tmpRowsBuf[0], i32, mspace, fspace );
        }
    }
    /*


public: // main interaction methods

    // should only call this once, with the initial visual hull points
    void
    Insert_Points(
        const Numpy& points,        // N rows x 3 cols
        const Numpy& indices,       // N rows x 1 col
        const Numpy& gradients      // N rows x 3 cols
    ) {
        int nPoints = points.height();
        _flann_datapoints.clear();
        _flann_datapoints.reserve( points.size() );
        copy( points.data(), points.data()+points.size(), &_flann_datapoints[0] );
        flann::Matrix<float> xyz( &_flann_datapoints[0], nPoints, 3 );
        // KDTreeSingleIndex(1) faster than 10 or 16, or Linear or Autotuned(0.2%)
        _flann.reset( new Flann_index_base(xyz, flann::KDTreeSingleIndexParams(1)) );
        _flann->buildIndex();

        // Store complete internal list of points
        _voxels.clear();
        for( int i=0; i<nPoints; i++ )
        {
            const float x  = xyz[i][0];
            const float y  = xyz[i][1];
            const float z  = xyz[i][2];
            const float r  = indices(i);
            const float gx = gradients(0,i);
            const float gy = gradients(1,i);
            const float gz = gradients(2,i);
            _voxels.push_back( Voxel(i,x,y,z,r,gx,gy,gz) );
        }

        Construct_AABB();
    }


    // the ODE is autonomous (independent of time) but we need to pass in a dummy t
    // parameter (and extra params pointer) to satisfy the GSL function signature
    // this is a constant spatial step size parametrisation. See Art Tevs' MSc thesis
    // for an alternate, constant temporal step size, parametrisation
    int
    Y_Prime(
              double    t,
        const double    y[],
              double    dydt[],
              void*     params
    ) {
        // y holds current point and velocity
        static Voxel voxel;
        voxel.x = y[0];
        voxel.y = y[1];
        voxel.z = y[2];
        const float vx = y[3], vy = y[4], vz = y[5];

        // interpolate rindex and gradient at the current point
        Interpolate( voxel );

        // compute dydt
        // change in position := velocity scaled by rindex
        assert( voxel.rindex != 0 );
        float inv_denom = 1/voxel.rindex;
        dydt[0] = vx * inv_denom;
        dydt[1] = vy * inv_denom;
        dydt[2] = vz * inv_denom;
        // change in velocity := rindex gradient
        dydt[3] = voxel.gx;
        dydt[4] = voxel.gy;
        dydt[5] = voxel.gz;

        return GSL_SUCCESS;
    }


    int             // If in matrix-output mode (when Init_HDF5 as been called)
                    // return a "tuple" of two 16bit ints packed into a 32bit int
                    // the most significant is the number of (uniformly resampled)
                    // steps in the path. The least significant is the number of
                    // nonzeros in the matrix row
                    // Otherwise just return 0 if no voxels were hit and 1 if some
                    // voxels were hit
                    // This function also writes to the _exit_dir member
    Solve_ODE()
    {
        // caller should write to the _ray member before calling Solve_ODE and
        // should ensure that it starts somewhere before bounding box,
        // pointing towards it
        Ray ray(_ray);

        // for recording (spatial) step positions. These will be spaced unevenly
        // due to adaptive step sizes. They will be evenly resampled later.
        _path_t.clear();
        _path_x.clear();
        _path_y.clear();
        _path_z.clear();

        // initialise ray exit direction to ray entry direction
        _exit_dir(0) = ray.dx;
        _exit_dir(1) = ray.dy;
        _exit_dir(2) = ray.dz;

        // move ray to first intersection with bounding box
        double tmax;
        if( !_bbox.Advance_to_Intersection(ray,&tmax) )
        {
            // if ray missed bbox completely, append null row to matrix
            Append_HDF5( 0 );
            return 0;
        }

        // stepping parameters
        double t = 0;               // time
        double h = _hmin;           // step size
        // termination criteria: allow ray to travel further inside than the
        // straight line distance between bbox entry/exit points due to bending
        tmax *= _tpad;

        // initial value (GSL uses doubles, Ray uses floats, so must cast here)
        double y[6] = { (double)ray.x,  (double)ray.y,  (double)ray.z,
                        (double)ray.dx, (double)ray.dy, (double)ray.dz };

        // store initial ray position (on bbox boundary)
        _path_t.push_back( t    );
        _path_x.push_back( y[0] );
        _path_y.push_back( y[1] );
        _path_z.push_back( y[2] );

        // tracing
        bool inside = true;
        while( (t<tmax) && inside && (_path_t.size()<_path.height()) )
        {
            int status = gsl_odeiv_evolve_apply(
                    _ode_e, _ode_c, _ode_s, &_ode_sys,
                    &t, tmax, &h, y
            );
            assert( status == GSL_SUCCESS );

            // note that we do NOT renormalise the velocity vector here
            // see http://www.utdallas.edu/~cantrell/ee6334/lectures.html

            // limit step size to sensible range, to prevent it jumping right over
            // interesting regions, or shrinking to zero at sharp discontinuities.
            // Note that this does not actually prevent the solver from choosing a
            // smaller/larger step size if necessary, it just prevents extreme
            // sizes from being used as initial sizes for the next iteration
            h = GSL_MIN_DBL( GSL_MAX_DBL(h,_hmin), _hmax );

            // here is where we check for total internal reflection to terminate
            // ray early. See simulation cython code for hints. Currently it
            // doesn't work reliably enough to include, needs more testing
            // TODO:...

            // stop tracing once we leave the bounding box, clamp to it
            inside = _bbox.Inside( y[0], y[1], y[2] );
            if( !inside )
            {
                // store exit direction
                _exit_dir(0) = y[3];
                _exit_dir(1) = y[4];
                _exit_dir(2) = y[5];
                // turn around and back up to the bbox
                ray.x  =  y[0];   ray.y  =  y[1];    ray.z  =  y[2];
                ray.dx = -y[3];   ray.dy = -y[4];    ray.dz = -y[5];
                t -= _bbox.Advance_to_Intersection( ray );
                y[0] = ray.x;     y[1] = ray.y;      y[2] = ray.z;
                // t is advanced an extra epsilon towards bbox, but we must ensure
                // that _path_t remains monotonically increasing otherwise spline
                // initialisation will fail
                if( t <= _path_t.back() )
                {
                    t = _path_t.back() + 0.00001;
                }
            }

            // store history of ray states
            _path_t.push_back( t    );
            _path_x.push_back( y[0] );
            _path_y.push_back( y[1] );
            _path_z.push_back( y[2] );
        }

        // return true if ray successfully traced and exited bbox in reasonable time
        // returns number of [evenly resampled] steps if successful, 0 otherwise
        if( (_path_t.size()>1) && (t<tmax) && (!inside) )
        {
            // resample path and build matrix row if we're outputting the matrix
            if( _h5.get() )
            {
                const vector<double>* path_xyz[] = {&_path_x, &_path_y, &_path_z};
                Spline spline( _path_t );
                double tStep = _rbf_radius;
                int idx;
                // resample x then y then z to constant spatial step size
                for( int xyz=0; xyz<3; xyz++ )
                {
                    idx = 0;
                    spline.init( *path_xyz[xyz] );
                    for( double ti = _path_t.front(); ti < _path_t.back(); ti += tStep )
                    {
                        _path(xyz,idx++) = spline[ti];
                    }
                    // one more for the final t value, which may be less that a full
                    // tStep away from the penultimate one
                    _path(xyz,idx++) = spline[_path_t.back()];
                }

                // export to disk (if Init_HDF5 has been called)
                int nnzRow = Build_Matrix_Row( _path, idx );
                Append_HDF5( nnzRow );

                // return 2 integers from a single function by packing them into the
                // low and high 16bit words of a 32bit integer
                return (idx << 16) + nnzRow;
            }
            // otherwise skip that extra work and just return true
            // The exit direction will already have been written to the _exit_dir array
            else
            {
                return 1;
            }
        }
        else
        {
            // append null row to matrix and exit direction
            Append_HDF5( 0 );
            return 0;
        }
    }


    int                             // return num nonzeros in row
                                    // expect about 9*N nonzeros for NxNxN full voxel grid
                                    // or about 9*N^(1/3) for N voxels (with 1-ring RBF)
    Build_Matrix_Row(
        const Numpy& path,          // mPoints rows x 3 cols, worldspace
        const int    m
    ) {
        assert( m > 1 );
        float invRadius = 1 / _rbf_radius;
        vector<Pair> row;
        row.reserve(1000);
        Voxel query;
        float curDir[3], neighbourToPoint[3], u[3];
        // for each of the uniformly-spaced points in 'path'
        for( int idx=0; idx<m; idx++ )
        {
            // current point along path
            query.x = path(0,idx);
            query.y = path(1,idx);
            query.z = path(2,idx);
            // central difference approx to current ray direction
            int prevIdx = (idx > 0  ) ? idx-1 : idx;
            int nextIdx = (idx < m-1) ? idx+1 : idx;
            curDir[0] = path(0,nextIdx) - path(0,prevIdx);
            curDir[1] = path(1,nextIdx) - path(1,prevIdx);
            curDir[2] = path(2,nextIdx) - path(2,prevIdx);
            normalise( curDir );
            // find all voxels within range of this point
            int n = Range_Search( query );
            for( int nIdx=0; nIdx<n; nIdx++ )
            {
                // take one of those nearby voxels
                int label = _flann_indices.data[nIdx];
                const Voxel& neighbour = _voxels[label];
                // find perpendicular distance from neighbour voxel to the
                // line through curPoint, parallel to curDir
                // http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
                neighbourToPoint[0] = query.x - neighbour.x;
                neighbourToPoint[1] = query.y - neighbour.y;
                neighbourToPoint[2] = query.z - neighbour.z;
                cross( u, curDir, neighbourToPoint );
                float d = norm( u ) * invRadius; // [0..1]
                if( d > 1 )
                {
                    cerr << "Warning: d > 1 in Build_Matrix_Row" << endl;
                }
                // lookup pre-integrated weight
                int lutBin = static_cast<int>( d * _rbf_integral_lut.size() );
                float weight = _rbf_integral_lut(lutBin);
                // accumulate into row buffer (with duplicates)
                row.push_back( Pair(label,weight) );
            }
        }
        // remove duplicates in row buffer by coalescing entries with same label
        sort( row.begin(), row.end() );
        int outIdx = 0;
        float totalWeight = 0;
        for( int i=0,j=0; i<row.size(); i=j,outIdx++ )
        {
            // accumulate weights for all entries with same label
            float labelSum = row[i].weight;
            j = i + 1;
            while( (row[j].label==row[i].label) && (j<row.size()) )
            {
                labelSum += row[j++].weight;
            }
            // normalise single voxel
            labelSum    /= (j-i);
            totalWeight += labelSum;
            _row_indices(outIdx) = row[i].label;
            _row_weights(outIdx) = labelSum;
        }
        // normalise across entire row
        // this will set all rows to have weights that sum to 1.0
        // for absorption tomography, Anderson & Kak "SART" [1984] mention that
        // the sum of weights should equal the physical length of the ray. That
        // is ignored here because we have rays that don't refract at all
        // outside the reconstruction volume.
        // NB. Actually, check this, is it correct to ignore length of ray???
        float invTotalWeight = 1 / totalWeight;
        for( int i=0; i<outIdx; i++ )
        {
            _row_weights(i) *= invTotalWeight;
        }
        return outIdx;
    }


    void
    Init_HDF5(
        const char* filename
    ) {
        // close previous file
        Close_HDF5();

        int rank          =   1;
        hsize_t dims[]    = { 0             };
        hsize_t maxdims[] = { H5S_UNLIMITED };
        DataSpace space( rank, dims, maxdims );

        // overwrite if file exists
        _h5.reset( new H5File(filename,H5F_ACC_TRUNC) );

        // chunking mandatory for extendible datasets
        // don't bother setting fill value
        DSetCreatPropList cparms;
        hsize_t chunk_dims[] = { _h5_chunk_size };
        cparms.setChunk( rank, chunk_dims );
        cparms.setDeflate( 3 );
        cparms.setShuffle();

        // 1D vectors to store row/col/val tuples
        PredType f32 = PredType::NATIVE_FLOAT;
        PredType i32 = PredType::NATIVE_INT32;
        _h5->createDataSet( "data", f32, space, cparms );
        _h5->createDataSet( "cols", i32, space, cparms );
        _h5->createDataSet( "rows", i32, space, cparms );

        // exit direction vectors
        rank                  =   2;
        hsize_t dimsb[]       = { 0,                3 };
        hsize_t maxdimsb[]    = { H5S_UNLIMITED,    3 };
        hsize_t chunk_dimsb[] = { _h5_chunk_size/3, 3 };
        space = DataSpace( rank, dimsb, maxdimsb );
        cparms.setChunk( rank, chunk_dimsb );
        _h5->createDataSet( "exit", f32, space, cparms );

        // matrix dimension attributes
        Group root = _h5->openGroup( "/" );
        hsize_t dims2[] = { 2 };
        DataSpace vec2( 1, dims2 );
        IntType inttype( i32 );
        int zeros[] = { 0, 0 };
        Attribute size = root.createAttribute( "size", inttype, vec2 );
        size.write( i32, zeros );

        // the destructor will close the file if Close_HDF5 isn't called
        // but be sure not to try to modify the file in another process
        // (i.e. from python) while this class still holds an open handle
    }
    */

    void
    Close_HDF5()
    {
        if( _h5.get() )
        {
            // flush buffers to disk
            Write_ExitBuf_to_HDF5( _h5_exit_buf.size() );
            Write_Buffers_to_HDF5( _h5_data_buf.size() );

            // invalidate file parameters
            _h5_curRow = 0;
            _h5_nRows  = 0;
            _h5_nnz    = 0;
            _h5_data_buf.clear();
            _h5_cols_buf.clear();
            _h5_rows_buf.clear();
            _h5_exit_buf.clear();

            // clear pointer
            _h5.reset();
        }
    }

    /*
    // call this after tracing each ray. It maintains a 1:1 correspondence
    // between matrix rows and rays, even if the ray hits no voxels (nnz=0)
    // This function assumes that _exit_dir vector has already been set
    void
    Append_HDF5(
        const int nnz               // number nonzeros in current row buffer
    ) {
        // don't do anything if Init_HDF5 hasn't been called, ie if we're not
        // set to output the matrix to disk
        if( _h5.get() )
        {
            // append to buffers (do nothing if nnz=0, ie ray missed all voxels)
            float* ptrf = _row_weights.data();
            int*   ptri = _row_indices.data();
            _h5_data_buf.insert( _h5_data_buf.end(), ptrf, ptrf+nnz     );
            _h5_cols_buf.insert( _h5_cols_buf.end(), ptri, ptri+nnz     );
            _h5_rows_buf.insert( _h5_rows_buf.end(), nnz,  _h5_curRow++ );

            // copy exit direction from _exit_dir vector
            normalise( _exit_dir.data() );
            ptrf = _exit_dir.data();
            _h5_exit_buf.insert( _h5_exit_buf.end(), ptrf, ptrf+3       );

            // write to disk
            if( _h5_exit_buf.size() >= _h5_chunk_size )
            {
                Write_ExitBuf_to_HDF5( _h5_exit_buf.size() );
            }
            if( _h5_data_buf.size() >= _h5_chunk_size )
            {
                Write_Buffers_to_HDF5( _h5_data_buf.size() );
            }
        }
    }
    */

};



// interface to expose class methods to Python via C API (using ctypes)
RasterSystem* asrt(void* ptr)
{
    return static_cast<RasterSystem*>(ptr);
}
extern "C"
{
    void* RasterSystem_Get_Instance(void)
    {
        return &(RasterSystem::Get_Instance());
    }

    void* RasterSystem_Calculate_PHash_Digest(void* ptr,
        const unsigned char* numpyImg,
        double sigma, double gamma,
        Digest digest, int N
    ) {
        const Numpyi tmp( numpyImg, 3, 1, 1, 1, true );
        asrt(ptr)->Calculate_PHash_Digest( tmp, sigma, gamma, digest, N );
    }
    void* RasterSystem_Calculate_PHash_DCT(void* ptr,
        const unsigned char* numpySrcImg,
        unsigned long long hash
    ) {
        const Numpyi tmp( numpySrcImg, 3, 1, 1, 1, true );
        asrt(ptr)->Calculate_PHash_DCT( tmp, hash );
    }

    /*
    void RasterSystem_Bind_Path(void* ptr,
        const float* numpyPath, const int maxSteps   // maxSteps rows x 3 cols
    ) {
        const Numpy tmp( numpyPath, 3, maxSteps, 1, 1, true );
        asrt(ptr)->Set_path( tmp );
    }

    void RasterSystem_Bind_Ray(void* ptr,
        const float* ray                            // 6 element vector
    ) {
        const Numpy tmp( ray, 6, 1, 1, 1, true );
        asrt(ptr)->Set_ray( tmp );
    }

    void RasterSystem_Bind_Exit_Dir(void* ptr,
        const float* numpyDir                       // 3 element vector
    ) {
        const Numpy tmp( numpyDir, 3, 1, 1, 1, true );
        asrt(ptr)->Set_exit_dir( tmp );
    }

    void RasterSystem_Bind_Row_Indices(void* ptr,
        const int* numpyIndices, const int n        // n element vector
    ) {
        const Numpyi tmp( numpyIndices, n, 1, 1, 1, true );
        asrt(ptr)->Set_row_indices( tmp );
    }

    void RasterSystem_Bind_Row_Weights(void* ptr,
        const float* numpyWeights, const int n      // n element vector
    ) {
        const Numpy tmp( numpyWeights, n, 1, 1, 1, true );
        asrt(ptr)->Set_row_weights( tmp );
    }

    void RasterSystem_Bind_RBF_LUT(void* ptr,
        const float* numpyLUT, const int n          // n element vector
    ) {
        const Numpy tmp( numpyLUT, n, 1, 1, 1, true );
        asrt(ptr)->Set_rbf_lut( tmp );
    }

    void RasterSystem_Bind_RBF_Integral_LUT(void* ptr,
        const float* numpyLUT, const int n          // n element vector
    ) {
        const Numpy tmp( numpyLUT, n, 1, 1, 1, true );
        asrt(ptr)->Set_rbf_integral_lut( tmp );
    }

    // in world units (ie inches)
    void RasterSystem_Set_RBF_Radius(void* ptr,
        const float radius                          // note: NOT squared
    ) {
        asrt(ptr)->Set_rbf_radius( radius );
        asrt(ptr)->Set_hmin( radius / 3 );
        asrt(ptr)->Set_hmax( radius * 3 );
    }

    void RasterSystem_Init_HDF5(void* ptr,
        const char* filename                        // fully qualified path
    ) {
        asrt(ptr)->Init_HDF5( filename );
    }
    */

    void RasterSystem_Close_HDF5(void* ptr)
    {
        asrt(ptr)->Close_HDF5();
    }

    /*
    void RasterSystem_Insert_Points(void* ptr,
        const float* numpyPoints, const int nPoints,// nPoints rows x 3 cols
        const float* rindices,                      // nPoints rows x 1 col
        const float* gradients                      // nPoints rows x 3 cols
    ) {
        const Numpy pts(  numpyPoints, 3, nPoints,  1, 1, true );
        const Numpy ind(  rindices,    nPoints, 1,  1, 1, true );
        const Numpy grad( gradients,   3, nPoints,  1, 1, true );
        asrt(ptr)->Insert_Points( pts, ind, grad );
    }

    int RasterSystem_Solve_ODE(void* ptr
    ) {
        return asrt(ptr)->Solve_ODE();
    }

    // these are not part of the public API
    // just here for debugging

    int RasterSystem_Build_Matrix_Row(void* ptr,
        const float* numpyPoints, const int mPoints
    ) {
        const Numpy tmp( numpyPoints, 3, mPoints, 1, 1, true );
        return asrt(ptr)->Build_Matrix_Row( tmp, mPoints );
    }

    int RasterSystem_Range_Search(void* ptr,
        const float x, const float y, const float z
    ) {
        const Voxel voxel( 0, x,y,z, 0, 0,0,0 );
        return asrt(ptr)->Range_Search( voxel );
    }

    float RasterSystem_Interpolate(void* ptr,
        const float x, const float y, const float z,
        float* gradient
    ) {
        Voxel voxel( 0, x,y,z, 0, 0,0,0 );
        asrt(ptr)->Interpolate( voxel );
        gradient[0] = voxel.gx;
        gradient[1] = voxel.gy;
        gradient[2] = voxel.gz;
        return voxel.rindex;
    }

    float RasterSystem_Get_RBF_Radius(void* ptr)
    {
        return asrt(ptr)->Get_rbf_radius();
    }

    float RasterSystem_Get_HMin(void* ptr)
    {
        return asrt(ptr)->Get_hmin();
    }

    float RasterSystem_Get_HMax(void* ptr)
    {
        return asrt(ptr)->Get_hmax();
    }
    */
}
