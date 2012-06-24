#!/usr/bin/env python
#
# NB: see warning regarding matrix row/column ordering in raytracer.cpp file header
#
# The matrix is stored in sparse format.
# For each row we have an (index,value) pair.



import numpy as np
from os.path import isdir, join, basename, dirname, abspath, splitext
from os import getcwd, chdir
from scipy.integrate import quadrature
from scipy import special
from h5py import File

from datetime import datetime
from glob import glob
from optparse import OptionParser
import sys


ipython = hasattr(sys, 'ipcompleter')
if ipython:
    from IPython.Debugger import Tracer
    dbg = Tracer()
else:
    def dbg():
        print("Run inside ipython to activate debugging")

def log_verbose(msg, force=True):
    print( "{} {}".format(datetime.now(),msg) )
def log_quiet(msg, force=False):
    if force:
        log_verbose( msg )
    pass
def logerr(msg):
    log_verbose( "ERROR: " + msg )
def logwarn(msg):
    log_verbose( "WARNING: " + msg )
log = log_quiet




#############################################################################################
# C++ interface glue
#############################################################################################
import ctypes
from ctypes import c_void_p as cvoidp, c_int as cint, c_float as cfloat
from ctypes import c_char_p as ccharp
#from ctypes import c_bool as cbool
#from ctypes import c_uint64 as cuint64
from ctypes import c_uint8 as cuint8
from ctypes import c_double as cdouble
from ctypes import POINTER
from functools import partial

ptrUInt8_3D = np.ctypeslib.ndpointer(dtype=np.uint8,   ndim=3, flags='CONTIGUOUS')
ptrUInt8_2D = np.ctypeslib.ndpointer(dtype=np.uint8,   ndim=2, flags='CONTIGUOUS')
ptrUInt8    = np.ctypeslib.ndpointer(dtype=np.uint8,   ndim=1, flags='CONTIGUOUS')
ptrInt32_3D = np.ctypeslib.ndpointer(dtype=np.int32,   ndim=3, flags='CONTIGUOUS')
ptrInt32_2D = np.ctypeslib.ndpointer(dtype=np.int32,   ndim=2, flags='CONTIGUOUS')
ptrInt32    = np.ctypeslib.ndpointer(dtype=np.int32,   ndim=1, flags='CONTIGUOUS')
ptrFloat_3D = np.ctypeslib.ndpointer(dtype=np.float32, ndim=3, flags='CONTIGUOUS')
ptrFloat_2D = np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='CONTIGUOUS')
ptrFloat    = np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='CONTIGUOUS')

ptrUInt64_3D = np.ctypeslib.ndpointer(dtype=np.uint64,   ndim=3, flags='CONTIGUOUS')
ptrUInt64_2D = np.ctypeslib.ndpointer(dtype=np.uint64,   ndim=2, flags='CONTIGUOUS')
ptrUInt64    = np.ctypeslib.ndpointer(dtype=np.uint64,   ndim=1, flags='CONTIGUOUS')

class CPlusPlusInterface(object):
    
    def bindfunc(self, name, restype=None, argtypes=[]):
        mangledName = self.__class__.__name__ + '_' + name
        self._dll.__setattr__(mangledName, self._dll.__getattr__(mangledName))
        
        func = self._dll.__getattr__(mangledName)
        func.restype = restype
        func.argtypes = [cvoidp] + argtypes
        
        self.__setattr__(name, partial(func, self.obj))
        self._name = {
            'name': name,
            'mangledName': mangledName}



#############################################################################################
# Example external C++ library interface
#############################################################################################

'''

class ExampleWrapper(CPlusPlusInterface):

    def __init__(self):
        self._dll = ctypes.cdll.LoadLibrary('./build/libexample.dylib')
        self._dll.Example_New.restype = cvoidp
        self.obj = self._dll.Example_New()

        self.bindfunc('Release')
        self.bindfunc('Foo', None, [ccharp])
        self.bindfunc('Sum_Array', cfloat, [ptrFloat_3D, cint, cint, cint])
        self.bindfunc('Get_Cached_Sum', cfloat)

        self._Sum_Array = self.Sum_Array
        self.Sum_Array = lambda X: self._Sum_Array(np.ascontiguousarray(X), *X.shape)

'''

try:
    if __file__:
        pass
except:
    from os import environ
    __file__ = environ.get('PWD')


#############################################################################################
# External C++ library interface
#############################################################################################
RS = None
class RasterSystem(CPlusPlusInterface):

    class Digest(ctypes.Structure):
        _fields_ = [("id",      POINTER(cint)),
                    ("coeffs",  POINTER(cuint8)),
                    ("size",    cint)]

    def __init__(self, dataDir):
        # currently need to be in the raytracer dir for this to work
        chdir(abspath(dirname(__file__)))
        self._dll = ctypes.cdll.LoadLibrary(join(getcwd(), 'ext/librastersystem.so'))
        self._dll.RasterSystem_Get_Instance.restype = cvoidp
        self.obj = self._dll.RasterSystem_Get_Instance()

        # setters
        self.bindfunc( 'Bind_Path',             None,
                        [ptrFloat_2D, cint] )
        self.bindfunc( 'Bind_Ray',              None,
                        [ptrFloat] )
        self.bindfunc( 'Bind_Exit_Dir',         None,
                        [ptrFloat] )
        self.bindfunc( 'Bind_Row_Indices',      None,
                        [ptrInt32, cint] )
        self.bindfunc( 'Bind_Row_Weights',      None,
                        [ptrFloat, cint] )
        self.bindfunc( 'Bind_RBF_LUT',          None,
                        [ptrFloat, cint] )
        self.bindfunc( 'Bind_RBF_Integral_LUT', None,
                        [ptrFloat, cint] )
        self.bindfunc( 'Set_RBF_Radius',        None,
                        [cfloat] )

        # main interaction functions
        self.bindfunc( 'Init_HDF5',             None,
                        [ccharp] )
        self.bindfunc( 'Close_HDF5',            None,
                        [] )
        self.bindfunc( 'Insert_Points',         None,
                        [ptrFloat_2D, cint, ptrFloat, ptrFloat_2D] )
        self.bindfunc( 'Solve_ODE',             cint,
                        [] )

        # MY SHIT
        self.bindfunc( 'Calculate_PHash_Digest',    None,
                        [ptrUInt8_3D, cdouble, cdouble, self.Digest(), cint] )
        self.bindfunc( 'Calculate_PHash_DCT',       None,
                        [ptrUInt8_3D, ] )

        # for debugging
        self.bindfunc( 'Build_Matrix_Row',      cint,
                        [ptrFloat_2D, cint] )
        self.bindfunc( 'Range_Search',          cint,
                        [cfloat, cfloat, cfloat] )
        self.bindfunc( 'Interpolate',           cfloat,
                        [cfloat, cfloat, cfloat, ptrFloat] )
        self.bindfunc( 'Get_RBF_Radius',        cfloat,
                        [] )
        self.bindfunc( 'Get_HMin',              cfloat,
                        [] )
        self.bindfunc( 'Get_HMax',              cfloat,
                        [] )

        # init ray in/out buffers
        self.ray = np.zeros( 6, dtype=np.float32 )
        self.Bind_Ray( self.ray )
        self.exitDir = np.zeros( 3, dtype=np.float32 )
        self.Bind_Exit_Dir( self.exitDir )

        # init path
        # bad things will happen if these buffers overflow in the C++ code
        BUFSIZE = 4096
        self.path = np.zeros( (BUFSIZE,3), dtype=np.float32 )
        self.Bind_Path( self.path, self.path.shape[0] )
        self.rowIndices = np.zeros( BUFSIZE, dtype=np.int32 )
        self.Bind_Row_Indices( self.rowIndices, self.rowIndices.shape[0] )
        self.rowWeights = np.zeros( BUFSIZE, dtype=np.float32 )
        self.Bind_Row_Weights( self.rowWeights, self.rowWeights.shape[0] )

        # init voxels and rays
        assert isdir( dataDir ), "Invalid data directory"
        self.dataDir = dataDir
        voxelFile  = join( dataDir, 'voxels.h5' )
        #rayFile    = join( dataDir, 'rays.h5' )

        # get interVoxelSpace of finest resolution grid
        log( "Loading " + voxelFile )
        with File( voxelFile, 'r' ) as dat:
            self.interVoxelSpace = sorted( map(float,dat.keys()) )[0]
            group = dat[str(self.interVoxelSpace)]
            try:
                pos       = group['pos'][:]
                rindex    = group['rindex'][:]
                gradients = group['gradients'][:]
                if pos.shape != gradients.shape or len(pos) != len(rindex):
                    logerr( "Mismatching voxel dataset sizes" )
                    return
            except KeyError, err:
                logerr( "{} in group '{}'".format(err.message,group.name) )
                return
            log( "Initialising voxels for raytracer" )
            self.insert_points( pos, rindex, gradients )

        # init RBF (radial basis function / filter kernel)
        nBins = 256
        self.radius = self.select_rbf_radius()
        (self.rbf_lut, self.rbf_integral_lut) = self.get_rbf( nBins )
        self.Bind_RBF_LUT( self.rbf_lut, nBins )
        self.Bind_RBF_Integral_LUT( self.rbf_integral_lut, nBins )
        self.Set_RBF_Radius( self.radius )


    def init_traced_file(self, cam):
        assert '{:03d}'.format(int(cam)) == cam
        self.tracedFile = join( self.dataDir, 'traced_{}.h5'.format(cam) )
        self.Init_HDF5( self.tracedFile )


    def kaiser_bessel_filter(self, resolution, alpha=2.0):
        """Tabled 1D function from 1 to [slightly above] 0 in 'resolution' steps"""
        # let radius = 1 unit
        samples = np.linspace( 0, 1, resolution )
        def kb(d):
            """Value of KaiserBessel filter of radius 1 at distance d from centre"""
            bessel = lambda y: special.iv(0,y)
            return bessel( np.pi*alpha*np.sqrt(1-d**2) ) / bessel( np.pi*alpha )
        rbf_lut = kb( samples ).astype(np.float32)
        def integrand(d, y):
            """kb filter value at distance y away from midpoint of a chord
            d world units away (orthogonally) from basis function centre"""
            hypot = np.sqrt( d**2 + y**2 )
            return kb( hypot )
        rbf_integral_lut = np.empty( resolution, dtype=np.float32 )
        for (i,d) in enumerate( samples ):
            ymax = np.sqrt( 1 - d**2 )
            f = lambda y: integrand( d, y )
            rbf_integral_lut[i] = quadrature( f, -ymax,ymax, vec_func=False )[0]
        return (rbf_lut, rbf_integral_lut/rbf_integral_lut[0])


    def get_rbf(self, resolution, alpha=2.0):
        """Load existing RBF kernel from disk, otherwise generate it"""
        lutFilename = join( self.dataDir, 'rbf.h5' )
        try:
            # load existing kernel from disk
            with File( lutFilename, 'r' ) as dat:
                log( "Loading RBF kernel" )
                rbf_lut = dat['rbf_lut'][:]
                rbf_integral_lut = dat['rbf_integral_lut'][:]
                if rbf_lut.shape != rbf_integral_lut.shape:
                    raise ValueError( "LUT resolution mismatch" )
                if rbf_lut.size != resolution:
                    raise ValueError( "LUT resolution mismatch" )
        except (IOError, ValueError):
            # if it didn't already exist in correct size, create data and save
            (rbf_lut, rbf_integral_lut) = self.kaiser_bessel_filter( resolution, alpha )
            with File( lutFilename, 'w' ) as dat:
                log( "Precomputing RBF kernel" )
                def insert(name, data):
                    dat.create_dataset( name, data=data, compression=3, shuffle=True )
                insert( 'rbf_lut', rbf_lut )
                insert( 'rbf_integral_lut', rbf_integral_lut )
                dat.attrs['timestamp'] = str( datetime.now() )
        return (rbf_lut, rbf_integral_lut)


    def select_rbf_radius(self):
        """Find reasonable size for RBF (width in inches)"""
        # want to reach the voxel diagonally adjacent to a voxel, but not overlap
        # the one two voxels off along a single dimension
        minFactor = np.sqrt(3)
        maxFactor = 2
        factor = (minFactor + maxFactor) / 2
        log( "InterVoxel:{}, RBF_Radius:{}".format(self.interVoxelSpace, self.interVoxelSpace*factor) )
        return self.interVoxelSpace * factor


    def insert_points(self, points, rindices, gradients):
        """points = Nx3, will be labelled sequentially"""
        ascont = lambda arr: np.ascontiguousarray( arr, dtype=np.float32 )
        self.Insert_Points( ascont(points), points.shape[0],
                            ascont(rindices), ascont(gradients) )


    def solve_ode(self, ray):
        """ray = 1x6 vector (x,y,z,dx,dy,dz), unit length dir"""
        # passing ctypes arrays to C++ funcs incurs some nontrivial overhead
        # so rather than passing ray as an argument, write its values directly
        # to memory shared by Python/C++ and then call Solve_ODE without args
        self.ray[:] = ray[:]
        pair = self.Solve_ODE()
        # decode the 'tuple' of two 16bit ints packed into single 32bit int
        pathLength = pair >> 16
        #nnz = pair - (pathLength << 16)
        if pathLength > 0:
            # return angle difference between ingoing and exiting directions (radians)
            dot = np.dot( ray[3:], self.exitDir )
            return np.arccos( dot ) if abs(dot)<=1 else 0
            # NOTE there may be something wrong with computed delta.
            # they are always zero or occassionally 0.0197823 (homogeneous voxels)
        else:
            return 0


    def build_matrix_row(self, path):
        """path = mx3 matrix (x,y,z)
        m = num steps inside bbox after resampling
        After this, self.path and self.indices are sparse vector of points"""
        path = np.ascontiguousarray( path, dtype=np.float32 )
        return self.Build_Matrix_Row( path, path.shape[0] )


    # for debugging
    def range_search(self, x, y, z):
        return self.Range_Search( x, y, z )


    # for debugging
    def interpolate(self, x, y, z):
        gradient = np.empty( 3, dtype=np.float32 )
        rindex = self.Interpolate( x,y,z, gradient )
        return (rindex, gradient)






def main(cam):
    """Trace and output results for all rays belonging to this camera"""
    # note that we intentionally output traced rays to multiple files rather
    # then keeping them more neatly as groups within a single HDF5 file. This
    # is because raytracing is the slowest operation and we may want to
    # perform it on multiple computers in parallel, so they cannot all be
    # writing to the same file simultaneously.
    RS.init_traced_file(cam)
    # now trace each camera
    rayFile = join( RS.dataDir, 'rays.h5' )
    log( "Loading " + rayFile )
    with File( rayFile, 'r' ) as dat:
        if cam not in dat.keys():
            logerr( "Missing '{}' in '{}'".format(cam,dat.name) )
            return
        group = dat[cam]
        log( "Tracing camera " + cam )
        origin    = group['origin']
        direction = group['direction']
        rays = np.hstack( [origin, direction] )
        delta = np.rad2deg( np.array(map(RS.solve_ode,rays),dtype=np.float32) )

    # now write output to disk
    # don't open file with 'w' otherwise you clobber the existing data
    RS.Close_HDF5()
    log( "Saving data to " + RS.tracedFile )
    with File( RS.tracedFile, 'r+' ) as dat:
        dat.attrs['timestamp'] = str( datetime.now() )
        def insert(name, data):
            dat.create_dataset( name, data=data, compression=3, shuffle=True )
        insert( 'degreesTraced', delta )

    log( "Done tracing " + cam )



if __name__ == '__main__':

    cmd = OptionParser()
    cmd.add_option("--dataDir", default=None,
                   type="string", metavar="DIR",
                   help="location of data files to process [%default]")
    cmd.add_option("--camera", default=-1,
                   type="int", metavar="N",
                   help="camera angle to trace (-1 for all) [%default]")
    cmd.add_option("-v", "--verbose", default=False,
                   action="store_true",
                   help="output additional progress information [%default]")
    (opt,args) = cmd.parse_args()


    if opt.verbose:
        log = log_verbose

    RS = RasterSystem( opt.dataDir )

    if opt.camera == -1:
        cameras = [splitext( basename(c) )[0].split( '_' )[-1]
                   for c in glob( join(opt.dataDir,'data_far_???.h5') )]
        for cam in cameras:
            main( cam )
    else:
        cam = '{:03d}'.format( opt.camera )
        main( cam )
