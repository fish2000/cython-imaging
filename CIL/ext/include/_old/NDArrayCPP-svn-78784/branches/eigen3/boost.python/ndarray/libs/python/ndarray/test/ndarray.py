import ndarray_mod
import unittest
import numpy

class TestNdArrayWrappers(unittest.TestCase):

    def testReturnArray(self):
        for suffix in ("11","10","22","21","20","33","32","31","30"):
            func = getattr(ndarray_mod, "returnArray_d%s" % suffix)
            cfunc = getattr(ndarray_mod, "returnArray_dc%s" % suffix)
            array = func()
            carray = cfunc()
            self.assert_((array == numpy.arange(0, array.size).reshape(array.shape)).all())
            self.assert_((carray == numpy.arange(0, carray.size).reshape(carray.shape)).all())
            self.assert_(array.flags["WRITEABLE"])
            self.assert_(not carray.flags["WRITEABLE"])

    def testAcceptArray(self):
        invalid_examples = {
            "11": (numpy.zeros((12,), dtype=float)[::2],
                   ),
            "10": (),
            "22": (numpy.zeros((6, 5), dtype=float)[::2,:], 
                   numpy.zeros((6, 8), dtype=float)[::2,::2], 
                   numpy.zeros((5, 8), dtype=float)[:,::2], 
                   ),
            "21": (numpy.zeros((6, 8), dtype=float)[::2,::2], 
                   numpy.zeros((5, 8), dtype=float)[:,::2], 
                   ),
            "20": (),
            "33": (), "32": (), "31": (), "30": (),
            }
        valid_examples = {
            "10": (numpy.zeros((12,), dtype=float)[::2],
                   ),
            "11": (),
            "20": (numpy.zeros((6, 5), dtype=float)[::2,:], 
                   numpy.zeros((6, 8), dtype=float)[::2,::2], 
                   numpy.zeros((5, 8), dtype=float)[:,::2], 
                   ),
            "21": (numpy.zeros((6, 5), dtype=float)[::2,:], 
                   ),
            "22": (),
            "33": (), "32": (), "31": (), "30": (),
            }
        for suffix in ("11","10","22","21","20","33","32","31","30"):
            func = getattr(ndarray_mod, "acceptArray_d%s" % suffix)
            cfunc = getattr(ndarray_mod, "acceptArray_dc%s" % suffix)
            nd = int(suffix[0])
            shape = tuple(int(i) for i in numpy.random.randint(low=2,high=5,size=nd))
            array = numpy.zeros(shape, dtype=float)
            array[:] = numpy.arange(0, array.size).reshape(array.shape)
            self.assert_(func(array))
            self.assert_(cfunc(array))
            for array in invalid_examples[suffix]:
                array[:] = numpy.arange(array.size, dtype=float).reshape(*array.shape)
                self.assertRaises(ValueError, func, array)
                self.assert_(cfunc(array))
            for array in valid_examples[suffix]:
                array[:] = numpy.arange(array.size, dtype=float).reshape(*array.shape)
                self.assert_(func(array))
                self.assert_(cfunc(array))
            

    def testAcceptArrayVal(self):
        for suffix in ("11","10","22","21","20","33","32","31","30"):
            func = getattr(ndarray_mod, "acceptArrayVal_d%s" % suffix)
            cfunc = getattr(ndarray_mod, "acceptArrayVal_dc%s" % suffix)
            nd = int(suffix[0])
            shape = tuple(int(i) for i in numpy.random.randint(low=2,high=5,size=nd))
            array = numpy.zeros(shape, dtype=float)
            array[:] = numpy.arange(0, array.size).reshape(array.shape)
            self.assert_(func(array))
            self.assert_(cfunc(array))

    def testExtractArray(self):
        for suffix in ("11","10","22","21","20","33","32","31","30"):
            func = getattr(ndarray_mod, "extractArray_d%s" % suffix)
            cfunc = getattr(ndarray_mod, "extractArray_dc%s" % suffix)
            nd = int(suffix[0])
            shape = tuple(int(i) for i in numpy.random.randint(low=2,high=5,size=nd))
            array = numpy.zeros(shape, dtype=float)
            array[:] = numpy.arange(0, array.size).reshape(array.shape)
            self.assert_(func(array))
            self.assert_(cfunc(array))

    def testReturnVector(self):
        for suffix in ("2","3"):
            func = getattr(ndarray_mod, "returnVector_d%s" % suffix)
            vector = func()
            self.assertEqual(vector, tuple(numpy.arange(len(vector), dtype=float)))

    def testAcceptVector(self):
        for suffix in ("2","3"):
            func = getattr(ndarray_mod, "acceptVector_d%s" % suffix)
            nd = int(suffix[0])
            vector = tuple(numpy.arange(nd, dtype=float))
            self.assert_(func(vector))

    def _testMemory(self):
        shape = (400, 400, 10)
        for n in range(1000000):
            a = ndarray_mod.makeArray_d33(shape)

if __name__=="__main__":
    unittest.main()
