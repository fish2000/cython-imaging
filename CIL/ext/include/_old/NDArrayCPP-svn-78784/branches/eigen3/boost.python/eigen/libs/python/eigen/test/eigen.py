import eigen_mod
import unittest
import numpy

class TestEigenWrappers(unittest.TestCase):

    def setUp(self):
        self.matrix_i = numpy.matrix([[1,2,3],[4,5,6]], dtype=int)
        self.matrix_d = numpy.matrix([[1,2,3],[4,5,6]], dtype=float)
        self.vector_i = numpy.array([1,2,3,4], dtype=int)
        self.vector_d = numpy.array([1,2,3,4], dtype=float)

    def testAcceptMatrix(self):
        self.assert_(eigen_mod.acceptMatrix_23d_cref(self.matrix_d))
        self.assert_(eigen_mod.acceptMatrix_X3d_cref(self.matrix_d))
        self.assert_(eigen_mod.acceptMatrix_2Xd_cref(self.matrix_d))
        self.assert_(eigen_mod.acceptMatrix_XXd_cref(self.matrix_d))

    def testAcceptVector(self):
        self.assert_(eigen_mod.acceptVector_41d_cref(self.vector_d))
        self.assert_(eigen_mod.acceptVector_X1d_cref(self.vector_d))
        self.assert_(eigen_mod.acceptVector_14d_cref(self.vector_d))
        self.assert_(eigen_mod.acceptVector_1Xd_cref(self.vector_d))

    def testReturnMatrix(self):
        self.assert_((eigen_mod.returnMatrix_23d() == self.matrix_d).all())
        self.assert_((eigen_mod.returnMatrix_X3d() == self.matrix_d).all())
        self.assert_((eigen_mod.returnMatrix_2Xd() == self.matrix_d).all())
        self.assert_((eigen_mod.returnMatrix_XXd() == self.matrix_d).all())
        self.assert_((eigen_mod.returnMatrix_23d_c() == self.matrix_d).all())
        self.assert_((eigen_mod.returnMatrix_X3d_c() == self.matrix_d).all())
        self.assert_((eigen_mod.returnMatrix_2Xd_c() == self.matrix_d).all())
        self.assert_((eigen_mod.returnMatrix_XXd_c() == self.matrix_d).all())

    def testReturnObject(self):
        self.assert_((eigen_mod.returnObject_23d() == self.matrix_d).all())
        self.assert_((eigen_mod.returnObject_X3d() == self.matrix_d).all())
        self.assert_((eigen_mod.returnObject_2Xd() == self.matrix_d).all())
        self.assert_((eigen_mod.returnObject_XXd() == self.matrix_d).all())

    def testMatrixOwner(self):
        for suffix in ("23d", "X3d", "2Xd", "XXd"):
            cls = getattr(eigen_mod, "MatrixOwner_%s" % suffix)
            obj = cls()
            self.assert_(obj.getMatrix_ref().base is obj)
            self.assert_(obj.getMatrix_cref().base is obj)
            self.assert_((obj.getMatrix_ref() == self.matrix_d).all())
            self.assert_((obj.getMatrix_cref() == self.matrix_d).all())
            self.assert_(obj.getMatrix_ref().flags["WRITEABLE"])
            self.assert_(not obj.getMatrix_cref().flags["WRITEABLE"])
            self.assert_(not obj.getMatrix_ref().flags["OWNDATA"])
            self.assert_(not obj.getMatrix_cref().flags["OWNDATA"])

    def testSqueeze(self):
        for suffix in ("41d", "14d", "X1d", "1Xd"):
            func1 = getattr(eigen_mod, "returnVector_%s_sq" % suffix)
            self.assert_((func1() == self.vector_d).all())
            self.assertEqual(type(func1()), numpy.ndarray)
            self.assertEqual(func1().ndim, 1)
            func2 = getattr(eigen_mod, "returnVector_%s" % suffix)
            self.assert_((func1() == self.vector_d).all())
            self.assertEqual(type(func2()), numpy.matrix)
            self.assertEqual(func2().ndim, 2)

if __name__=="__main__":
    unittest.main()
