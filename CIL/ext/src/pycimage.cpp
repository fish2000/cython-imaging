
#include <Python.h>
#include "pycimage.h"

typedef struct {
    PyObject_HEAD

} pycimage_PyCImageObject;



PyMODINIT_FUNC
initpHash(void) {
    PyObject *m;

    m = Py_InitModule("pHash", pHashMethods);
    if(m == NULL)
	return;

    pHashError = PyErr_NewException("pHash.error", NULL, NULL);
    Py_INCREF(pHashError);
    PyModule_AddObject(m, "error", pHashError);
}

static PyObject *
phash_imagehash(PyObject *self, PyObject *args) {
    const char *filename;
    ulong64 hash = 0;

    if(!PyArg_ParseTuple(args, "s", &filename))
	return NULL;
    ph_dct_imagehash(filename, hash);
    return PyLong_FromUnsignedLongLong(hash);
}

static PyObject *
phash_hamming_distance(PyObject *self, PyObject *args) {
    ulong64 hash1, hash2;
    PyObject *py_hash1, *py_hash2;
    int ret;

    if(!PyArg_ParseTuple(args, "OO", &py_hash1, &py_hash2))
	    return NULL;
    hash1 = PyLong_AsUnsignedLongLong(py_hash1);
    hash2 = PyLong_AsUnsignedLongLong(py_hash2);
    ret = ph_hamming_distance(hash1, hash2);
    return Py_BuildValue("i", ret);
}





static PyTypeObject pycimage_PyCImageType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "pycimage.PyCImage",             /*tp_name*/
    sizeof(pycimage_PyCImageObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    0,                         /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    0,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,        /*tp_flags*/
    "Noddy objects",           /* tp_doc */
};

static PyMethodDef noddy_methods[] = {
    {NULL}  /* Sentinel */
};
