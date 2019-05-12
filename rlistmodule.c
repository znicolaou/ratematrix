#define PY_SSIZE_T_CLEAN  /* Make "s#" use Py_ssize_t rather than int. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/arrayobject.h"
//#include <omp.h>

static PyObject *rlist_list(PyObject *self, PyObject* args){
    PyArrayObject *arr1=NULL;
    PyArrayObject *arr2=NULL;
    int nd1,nd2;
    npy_intp *dims1, *dims2;
    double* dptr1, *dptr2;

    if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &arr1, &PyArray_Type, &arr2)) return NULL;

    nd1 = PyArray_NDIM(arr1);
    dims1 = PyArray_DIMS(arr1);
    dptr1 = (double *)PyArray_DATA(arr1);
    nd2 = PyArray_NDIM(arr2);
    dims2 = PyArray_DIMS(arr2);
    dptr2 = (double *)PyArray_DATA(arr2);
    for(int i =0; i<dims1[0]; i++)
      dptr2[i]=dptr1[i]+1;

    Py_INCREF(arr2);
    return (PyObject *) arr2;

}

static PyMethodDef rlist_methods[] = {
    {"list",             rlist_list,      METH_VARARGS | METH_KEYWORDS,
     "Return the recursive list."},
    {NULL, NULL, 0, NULL}           /* sentinel */
};

static struct PyModuleDef rlistmodule = {
    PyModuleDef_HEAD_INIT,
    "rlist",
    NULL,
    -1,
    rlist_methods
};

PyMODINIT_FUNC PyInit_rlist(void){
    import_array();
    return PyModule_Create(&rlistmodule);
}

int main(int argc, char *argv[]){
    wchar_t *program = Py_DecodeLocale(argv[0], NULL);
    if (program == NULL) {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        exit(1);
    }
    PyImport_AppendInittab("rlist", PyInit_rlist);
    Py_SetProgramName(program);
    Py_Initialize();
    PyImport_ImportModule("rlist");
    PyMem_RawFree(program);
    return 0;
}
