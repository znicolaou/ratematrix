#define PY_SSIZE_T_CLEAN  /* Make "s#" use Py_ssize_t rather than int. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/arrayobject.h"
#include <omp.h>
using namespace std;

static PyObject *tlist_list(PyObject *self, PyObject* args){
  PyArrayObject *array;
  int nd;
  npy_intp *dims;
  long int *aptr;
  if (!PyArg_ParseTuple(args, "O!", &PyArray_Type, &array)) return NULL;
  nd = PyArray_NDIM(array);
  dims = PyArray_DIMS(array);
  aptr = (long int *)PyArray_DATA(array);
  long int n=dims[0]; //num species
  for(int i=0; i<n; i++){
    printf("%li ", aptr[i]);
  }
  printf("\n");

  int maxmax=3;
  double *maxes=new double[maxmax+1];
  #pragma omp parallel for
  for(int i=0; i<=maxmax; i++)
  {
    printf("%i\n",i);
    double *temp = &(maxes[i]);
    *(temp) = 0;
    fflush(stdout);
    for(long int j=0; j<1e10; j++){
      *(temp) += 1;
    }
  }

  for(int i =0; i<=maxmax; i++){
    printf("%i %f\n", i, maxes[i]);
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static PyMethodDef tlist_methods[] = {{"list", (PyCFunction)tlist_list, METH_VARARGS, "Return the recursive list."}, {NULL, NULL}};

static struct PyModuleDef tlistmodule = {PyModuleDef_HEAD_INIT, "tlist", NULL, -1, tlist_methods};

PyMODINIT_FUNC PyInit_tlist(void){
    import_array();
    return PyModule_Create(&tlistmodule);
}

int main(int argc, char *argv[]){
    wchar_t *program = Py_DecodeLocale(argv[0], NULL);
    if (program == NULL) {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        exit(1);
    }
    PyImport_AppendInittab("tlist", PyInit_tlist);
    Py_SetProgramName(program);
    Py_Initialize();
    PyImport_ImportModule("tlist");
    PyMem_RawFree(program);
    return 0;
}
