#define PY_SSIZE_T_CLEAN  /* Make "s#" use Py_ssize_t rather than int. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/arrayobject.h"
//#include <omp.h>


//Translate this to C code, and call it from rlist_list with data from the parsed numpy arrays
//We may want to use stacks, lists, queues, or some other C data structures
// def recursive_list(remaining_atoms, multiindex, last_avail, previously_enumerated=[],level=0):
//     #Add current multiindex to previously enumerated list so it is not repeated
//     previously_enumerated.append(multiindex.copy())
//
//     #Find available species to add out of last available set
//     avail=[[],[]]
//     for i in range(len(last_avail[0])):
//         if (np.all(remaining_atoms-last_avail[1][i]>=0)):
//             avail[0].append(last_avail[0][i])
//             avail[1].append(last_avail[1][i])
//
//     #Recurse for each new multiindex that has not been previously enumerated and return list of multiindices
//     if(avail!=[[],[]]):
//         ret_lists=[]
//         ret_counts=0
//         max_level=0
//         for i in range(len(avail[0])):
//             multiindex[avail[0][i]]+=1
//             if not (multiindex in previously_enumerated):
//                 returned_list,returned_count,returned_level=recursive_list(remaining_atoms-avail[1][i], multiindex, avail, previously_enumerated, level+1)
//                 ret_lists+=returned_list
//                 ret_counts+=returned_count
//                 if returned_level>max_level:
//                     max_level=returned_level
//             multiindex[avail[0][i]]-=1
//         return ret_lists,1+ret_counts,max_level
//
//     #Base case if no new species can be added
//     else:
//         #Return list with current multiindex if all atoms exausted
//         if np.all(remaining_atoms == np.zeros(len(elements))):
//             return [multiindex.copy()],1,level
//         #Could not exaust atoms with this branch; return nothing
//         else:
//             return [],1,level

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
