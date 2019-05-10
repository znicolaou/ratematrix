#include "Python.h"
#include <omp.h>

void rlist(void);

int main(int argc, char **argv)
{
    Py_SetProgramName(argv[0]);
    Py_Initialize();
    rlist();
    printf("Hello, brave new world\n\n");
    printf("\nGoodbye, cruel world\n");
    Py_Exit(0);
}

static PyObject *
rlist_list(PyObject *self, PyObject* args)
{
    static char *kwlist[] = {"atoms", "species", "remaining_atoms", NULL};
    int* atoms;
    int* species;
    int* remaining_atoms;

    return PyInt_FromLong(42L);
}

static PyMethodDef rlist_methods[] = {
    {"list",             rlist_list,      METH_VARARGS | METH_KEYWORDS,
     "Return the recursive list."},
    {NULL, NULL, 0, NULL}           /* sentinel */
};

void rlist(void)
{
    PyImport_AddModule("rlist");
    Py_InitModule("rlist", rlist_methods);
}
