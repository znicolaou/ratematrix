#include"rlistmodule.h"
#include "Python.h"

static void printMultiindex (Multiindex mult) {
  printf("atoms: ");
  for (int i = 0; i<mult.na; i++){
    printf("%hi ", mult.atoms[i]);
  }
  printf("species: ");
  for (int i=0; i<mult.ns; i++) {
    printf("%hi ", mult.species[i]);
  }
  printf("\n");
}

// static void recursive_list(Multiindex current, MultiindexVector reactions, MultiindexSet &previously_enumerated){
static void recursive_list(Multiindex current, long int *reactionsptr, int ns, int nr, MultiindexSet &previously_enumerated){
  MultiindexVector stack;
  stack.push_back(current);
  previously_enumerated.insert(current);
  int remove;

  while(!stack.empty()) {
    current=stack.back();
    // printf("%i\t%i\t\r",stack.size(),previously_enumerated.size());
    stack.pop_back();
    //since reactions has fixed length, we could do a pragma omp parallel here
    // for (auto itr = reactions.begin(); itr != reactions.end(); ++itr) {
    for (int j =0; j < nr; j++) {
      Multiindex next(current.na, current.ns);
      remove=0;
      for(int i=0; i<next.ns; i++){
        // next.species[i]=current.species[i]+itr->species[i];
        next.species[i]=current.species[i]+reactionsptr[ns*i+j];
        if(next.species[i] < 0){
          remove=1;
          break;
        }
      }
      // for(int i=0; i<next.na; i++)
      //   next.atoms[i]=current.atoms[i];

      if(remove==0 && previously_enumerated.find(next)==previously_enumerated.end()) {
        // printMultiindex(next);
        stack.push_back(next);
        previously_enumerated.insert(next);
      }
    }
  }
}


static PyObject *rlist_list(PyObject *self, PyObject* args){
    //Parse arguments and store in C arrays
    PyArrayObject *species, *reactions;
    int ndspecies, ndreactions;
    npy_intp *dspecies, *dreactions;
    long int *speciesptr, *reactionsptr;
    if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &species, &PyArray_Type, &reactions)) return NULL;
    ndspecies = PyArray_NDIM(species);
    ndreactions = PyArray_NDIM(reactions);
    dspecies = PyArray_DIMS(species);
    dreactions = PyArray_DIMS(reactions);
    if(!(ndspecies==1)){
      printf("Bad input; atoms must be one dimensional, exclude must be 1 dimensional, spatoms must be 2 dimensional, and the size of the second dimension of spatoms must equal the size of atoms.");
      return NULL;
    }
    speciesptr = (long int *)PyArray_DATA(species);
    reactionsptr = (long int *)PyArray_DATA(reactions);
    long int ns=dreactions[1]; //num species
    long int na=0; //num atoms
    long int nr=dreactions[0]; //num atoms
    MultiindexSet previously_enumerated;
    MultiindexVector allreac;
    //initialize the current multiindex with zeros and the number of atoms
    Multiindex current(na, ns);
    for (int i=0; i<na; i++)
      current.atoms[i]=0;
    for (int i=0; i<ns; i++){
      current.species[i]=speciesptr[i];
      for (int j=0; j<na; j++)
        current.atoms[j]+=0;
    }
    //populate reactions and avail with available reactions
    for (int i =0; i<nr; i++){
      Multiindex add(na, ns);
      for(int j =0; j<na; j++)
        add.atoms[j]=0;
      for(int j=0; j<ns; j++){
        // add.species[j]=reactionsptr[ns*i+j];
        add.species[j]=reactionsptr[ns*i+j];
      }
      allreac.push_back(add);
      // printMultiindex(add);
    }

    // recursive_list(current, allreac, previously_enumerated);
    recursive_list(current, reactionptr, ns, nr, previously_enumerated);

    //Populate the numpy output arrays with the information to return to Python
    int dim=previously_enumerated.size();
    long int* data = (long int*)malloc(dim*ns*sizeof(long int));
    int i=0;
    for (auto itr = previously_enumerated.begin(); itr != previously_enumerated.end(); ++itr) {
      for(int j=0; j<ns; j++){
        data[ns*i+j] = itr->species[j];
      }
      i++;
    }
    long int retdims[2];
    retdims[0]=dim;
    retdims[1]=ns;
    PyObject *ret0 = PyArray_SimpleNewFromData (2,retdims,NPY_LONG,data);
    Py_INCREF(ret0);
    return ret0;
}

static PyMethodDef rlist_methods[] = {{"list", (PyCFunction)rlist_list, METH_VARARGS, "Return the recursive list."}, {NULL, NULL}};

static struct PyModuleDef rlistmodule = {PyModuleDef_HEAD_INIT, "rlist", NULL, -1, rlist_methods};

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
