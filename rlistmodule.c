#include"rlistmodule.h"
#include "Python.h"

// static void printMultiindex (Multiindex mult) {
//   printf("atoms: ");
//   for (int i = 0; i<mult.na; i++){
//     printf("%hi ", mult.atoms[i]);
//   }
//   printf("species: ");
//   for (int i=0; i<mult.ns; i++) {
//     printf("%hi", mult.species[i]);
//   }
//   printf("\n");
// }

static void recursive_list(Multiindex current, MultiindexVector last_avail, MultiindexSet &previously_enumerated, MultiindexVector &ret, long int &count, long int &maxlevel, long int level){
  previously_enumerated.insert(current);
  count++;
  //Base case - check is all atoms have been exhausted, if so add current to ret
  int exhausted=1;
  for(int i =0; i < current.na; i++)
    if(current.atoms[i] != 0)
      exhausted=0;
  if(exhausted==1)
    ret.push_back(current);
  //recursive cases
  else{
    //Check from the previous set of available species which species are still available
    MultiindexVector avail;
    avail.reserve(last_avail.size());
    for (auto itr = last_avail.begin(); itr != last_avail.end(); ++itr) {
      int remove=0;
      for(int i =0; i < current.na; i++)
        if(current.atoms[i] - itr->atoms[i] < 0)
          remove=1;
      if(remove==0)
        avail.push_back(*itr);
    }
    //Loop through available species, create the next multiindex, and recurse if next has not been previously enumerated
    for (auto itr = avail.begin(); itr != avail.end(); ++itr) {
      Multiindex next(current.na, current.ns);
      for(int i=0; i<next.ns; i++)
        next.species[i]=current.species[i]+itr->species[i];
      for(int i=0; i<next.na; i++)
        next.atoms[i]=current.atoms[i]-itr->atoms[i];
      if(previously_enumerated.find(next)==previously_enumerated.end()){
        if(level+1>maxlevel)
          maxlevel=level+1;
        recursive_list(next,avail,previously_enumerated,ret,count,maxlevel,level+1);
      }
    }
  }
}

static PyObject *rlist_list(PyObject *self, PyObject* args){
    //Parse arguments and store in C arrays
    PyArrayObject *atoms, *spatoms, *exclude;
    int ndatoms, ndspatoms, ndexclude;
    npy_intp *datoms, *dspatoms, *dexclude;
    long int *atomsptr, *spatomsptr, *excludeptr;
    if (!PyArg_ParseTuple(args, "O!O!O!", &PyArray_Type, &atoms, &PyArray_Type, &spatoms,  &PyArray_Type, &exclude)) return NULL;
    ndatoms = PyArray_NDIM(atoms);
    ndspatoms = PyArray_NDIM(spatoms);
    ndexclude = PyArray_NDIM(exclude);
    datoms = PyArray_DIMS(atoms);
    dspatoms = PyArray_DIMS(spatoms);
    dexclude = PyArray_DIMS(exclude);
    if(!(ndatoms==1 && ndspatoms==2 && dspatoms[1] == datoms[0] && ndexclude == 1)){
      printf("Bad input; atoms must be one dimensional, exclude must be 1 dimensional, spatoms must be 2 dimensional, and the size of the second dimension of spatoms must equal the size of atoms.");
      return NULL;
    }
    atomsptr = (long int *)PyArray_DATA(atoms);
    spatomsptr = (long int *)PyArray_DATA(spatoms);
    excludeptr = (long int *)PyArray_DATA(exclude);
    long int ns=dspatoms[0]; //num species
    long int na=datoms[0]; //num atoms
    long int ne=dexclude[0]; //num atoms

    MultiindexSet previously_enumerated;
    MultiindexVector avail, list;
    //initialize the current multiindex with zeros and the number of atoms
    Multiindex current(na, ns);
    for (int i=0; i<na; i++)
      current.atoms[i]=atomsptr[i];
    for (int i=0; i<ns; i++)
      current.species[i]=0;

    // printMultiindex(current);
    //initialize avail with each available species' atoms and indices
    for (int i =0; i<ns; i++){
      Multiindex add(na, ns);
      for(int j =0; j<na; j++)
        add.atoms[j]=spatomsptr[na*i+j];
      for(int j=0; j<ns; j++)
        add.species[j]=0;
      add.species[i]=1;
      int remove=0;
      for(int j =0; j < na; j++)
        if(current.atoms[j] - add.atoms[j] < 0)
          remove=1;
      for(int j = 0; j < ne; j++){
        if(i == excludeptr[j])
          remove=1;
      }
      if(remove!=1){
        avail.push_back(add);
      }
    }

    long int count=0;
    long int maxlevel=0;
    recursive_list(current, avail, previously_enumerated, list, count, maxlevel, 0);

    //Populate the numpy output arrays with the information to return to Python
    // long int* data = new long int[list.size()*ns];
    long int* data = (long int*)malloc(list.size()*ns*sizeof(long int));
    int i=0;
    for (auto itr = list.begin(); itr != list.end(); ++itr) {
      for(int j=0; j<ns; j++){
        data[ns*i+j] = itr->species[j];
      }
      i++;
    }
    long int retdims[2];
    retdims[0]=list.size();
    retdims[1]=ns;
    PyObject *ret0 = PyArray_SimpleNewFromData (2,retdims,NPY_LONG,data);
    PyObject *ret=Py_BuildValue("(Oii)",ret0,count,maxlevel);
    // Py_INCREF(ret0);
    Py_INCREF(ret);
    return ret;
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
