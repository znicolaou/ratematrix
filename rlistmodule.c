#include"rlistmodule.h"
#include "Python.h"

// static void printMultiindex (Multiindex mult) {
//   printf("atoms: ");
//   for (int i = 0; i<mult.na; i++){
//     printf("%li ", mult.atoms[i]);
//   }
//   printf("species: ");
//   for (int i=0; i<mult.ns; i++) {
//     printf("%li", mult.species[i]);
//   }
//   printf("\n");
// }

static void recursive_list(Multiindex *current, MultiindexVector &last_avail, MultiindexSet &previously_enumerated, MultiindexVector &ret, long int &count, long int &maxlevel, long int level){
  previously_enumerated.insert(*current);
  count++;
  //Base case - check is all atoms have been exhausted, if so add current to ret
  int exhausted=1;
  for(int i =0; i < current->na; i++)
    if(current->atoms[i] != 0)
      exhausted=0;
  if(exhausted==1)
    ret.push_back(*current);
  //recursive cases
  else{
    //Check from the previous set of available species which species are still available
    MultiindexVector avail;
    avail.reserve(last_avail.size());
    for (auto itr = last_avail.begin(); itr != last_avail.end(); ++itr) {
      int remove=0;
      for(int i =0; i < current->na; i++)
        if(current->atoms[i] - itr->atoms[i] < 0)
          remove=1;
      if(remove==0)
        avail.push_back(*itr);
    }
    //Loop through available species, create the next multiindex, and recurse if next has not been previously enumerated
    for (auto itr = avail.begin(); itr != avail.end(); ++itr) {
      Multiindex *next = new Multiindex();
      next->na=current->na;
      next->ns=current->ns;
      next->atoms=new long int[current->na];
      next->species=new long int[current->ns];
      for(int i=0; i<next->ns; i++)
        next->species[i]=current->species[i]+itr->species[i];
      for(int i=0; i<next->na; i++)
        next->atoms[i]=current->atoms[i]-itr->atoms[i];
      if(previously_enumerated.find(*next)==previously_enumerated.end()){
        if(level+1>maxlevel)
          maxlevel=level+1;
        recursive_list(next,avail,previously_enumerated,ret,count,maxlevel,level+1);
      }
    }
  }
}



static PyObject *rlist_list(PyObject *self, PyObject* args){
    //Parse arguments and store in C arrays
    PyArrayObject *atoms, *spatoms;
    int ndatoms, ndspatoms;
    npy_intp *datoms, *dspatoms;
    long int *atomsptr, *spatomsptr;
    if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &atoms, &PyArray_Type, &spatoms)) return NULL;
    ndatoms = PyArray_NDIM(atoms);
    ndspatoms = PyArray_NDIM(spatoms);
    datoms = PyArray_DIMS(atoms);
    dspatoms = PyArray_DIMS(spatoms);
    if(!(ndatoms==1 && ndspatoms==2 && dspatoms[1] == datoms[0])){
      printf("Bad input; atoms must be one dimensional, spatoms must be 2 dimensional, and the size of the second dimension of spatoms must equal the size of atoms.");
      return NULL;
    }
    atomsptr = (long int *)PyArray_DATA(atoms);
    spatomsptr = (long int *)PyArray_DATA(spatoms);
    long int ns=dspatoms[0]; //num species
    long int na=datoms[0]; //num atoms


    //Find the species that can be added the most times
    //We will parallelize over the number of this species that we include
    long int *maxes = new long int[ns];
    for (int i=0; i<ns; i++){
      long int max=0;
      long int stop=0;
      while(stop!=1){
        max++;
        for(int j=0; j<na; j++){
          if(atomsptr[j] - max*spatomsptr[na*i+j] < 0){
            stop=1;
          }
        }
      }
      maxes[i]=max;
    }
    int maxind=0;
    int maxmax=0;
    for(int i = 0; i<ns; i++){
      if(maxes[i]>maxmax){
        maxmax=maxes[i];
        maxind=i;
      }
    }

    MultiindexVector *lists = new MultiindexVector[maxmax];
    long int *counts=new long int[maxmax];
    long int *maxlevels=new long int[maxmax];
    //Run the recursive listing in parallel
    int k=0;
    int completed=0;

    #pragma omp parallel for shared(completed)
    for(k=0; k<=maxmax; k++)
    {
      printf("%.3f  \r",completed*1.0/maxmax);
      fflush(stdout);
      Multiindex current;
      // Multiindex *current = new Multiindex();
      MultiindexSet previously_enumerated;
      MultiindexVector avail;
      //initialize the current multiindex with zeros and the number of atoms
      // current->na=na;
      // current->ns=ns;
      // current->atoms=new long int[na];
      // current->species=new long int[ns];
      // for (int i=0; i<na; i++)
      //   current->atoms[i]=atomsptr[i]-k*spatomsptr[na*maxind+i];
      // for (int i=0; i<ns; i++)
      //   current->species[i]=0;
      // current->species[maxind]=k;
      current.na=na;
      current.ns=ns;
      current.atoms=new long int[na];
      current.species=new long int[ns];
      for (int i=0; i<na; i++)
        current.atoms[i]=atomsptr[i]-k*spatomsptr[na*maxind+i];
      for (int i=0; i<ns; i++)
        current.species[i]=0;
      current.species[maxind]=k;

      //initialize avail with each available species' atoms and indices
      for (int i =0; i<ns; i++){
        //Don't add maxind species to avail - it is fixed.
        if(i!=maxind){

          Multiindex add;
          add.na=na;
          add.ns=ns;
          add.atoms=new long int[na];
          add.species=new long int[ns];
          for(int j =0; j<na; j++){
            add.atoms[j]=spatomsptr[na*i+j];
          }
          for(int j=0; j<ns; j++){
            add.species[j]=0;
          }
          add.species[i]=1;
          int remove=0;
          for(int i =0; i < current.na; i++)
            if(current.atoms[i] - add.atoms[i] >= 0)
              remove=1;
          if(remove==1){
            avail.push_back(add);
          }
        }
      }
      //Run the recursive listing
      counts[k]=0;
      maxlevels[k]=0;
      recursive_list(&current, avail, previously_enumerated, lists[k], counts[k], maxlevels[k], 0);
      completed++;
    }

    int count=0;
    int maxlevel=0;
    MultiindexVector list;

    for(int k=0; k<maxmax; k++){
      // printf("%i %i %i %i\n",k, counts[k], maxlevels[k], lists[k].size());
      count+=counts[k];
      if(maxlevels[k]>maxlevel){
        maxlevel=k;
      }
      for(unsigned int l=0; l<lists[k].size(); l++){
        list.push_back(lists[k][l]);
        // printMultiindex(lists[k][l]);
      }
    }
    //We should technically free all the Multiindices allocated and stored in previously_enumerated that were not returned
    //Memory is not the biggest constraint though at the moment


    //Populate the numpy output arrays with the information to return to Python
    long int* data = new long int[list.size()*ns];
    int i=0;
    for (auto itr = list.begin(); itr != list.end(); ++itr) {
      for(int j=0; j<ns; j++){
        data[ns*i+j] = itr->species[j];
      }
      i++;
    }

    long int *retdims = new long int[2];
    retdims[0]=list.size();
    retdims[1]=ns;
    PyObject *ret0 = PyArray_SimpleNewFromData (2,retdims,NPY_LONG,data);
    PyObject *ret=Py_BuildValue("(Oii)",ret0,count,maxlevel);
    Py_INCREF(ret0);
    Py_INCREF(ret);

    for (auto itr = list.begin(); itr != list.end(); ++itr) {
      delete[] itr->species;
      delete[] itr->atoms;
    }
    delete[] lists;
    delete[] counts;
    delete[] maxlevels;
    delete[] retdims;


    return ret;
}

static PyObject *rlist_list2(PyObject *self, PyObject* args){
    //Parse arguments and store in C arrays
    PyArrayObject *atoms, *spatoms;
    int ndatoms, ndspatoms;
    npy_intp *datoms, *dspatoms;
    long int *atomsptr, *spatomsptr;
    if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &atoms, &PyArray_Type, &spatoms)) return NULL;
    ndatoms = PyArray_NDIM(atoms);
    ndspatoms = PyArray_NDIM(spatoms);
    datoms = PyArray_DIMS(atoms);
    dspatoms = PyArray_DIMS(spatoms);
    if(!(ndatoms==1 && ndspatoms==2 && dspatoms[1] == datoms[0])){
      printf("Bad input; atoms must be one dimensional, spatoms must be 2 dimensional, and the size of the second dimension of spatoms must equal the size of atoms.");
      return NULL;
    }
    atomsptr = (long int *)PyArray_DATA(atoms);
    spatomsptr = (long int *)PyArray_DATA(spatoms);
    long int ns=dspatoms[0]; //num species
    long int na=datoms[0]; //num atoms

    Multiindex *current = new Multiindex();
    MultiindexSet previously_enumerated;
    MultiindexVector avail, list;
    //initialize the current multiindex with zeros and the number of atoms
    current->na=na;
    current->ns=ns;
    current->atoms=new long int[na];
    current->species=new long int[ns];
    for (int i=0; i<na; i++)
      current->atoms[i]=atomsptr[i];
    for (int i=0; i<ns; i++)
      current->species[i]=0;

    //initialize avail with each available species' atoms and indices
    for (int i =0; i<ns; i++){
      Multiindex *add = new Multiindex();
      add->na=na;
      add->ns=ns;
      add->atoms=new long int[na];
      add->species=new long int[ns];
      for(int j =0; j<na; j++){
        add->atoms[j]=spatomsptr[na*i+j];
      }
      for(int j=0; j<ns; j++){
        add->species[j]=0;
      }
      add->species[i]=1;
      int remove=0;
      for(int i =0; i < current->na; i++)
        if(current->atoms[i] - add->atoms[i] >= 0)
          remove=1;
      if(remove==1){
        avail.push_back(*add);
      }
    }
    //Run the recursive listing
    long int count=0;
    long int maxlevel=0;
    recursive_list(current, avail, previously_enumerated, list, count, maxlevel, 0);
    //We should technically free all the Multiindices allocated and stored in previously_enumerated that were not returned
    //Memory is not the biggest constraint though at the moment

    //Populate the numpy output arrays with the information to return to Python
    long int* data = new long int[list.size()*ns];
    int i=0;
    for (auto itr = list.begin(); itr != list.end(); ++itr) {
      for(int j=0; j<ns; j++){
        data[ns*i+j] = itr->species[j];
      }
      i++;
    }
    long int *retdims = new long int[2];
    retdims[0]=list.size();
    retdims[1]=ns;
    PyObject *ret0 = PyArray_SimpleNewFromData (2,retdims,NPY_LONG,data);
    PyObject *ret=Py_BuildValue("(Oii)",ret0,count,maxlevel);
    Py_INCREF(ret0);
    Py_INCREF(ret);
    return ret;
}


static PyMethodDef rlist_methods[] = {{"list", (PyCFunction)rlist_list, METH_VARARGS, "Return the recursive list."}, {"list2", (PyCFunction)rlist_list2, METH_VARARGS, "Return the recursive list."}, {NULL, NULL}};

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
