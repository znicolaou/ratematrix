#define PY_SSIZE_T_CLEAN  /* Make "s#" use Py_ssize_t rather than int. */
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include "Python.h"
#include "numpy/arrayobject.h"
#include <unordered_set>
#include <vector>
#include <omp.h>
using namespace std;

typedef struct Multiindex{
	short int na, ns;
	short int *species, *atoms;
	Multiindex(long int Na, long int Ns){na=Na; ns=Ns; species=(short int*)malloc(ns*sizeof(short int)); atoms=(short int*)malloc(na*sizeof(short int));};
} Multiindex;

struct MultiindexHash {
	public:
		size_t operator()(const Multiindex & ind) const {
    int ret=0;
    for(int i = 0; i<(ind.ns); i++)
      ret += pow(ind.ns, i)*ind.species[i];
		return std::hash<int>()(ret);
	}
};

struct EqualMultiindices {
	public:
		bool operator()(const Multiindex & ind1, const Multiindex & ind2) const {
    bool ret=true;
    for(int i = 0; i < (ind1.ns); i++)
      if(ind1.species[i] != ind2.species[i])
        ret=false;
		return ret;
	}
};

typedef unordered_set<Multiindex, MultiindexHash, EqualMultiindices> MultiindexSet;
typedef vector<Multiindex> MultiindexVector;
