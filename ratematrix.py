#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import cantera as ct
import timeit
import os
import re
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigs
import argparse

#Command-line arguments
parser = argparse.ArgumentParser(description='Generate a sparse rate matrix from cantera model.')
parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for npy file output. Three files will be created for each reaction, storing rates, row indices, and column indices.')
parser.add_argument("--accumulate", type=int, default=0, choices=[0,1], help='If 1, search filebase directory for npy files, generate sparse matrix, and plot eigenvalues. Default 0.')
parser.add_argument("--mechanism", type=str, required=False, default='mechanisms/h2o2.cti', dest='mechanism', help='Mechanism cti file. Default mechanisms/h2o2.cti.')
parser.add_argument("--reaction", type=int, required=False, default=None, dest='rind', help='Reaction index, provided for parallelization. If none is specified, the program will loop through all reactions in the model in sequence. Default None.')
parser.add_argument("--Nmax", type=int, required=False, default=5, dest='Nmax', help='Maximum number of molecules for each species. Default 5.')
parser.add_argument("--Nvals", type=int, required=False, default=1000, dest='Nvals', help='Number of eigenvalues to calculate, when --accumulate 1 is set. Default 1000')
parser.add_argument("--progress", type=int, required=False, default=1, choices=[0,1], help='Print progress during calculation. Default 1.')
args = parser.parse_args()

#Functions for relating multiindices to matrix indices
def get_digit(number, n, base):
    if number - base**n < 0:
        return 0
    return number // base**n % base
def get_multiindex(index):
    return np.array([get_digit(index, n, Nmax) for n in range(ns)])
def get_index(nums):
    return int(np.sum(nums*Nmax**np.arange(ns)))
#Function to get the transition rate given concentrations and rate constant
def get_rate (multiindex, rstoi, k, reaction):
    if reaction.reaction_type == 1: #bimolecular
        return k*np.product(multiindex**rstoi)
    else: #three body reactions
        ret=0
        for third_body in range(ns):
            rstoi2=rstoi
            rstoi2[third_body]+=1
            efficiency=reaction.default_efficiency
            if species[third_body] in reaction.efficiencies.keys():
                efficiency=reaction.efficiencies[species[third_body]]
            ret+=efficiency*k*np.product(multiindex**rstoi2)
        return ret
def calculate_sparse_elements(rind, filebase):
    #Main loop over rows to enumerate sparse data
    data=[]
    rows=[]
    columns=[]
    if(progress==1):
        print("Percent\t\tElapsed\t\tRemaining\t\t")
    stop=timeit.default_timer()
    reaction=gas.reactions()[rind]
    rstoi=np.array([reaction.reactants[x] if x in reaction.reactants.keys() else 0 for x in species])
    pstoi=np.array([reaction.products[x] if x in reaction.products.keys() else 0 for x in species])
    start2=timeit.default_timer()
    for  i in range(Nmax**ns):
        stop=timeit.default_timer()
        if progress==1:
            print("%.5f\t\t%.1fs\t\t%.1fs\t\t"%(i*1.0/Nmax**ns, stop-start2,(stop-start2)*(Nmax**ns-i-1)/(i+1)),end='\t\r')
        multiindex=get_multiindex(i)
        gas.X=multiindex
        #forward reaction
        k=gas.forward_rate_constants[rind]
        multiindex2=multiindex-rstoi+pstoi
        if np.all((multiindex-rstoi)[np.where(rstoi>0)]>=0) and np.all(multiindex2<Nmax) and k>0:
            j=get_index(multiindex2)
            data.append(get_rate(multiindex,rstoi,k,reaction))
            rows.append(i)
            columns.append(j)
            data.append(-get_rate(multiindex,rstoi,k,reaction))
            rows.append(i)
            columns.append(i)
        #reverse reaction
        k=gas.reverse_rate_constants[rind]
        multiindex2=multiindex+rstoi-pstoi
        if np.all((multiindex-pstoi)[np.where(pstoi>0)]>=0) and np.all(multiindex2<Nmax) and k>0:
            j=get_index(multiindex2)
            data.append(get_rate(multiindex,pstoi,k,reaction))
            rows.append(i)
            columns.append(j)
            data.append(-get_rate(multiindex,pstoi,k,reaction))
            rows.append(i)
            columns.append(i)
    #save to files
    np.save("%s_data"%(filebase),data)
    np.save("%s_rows"%(filebase),rows)
    np.save("%s_columns"%(filebase),columns)

#Main
filebase=args.filebase
mechanism=args.mechanism
Nmax=args.Nmax
Nvals=args.Nvals
rateindex=args.rind
progress=args.progress
accumulate=args.accumulate
gas=ct.Solution(mechanism)
ns=gas.n_species
nr=gas.n_reactions
species=gas.species_names
start=timeit.default_timer()

#Accumulate sparse data and plot
if(accumulate==1):
    p = re.compile("(.*)_data.npy")
    filebases=[item for sublist in np.unique([p.findall(name) for name in np.sort(os.listdir(filebase))]) for item in sublist]
    data=np.array([])
    columns=np.array([])
    rows=np.array([])
    for name in filebases:
        data=np.append(data,np.load(filebase+"/"+name+"_data.npy"))
        columns=np.append(columns,np.load(filebase+"/"+name+"_columns.npy"))
        rows=np.append(rows,np.load(filebase+"/"+name+"_rows.npy"))
    print("Sparsity: %f"%(1-len(data)*1.0/(Nmax**(2*ns))))
    print("Average non-zero entry: %f"%(np.linalg.norm(data)/len(data)))
    ratematrix=coo_matrix((np.array(data),(np.array(rows),np.array(columns))),(Nmax**ns,Nmax**ns))
    eigenvalues,eigenvectors=eigs(ratematrix,k=Nvals)
    np.save(filebase+"/eigenvalues",eigenvalues)
    np.save(filebase+"/eigenvectors",eigenvectors)
    fig=plt.figure()
    fig.gca().set_ylabel(r'$\mathrm{Im}\left(\lambda\right)$')
    fig.gca().set_xlabel(r'$\mathrm{Re}\left(\lambda\right)$')
    plt.plot(np.real(eigenvalues),np.imag(eigenvalues), 'bx', markersize=3.0)
    plt.tight_layout()
    plt.show()
    fig.savefig(filebase+"/eigenvalues.pdf", bbox_inches='tight')
    print("\nRuntime: %.1f s"%(timeit.default_timer()-start))
elif rateindex == None:
    for rind in range(nr):
        print('Reaction %i'%(rind))
        calculate_sparse_elements(rind,filebase+"/%i"%(rind))
        if(progress==1):
            print('')
    print("\nRuntime: %.1fs"%(timeit.default_timer()-start))
else:
    print('Reaction %i'%(rateindex))
    calculate_sparse_elements(rateindex, filebase)
    print("\nRuntime: %.1fs"%(timeit.default_timer()-start))
