#!/usr/bin/env python
import matplotlib
matplotlib.use("TKAgg")
import matplotlib.pyplot as plt
import numpy as np
import cantera as ct
import timeit
import argparse
from scipy.sparse import coo_matrix
from scipy.linalg import eig
from scipy.special import factorial
import sys
import rlist

#Command-line arguments
parser = argparse.ArgumentParser(description='Generate a sparse rate matrix from cantera model.')
parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for npy file output. Three files will be created for each reaction, storing rates, row indices, and column indices.')
parser.add_argument("--mechanism", type=str, required=False, default='mechanisms/h2o2.cti', dest='mechanism', help='Mechanism cti file. Default mechanisms/h2o2.cti.')
parser.add_argument("--calculate", type=int, required=False, default=1, choices=[0,1], help='Flag to calculate rate matrix and eigenvalues. If 1, calculate space, then print args.temperature, args.pressure, total atoms, dimension, runtime, recursive calls, recursive levels, sparsity, and Nvals largest eigenvalues, save rate matrix, eigenvalues, and eigenvectors, then quit. If 0, only calculate space, then print total atoms, dimension, runtime, recursive calls, and recursive levels, then quit. Default 1.')
parser.add_argument("--Nvals", type=int, required=False, default=-1, help='Number of eigenvalues to print; 0 for all. Default 0.')
parser.add_argument("--plot", type=int, required=False, default=1, choices=[0,1], help='Flag to plot eigenvalues and save filebase/eigenvales.pdf. Default 1.')
parser.add_argument("--save", type=int, required=False, default=1, choices=[0,1], help='Flag to save the results to filebaserate.npy, filebaseeigenvalues.npy, and filebaseeigenvectors.npy. Default 1.')
parser.add_argument("--temperature", type=float, required=False, default=1500, help='Temperature in Kelvin. Default 1500.')
parser.add_argument("--pressure", type=float, required=False, default=1, help='Pressure in atm. Default 1.')
parser.add_argument("--atoms", nargs='+', type=int, required=False, default=[3, 3, 3], help='Number of each atom, in order of their appearance in the .cti file. If number of values is not number of atoms, print the atoms. Default 3 3 3.')
parser.add_argument("--fix", nargs='+', type=int, required=False, default=[], help='Fix species numbers. Include each species index followed by the number of molecules to fix.')
args = parser.parse_args()

#Functions for relating multiindices to matrix indices
def get_multiindex(index):
    return multiindices[index]
def get_index(multiindex):
    return np.where(np.all(multiindices==multiindex,axis=1))[0][0]

#Function to get the transition rate given concentrations and rate constant
def get_rate (multiindex, rstoi, k, reaction):
    if reaction.reaction_type == 1: #bimolecular
        # return k*np.product(multiindex**rstoi)
        if np.all(multiindex>=rstoi):
            # return k*np.product(multiindex**rstoi)
            return k*np.product(factorial(multiindex)/factorial(multiindex-rstoi))
        else:
            return 0.
    else: #three body reactions
        ret=0
        for third_body in range(ns):
            rstoi[third_body]+=1
            if np.all(multiindex>=rstoi):
                efficiency=reaction.default_efficiency
                if species[third_body] in reaction.efficiencies.keys():
                    efficiency=reaction.efficiencies[species[third_body]]
                # ret+=efficiency*k*np.product(multiindex**rstoi)
                ret+=efficiency*k*np.product(factorial(multiindex)/factorial(multiindex-rstoi))
            rstoi[third_body]-=1
        return ret
#Main loop over rows to enumerate sparse data
def calculate_sparse_elements(rind):
    data=[]
    rows=[]
    columns=[]
    reaction=gas.reactions()[rind]
    rstoi=np.array([reaction.reactants[x] if x in reaction.reactants.keys() else 0 for x in species])
    pstoi=np.array([reaction.products[x] if x in reaction.products.keys() else 0 for x in species])
    start2=timeit.default_timer()
    for  i in range(len(multiindices)):
        stop=timeit.default_timer()
        multiindex=get_multiindex(i)
        gas.X=multiindex
        #forward reaction
        k=gas.forward_rate_constants[rind]
        multiindex2=multiindex-rstoi+pstoi
        if np.all(multiindex2>=0.) and not np.isnan(k):
            rate=get_rate(multiindex,rstoi,k,reaction)
            j=get_index(multiindex2)
            data.append(rate)
            rows.append(i)
            columns.append(j)
            data.append(-rate)
            rows.append(i)
            columns.append(i)
        #reverse reaction
        k=gas.reverse_rate_constants[rind]
        multiindex2=multiindex+rstoi-pstoi
        if np.all(multiindex2>=0) and not np.isnan(k):
            rate=get_rate(multiindex,pstoi,k,reaction)
            j=get_index(multiindex2)
            data.append(rate)
            rows.append(i)
            columns.append(j)
            data.append(-rate)
            rows.append(i)
            columns.append(i)
    return data,rows,columns

#Main
filebase=args.filebase
mechanism=args.mechanism
gas=ct.Solution(mechanism)
gas.TP=args.temperature,args.pressure*ct.one_atm
ns=gas.n_species
nr=gas.n_reactions
species=gas.species_names
elements=gas.element_names
atoms=np.array(args.atoms)
fixed=np.array(args.fix)
if(len(atoms)!=len(elements)):
    print("elements are: ", *elements)
    quit()
start=timeit.default_timer()

#Calculate the space of possible states
atoms=np.array(atoms)
multiindex=[]
last_avail=[[],[]]
for i in range(ns):
    multiindex.append(0)
    last_avail[0].append(i)
    last_avail[1].append(np.array([int(gas.species()[i].composition[el] if el in gas.species()[i].composition.keys() else 0) for el in elements]))

sp_atoms=[]
for i in range(ns):
    multiindex.append(0)
    sp_atoms.append(np.array([int(gas.species()[i].composition[el] if el in gas.species()[i].composition.keys() else 0) for el in elements]))
sp_atoms=np.array(sp_atoms)

remove_atoms = np.zeros(len(atoms));
for i in range(0,len(fixed),2):
    remove_atoms += fixed[i+1]*sp_atoms[fixed[i]]

multiindices,count,level=rlist.list(atoms-remove_atoms.astype(int), sp_atoms, fixed[::2].astype(int))

dim=len(multiindices)
if args.calculate==0:
    #Print total atoms, dimension, runtime, recursive calls, and recursive levels
    print(np.sum(atoms), dim, timeit.default_timer()-start, count, level)
    out=open(filebase+"out.dat","a+")
    print(np.sum(atoms), dim, timeit.default_timer()-start, count, level,file=out)
    out.close()
    sys.stdout.flush()
    quit()

#Loop through each reaction index and calculate spase elements
data=[]
rows=[]
columns=[]
for rind in range(nr):
    reac_data,reac_rows,reac_columns=calculate_sparse_elements(rind)
    data+=reac_data
    rows+=reac_rows
    columns+=reac_columns
nonzero=np.array([rows,columns])
ratematrix=coo_matrix((np.array(data),(np.array(rows),np.array(columns))),(dim,dim))
#The dimension is not that big - don't use spase matrix algorithms for eigenvalues
eigenvalues,eigenvectors=eig(np.transpose(ratematrix.toarray()))
sorted=np.argsort(eigenvalues)
if args.save==1:
    np.save(filebase+"ratematrix.npy",ratematrix.toarray())
    np.save(filebase+"eigenvalues.npy",eigenvalues.astype(complex)[sorted])
    np.save(filebase+"eigenvectors.npy",eigenvectors.astype(complex)[:,sorted])
    np.save(filebase+"multiindices.npy",multiindices)
    np.save(filebase+"spatoms.npy",sp_atoms)

#Print dimension, runtime, sparsity, and three smallest eigenvalues
print(args.temperature, args.pressure, np.sum(atoms), dim, timeit.default_timer()-start, count, level, (1-len(np.unique(nonzero))*1.0/(dim)**2), *np.sort(np.real(eigenvalues[-args.Nvals:])))
out=open(filebase+"out.dat","w")
print(args.temperature, args.pressure, np.sum(atoms), dim, timeit.default_timer()-start, count, level, (1-len(np.unique(nonzero))*1.0/(dim)**2), *np.sort(np.real(eigenvalues[-args.Nvals:])),file=out)
print(*elements, file=out)
out.close()
sys.stdout.flush()

if args.plot==1:
    fig=plt.figure()
    fig.gca().set_ylabel(r'$\mathrm{Im}\left(\lambda\right)$')
    fig.gca().set_xlabel(r'$\mathrm{Re}\left(\lambda\right)$')
    plt.plot(np.real(eigenvalues),np.imag(eigenvalues), 'bx', markersize=3.0)
    plt.tight_layout()
    plt.show()
    fig.savefig(filebase+"eigenvalues.pdf", bbox_inches='tight')
