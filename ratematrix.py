#!/usr/bin/env python
from __future__ import print_function
import matplotlib
matplotlib.use("TKAgg")
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
parser.add_argument("--mechanism", type=str, required=False, default='mechanisms/h2o2.cti', dest='mechanism', help='Mechanism cti file. Default mechanisms/h2o2.cti.')
parser.add_argument("--Nvals", type=int, required=False, default=10, dest='Nvals', help='Number of eigenvalues to calculate, when --accumulate 1 is set. Default 1000')
parser.add_argument("--calculate", type=int, required=False, default=1, choices=[0,1], help='Flag to calculate rate matrix and eigenvalues. If 0, only calculate dimension and quit. Default 1.')
parser.add_argument("--temperature", type=float, required=False, default=1500, help='Temperature in Kelvin. Default 1500.')
parser.add_argument("--pressure", type=float, required=False, default=1, help='Pressure in atm. Default 1')
parser.add_argument("--atoms", nargs='+', type=int, required=False, default=[3, 3, 3], help='Number of each atom, in order of their appearance in the .cti file.')
args = parser.parse_args()

#Functions for relating multiindices to matrix indices
def get_multiindex(index):
    return multiindices[index]
def get_index(multiindex):
    return np.where(multiindices==multiindex)[0][0]
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
#Main loop over rows to enumerate sparse data
def calculate_sparse_elements(rind, filebase):
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

#Find all list of species that have specified number of atoms
def recursive_list(remaining_atoms, multiindex, last_avail, previously_enumerated=[],level=0):
    #Add current multiindex to previously enumerated list so it is not repeated
    previously_enumerated.append(multiindex.copy())

    #Find available species to add out of last available set
    avail=[[],[]]
    for i in range(len(last_avail[0])):
        if (np.all(remaining_atoms-last_avail[1][i]>=0)):
            avail[0].append(last_avail[0][i])
            avail[1].append(last_avail[1][i])

    #Recurse for each new multiindex that has not been previously enumerated and return list of multiindices
    if(avail!=[[],[]]):
        ret_lists=[]
        ret_counts=0
        max_level=0
        for i in range(len(avail[0])):
            multiindex[avail[0][i]]+=1
            if not (multiindex in previously_enumerated):
                returned_list,returned_count,returned_level=recursive_list(remaining_atoms-avail[1][i], multiindex, avail, previously_enumerated, level+1)
                ret_lists+=returned_list
                ret_counts+=returned_count
                if returned_level>max_level:
                    max_level=returned_level
            multiindex[avail[0][i]]-=1
        return ret_lists,1+ret_counts,max_level

    #Base case if no new species can be added
    else:
        #Return list with current multiindex if all atoms exausted
        if np.all(remaining_atoms == np.zeros(len(elements))):
            return [multiindex.copy()],1,level
        #Could not exaust atoms with this branch; return nothing
        else:
            return [],1,level
#Main
filebase=args.filebase
if not os.path.isdir(filebase):
    os.mkdir(filebase)
mechanism=args.mechanism
Nvals=args.Nvals
gas=ct.Solution(mechanism)
gas.TP=args.temperature,args.pressure*ct.one_atm
ns=gas.n_species
nr=gas.n_reactions
species=gas.species_names
elements=gas.element_names
atoms=np.array(args.atoms)
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
    last_avail[1].append(np.array([gas.species()[i].composition[el] if el in gas.species()[i].composition.keys() else 0 for el in elements]))

multiindices,count,level=recursive_list(atoms,multiindex,last_avail)
dim=len(multiindices)
print(np.sum(atoms), dim, timeit.default_timer()-start, count, level)
if args.calculate==0:
    quit()

#Loop through each reaction index and calculate spase elements
data=[]
rows=[]
columns=[]
for rind in range(nr):
    reac_data,reac_rows,reac_columns=calculate_sparse_elements(rind,filebase+"/%i"%(rind))
    data+=reac_data
    rows+=reac_rows
    columns+=reac_columns
nonzero=np.array([rows,columns])
print("Sparsity: %f"%(1-len(np.unique(nonzero))*1.0/(dim)**2))
print("Average non-zero entry: %f"%(np.linalg.norm(data)/len(data)))
ratematrix=coo_matrix((np.array(data),(np.array(rows),np.array(columns))),(dim,dim))
np.save(filebase+"/rate.npy",ratematrix)
eigenvalues,eigenvectors=eigs(ratematrix,k=Nvals,which='LR')
np.save(filebase+"/eigenvalues",eigenvalues)
np.save(filebase+"/eigenvectors",eigenvectors)
rmax=np.abs(np.max(np.real(eigenvalues)))
imax=np.max(np.abs(np.imag(eigenvalues[np.where(np.real(eigenvalues)>-2*rmax)[0]])))
fig=plt.figure()
fig.gca().set_ylabel(r'$\mathrm{Im}\left(\lambda\right)$')
fig.gca().set_xlabel(r'$\mathrm{Re}\left(\lambda\right)$')
fig.gca().set_xlim(-2*rmax,2*rmax)
fig.gca().set_ylim(-2*imax-0.1,2*imax+0.1)
print(eigenvalues)
plt.plot(np.real(eigenvalues),np.imag(eigenvalues), 'bx', markersize=3.0)
plt.tight_layout()
# plt.show()
fig.savefig(filebase+"/eigenvalues.pdf", bbox_inches='tight')
print("\nRuntime: %.1f s"%(timeit.default_timer()-start))
