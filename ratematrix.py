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
parser.add_argument("--accumulate", type=int, default=0, choices=[0,1], help='If 1, search filebase directory for npy files, generate sparse matrix, and plot eigenvalues. Default 0.')
parser.add_argument("--mechanism", type=str, required=False, default='mechanisms/h2o2.cti', dest='mechanism', help='Mechanism cti file. Default mechanisms/h2o2.cti.')
parser.add_argument("--reaction", type=int, required=False, default=None, dest='rind', help='Reaction index, provided for parallelization. If none is specified, the program will loop through all reactions in the model in sequence. Default None.')
parser.add_argument("--Nmax", type=int, required=False, default=3, dest='Nmax', help='Maximum number of molecules for each species. Default 3.')
parser.add_argument("--Nvals", type=int, required=False, default=100, dest='Nvals', help='Number of eigenvalues to calculate, when --accumulate 1 is set. Default 1000')
parser.add_argument("--progress", type=int, required=False, default=1, choices=[0,1], help='Print progress during calculation. Default 1.')
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
        if np.all(multiindex2>=0) and np.all(multiindex2<Nmax) and not np.isnan(k):
            rate=get_rate(multiindex,pstoi,k,reaction)
            j=get_index(multiindex2)
            data.append(rate)
            rows.append(i)
            columns.append(j)
            data.append(-rate)
            rows.append(i)
            columns.append(i)
    #save to files
    np.save("%s_data"%(filebase),data)
    np.save("%s_rows"%(filebase),rows)
    np.save("%s_columns"%(filebase),columns)

#Find all list of species that have specified number of atoms
def recursive_list(remaining_atoms, list=[], previously_enumerated=[]):
    avail=available_species(remaining_atoms)
    previously_enumerated.append(list_to_multiindex(list)
    if(avail==[]):
        if np.all(remaining_atoms == np.zeros(len(elements))):
            #Return the multiindex of the state given by the list of species
            return [list_to_multiindex(list)]
        else:
            #could not exaust atoms with this branch; return nothing
            return []
    else:
        #call recursive_list for each available species
        ret=[]
        for spec in avail:
            #Check that the new list of species has not already been enumerated
            if not (list_to_multiindex(list+[spec]) in previously_enumerated):
                ret+=recursive_list(remaining_atoms-np.array([gas.species()[spec].composition[el] if el in gas.species()[spec].composition.keys() else 0 for el in elements]),list+[spec],previously_enumerated)
        return ret

#Return indices of species that can be added given remaining atoms
def available_species(remaining_atoms):
    avail=[]
    for i in range(ns):
        datoms=remaining_atoms-np.array([gas.species()[i].composition[el] if el in gas.species()[i].composition.keys() else 0 for el in elements])
        if (np.all(datoms>=0)):
            avail.append(i)
    return avail

#Convert the list of species added sequentially to a multiindex
def list_to_multiindex(list):
    unique,counts=np.unique(list,return_counts=True)
    return [counts[np.where(unique==index)[0][0]] if index in unique else 0 for index in range(ns)]

#Main
filebase=args.filebase
if not os.path.isdir(filebase):
    os.mkdir(filebase)
mechanism=args.mechanism
Nmax=args.Nmax
Nvals=args.Nvals
rateindex=args.rind
progress=args.progress
accumulate=args.accumulate
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

if(accumulate==1):
    #Accumulate sparse data and plot
    p = re.compile("(.*)_data.npy")
    filebases=[item for sublist in np.unique([p.findall(name) for name in np.sort(os.listdir(filebase))]) for item in sublist]
    data=np.array([])
    columns=np.array([])
    rows=np.array([])
    for name in filebases:
        data=np.append(data,np.load(filebase+"/"+name+"_data.npy"))
        columns=np.append(columns,np.load(filebase+"/"+name+"_columns.npy"))
        rows=np.append(rows,np.load(filebase+"/"+name+"_rows.npy"))
    dim=int(np.max(rows)+1)
    print("Sparsity: %f"%(1-len(data)*1.0/(dim)**2))
    print("Average non-zero entry: %f"%(np.linalg.norm(data)/len(data)))
    ratematrix=coo_matrix((np.array(data),(np.array(rows),np.array(columns))),(dim,dim))
    eigenvalues,eigenvectors=eigs(ratematrix,k=Nvals,which='LR')
    np.save(filebase+"/eigenvalues",eigenvalues)
    np.save(filebase+"/eigenvectors",eigenvectors)
    rmax=np.abs(np.max(np.real(eigenvalues)))
    imax=np.max(np.abs(np.imag(eigenvalues[np.where(np.real(eigenvalues)>-2*rmax)[0]])))
    fig=plt.figure()
    fig.gca().set_ylabel(r'$\mathrm{Im}\left(\lambda\right)$')
    fig.gca().set_xlabel(r'$\mathrm{Re}\left(\lambda\right)$')
    fig.gca().set_xlim(-2*rmax,2*rmax)
    fig.gca().set_ylim(-2*imax,2*imax)
    plt.plot(np.real(eigenvalues),np.imag(eigenvalues), 'bx', markersize=3.0)
    plt.tight_layout()
    plt.show()
    fig.savefig(filebase+"/eigenvalues.pdf", bbox_inches='tight')
    print("\nRuntime: %.1f s"%(timeit.default_timer()-start))
elif rateindex == None:
    #Loop through each reaction index and calculate spase elements
    atoms=np.array(atoms)
    multiindices=recursive_list(atoms)
    print("\nGenerated %i-dimensional space in %.1f s"%(len(multiindices), timeit.default_timer()-start))

    for rind in range(nr):
        # print('Reaction %i'%(rind))
        calculate_sparse_elements(rind,filebase+"/%i"%(rind))
        # if(progress==1):
            # print('')
    print("\nRuntime: %.1fs"%(timeit.default_timer()-start))
else:
    #Caclulate sparse elements for specified reaction index
    atoms=np.array(atoms)
    multiindices=get_states(atoms)
    print("\nGenerated %i-dimensional space in %.1f s"%(len(multiindices), timeit.default_timer()-start))
    #print('Reaction %i'%(rateindex))
    calculate_sparse_elements(rateindex, filebase)
    print("\nRuntime: %.1fs"%(timeit.default_timer()-start))
