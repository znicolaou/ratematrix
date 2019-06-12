#!/usr/bin/env python
import sys
import numpy as np
import cantera as ct
import timeit
import os
import argparse
from scipy.sparse import coo_matrix
from scipy.linalg import eig


def runsim_nosense ( tmax, temp, pres, ics ):
    dt = tmax/Npoints

    gas.TPX=temp,pres*ct.one_atm,ics
    print(gas())


    r = ct.IdealGasReactor(gas, name='R1')
    sim = ct.ReactorNet([r])

    sim.rtol = 1.0e-6
    sim.atol = 1.0e-14

    states = ct.SolutionArray(gas, extra=['t', 'A', 'evals', 'evecs', 'frates', 'rrates'])
    for i in range(0, int(Npoints)):
        t=i*dt
        sim.advance(t)
        print("%5f"%(t/tmax), end="\n\n")
        data=[]
        rows=[]
        columns=[]
        for n in range(gas.n_species):
            for rind in range(gas.n_reactions):
                reaction=gas.reactions()[rind]
                rstoi=np.array([reaction.reactants[x] if x in reaction.reactants.keys() else 0 for x in species])
                pstoi=np.array([reaction.products[x] if x in reaction.products.keys() else 0 for x in species])
                if reaction.reaction_type == 1: #bimolecular
                    if species[n] in reaction.reactants.keys():
                        #forward reaction with species[n] a reactant
                        for m in np.where(rstoi>0)[0]:
                            order=rstoi[m]
                            rstoi[m]-=1 #remaining reactants after removing species[m]
                            data.append(-order*0.5*gas.forward_rate_constants[rind]*np.product(gas.X**(rstoi)))
                            rows.append(n)
                            columns.append(m)
                            rstoi[m]+=1
                        #reverse reaction with species[n] a product
                        for m in np.where(pstoi>0)[0]:
                            order=pstoi[m]
                            pstoi[m]-=1  #remaining reactants after removing species[m]
                            data.append(order*0.5*gas.reverse_rate_constants[rind]*np.product(gas.X**(pstoi)))
                            rows.append(n)
                            columns.append(m)
                            pstoi[m]+=1
                    if species[n] in reaction.products.keys():
                        #forward reaction with species[n] a product
                        for m in np.where(rstoi>0)[0]:
                            order=rstoi[m]
                            rstoi[m]-=1 #remaining reactants after removing species[m]
                            data.append(order*0.5*gas.forward_rate_constants[rind]*np.product(gas.X**(rstoi)))
                            rows.append(n)
                            columns.append(m)
                            rstoi[m]+=1
                        #reverse reaction with species[n] a reactant
                        for m in np.where(pstoi>0)[0]:
                            order=pstoi[m]
                            pstoi[m]-=1 #remaining reactants after removing species[m]
                            data.append(-order*0.5*gas.reverse_rate_constants[rind]*np.product(gas.X**(pstoi)))
                            rows.append(n)
                            columns.append(m)
                            pstoi[m]+=1
                else: #three body reaction
                # We need to do this third body loop in each part below.
                # for third_body in range(ns):
                #     rstoi[third_body]+=1
                #     if np.all(multiindex>=rstoi):
                #         efficiency=reaction.default_efficiency
                #         if species[third_body] in reaction.efficiencies.keys():
                #             efficiency=reaction.efficiencies[species[third_body]]
                #         # ret+=efficiency*k*np.product(factorial(multiindex)/factorial(multiindex-rstoi))/(gas.volume_mole*np.sum(multiindex))**2
                #         ret+=efficiency*k*np.product(factorial(multiindex)/factorial(multiindex-rstoi))/(ct.avogadro*vol)**2
                #     rstoi[third_body]-=1
                    if species[n] in reaction.reactants.keys():
                        #forward reaction with species[n] a reactant
                        for m in np.where(rstoi>0)[0]:
                            order=rstoi[m]
                            rstoi[m]-=1 #remaining reactants after removing species[m]
                            data.append(-order*0.5*gas.forward_rate_constants[rind]*np.product(gas.X**(rstoi)))
                            rows.append(n)
                            columns.append(m)
                            rstoi[m]+=1
                        #reverse reaction with species[n] a product
                        for m in np.where(pstoi>0)[0]:
                            order=pstoi[m]
                            pstoi[m]-=1  #remaining reactants after removing species[m]
                            data.append(order*0.5*gas.reverse_rate_constants[rind]*np.product(gas.X**(pstoi)))
                            rows.append(n)
                            columns.append(m)
                            pstoi[m]+=1
                    if species[n] in reaction.products.keys():
                        #forward reaction with species[n] a product
                        for m in np.where(rstoi>0)[0]:
                            order=rstoi[m]
                            rstoi[m]-=1 #remaining reactants after removing species[m]
                            data.append(order*0.5*gas.forward_rate_constants[rind]*np.product(gas.X**(rstoi)))
                            rows.append(n)
                            columns.append(m)
                            rstoi[m]+=1
                        #reverse reaction with species[n] a reactant
                        for m in np.where(pstoi>0)[0]:
                            order=pstoi[m]
                            pstoi[m]-=1 #remaining reactants after removing species[m]
                            data.append(-order*0.5*gas.reverse_rate_constants[rind]*np.product(gas.X**(pstoi)))
                            rows.append(n)
                            columns.append(m)
                            pstoi[m]+=1
        quit()

        mat =coo_matrix((np.array(data),(np.array(rows),np.array(columns))),(int(gas.n_species),int(gas.n_species))).toarray()
        eigenvalues,eigenvectors=eig(mat)
        states.append(r.thermo.state, t=t, A=mat, evals=eigenvalues, evecs=eigenvectors, frates=gas.forward_rate_constants, rrates=gas.reverse_rate_constants)

    return states

#Command-line arguments
parser = argparse.ArgumentParser(description='Generate a sparse rate matrix from cantera model.')
parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for npy file output. Three files will be created for each reaction, storing rates, row indices, and column indices.')
parser.add_argument("--mechanism", type=str, required=False, default='mechanisms/gri30.cti', dest='mechanism')
parser.add_argument("--Npoints", type=float, required=False, default=1e2, dest='Npoints')
parser.add_argument("--tmax", type=float, required=False, default=1e-3, dest='tmax')
parser.add_argument("--temperature", type=float, required=False, default=1500, dest='temperature')
parser.add_argument("--pressure", type=float, required=False, default=1, dest='pressure')
parser.add_argument("--ics", type=float, nargs='+', required=False, default=[13, 1, 3, 2, 47, 8], help='Initial mole fractions. Default 13 1 3 2 47 8 stoichiometric methane in air.')

args = parser.parse_args()

start=timeit.default_timer()
gas = ct.Solution(args.mechanism)
reactions=gas.reactions()
Npoints = args.Npoints
filebase = args.filebase
temp = args.temperature
pres = args.pressure
ics = args.ics
tmax=args.tmax
ns=gas.n_total_species
nr=gas.n_reactions

multiindex=np.zeros(ns)
for i in range(0,len(args.ics),2):
    multiindex[int(args.ics[i])]=args.ics[i+1]

species=species=gas.species_names
reactions = gas.reactions()

observation=runsim_nosense(tmax, temp, pres, multiindex)

f=open('%sout.dat'%(filebase),'w')
print(Npoints, ns, nr, temp, pres, tmax, file=f)
print(*species, file=f)
print(*(observation.t), file=f)
print(*(observation.T), file=f)
print(*(observation.P), file=f)
for j in range(0, ns):
    print(*([item for sublist in observation(gas.species_names[j]).X for item in sublist]), sep=' ', file=f)
f.close()

np.save(filebase+"matrices.npy",observation.A)
np.save(filebase+"eigenvalues.npy",observation.evals)
np.save(filebase+"eigenvectors.npy",observation.evecs)
np.save(filebase+"frates.npy",observation.frates)
np.save(filebase+"rrates.npy",observation.rrates)

stop=timeit.default_timer()
print ('runtime: %f'%(stop-start))
sys.stdout.flush()
