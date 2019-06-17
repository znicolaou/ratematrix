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

    r = ct.IdealGasReactor(gas, energy='on')
    # r = ct.IdealGasConstPressureReactor(gas, energy='off')
    sim = ct.ReactorNet([r])

    sim.rtol = 1.0e-8
    sim.atol = 1.0e-14
    output=False
    states = ct.SolutionArray(gas, extra=['t', 'A', 'evals', 'evecs', 'frates', 'rrates'])
    for i in range(0, int(Npoints)):
        t=i*dt
        sim.advance(t)
        print("%3f"%(t/tmax), end="\r")
        data=[]
        rows=[]
        columns=[]
        #Enumerate over reactions and add corresponding sparse terms
        concentrations=gas.X/(ct.gas_constant*gas.T/gas.P)
        for rind in range(len(gas.reactions())):
            reaction=gas.reactions()[rind]
            rstoi=np.array([reaction.reactants[x] if x in reaction.reactants.keys() else 0 for x in species])
            pstoi=np.array([reaction.products[x] if x in reaction.products.keys() else 0 for x in species])
            for n in np.where(rstoi>0)[0]: #row of the matrix, equation for species n
                #forward reaction with species[n] a reactant
                for m in np.where(rstoi>0)[0]:
                    order=rstoi[m]
                    num=np.sum(rstoi)
                    rstoi[m]-=1 #remaining reactants after factoring out species[m]
                    k=gas.forward_rate_constants[rind]
                    data.append(-order/num*k*np.product(concentrations**(rstoi)))
                    rows.append(n)
                    columns.append(m)
                    if output:
                        print("forward reactant", rind, order, num, k, -order/num*k*np.product(concentrations**(rstoi)), species[n], species[m])
                    rstoi[m]+=1
                #reverse reaction with species[n] a product
                for m in np.where(pstoi>0)[0]:
                    order=pstoi[m]
                    num=np.sum(pstoi)
                    pstoi[m]-=1  #remaining reactants after factoring out species[m]
                    k=gas.reverse_rate_constants[rind]
                    data.append(order/num*k*np.product(concentrations**(pstoi)))
                    rows.append(n)
                    columns.append(m)
                    if output:
                        print("reverse product", rind, order, num, k, rind, order/num*k*np.product(concentrations**(pstoi)), species[n], species[m])
                    pstoi[m]+=1
            for n in np.where(pstoi>0)[0]: #row of the matrix, equation for species n
                #forward reaction with species[n] a product
                for m in np.where(rstoi>0)[0]:
                    order=rstoi[m]
                    num=np.sum(rstoi)
                    rstoi[m]-=1 #remaining reactants after removing species[m]
                    k=gas.forward_rate_constants[rind]
                    data.append(order/num*k*np.product(concentrations**(rstoi)))
                    rows.append(n)
                    columns.append(m)
                    if output:
                        print("forward product", rind, order, num, k, order/num*k*np.product(concentrations**(rstoi)), species[n], species[m])
                    rstoi[m]+=1
                #reverse reaction with species[n] a reactant
                for m in np.where(pstoi>0)[0]:
                    order=pstoi[m]
                    num=np.sum(pstoi)
                    pstoi[m]-=1 #remaining reactants after removing species[m]
                    k=gas.reverse_rate_constants[rind]
                    data.append(-order/num*k*np.product(concentrations**(pstoi)))
                    rows.append(n)
                    columns.append(m)
                    if output:
                        print("reverse reactant", rind, order, num, k, -order/num*k*np.product(concentrations**(pstoi)), species[n], species[m])
                    pstoi[m]+=1

        mat =coo_matrix((np.array(data),(np.array(rows),np.array(columns))),(int(gas.n_species),int(gas.n_species))).toarray()
        eigenvalues,eigenvectors=eig(mat)
        sorted=np.argsort(np.abs(eigenvalues))
        states.append(r.thermo.state, t=t, A=mat, evals=eigenvalues[sorted], evecs=eigenvectors[sorted], frates=gas.forward_rate_constants, rrates=gas.reverse_rate_constants)

    return states

#Command-line arguments
parser = argparse.ArgumentParser(description='Generate a sparse rate matrix from cantera model.')
parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for npy file output. Three files will be created for each reaction, storing rates, row indices, and column indices.')
parser.add_argument("--mechanism", type=str, required=False, default='mechanisms/gri30.cti', dest='mechanism')
parser.add_argument("--Npoints", type=float, required=False, default=1e2, dest='Npoints')
parser.add_argument("--time", type=float, required=False, default=3e-3, dest='tmax')
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

species=gas.species_names
reactions = gas.reactions()

observation=runsim_nosense(tmax, temp, pres, multiindex)

f=open('%sout.dat'%(filebase),'w')
print(Npoints, ns, nr, temp, pres, tmax, file=f)
print(*species, file=f)
print(*(observation.t), file=f)
print(*(observation.T), file=f)
print(*(observation.P)/ct.one_atm, file=f)

concentrations=np.array([observation(gas.species_names[j]).X[:,0] for j in range(ns)])
np.save(filebase+"matrices.npy",observation.A)
np.save(filebase+"eigenvalues.npy",observation.evals)
np.save(filebase+"eigenvectors.npy",observation.evecs)
np.save(filebase+"concentrations.npy", concentrations)

stop=timeit.default_timer()
print ('runtime: %f'%(stop-start))
sys.stdout.flush()
