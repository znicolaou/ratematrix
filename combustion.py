#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import cantera as ct
import timeit
import argparse
from scipy.sparse import coo_matrix


def runsim (times):
    if(pressureflag==0):
        if(energyflag==1):
            r = ct.IdealGasReactor(gas, name='R1')
        else:
            r = ct.IdealGasReactor(gas, name='R1', energy='off')
    else:
        if(energyflag==1):
            r = ct.IdealGasConstPressureReactor(gas, name='R1')
        else:
            r = ct.IdealGasConstPressureReactor(gas, name='R1', energy='off')
    sim = ct.ReactorNet([r])

    sim.rtol = 1.0e-6
    sim.atol = 1.0e-14

    if(sensflag==1):
        nr=gas.n_reactions
        for i in range(0,nr):
            r.add_sensitivity_reaction(i)
        sim.rtol_sensitivity = 1.0e-6
        sim.atol_sensitivity = 1.0e-6

    states = ct.SolutionArray(gas, extra=['t','matrices', 'evals', 'evecs', 'sens', 'norm'])

    for t in times:
        sim.advance(t)

        sensitivities=[]
        data=[]
        rows=[]
        columns=[]
        mat=[]
        eigenvalues=[]
        eigenvectors=[]
        norm=0
        if(sensflag==1):
            for j in range(0, gas.n_total_species):
                for k in range(0, gas.n_reactions):
                    sensitivities.append(sim.sensitivity(gas.species_names[j], k))

        if(quasilinear==1):
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
                        if k>0 and np.isfinite(k):
                            data.append(-order/num*k*np.product(concentrations**(rstoi)))
                            rows.append(n)
                            columns.append(m)
                        rstoi[m]+=1
                    #reverse reaction with species[n] a product
                    for m in np.where(pstoi>0)[0]:
                        order=pstoi[m]
                        num=np.sum(pstoi)
                        pstoi[m]-=1  #remaining reactants after factoring out species[m]
                        k=gas.reverse_rate_constants[rind]
                        if k>0 and np.isfinite(k):
                            data.append(order/num*k*np.product(concentrations**(pstoi)))
                            rows.append(n)
                            columns.append(m)
                        pstoi[m]+=1
                for n in np.where(pstoi>0)[0]: #row of the matrix, equation for species n
                    #forward reaction with species[n] a product
                    for m in np.where(rstoi>0)[0]:
                        order=rstoi[m]
                        num=np.sum(rstoi)
                        rstoi[m]-=1 #remaining reactants after removing species[m]
                        k=gas.forward_rate_constants[rind]
                        if k>0 and np.isfinite(k):
                            data.append(order/num*k*np.product(concentrations**(rstoi)))
                            rows.append(n)
                            columns.append(m)
                        rstoi[m]+=1
                    #reverse reaction with species[n] a reactant
                    for m in np.where(pstoi>0)[0]:
                        order=pstoi[m]
                        num=np.sum(pstoi)
                        pstoi[m]-=1 #remaining reactants after removing species[m]
                        k=gas.reverse_rate_constants[rind]
                        if k>0 and np.isfinite(k):
                            data.append(-order/num*k*np.product(concentrations**(pstoi)))
                            rows.append(n)
                            columns.append(m)
                        pstoi[m]+=1

            mat =coo_matrix((np.array(data),(np.array(rows),np.array(columns))),(int(gas.n_species),int(gas.n_species))).toarray()
            eigenvalues,eigenvectors=np.linalg.eig(mat)
            u, singularvalues, vh = np.linalg.svd(mat)
            norm=(np.sum(np.abs(eigenvalues)**2)/np.sum(np.abs(singularvalues)**2))**(0.5)
            sorted=np.argsort(np.abs(eigenvalues))
            eigenvalues=eigenvalues[sorted]
            eigenvectors=eigenvectors[sorted]
        states.append(r.thermo.state, t=t, sens=sensitivities, matrices=mat, evals=eigenvalues, evecs=eigenvectors, norm=norm)

    return states

parser = argparse.ArgumentParser(description='Generate a sparse rate matrix from cantera model.')
parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for npy file output. Three files will be created for each reaction, storing rates, row indices, and column indices.')
parser.add_argument("--temperature", type=float, required=False, default=1000, help='Temperature in Kelvin. Default 1000.')
parser.add_argument("--adiabatic", type=int, choices=[0, 1, 2], required=False, default=1, help='Convert energy from reactions to heat. The values 0, 1, and 2 correspond to constant volume/temperature, constant volume/energy, and constant pressure/enthalpy, respectively. The temperature is specify the reference multiindix specified with --reference. ')
parser.add_argument("--reference", type=int, nargs='+', required=False, default=[0, 10, 3, 5, 8, 100], help='Reference multiindex for which the temperature and number of atoms are specified. Default 0 4 3 2 8 0.')
parser.add_argument("--pressure", type=float, required=False, default=1, help='Pressure in atm. Default 1.')
parser.add_argument("--t0", type=float, required=False, default=1e-8, help='Initial integration time for propogating.')
parser.add_argument("--tmax", type=float, required=False, default=1, help='Final integration time for propogating.')
parser.add_argument("--Npoints", type=int, required=False, default=100, help='Number of times to propogate.')
parser.add_argument("--sensitivities", type=int, required=False, choices=[0,1],default=0, help='Flag to calculate sensitivities.')
parser.add_argument("--quasilinear", type=int, required=False, choices=[0,1],default=0, help='Flag to calculate linearized matrix and eigenvalues.')
parser.add_argument("--mechanism", type=str, required=False,default='h2o2.cti', help='Cantera mechanism file.')
args = parser.parse_args()

start=timeit.default_timer()
Npoints = args.Npoints
t0=args.t0
tmax=args.tmax
temperature=args.temperature
pressure=args.pressure*ct.one_atm
sensflag=args.sensitivities
quasilinear=args.quasilinear

if args.adiabatic == 0:
    energyflag=0
    pressureflag=1
elif args.adiabatic == 1:
    energyflag=1
    pressureflag=0
else:
    energyflag=1
    pressureflag=1
out=args.filebase

#Change initial conditions here
gas = ct.Solution(args.mechanism)
species=gas.species_names
reactions = gas.reactions()
ns=gas.n_species
refmultiindex=np.zeros(ns,dtype=int)
for i in range(0,len(args.reference),2):
    refmultiindex[args.reference[i]]=args.reference[i+1]
gas.TPX = temperature, pressure, refmultiindex
times=[t0*(tmax/t0)**(n*1.0/Npoints) for n in range(Npoints)]

observation=runsim(times)


np.save(args.filebase+"times.npy", observation.t)
np.save(args.filebase+"concentrations.npy", observation.X)
np.save(args.filebase+"temperatures.npy", observation.T)
np.save(args.filebase+"pressures.npy", observation.P/ct.one_atm)
# if quasilinear==1:
np.save(args.filebase+"matrices.npy", observation.matrices)
np.save(args.filebase+"eigenvalues.npy", observation.evals)
np.save(args.filebase+"eigenvectors.npy", observation.evecs)
# if sensflag==1:
np.save(args.filebase+"sensitivities.npy", observation.sens)
np.save(args.filebase+"norms.npy", observation.norm)

stop=timeit.default_timer()
print ('runtime: %f'%(stop-start))
sys.stdout.flush()
