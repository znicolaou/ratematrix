#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import cantera as ct
import timeit
import argparse

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

    states = ct.SolutionArray(gas, extra=['t','sens'])

    for t in times:
        sim.advance(t)

        sensitivities=[]
        if(sensflag==1):
            for j in range(0, gas.n_total_species):
                for k in range(0, gas.n_reactions):
                    sensitivities.append(sim.sensitivity(gas.species_names[j], k))

        states.append(r.thermo.state, t=t, sens=sensitivities)
    return states

parser = argparse.ArgumentParser(description='Generate a sparse rate matrix from cantera model.')
parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for npy file output. Three files will be created for each reaction, storing rates, row indices, and column indices.')
parser.add_argument("--temperature", type=float, required=False, default=1500, help='Temperature in Kelvin. Default 1500.')
parser.add_argument("--adiabatic", type=int, choices=[0, 1, 2], required=False, default=0, help='Convert energy from reactions to heat. The values 0, 1, and 2 correspond to constant volume/temperature, constant volume/energy, and constant pressure/enthalpy, respectively. The temperature is specify the reference multiindix specified with --reference. ')
parser.add_argument("--reference", type=int, nargs='+', required=False, default=[0, 4, 3, 2, 8, 0], help='Reference multiindex for which the temperature and number of atoms are specified. Default 0 4 3 2 8 0.')
parser.add_argument("--pressure", type=float, required=False, default=1, help='Pressure in atm. Default 1.')
parser.add_argument("--t0", type=float, required=False, default=1e-8, help='Initial integration time for propogating.')
parser.add_argument("--tmax", type=float, required=False, default=1e-2, help='Final integration time for propogating.')
parser.add_argument("--Npoints", type=int, required=False, default=25, help='Number of times to propogate.')
parser.add_argument("--sensitivities", type=int, required=False, choices=[0,1],default=1, help='Flag to calculate sensitivities.')
args = parser.parse_args()

start=timeit.default_timer()
Npoints = args.Npoints
t0=args.t0
tmax=args.tmax
temperature=args.temperature
pressure=args.pressure*ct.one_atm
sensflag=args.sensitivities
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
gas = ct.Solution('h2o2.cti')
ns=gas.n_species
refmultiindex=np.zeros(ns,dtype=int)
for i in range(0,len(args.reference),2):
    refmultiindex[args.reference[i]]=args.reference[i+1]
gas.TPX = temperature, pressure, refmultiindex
times=[t0*(tmax/t0)**(n*1.0/Npoints) for n in range(Npoints)]

observation=runsim(times)

species = gas.species_names
reactions = gas.reactions()
f=open(out,'w')
outarray=[]
outarray.append(observation.t)
outarray.append(observation.T)
outarray.append((observation.P/ct.one_atm))
columns='t,T,P'
for j in range(0, gas.n_total_species):
    columns=columns + ',[%s]'%(gas.species_names[j])
    outarray.append([item for sublist in observation(gas.species_names[j]).X for item in sublist])
for j in range(0, gas.n_reactions):
    columns=columns + ',%s Flux'%(str(reactions[j].reactants).replace(" ","").replace(",","+").replace("{","").replace("}","")+"->"+str(reactions[j].products).replace(" ","").replace(",","+").replace("{","").replace("}",""))
    outarray.append(np.transpose(observation.net_rates_of_progress)[j])
if(sensflag==1):
    for j in range(0, gas.n_total_species):
        for k in range(0, gas.n_reactions):
            columns=columns + ',%s %s-Sensitivity'%(str(reactions[k].reactants).replace(" ","").replace(",","+").replace("{","").replace("}","")+"->"+str(reactions[k].products).replace(" ","").replace(",","+").replace("{","").replace("}",""), gas.species_names[j])
            outarray.append(np.array(observation.sens).transpose()[j*gas.n_reactions+k])

np.savetxt(out, np.transpose(outarray), delimiter=',', header=columns, comments='')

stop=timeit.default_timer()
print ('runtime: %f'%(stop-start))
sys.stdout.flush()
