#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import cantera as ct
import timeit

def runsim ( ):
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

    for t in np.arange(0, tmax, dt):
        sim.advance(t)

        sensitivities=[]
        if(sensflag==1):
            for j in range(0, gas.n_total_species):
                for k in range(0, gas.n_reactions):
                    sensitivities.append(sim.sensitivity(gas.species_names[j], k))

        states.append(r.thermo.state, t=t, sens=sensitivities)
    return states

if(len(sys.argv) != 12):
    print(len(sys.argv))
    print('usage: python2 experiment.py [Npoints] [tmax] [temperature] [pressure] [phi] [dilute] [sensitivities] [energy] [constantpressure] [methane] [out]')
    print('Npoints is the number of points between 0 and tmax to output')
    print('tmax is maximum time to integrate')
    print('temperature is initial temperature')
    print('pressure is initial pressure')
    print('phi is initial ratio of fuel/O2 (one OH per hundred O2 is used initially to start chain reaction)')
    print('dilute is initial ratio of Ar/O2')
    print('sensitivities is 1 to calculate sensitivities, 0 otherwise')
    print('energy is 1 for adiabatic, 0 for isothermal')
    print('constantpressure is 1 for constant pressure, 0 for constant volume')
    print('methane is 1 for methane combustion, 0 for hydrogen volume')
    print('out is output file')
    print('example: python combustion.py  10000 2e-6 2000 100 1 0 0 1 0 1 ch4test.csv')
    print('example: python combustion.py  10000 2e-8 2000 100 2 0 1 1 0 0 h2test.csv')
    exit()

start=timeit.default_timer()
Npoints = int(sys.argv[1])
tmax=float(sys.argv[2])
temperature=float(sys.argv[3])
pressure=float(sys.argv[4])*ct.one_atm
phi=float(sys.argv[5])
dilute=float(sys.argv[6])
sensflag=int(sys.argv[7])
energyflag=int(sys.argv[8])
pressureflag=int(sys.argv[9])
methaneflag=int(sys.argv[10])
out=sys.argv[11]

#Change initial conditions here
if (methaneflag==1):
    gas = ct.Solution('gri30.cti')
    gas.X = ('CH4:%f, O2:1, N2:4, OH:0.01, Ar:%f'%(phi, 0.05+dilute))
else:
    gas = ct.Solution('h2o2.cti')
    gas.X = ('H2:%f, O2:1, OH:0.0, Ar:%f'%(phi, dilute))

dt = tmax/Npoints
gas.TP = temperature, pressure

observation=runsim()

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
