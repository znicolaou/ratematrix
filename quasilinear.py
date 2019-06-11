#!/usr/bin/env python
import sys
import numpy as np
import cantera as ct
import timeit
import os
import argparse
from scipy.sparse import coo_matrix

def runsim_nosense ( tmax, temp, pres, phi, dilute ):
    dt = tmax/Npoints
    gas.X = ('CH4:%f, O2:1, N2:4, Ar:%f'%(phi, dilute))
    gas.TP = temp, pres
    r = ct.IdealGasReactor(gas, name='R1')
    sim = ct.ReactorNet([r])

    sim.rtol = 1.0e-8
    sim.atol = 1.0e-14

    states = ct.SolutionArray(gas, extra=['t', 'A'])
    for i in range(0, int(Npoints)):
        t=i*dt
        sim.advance(t)
        print("%5f"%(t/tmax), end="\r")
        mat=np.zeros((gas.n_species,gas.n_species))

    #We would go faster to do this in the sparse matrix fashion - only loop over the reactions.
    data=[]
    rows=[]
    columns=[]
    for n in range(gas.n_species):
        for rind in range(gas.n_reactions):
            reaction=gas.reactions()[rind]
            if species[n] in reaction.reactants.keys():
                rstoi=np.array([reaction.reactants[x] if x in reaction.reactants.keys() else 0 for x in species])
                num=len(np.where(rstoi>0)[0])
                for m in np.where(rstoi>0)[0]:
                    rstoi[m]-=1
                    #forward reaction
                    data.append(-1.0/num*gas.forward_rate_constants[rind]*np.product(gas.X**(rstoi)))
                    rows.append(n)
                    columns.append(m)
                    #reverse reaction
                    data.append(1.0/num*gas.reverse_rate_constants[rind]*np.product(gas.X**(rstoi)))
                    rows.append(m)
                    columns.append(n)
                    rstoi[m]+=1

            if species[n] in reaction.products.keys():
                rstoi=np.array([reaction.reactants[x] if x in reaction.reactants.keys() else 0 for x in species])
                num=len(np.where(rstoi>0)[0])
                for m in np.where(rstoi>0)[0]:
                    rstoi[m]-=1
                    #forward reaction
                    data.append(1.0/num*gas.forward_rate_constants[rind]*np.product(gas.X**(rstoi)))
                    rows.append(n)
                    columns.append(m)
                    #reverse reaction
                    data.append(-1.0/num*gas.reverse_rate_constants[rind]*np.product(gas.X**(rstoi)))
                    rows.append(m)
                    columns.append(n)
                    rstoi[m]+=1

            mat =coo_matrix((np.array(data),(np.array(rows),np.array(columns))),(int(gas.n_species),int(gas.n_species)))
            states.append(r.thermo.state, t=t, A=mat.toarray())

    return states

#Command-line arguments
parser = argparse.ArgumentParser(description='Generate a sparse rate matrix from cantera model.')
parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for npy file output. Three files will be created for each reaction, storing rates, row indices, and column indices.')
parser.add_argument("--mechanism", type=str, required=False, default='mechanisms/gri30.cti', dest='mechanism')
parser.add_argument("--Npoints", type=float, required=False, default=1e2, dest='Npoints')
parser.add_argument("--tmax", type=float, required=False, default=1e-3, dest='tmax')
parser.add_argument("--temperature", type=float, required=False, default=1500, dest='temperature')
parser.add_argument("--pressure", type=float, required=False, default=1, dest='pressure')
parser.add_argument("--dilute", type=float, required=False, default=0.1, dest='dilute')
parser.add_argument("--phi", type=float, required=False, default=2.0, dest='phi')

args = parser.parse_args()

start=timeit.default_timer()
gas = ct.Solution(args.mechanism)
reactions=gas.reactions()
Npoints = args.Npoints
filebase = args.filebase
temp = args.temperature
pres = args.pressure
phi = args.phi
dilute=args.dilute
tmax=args.tmax

species=species=gas.species_names
reactions = gas.reactions()

observation=runsim_nosense(tmax, temp, pres, phi, dilute)

f=open('%s.dat'%(filebase),'w')
print('species', *species, sep=' ', file=f)
print('t', *(observation.t), sep=' ', file=f)
print('T', *(observation.T), sep=' ', file=f)
print('P', *(observation.P), sep=' ', file=f)

rates=np.transpose(observation.net_rates_of_progress)
for j in range(0, gas.n_reactions):
    print('%s->%s'%(str(reactions[j].reactants).replace(" ","").replace("\"","").replace("'",""),str(reactions[j].products).replace(" ","").replace("\"","").replace("'","")), *(rates[j]), sep=' ', file=f)

for j in range(0, gas.n_total_species):
    print('[%s]'%(gas.species_names[j]), *([item for sublist in observation(gas.species_names[j]).X for item in sublist]), sep=' ', file=f) #flatten that list

f.close()

stop=timeit.default_timer()
print ('runtime: %f'%(stop-start))
sys.stdout.flush()
