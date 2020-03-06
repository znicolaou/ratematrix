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

    states = ct.SolutionArray(gas, extra=['t','matrices', 'evals', 'evecs'])

    for t in times:
        sim.advance(t)

        data=[]
        rows=[]
        columns=[]
        mat=[]
        eigenvalues=[]
        eigenvectors=[]

        states.append(r.thermo.state, t=t, matrices=mat)

    return states

parser = argparse.ArgumentParser(description='Generate a sparse rate matrix from cantera model.')
parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for npy file output. Three files will be created for each reaction, storing rates, row indices, and column indices.')
parser.add_argument("--temperature", type=float, required=False, default=1000, help='Temperature in Kelvin. Default 1000.')
parser.add_argument("--adiabatic", type=int, choices=[0, 1, 2], required=False, default=0, help='Convert energy from reactions to heat. The values 0, 1, and 2 correspond to constant pressure/temperature, constant volume/energy, and constant pressure/enthalpy, respectively. The temperature is specify the reference multiindix specified with --reference. ')
parser.add_argument("--reference", type=int, nargs='+', required=False, default=[0, 20, 3, 10, 4, 0, 8, 200], help='Reference multiindex for which the temperature and number of atoms are specified. Default 0 20 3 10 4 1 8 200.')
parser.add_argument("--pressure", type=float, required=False, default=1, help='Pressure in atm. Default 1.')
parser.add_argument("--t0", type=float, required=False, default=1e-8, help='Initial integration time for propogating.')
parser.add_argument("--tmax", type=float, required=False, default=1, help='Final integration time for propogating.')
parser.add_argument("--Npoints", type=int, required=False, default=100, help='Number of times to propogate.')
parser.add_argument("--mechanism", type=str, required=False,default='h2o2.cti', help='Cantera mechanism file.')
args = parser.parse_args()

start=timeit.default_timer()
Npoints = args.Npoints
t0=args.t0
tmax=args.tmax
temperature=args.temperature
pressure=args.pressure*ct.one_atm

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

ns=gas.n_species
nr=gas.n_reactions
species=gas.species_names
Ar=np.zeros((ns,ns))
Ap=np.zeros((ns,ns))

uninu=[]
binu=[]
count=0
totrstoi=np.zeros(ns,dtype=int)
totpstoi=np.zeros(ns,dtype=int)
for rind in range(nr):
    reaction=gas.reactions()[rind]
    rstoi=np.array([reaction.reactants[x] if x in reaction.reactants.keys() else 0 for x in species],dtype=int)
    pstoi=np.array([reaction.products[x] if x in reaction.products.keys() else 0 for x in species],dtype=int)

    totrstoi+=rstoi
    totpstoi+=pstoi
    if reaction.reversible:
        totrstoi+=pstoi
        totpstoi+=rstoi

    print(rind,rstoi,pstoi)
    print("reaction",rind,"type",gas.reaction_type(rind))


    if(gas.reaction_type(rind) == 1):
        rinds=np.where(rstoi>0)[0]
        pinds=np.where(pstoi>0)[0]
        #forward reactions
        print("A=",reaction.rate.pre_exponential_factor,"b=",reaction.rate.temperature_exponent,"E=",reaction.rate.activation_energy)
        k=reaction.rate.pre_exponential_factor*temperature**reaction.rate.temperature_exponent*np.exp(-reaction.rate.activation_energy/(temperature*ct.gas_constant))
        print("k=",k)

        if(len(rinds)==1) :
            print("unispecies indices", rinds)
            uninu.append([rinds,rstoi[rinds]])
            Ar[rinds[0],rinds[0]]=k*rstoi[rinds[0]]
        elif(len(rinds)==2) :
            print("bispecies indices", rinds)
            binu.append([rinds,rstoi[rinds]])
            Ar[rinds[0],rinds[0]]=k*rstoi[rinds[0]]
        else:
            count=count+1


        #reverse reactions
        if reaction.reversible:
            if(len(pinds)==1) :
                print("unispecies indices", pinds)
                uninu.append([pinds,pstoi[pinds]])
            elif(len(pinds)==2) :
                print("bispecies indices", pinds)
                binu.append([pinds,pstoi[pinds]])
            else:
                count=count+1

            k=gas.equilibrium_constants[i]/k
            print("k=",k)
    else:
        count=count+1
        if reaction.reversible:
            count=count+1

    # mat =coo_matrix((np.array(data),(np.array(rows),np.array(columns))),(int(gas.n_species),int(gas.n_species))).toarray()
    # eigenvalues,eigenvectors=np.linalg.eig(mat)
    # sorted=np.argsort(np.abs(eigenvalues))
    # eigenvalues=eigenvalues[sorted]
    # eigenvectors=eigenvectors[sorted]

# observation=runsim(times)


# np.save(args.filebase+"times.npy", observation.t)
# np.save(args.filebase+"concentrations.npy", observation.X)
# np.save(args.filebase+"temperatures.npy", observation.T)
# np.save(args.filebase+"pressures.npy", observation.P/ct.one_atm)
print("Times species appeared as reactants")
print(totrstoi)
print("Times species appeared as products")
print(totpstoi)
print("Total reactants and products")
print(np.sum(totrstoi),np.sum(totpstoi))

print(len(binu)+len(uninu)+count, " total reactions")
print(count, " reactions are not unispecies or bispecies")
unistois,unicounts=np.unique(np.array(uninu)[:,0],axis=0,return_counts=True)
bistois,bicounts=np.unique(np.array(binu)[:,0],axis=0,return_counts=True)
print("Bispecies constraint violating")
viol=0
for max in np.argsort(bicounts):
    pos=np.where(np.all(np.array(binu)[:,0]==bistois[max],axis=1))[0]
    if(len(pos)>1):
        uniq,counts=np.unique(np.array(binu)[pos,1], axis=0, return_counts=True)
        if(len(uniq)>1):
            print([species[i] for i in np.array(bistois[max],dtype=int)], bicounts[max])
            viol=viol+bicounts[max]
            print(uniq,counts)
print(viol," constraint violating reactions")

print("Unispecies constraint violating")
viol=0
for max in np.argsort(unicounts):
    pos=np.where(np.all(np.array(uninu)[:,0]==unistois[max],axis=1))[0]
    if(len(pos)>1):
        uniq,counts=np.unique(np.array(uninu)[pos,1], axis=0, return_counts=True)
        if(len(uniq)>1):
            print([species[i] for i in np.array(unistois[max],dtype=int)], unicounts[max])
            viol=viol+unicounts[max]
            print(uniq,counts)
print(viol," constraint violating reactions")

stop=timeit.default_timer()
print ('runtime: %f'%(stop-start))
sys.stdout.flush()
