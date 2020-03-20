#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np
import cantera as ct
import timeit
import argparse
from scipy.sparse import coo_matrix
# import progressbar



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

    states = ct.SolutionArray(gas, extra=['t'])

    for t in times:
        sim.advance(t)

        data=[]
        rows=[]
        columns=[]
        mat=[]
        eigenvalues=[]
        eigenvectors=[]

        states.append(r.thermo.state, t=t)

    return states

parser = argparse.ArgumentParser(description='Generate a sparse rate matrix from cantera model.')
parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for npy file output. Three files will be created for each reaction, storing rates, row indices, and column indices.')
parser.add_argument("--temperature", type=float, required=False, default=1000, help='Temperature in Kelvin. Default 1000.')
parser.add_argument("--adiabatic", type=int, choices=[0, 1, 2], required=False, default=0, help='Convert energy from reactions to heat. The values 0, 1, and 2 correspond to constant pressure/temperature, constant volume/energy, and constant pressure/enthalpy, respectively. The temperature is specify the reference multiindix specified with --reference. ')
parser.add_argument("--reference", type=int, nargs='+', required=False, default=[0, 2, 3, 1], help='Reference multiindex for which the temperature and number of atoms are specified. Default 0 20 3 10 4 1 8 200.')
parser.add_argument("--pressure", type=float, required=False, default=1, help='Pressure in atm. Default 1.')
parser.add_argument("--t0", type=float, required=False, default=1e-8, help='Initial integration time for propogating.')
parser.add_argument("--tmax", type=float, required=False, default=1, help='Final integration time for propogating.')
parser.add_argument("--Npoints", type=int, required=False, default=100, help='Number of times to propogate.')
parser.add_argument("--mechanism", type=str, required=False,default='h2o2.cti', help='Cantera mechanism file.')
parser.add_argument("--seed", type=int, required=False, default=0, help='Random seed')
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
np.random.seed(args.seed)
Ar=np.zeros((ns,ns))
Ap=np.zeros((ns,ns))


rstois=[]
pstois=[]
unifwd=[]
bifwd=[]
unirev=[]
birev=[]
rmfwd=[]
rmrev=[]

uninu=[]
binu=[]
count=0
totrstoi=np.zeros(ns,dtype=int)
totpstoi=np.zeros(ns,dtype=int)

reversible=0
for rind in range(nr):
    if gas.is_reversible(rind):
        reversible+=1
print(ns,"species",nr,"reactions,",reversible,"of which are reversible.")

# pbar=progressbar.ProgressBar(widgets=['Reactions: ', progressbar.Percentage(), progressbar.Bar(), ' ', progressbar.ETA()], maxval=nr)
# pbar.start()
#Find list of reactions that must be removed for constraints
for rind in range(nr):
    # pbar.update(nr)
    reaction=gas.reactions()[rind]
    rstoi=np.array([reaction.reactants[x] if x in reaction.reactants.keys() else 0 for x in species],dtype=int)
    pstoi=np.array([reaction.products[x] if x in reaction.products.keys() else 0 for x in species],dtype=int)
    rstois.append(rstoi)
    pstois.append(pstoi)

    if(gas.reaction_type(rind) == 1):
        rinds=np.where(rstoi>0)[0]
        pinds=np.where(pstoi>0)[0]
        #forward reactions
        if(len(rinds)==1):
            unifwd.append(rind)
        elif(len(rinds)==2):
            bifwd.append(rind)
        else:
            rmfwd.append(rind)
            count=count+1
        #reverse reactions
        if reaction.reversible:
            if(len(pinds)==1):
                unirev.append(rind)
            elif(len(pinds)==2):
                birev.append(rind)
            else:
                rmrev.append(rind)
                count=count+1
    else:
        rmfwd.append(rind)
        count=count+1
        if reaction.reversible:
            rmrev.append(rind)
            count=count+1

#Identify reactions that violate constraints
#Unispecies
unis=np.concatenate((unifwd,unirev))
stois=np.concatenate((np.array(rstois)[unifwd],np.array(pstois)[unirev]))
uniindices=np.random.permutation(np.arange(len(unifwd)+len(unirev)))
uniremoves=np.zeros(len(uniindices),dtype=int)
for ind in range(len(uniindices)):
    index=uniindices[ind]
    unispecies=np.array(np.where(stois[index] > 0)[0])
    for ind2 in range(ind):
        index2 = uniindices[ind2]
        unispecies2=np.array(np.where(stois[index2] > 0)[0])
        #remove if they involve the same species but with different stoichiometric coefficients
        if(np.all(unispecies==unispecies2)) and not np.all(stois[index] == stois[index2]):
            uniremoves[ind]=1
#Bispecies
bis=np.concatenate((bifwd,birev))
stois=np.concatenate((np.array(rstois)[bifwd],np.array(pstois)[birev]))
biindices=np.random.permutation(np.arange(len(bifwd)+len(birev)))
biremoves=np.zeros(len(biindices),dtype=int)
for ind in range(len(biindices)):
    index=biindices[ind]
    bispecies=np.array(np.where(stois[index] > 0)[0])
    for ind2 in range(ind):
        index2 = biindices[ind2]
        bispecies2=np.array(np.where(stois[index2] > 0)[0])
        #remove if they involve the same species but with different stoichiometric coefficients
        if(np.all(bispecies==bispecies2)) and not np.all(stois[index] == stois[index2]):
            biremoves[ind]=1
rm=np.unique(np.concatenate((rmfwd,rmrev,unis[uniindices[np.where(uniremoves>0)[0]]],bis[biindices[np.where(biremoves>0)[0]]])))
print(len(rm)*1.0/nr, " percent of reactions removed for constraints.")

#Remove those reactions and find the Ar and Ap matrices for the remainder
Sr=np.zeros((ns,ns,ns), dtype=int)
Sp=np.zeros((ns,ns,ns), dtype=int)
Ar=np.zeros((ns,ns))
Ap=np.zeros((ns,ns))

productsets=[[] for i in range(ns)]
for ind in range(nr):
    if ind not in rm:
        for prod in np.where(pstois[ind]>0)[0]:
            productsets[prod].append(np.where(rstois[ind]>0)[0])
        if gas.is_reversible(ind):
            for prod in np.where(rstois[ind]>0)[0]:
                productsets[prod].append(np.where(pstois[ind]>0)[0])

count=0
for i in range(ns):
    if(len(productsets[i]) > 0):
        if (len(productsets[i]) > len(np.unique(np.concatenate(productsets[i])))):
            count=count+1
print("Cannot generate the list L for", count, "of", ns, "species.")

#Calculate matrices
for rind in np.setdiff1d(np.arange(nr), rm):
    reaction=gas.reactions()[rind]
    rstoi=np.array([reaction.reactants[x] if x in reaction.reactants.keys() else 0 for x in species],dtype=int)
    pstoi=np.array([reaction.products[x] if x in reaction.products.keys() else 0 for x in species],dtype=int)
    rstois.append(rstoi)
    pstois.append(pstoi)
    rinds=np.where(rstoi>0)[0]
    pinds=np.where(pstoi>0)[0]

    #forward reactions
    k=gas.forward_rate_constants[rind]
    if(len(rinds)==1):
        #reactants
        Ar[rinds[0],rinds[0]]=k
        if not np.all(Sr[rinds[0],rinds[0]] == np.zeros(ns)) and not np.all(Sr[rinds[0],rinds[0]] == rstoi):
            print("error!")
        Sr[rinds[0],rinds[0]]=rstoi
        #products
        # for pind in pinds:
        #     Ap[pind, rinds[0]] = k
        #     if not np.all(Sp[pind, rinds[0]] == np.zeros(ns)) and not np.all(Sp[pind, rinds[0]] == rstoi):
        #         print("error!")
        #     Sp[pind, rinds[0]] = rstoi
    elif(len(rinds)==2):
        #reactants
        Ar[rinds[0],rinds[1]]=k
        Ar[rinds[1],rinds[0]]=k
        if not np.all(Sr[rinds[0],rinds[1]] == np.zeros(ns)) and not np.all(Sr[rinds[0],rinds[1]] == rstoi):
            print("error!")
        Sr[rinds[0],rinds[1]]=rstoi
        if not np.all(Sr[rinds[1],rinds[0]] == np.zeros(ns)) and not np.all(Sr[rinds[1],rinds[0]] == rstoi):
            print("error!")
        Sr[rinds[1],rinds[0]]=rstoi
        #products
        # for pind in pinds:
        #     rind=np.random.choice(rinds)
        #     Ap[pind, rind] = k
        #     if not np.all(Sp[pind, rind] == np.zeros(ns)) and not np.all(Sp[pind, rind] == pstoi):
        #         print("error reverse products!")
        #     Sp[pind, rind] = pstoi
    else:
        print("error!")

    #reverse reactions
    if gas.is_reversible(rind):
        k=gas.reverse_rate_constants[rind]
        if(len(pinds)==1):
            #reactants
            Ar[pinds[0],pinds[0]]=k
            if not np.all(Sr[pinds[0],pinds[0]] == np.zeros(ns)) and not np.all(Sr[pinds[0],pinds[0]] == pstoi):
                print("error!")
            Sr[pinds[0],pinds[0]]=pstoi
            #products
            # for rind in rinds:
            #     Ap[rind, pinds[0]] = k
            #     if not np.all(Sp[rind, pinds[0]] == np.zeros(ns)) and not np.all(Sp[rind, pinds[0]] == pstoi):
            #         print("error!")
            #     Sp[rind, pinds[0]] = pstoi
        elif(len(pinds)==2):
            #reactants
            Ar[pinds[0],pinds[1]]=k
            Ar[pinds[1],pinds[0]]=k
            if not np.all(Sr[pinds[0],pinds[1]] == np.zeros(ns)) and not np.all(Sr[pinds[0],pinds[1]] == pstoi):
                print("error!")
            Sr[pinds[0],pinds[1]]=pstoi
            if not np.all(Sr[pinds[1],pinds[0]] == np.zeros(ns)) and not np.all(Sr[pinds[1],pinds[0]] == pstoi):
                print("error!")
            Sr[pinds[1],pinds[0]]=pstoi
            #products
            # for rind in rinds:
            #     pind=np.random.choice(pinds)
            #     Ap[rind, pind] = k
            #     if not np.all(Sp[rind, pind] == np.zeros(ns)) and not np.all(Sp[rind, pind] == pstoi):
            #         print("error reverse products!")
            #     Sp[rind, pind] = pstoi
        else:
            print("error!")

np.save(args.filebase+"Ar", Ar)


for remove in rm:
    if(gas.reaction_type(remove) == 4):
        rate=gas.reactions()[remove]
        rate.high_rate = ct.Arrhenius(0, 0, 0)
        rate.low_rate = ct.Arrhenius(0, 0, 0)
        gas.modify_reaction(remove,rate)
    elif(gas.reaction_type(remove) == 5):
        rate=gas.reactions()[remove]
        for rt in rate.rates:
            rt = (rt[0],ct.Arrhenius(0, 0, 0))
        gas.modify_reaction(remove,rate)
    else:
        rate=gas.reactions()[remove]
        rate.rate = ct.Arrhenius(0, 0, 0)
        gas.modify_reaction(remove,rate)

# observation=runsim(times)
# np.save(args.filebase+"times.npy", observation.t)
# np.save(args.filebase+"concentrations.npy", observation.X)
# np.save(args.filebase+"temperatures.npy", observation.T)
# np.save(args.filebase+"pressures.npy", observation.P/ct.one_atm)

stop=timeit.default_timer()
print ('runtime: %f'%(stop-start))
sys.stdout.flush()
