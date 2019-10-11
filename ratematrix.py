#!/usr/bin/env python
import os
import numpy as np
import cantera as ct
import timeit
import argparse
import resource
from scipy.sparse import coo_matrix, save_npz, load_npz
from scipy.sparse.linalg import eigs, svds
from scipy.special import factorial, binom
from scipy.integrate import ode, solve_ivp
import sys
import rlist
import glob

#Command-line arguments
parser = argparse.ArgumentParser(description='Generate a sparse rate matrix from cantera model.')
parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for npy file output. Three files will be created for each reaction, storing rates, row indices, and column indices.')
parser.add_argument("--mechanism", type=str, required=False, default='mechanisms/h2o2.cti', dest='mechanism', help='Mechanism cti file. Default mechanisms/h2o2.cti.')
parser.add_argument("--calculate", nargs=2, type=int, required=False, default=[0,-1], help='Rows to calculate in rate matrix. Takes starting index and ending index and save rate matrix corresponding to all reactions over these rows. Default [0 -1].')
parser.add_argument("--eigenvalues", type=int, required=False, default=50, help='Flag to calculate  eigenvalues. If 1, then print args.temperature, args.pressure, total atoms, dimension, runtime, recursive calls, recursive levels, and save rate matrix, eigenvalues, and eigenvectors, then quit. Default 1.')
parser.add_argument("--temperature", type=float, required=False, default=1000, help='Temperature in Kelvin. Default 1500.')
parser.add_argument("--adiabatic", type=int, choices=[0, 1, 2], required=False, default=1, help='Convert energy from reactions to heat. The values 0, 1, and 2 correspond to constant volume/temperature, constant volume/energy, and constant pressure/enthalpy, respectively. The temperature is specify the reference multiindix specified with --reference. ')
parser.add_argument("--reference", type=int, nargs='+', required=False, default=[0, 4, 3, 2, 8, 0], help='Reference multiindex for which the temperature and number of atoms are specified. Default 0 4 3 2 8 0.')
parser.add_argument("--pressure", type=float, required=False, default=1, help='Pressure in atm. Default 1.')
parser.add_argument("--fix", nargs='+', type=int, required=False, default=[], help='Fix species numbers for parallelization. Include each species index followed by the number of molecules to fix.')
parser.add_argument("--accumulate", type=int, required=False, default=0, choices=[0,1], help='Flag to accumulate the multiindices from parallel runs.')
parser.add_argument("--propogate", type=int, required=False, default=0, choices=[0,1], help='Flag to propogate reference multiindex.')
parser.add_argument("--t0", type=float, required=False, default=1e-8, help='Initial integration time for propogating.')
parser.add_argument("--tmax", type=float, required=False, default=1, help='Final integration time for propogating.')
parser.add_argument("--Nt", type=int, required=False, default=100, help='Number of times to propogate.')
parser.add_argument("--print", type=int, required=False, default=0, choices=[0,1], help='Print runtimes.')
args = parser.parse_args()

#Functions for relating multiindices to matrix indices
def get_multiindex(index):
    return multiindices[index]
def get_index(multiindex):
    return np.where(np.all(multiindices==multiindex,axis=1))[0][0]

def get_rate (multiindex, stoi, k, vol):

    if np.all(multiindex>=stoi):
        return k*np.product(binom(multiindex, stoi)*factorial(stoi))/(ct.avogadro*vol)**(np.sum(stoi)-1)
        # return k*np.product(factorial(multiindex)/factorial(multiindex-stoi))/(ct.avogadro*vol)**(np.sum(stoi)-1)
    else:
        return 0.

#We can parallelize this by calculating for each row, which should go relatively fast.
def calculate_sparse_elements_row(rind,i):
    reaction=gas.reactions()[rind]
    rstoi=np.array([reaction.reactants[x] if x in reaction.reactants.keys() else 0 for x in species])
    pstoi=np.array([reaction.products[x] if x in reaction.products.keys() else 0 for x in species])
    data=[]
    rows=[]
    columns=[]
    multiindex=get_multiindex(i)
    if(args.adiabatic == 1):
        gas.UVX=refenergy/refmass,refvol/refmass,multiindex
    elif(args.adiabatic == 2):
        gas.HPX=(refenth+np.dot(multiindex-refmultiindex,refpotentials)/(1000*ct.avogadro))/refmass,args.pressure*ct.one_atm,multiindex
    else:
        gas.TPX=args.temperature,args.pressure*ct.one_atm,multiindex

    #forward reaction
    k=gas.forward_rate_constants[rind]
    multiindex2=multiindex-rstoi+pstoi
    if np.all(multiindex2>=0.) and not np.isnan(k) and (np.any([np.all(multiindex2==multiindex) for multiindex in multiindices])):
        rate=get_rate(multiindex,rstoi,k,refvol)
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
    if np.all(multiindex2>=0) and not np.isnan(k) and (np.any([np.all(multiindex2==multiindex) for multiindex in multiindices])):
        rate=get_rate(multiindex,pstoi,k,refvol)
        j=get_index(multiindex2)
        data.append(rate)
        rows.append(i)
        columns.append(j)
        data.append(-rate)
        rows.append(i)
        columns.append(i)
    return data,rows,columns

#Main loop over rows to enumerate sparse data
def calculate_sparse_elements(rind, startrow, endrow):
    data=[]
    rows=[]
    columns=[]
    for  i in range(startrow,endrow):
        rowdata, rowrows, rowcolumns = calculate_sparse_elements_row(rind,i)
        data+=rowdata
        rows+=rowrows
        columns+=rowcolumns
    return data,rows,columns

#Main
filebase=args.filebase
mechanism=args.mechanism
gas=ct.Solution(mechanism)
gas.TP=args.temperature,args.pressure*ct.one_atm
ns=gas.n_species
nr=gas.n_reactions
species=gas.species_names
elements=gas.element_names
refmultiindex=np.zeros(ns,dtype=int)
for i in range(0,len(args.reference),2):
    refmultiindex[args.reference[i]]=args.reference[i+1]
sp_atoms=[]
for i in range(ns):
    sp_atoms.append(np.array([int(gas.species()[i].composition[el] if el in gas.species()[i].composition.keys() else 0) for el in elements]))
sp_atoms=np.array(sp_atoms)
atoms=np.zeros(len(elements),dtype=int)
for i in range(len(refmultiindex)):
    atoms += refmultiindex[i]*sp_atoms[i]
fixed=np.array(args.fix)

atoms=np.array(atoms)
gas.TPX=args.temperature,args.pressure*ct.one_atm,refmultiindex
refquant=ct.Quantity(gas,moles=np.sum(refmultiindex)/ct.avogadro)
refenth=refquant.enthalpy
refentr=refquant.entropy
refenergy=refquant.int_energy
refmass=refquant.mass
refvol=refquant.volume
refpotentials=refquant.chemical_potentials

#Calculate the space of possible states
if args.accumulate==0:
    start=timeit.default_timer()

    multiindex=[]
    last_avail=[[],[]]
    for i in range(ns):
        multiindex.append(0)
        last_avail[0].append(i)
        last_avail[1].append(np.array([int(gas.species()[i].composition[el] if el in gas.species()[i].composition.keys() else 0) for el in elements]))

    remove_atoms = np.zeros(len(atoms));
    for i in range(0,len(fixed),2):
        remove_atoms += fixed[i+1]*sp_atoms[fixed[i]]
    frvecs=np.transpose(gas.reactant_stoich_coeffs()-gas.product_stoich_coeffs())
    rvecs=np.zeros((2*nr,ns),dtype=int)
    for ind in range(nr):
        rvecs[2*ind]=frvecs[ind]
        rvecs[2*ind+1]=-frvecs[ind]
    multiindices=rlist.list(refmultiindex, rvecs)
    sorted=np.argsort(np.sum((ns**np.arange(ns)*multiindices),axis=1))
    multiindices=multiindices[sorted]


    for i in range(0,len(fixed),2):
        multiindices[:,fixed[i]]=fixed[i+1]

    accessible=[]
    temperatures=[]
    pressures=[]
    inaccessible=[]
    gas=ct.Solution(mechanism)

    if(args.adiabatic == 1):
        for multiindex in multiindices:
            try:
                gas.UVX=refenergy/refmass,refvol/refmass,multiindex
                quant=ct.Quantity(gas, moles=np.sum(multiindex)/ct.avogadro)
                if(quant.T > 100): #far outside the range where these constants are valid...
                    temperatures.append(quant.T)
                    pressures.append(quant.P/ct.one_atm)
                    accessible.append(multiindex)
                else:
                    inaccessible.append(multiindex)
            except:
                inaccessible.append(multiindex)
    elif(args.adiabatic == 2):
        for multiindex in multiindices:
            try:
                gas.HPX=(refenth+np.dot(multiindex-refmultiindex,refpotentials)/(1000*ct.avogadro))/refmass,args.pressure*ct.one_atm,multiindex
                quant=ct.Quantity(gas, moles=np.sum(multiindex)/ct.avogadro)
                if(quant.T > 100): #far outside the range where these constants are valid...
                    temperatures.append(quant.T)
                    pressures.append(quant.P/ct.one_atm)
                    accessible.append(multiindex)
                else:
                    inaccessible.append(multiindex)
            except:
                inaccessible.append(multiindex)
    else:
        accessible=multiindices
        temperatures=np.zeros(len(multiindices))+args.temperature
        pressures=np.zeros(len(multiindices))+args.pressure

    multiindices=np.array(accessible)
    dim=len(multiindices)

    atot=np.sum(atoms)
    runtime=timeit.default_timer()-start


    np.save(filebase+"multiindices.npy",multiindices)
    np.save(filebase+"temperatures.npy",temperatures)
    np.save(filebase+"pressures.npy",pressures)
    out=open(filebase+"out.dat","w")

    print(dim, ns, runtime, file=out)
    print(*refmultiindex, file=out)
    print(*elements, file=out)
    out.close()
    if(args.print == 1):
        print("state space dimension: ", dim)
        print("state space runtime: ", runtime)
else:
    if os.path.isfile(filebase+"multiindices.npy"):
        multiindices=np.load(filebase+"multiindices.npy")
    else:
        mfiles=sorted(glob.glob('%s*multiindices.npy'%filebase))
        tfiles=sorted(glob.glob('%s*temperatures.npy'%filebase))
        pfiles=sorted(glob.glob('%s*pressures.npy'%filebase))

        multiindices=[]
        temperatures=[]
        pressures=[]
        for file in mfiles:
            multiindices+=np.load(file).tolist()
        for file in tfiles:
            temperatures+=np.load(file).tolist()
        for file in pfiles:
            pressures+=np.load(file).tolist()
        multiindices=np.array(multiindices)
        temperatures=np.array(temperatures)
        pressures=np.array(pressures)
        np.save(filebase+"multiindices.npy",multiindices)
        np.save(filebase+"temperatures.npy",temperatures)
        np.save(filebase+"pressures.npy",pressures)

    file=open(filebase+"out.dat","r")
    instrings=file.readline().split()
    atot=int(instrings[0])
    dim=int(instrings[1])
    runtime=float(instrings[2])
    count=int(instrings[3])
    level=int(instrings[4])
    instrings=file.readline().split()
    refmultiindex=np.array([int(float(i)) for i in instrings])

startrow, endrow=args.calculate
if endrow==-1 or endrow > dim:
    endrow=dim
if endrow>startrow:
    #Loop through each reaction index and calculate spase elements
    start=timeit.default_timer()

    gas=ct.Solution(mechanism)
    gas.TPX=args.temperature,args.pressure*ct.one_atm,refmultiindex
    refquant=ct.Quantity(gas,moles=np.sum(refmultiindex)/ct.avogadro)
    refenth=refquant.enthalpy
    refmass=refquant.mass
    refvol=refquant.volume
    quant=ct.Quantity(gas, moles=np.sum(refmultiindex)/ct.avogadro)

    data=[]
    rows=[]
    columns=[]
    for rind in range(nr):
        reac_data,reac_rows,reac_columns=calculate_sparse_elements(rind, startrow, endrow)
        data+=reac_data
        rows+=reac_rows
        columns+=reac_columns
    np.save(filebase+"spatoms.npy",sp_atoms)
    if not os.path.exists(filebase+"rows"):
        os.mkdir(filebase+"rows")
    if not os.path.exists(filebase+"columns"):
        os.mkdir(filebase+"columns")
    if not os.path.exists(filebase+"data"):
        os.mkdir(filebase+"data")
    np.save(filebase+"rows/%i_%i.npy"%(startrow,endrow),rows)
    np.save(filebase+"columns/%i_%i.npy"%(startrow,endrow),columns)
    np.save(filebase+"data/%i_%i.npy"%(startrow,endrow),data)
    runtime=timeit.default_timer()-start
    out=open(filebase+"_%icout.dat"%(start),"w")
    print("sparsity ", len(rows)/(dim*dim), runtime, file=out)
    out.close()
    if(args.print == 1):
        print("calculate runtime:", runtime)

#Calculate eigenvalues
if args.eigenvalues==-1:
    args.eigenvalues=dim
if args.eigenvalues>0:
    #accumulate rows, data, and columns
    start=timeit.default_timer()

    if args.accumulate==1:
        if  os.path.isfile(args.filebase+"rows.npy") and os.path.isfile(args.filebase+"columns.npy") and os.path.isfile(args.filebase+"data.npy"):
            rows=np.load(args.filebase+"rows.npy")
            columns=np.load(args.filebase+"columns.npy")
            data=np.load(args.filebase+"data.npy")
        else:
            files=glob.glob(args.filebase+"rows/*.npy")
            rows=[]
            for file in files:
                row=np.load(file).tolist()
                rows+=row
            files=glob.glob(args.filebase+"columns/*.npy")
            columns=[]
            for file in files:
                column=np.load(file).tolist()
                columns+=column
            files=glob.glob(args.filebase+"data/*.npy")
            data=[]
            for file in files:
                dat=np.load(file).tolist()
                data+=dat
    np.save(filebase+"rows.npy",rows)
    np.save(filebase+"columns.npy",columns)
    np.save(filebase+"data.npy",data)

    ratematrix=coo_matrix((np.array(data),(np.array(columns),np.array(rows))),(int(dim),int(dim)))

    if args.eigenvalues < ratematrix.shape[0]:
        eigenvalues,eigenvectors=eigs(ratematrix, args.eigenvalues, sigma=-1e-1, which='LM')
        # svals=svds(ratematrix, args.eigenvalues, which='SM', return_singular_vectors=False)
        svals=np.array(dim)
    else:
        eigenvalues,eigenvectors=np.linalg.eig(ratematrix.toarray())
        svals=np.linalg.svd(ratematrix.toarray(),compute_uv=False)
    sorted=np.argsort(eigenvalues)
    np.save(filebase+"eigenvalues.npy",eigenvalues.astype(complex)[sorted])
    np.save(filebase+"singularvalues.npy",svals)
    np.save(filebase+"eigenvectors.npy",eigenvectors.astype(complex)[:,sorted])

    runtime=timeit.default_timer()-start
    out=open(filebase+"eout.dat","w")
    print(runtime, file=out)
    print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, file=out)
    print(dim, file=out)
    out.close()
    if(args.print == 1):
        print("eigenvalues runtime:", runtime)

if args.propogate == 1:
    start=timeit.default_timer()

    if os.path.isfile(args.filebase+"rows.npy") and os.path.isfile(args.filebase+"columns.npy") and os.path.isfile(args.filebase+"data.npy"):
        rows=np.load(args.filebase+"rows.npy")
        columns=np.load(args.filebase+"columns.npy")
        data=np.load(args.filebase+"data.npy")
    else:
        files=glob.glob(args.filebase+"rows/*.npy")
        rows=[]
        for file in files:
            row=np.load(file).tolist()
            rows+=row
        files=glob.glob(args.filebase+"columns/*.npy")
        columns=[]
        for file in files:
            column=np.load(file).tolist()
            columns+=column
        files=glob.glob(args.filebase+"data/*.npy")
        data=[]
        for file in files:
            dat=np.load(file).tolist()
            data+=dat
        np.save(filebase+"rows.npy",rows)
        np.save(filebase+"columns.npy",columns)
        np.save(filebase+"data.npy",data)

    ratematrix=coo_matrix((np.array(data),(np.array(columns),np.array(rows))),(int(dim),int(dim)))

    def func(t,y):
        return ratematrix.dot(y)
    #
    # jacarray=ratematrix.toarray()
    # def jac(t, y):
    #     return jacarray

    y0=np.zeros(dim)
    y0[get_index(refmultiindex)]=1.0
    times=[args.t0*(args.tmax/args.t0)**(n*1.0/args.Nt) for n in range(args.Nt)]
    sol=solve_ivp(func,[0,args.tmax], y0, method='BDF', t_eval=times, atol=1e-8, rtol=1e-6, first_step=args.t0/100, jac=ratematrix)
    # sol=solve_ivp(func,[0,args.tmax], y0, method='LSODA', t_eval=times, atol=1e-8, rtol=1e-6, first_step=args.t0/100, jac=jac)
    vals=np.transpose(np.array(sol.y)).tolist()

    np.save(filebase+"times.npy",np.array(times))
    np.save(filebase+"propogate.npy",vals)

    runtime=timeit.default_timer()-start
    out=open(filebase+"rout.dat","w")
    print(runtime, file=out)
    print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, file=out)
    print(dim, file=out)
    out.close()
    if(args.print == 1):
        print("propogate success:", sol.success)
        print("propogate runtime:", runtime)
if(args.print == 1):
    print("memory:", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
