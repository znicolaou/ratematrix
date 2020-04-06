# Files in the ratematrix repository
The file ratematrix.py contains Python code to calculate the sparse matrix elements of the transition rate matrix for a given Cantera .cti model and find eigenvalues and eigenvectors. The plot.nb file is a Mathematica notebook for plotting results. The folder mechanisms contains a few .cti models.

# System requirements
The python code has been run with anaconda, which can be downloaded here: https://www.anaconda.com/distribution/. The script requires packes numpy, scipy, cantera, and matplotlib, which can be installed in a new environment after installing anaconda with the shell command  
`conda create -n cantera_env -c cantera numpy scipy cantera matplotlib`  
Next, activate the environment with  
`source activate cantera_env`  

# Usage
Running `./ratematrix.py -h` produces the following usage message:
```
usage: ratematrix.py [-h] --filebase FILEBASE [--mechanism MECHANISM]
                     [--temperature TEMPERATURE] [--adiabatic {0,1,2}]
                     [--refspecies REFSPECIES [REFSPECIES ...]]
                     [--refcounts REFCOUNTS [REFCOUNTS ...]]
                     [--pressure PRESSURE] [--calculate {0,1}]
                     [--eigenvalues EIGENVALUES] [--propogate {0,1}]
                     [--thrs THRS] [--tau TAU] [--t0 T0] [--tmax TMAX]
                     [--Nt NT] [--print {0,1}] [--csv {0,1}]

Generate a sparse rate matrix from cantera model.

optional arguments:
  -h, --help            show this help message and exit
  --filebase FILEBASE   Base string for npy file output. Three files will be
                        created for each reaction, storing rates, row indices,
                        and column indices.
  --mechanism MECHANISM
                        Mechanism cti file. Default mechanisms/h2o2.cti.
  --temperature TEMPERATURE
                        Temperature in Kelvin. Default 1000.
  --adiabatic {0,1,2}   Convert energy from reactions to heat. The values 0,
                        1, and 2 correspond to constant volume/temperature,
                        constant volume/energy, and constant
                        pressure/enthalpy, respectively. The temperature is
                        specify the reference multiindix specified with
                        --reference.
  --refspecies REFSPECIES [REFSPECIES ...]
                        Reference multiindex for which the temperature and
                        number of atoms are specified. Default ['H2', 'O2',
                        'OH', 'AR'].
  --refcounts REFCOUNTS [REFCOUNTS ...]
                        Reference multiindex for which the temperature and
                        number of atoms are specified. Default [8, 4, 1, 80].
  --pressure PRESSURE   Pressure in atm. Default 1.
  --calculate {0,1}     Flag to calculate rate matrix. Default 1.
  --eigenvalues EIGENVALUES
                        Flag to calculate eigenvalues. If 1, then print
                        args.temperature, args.pressure, total atoms,
                        dimension, runtime, recursive calls, recursive levels,
                        and save rate matrix, eigenvalues, and eigenvectors,
                        then quit. Default 1.
  --propogate {0,1}     Flag to propogate reference multiindex.
  --thrs THRS           Threshold for including reactions.
  --tau TAU             Time scale for shift invert.
  --t0 T0               Initial integration time for propogating.
  --tmax TMAX           Final integration time for propogating.
  --Nt NT               Number of times to propogate.
  --print {0,1}         Print runtimes.
  --csv {0,1}           Save files to csv format.
  ```
  -----------
# Example
To find the sparse elements for the minimal mechanisms/h2o2.cti file starting from the reference state 8H2+4O2+80AR+1OH without calculate eigenvalues and saving output to data/test*, run
`./ratematrix.py --refspecies H2 O2 AR OH --refcounts 8 4 80 1 --eigenvalues 0 --filebase data/test`  
