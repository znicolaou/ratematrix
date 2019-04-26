# Files in the ratematrix repository
The file ratematrix.py contains python code to calculate the sparse matrix elements of the transition rate matrix for a given Cantera .cti model and find eigenvalues and eigenvectors. The folder mechanisms contains a few .cti models.

# System requirements
The python code has been run with anaconda, which can be downloaded here: https://www.anaconda.com/distribution/. The script requires packes numpy, scipy, cantera, and matplotlib, which can be installed in a new environment after installing anaconda with the shell command  
`conda create -n cantera_env -c cantera numpy scipy cantera matplotlib`  
To run the script, activate the environment with  
`source activate cantera_env` 

# Usage
Running `./ratematrix.py -h` produces the following usage message:
```
usage: ratematrix.py [-h] --filebase FILEBASE [--accumulate {0,1}]
                     [--mechanism MECHANISM] [--reaction RIND] [--Nmax NMAX]
                     [--Nvals NVALS] [--progress {0,1}]
                     [--temperature TEMPERATURE] [--pressure PRESSURE]

Generate a sparse rate matrix from cantera model.

optional arguments:
  -h, --help            show this help message and exit
  --filebase FILEBASE   Base string for npy file output. Three files will be
                        created for each reaction, storing rates, row indices,
                        and column indices.
  --accumulate {0,1}    If 1, search filebase directory for npy files,
                        generate sparse matrix, and plot eigenvalues. Default
                        0.
  --mechanism MECHANISM
                        Mechanism cti file. Default mechanisms/h2o2.cti.
  --reaction RIND       Reaction index, provided for parallelization. If none
                        is specified, the program will loop through all
                        reactions in the model in sequence. Default None.
  --Nmax NMAX           Maximum number of molecules for each species. Default
                        3.
  --Nvals NVALS         Number of eigenvalues to calculate, when --accumulate
                        1 is set. Default 1000
  --progress {0,1}      Print progress during calculation. Default 1.
  --temperature TEMPERATURE
                        Temperature. Default 1500K.
  --pressure PRESSURE   Pressure. Default 1atm
  ```
  -----------
# Examples
To find the sparse elements for the minimal mechanisms/h2o2.cti file with at most three molecules of each species on a single core, run  
`mkdir -p h2o2`  
`./ratematrix.py --filebase h2o2 --Nmax 3 --mechanism mechanisms/h2o2.cti`  
To calculate 100 eigenvalues and eigenvectors, and plot and store them, run  
`./ratematrix.py --accumulate 1 --filebase h2o2 --Nmax 3 --Nvals 100`  
To do this in parallel over the 28 reactions in the mechanism for 5 molecules per species, run  
`mkdir -p h2o2_2`  
`mkdir -p outs`  
`for i in {0..27}; do ./ratematrix.py --filebase h2o2_2/$i --reaction $i --progress 0 --Nmax 5 --mechanism mechanisms/h2o2.cti &> outs/$i.out & done`  
To plot and store these eigenvalues and eigenvectors, run   
`./ratematrix.py --accumulate 1 --filebase h2o2_2 --Nmax 5 --Nvals 100`  
This would take a very long time for this large matrix.
