# Files in the ratematrix repository
The file ratematrix.py contains Python code to calculate the sparse matrix elements of the transition rate matrix for a given Cantera .cti model and find eigenvalues and eigenvectors. The files rlistmodule.c, rlistmodule.h, and setup.py contain a C++ extension module for Python to recursively enumerate the space. The plot.nb file is a Mathematica notebook for plotting results. The folder mechanisms contains a few .cti models.

# System requirements
The python code has been run with anaconda, which can be downloaded here: https://www.anaconda.com/distribution/. The script requires packes numpy, scipy, cantera, and matplotlib, which can be installed in a new environment after installing anaconda with the shell command  
`conda create -n cantera_env -c cantera numpy scipy cantera matplotlib`  
To run the script, activate the environment with  
`source activate cantera_env`
Finally, compile and install the extension module with 
`python setup.py install` 

# Usage
Running `./ratematrix.py -h` produces the following usage message:
```
usage: ratematrix.py [-h] --filebase FILEBASE [--mechanism MECHANISM]
                     [--calculate {0,1}] [--Nvals NVALS] [--plot {0,1}]
                     [--save {0,1}] [--temperature TEMPERATURE]
                     [--pressure PRESSURE] [--atoms ATOMS [ATOMS ...]]

Generate a sparse rate matrix from cantera model.

optional arguments:
  -h, --help            show this help message and exit
  --filebase FILEBASE   Base string for npy file output. Three files will be
                        created for each reaction, storing rates, row indices,
                        and column indices.
  --mechanism MECHANISM
                        Mechanism cti file. Default mechanisms/h2o2.cti.
  --calculate {0,1}     Flag to calculate rate matrix and eigenvalues. If 1,
                        calculate space, then print args.temperature,
                        args.pressure, total atoms, dimension, runtime,
                        recursive calls, recursive levels, sparsity, and Nvals
                        largest eigenvalues, save rate matrix, eigenvalues,
                        and eigenvectors, then quit. If 0, only calculate
                        space, then print total atoms, dimension, runtime,
                        recursive calls, and recursive levels, then quit.
                        Default 1.
  --Nvals NVALS         Number of eigenvalues to print; 0 for all. Default 0.
  --plot {0,1}          Flag to plot eigenvalues and save
                        filebase/eigenvales.pdf. Default 1.
  --save {0,1}          Flag to save the results to filebaserate.npy,
                        filebaseeigenvalues.npy, and filebaseeigenvectors.npy.
                        Default 1.
  --temperature TEMPERATURE
                        Temperature in Kelvin. Default 1500.
  --pressure PRESSURE   Pressure in atm. Default 1.
  --atoms ATOMS [ATOMS ...]
                        Number of each atom, in order of their appearance in
                        the .cti file. If number of values is not number of
                        atoms, print the atoms. Default 3 3 3.
  ```
  -----------
# Examples
To find the sparse elements for the minimal mechanisms/h2o2.cti file with 7 oxygen atoms, 14 hydrogen atoms, and 5 argon atoms, with default temperature and pressure run  
`./ratematrix.py --filebase data/h2o2test --atoms 7 14 5`

