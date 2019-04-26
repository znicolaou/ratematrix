# ratematrix
Running `./ratematrix.py -h` produces the following usage message:
```
usage: ratematrix.py [-h] --filebase FILEBASE [--accumulate {0,1}]
                     [--mechanism MECHANISM] [--rind RIND] [--Nmax NMAX]
                     [--Nvals NVALS]

Generate a sparse rate matrix from cantera model.

optional arguments:
  -h, --help            show this help message and exit
  --filebase FILEBASE   Base string for file output
  --accumulate {0,1}    If 1, search filebase directory for npy files,
                        generate sparse matrix, and plot eigenvalues
  --mechanism MECHANISM
                        Mechanism cti file
  --rind RIND           Reaction index
  --Nmax NMAX           Maximum number of molecules
  --Nvals NVALS         Number of eigenvalues
  ```
