#!/usr/bin/env python
import numpy as np
import argparse
from scipy.sparse import coo_matrix
from scipy.sparse import load_npz
from scipy.linalg import eig
from scipy.special import factorial
import sys
import glob
parser = argparse.ArgumentParser(description='Generate a sparse rate matrix from cantera model.')
parser.add_argument("--filebase", type=str, required=True, dest='filebase', help='Base string for npy file output.')
parser.add_argument("--dimension", type=int, required=True, dest='dim', help='Dimension.')
args = parser.parse_args()
print(args.dim)
ratematrix=coo_matrix((int(args.dim),int(args.dim)))
files=glob.glob(args.filebase+"ratematrix*.npz")
for file in files:
    mat=load_npz(file)
    ratematrix+=mat
np.save(args.filebase+"ratematrix.npy",ratematrix.toarray())
