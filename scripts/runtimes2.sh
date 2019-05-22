#!/bin/bash
Natoms0=3
Natoms1=12
filebase=data/runtimes2
for i in `seq $Natoms0 $Natoms1`; do
  echo $i
  echo "./ratematrix.py --atoms $i $i $i --filebase $filebase --plot 0 --calculate 0 "
  ./ratematrix.py --atoms $i $i $i --filebase $filebase --plot 0 --calculate 0
done
