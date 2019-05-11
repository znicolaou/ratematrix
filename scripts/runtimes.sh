#!/bin/bash
Natoms0=3
Natoms1=9
filebase=data/runtimes
for i in `seq $Natoms0 $Natoms1`; do
  echo $i
  echo "./ratematrix.py --atoms $i $i $i --filebase $filebase --plot 0 --calculate 0 "
  ./ratematrix.py --atoms $i $i $i --filebase $filebase --plot 0 --calculate 0
done
