#!/bin/bash
Natoms0=3
Natoms1=20
filebase=data/runtimes2
for i in `seq $Natoms0 $Natoms1`; do
  echo $i
  ./ratematrix.py --atoms $i $i $i --filebase ${filebase}_${i} --plot 0 --calculate 0 &
done
wait

for i in `seq $Natoms0 $Natoms1`; do
cat ${filebase}_${i}out.dat >> ${filebase}out.dat
rm ${filebase}_${i}out.dat
done
