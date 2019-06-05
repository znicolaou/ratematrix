#!/bin/bash
N_temp=25
N_press=25
press0=0.5
press1=15
temp0=1000
temp1=2500
atoms="12 24 5"
filebase0=data/asweep
save=0
export OMP_NUM_THREADS=1

rm -r ${filebase0}*
# ./ratematrix.py --atoms $atoms --filebase ${filebase0} --plot 0 --calculate 0

mkdir -p $filebase0
for i in `seq 0 $N_temp`; do
  for j in `seq 0 $N_press`; do
    temp=`bc -l <<<"${temp0}+(${temp1}-${temp0})*${i}*1.0/${N_temp}"`
    press=`bc -l <<<"${press0}+(${press1}-${press0})*${j}*1.0/${N_press}"`
    echo $temp $press
    filebase=$filebase0/${i}_${j}
    # cp ${filebase0}multiindices.npy ${filebase}multiindices.npy
    # cp ${filebase0}out.dat ${filebase}out.dat
    ./ratematrix.py --temperature $temp --pressure $press --atoms $atoms --filebase ${filebase} --plot 0 --save 1 --calculate 1 --fix 8 5 --adiabatic 1 --reference 12 0 0 6 0 0 0 0 5 &
  done
  wait
done

for i in `seq 0 $N_temp`; do
  for j in `seq 0 $N_press`; do
    filebase=$filebase0/${i}_${j}
    # rm ${filebase}multiindices.npy
    tail -n 1 ${filebase}out.dat >> ${filebase0}out.dat
    # rm ${filebase}out.dat
  done
done
