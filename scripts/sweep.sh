#!/bin/bash
N_temp=25
N_press=25
press0=0.5
press1=10
temp0=1000
temp1=2500
atoms="5 10 5"
filebase0=data/sweep
save=0

rm -r ${filebase0}*
./ratematrix.py --atoms $atoms --filebase ${filebase0} --plot 0 --calculate 0

mkdir -p $filebase0
for i in `seq 0 $N_temp`; do
  for j in `seq 0 $N_press`; do
    temp=`bc -l <<<"${temp0}+(${temp1}-${temp0})*${i}*1.0/${N_temp}"`
    press=`bc -l <<<"${press0}+(${press1}-${press0})*${j}*1.0/${N_press}"`
    echo $temp $press
    filebase=$filebase0/${i}_${j}
    cp ${filebase0}multiindices.npy ${filebase}multiindices.npy
    cp ${filebase0}out.dat ${filebase}out.dat
    ./ratematrix.py --temperature $temp --pressure $press --atoms $atoms --filebase ${filebase} --accumulate 1 --plot 0 --save 0 --calculate 1 &
  done
  wait
done

for i in `seq 0 $N_temp`; do
  for j in `seq 0 $N_press`; do
    filebase=$filebase0/${i}_${j}
    rm ${filebase}multiindices.npy
    tail -n 1 ${filebase}out.dat >> ${filebase0}out.dat
    rm ${filebase}out.dat
  done
done
