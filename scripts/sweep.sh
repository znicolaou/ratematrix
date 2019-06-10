#!/bin/bash
N_temp=50
N_press=50
press0=0.01
press1=100
temp0=600
temp1=1000
atoms="4 8 5"
procs=100
filebase0=data/sweep
export OMP_NUM_THREADS=1

rm -r ${filebase0}*
./ratematrix.py --atoms $atoms --filebase ${filebase0} --plot 0 --calculate 0

mkdir -p $filebase0
for i in `seq 0 $N_temp`; do
  for j in `seq 0 $N_press`; do
    js=`jobs | wc -l`
    while [ $js -ge $procs ]; do
      js=`jobs | wc -l`
      sleep 1
    done
    temp=`bc -l <<<"${temp0}+(${temp1}-${temp0})*${j}*1.0/${N_temp}"`
    press=`bc -l <<<"${press0}*e(${i}*1.0/${N_press}*l(${press1}/${press0}))"`
    echo $temp $press
    filebase=$filebase0/${i}_${j}
    cp ${filebase0}multiindices.npy ${filebase}multiindices.npy
    cp ${filebase0}out.dat ${filebase}out.dat
    ./ratematrix.py --temperature $temp --pressure $press --atoms $atoms --filebase ${filebase} --accumulate 1 --plot 0 --save 1 --calculate 1 --fix 8 5 &
  done
done
wait

rm ${filebase0}out.dat

for i in `seq 0 $N_temp`; do
  for j in `seq 0 $N_press`; do
    filebase=$filebase0/${i}_${j}
    rm ${filebase}multiindices.npy
    tail -n 1 ${filebase}out.dat >> ${filebase0}out.dat
    rm ${filebase}out.dat
  done
done
rm -r $filebase0
