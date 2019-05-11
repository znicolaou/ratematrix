#!/bin/bash
N_temp=10
press=1
atoms0=3
atoms1=7
temp0=1000
temp1=2500
filebase0=data/sweep2
save=0

if [ -f ${filebase0}out.dat ]; then rm ${filebase0}out.dat; fi
for i in `seq 0 $N_temp`; do
  for j in `seq $atoms0 $atoms1`; do
    atoms="$j $j $j"
    temp=`bc -l <<<"${temp0}+(${temp1}-${temp0})*${i}*1.0/${N_temp}"`
    filebase=${filebase0}/${i}/${j}/
    mkdir -p $filebase
    echo $temp $press
    echo "./ratematrix.py --temperature $temp --pressure $press --atoms $atoms --filebase ${filebase} --plot 0 --save $save"
    ./ratematrix.py --temperature $temp --pressure $press --atoms $atoms --filebase ${filebase} --plot 0 --save $save &
  done
  wait

  for j in `seq $atoms0 $atoms1`; do
    filebase=${filebase0}/${i}/${j}/
    cat ${filebase}out.dat >> ${filebase0}out.dat
    rm ${filebase}out.dat
    rmdir $filebase
  done
  rmdir ${filebase0}/${i}

done
rmdir ${filebase0}
