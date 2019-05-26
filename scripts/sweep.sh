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

if [ -f ${filebase0}out.dat ]; then rm ${filebase0}out.dat; fi
for i in `seq 0 $N_temp`; do
  for j in `seq 0 $N_press`; do
    temp=`bc -l <<<"${temp0}+(${temp1}-${temp0})*${i}*1.0/${N_temp}"`
    press=`bc -l <<<"${press0}+(${press1}-${press0})*${j}*1.0/${N_press}"`
    filebase=${filebase0}/${i}/${j}/
    mkdir -p $filebase
    echo $temp $press
    echo "./ratematrix.py --temperature $temp --pressure $press --atoms $atoms --filebase ${filebase} --plot 0 --save $save"
    ./ratematrix.py --temperature $temp --pressure $press --atoms $atoms --filebase ${filebase} --plot 0 --save $save &
  done
  wait

  for j in `seq 0 $N_press`; do
    filebase=${filebase0}/${i}/${j}/
    cat ${filebase}out.dat >> ${filebase0}out.dat
    rm ${filebase}out.dat
    rmdir $filebase
  done
  rmdir ${filebase0}/${i}

done
rmdir ${filebase0}
