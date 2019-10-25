 #!/bin/bash
N_temp=50
N_press=50
press0=0.01
press1=10.0
temp0=400
temp1=1000
filebase0=data/sweep2
procs=32

num=5
adiabatic=1
AR=$((20*num))
H2=$((2*num))
mkdir -p ${filebase0}

mkdir -p ${filebase0}
for i in `seq 0 $N_temp`; do
  for j in `seq 0 $N_press`; do
    js=`jobs | wc -l`
    while [ $js -ge $procs ]; do
      js=`jobs | wc -l`
      sleep 1
    done
    temp=`bc -l <<<"${temp0}+(${temp1}-${temp0})*${i}*1.0/${N_temp}"`
    press=`bc -l <<<"${press0}*e(${j}*1.0/${N_press}*l(${press1}/${press0}))"`
    echo $temp $press
    filebase=$filebase0/${i}_${j}
    ./ratematrix.py --filebase ${filebase0}/${i}_${j} --reference 0 $H2 3 $num 4 1 8 $AR --calculate 0 -1 --eigenvalues -1 --propogate 0 --adiabatic $adiabatic --temperature $temp --pressure $press --print 1 &

  done
done
wait
