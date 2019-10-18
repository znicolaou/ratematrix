#!/bin/bash

threads=72 #number of jobs to run concurrently
mem=150

for num in `seq 3 10`; do
mkdir -p data/h2o2

echo $num
filebase0=data/h2o2/$num/
adiabatic=1
temperature=1000
AR=$((20*num))
H2=$((2*num))

#Calculate state space
start=`date +%s%N`
./ratematrix.py --filebase ${filebase0} --reference 0 $H2 3 $num 4 1 8 $AR --calculate 0 0 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature
end=`date +%s%N`
dim=`head -n 1 ${filebase0}out.dat | awk '{print $1}'`
cputime=`head -n 1 ${filebase0}out.dat | awk '{print $3}'`
runtime=`bc -l <<< "($end-$start)*0.000000001"`
echo "state space dimension: $dim"
echo "state space cputime: $cputime"
echo "state space runtime: $runtime"

#Calculate matrix entries
starttime=`date +%s%N`

dim=`head -n 1 ${filebase0}out.dat | awk '{print $2}'`

for start in `seq 0 100 $((dim+100))`; do
js=`jobs | wc -l`
while [ $js -ge $threads ]; do
  sleep 1
  js=`jobs | wc -l`
done
./ratematrix.py --filebase ${filebase0} --reference 0 $H2 3 $num 4 1 8 $AR --calculate $start $((start+100)) --accumulate 1 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature &
done
wait

cat ${filebase0}_*cout.dat >> ${filebase0}cout.dat
rm ${filebase0}_*cout.dat

cputime=`awk '{t+=$3}END{print(t)}' ${filebase0}cout.dat`
end=`date +%s%N`
runtime=`bc -l <<< "($end-$starttime)*0.000000001"`
echo "calculate cputime: $cputime"
echo "calculate runtime: $runtime"

evals=`bc <<< "$num*10"`
evals=-1

sleep 5
./ratematrix.py --filebase ${filebase0} --reference 0 $H2 3 $num 4 1 8 $AR --calculate 0 0 --accumulate 1 --eigenvalues $evals --propogate 1 --adiabatic $adiabatic --temperature $temperature
runtime=`awk '{print $1}' ${filebase0}rout.dat`
echo "propogate runtime: $runtime"

rm -r ${filebase0}rows
rm -r ${filebase0}columns
rm -r ${filebase0}data
done
