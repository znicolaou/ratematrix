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
mkdir -p ${filebase0}temp

#Calculate state space
start=`date +%s%N`
for i in `seq 0 $((4*num))`; do
  for j in `seq 0 $((2*num))`; do
    js=`jobs | wc -l`
    while [ $js -ge $threads ]; do
      sleep 0.01
      js=`jobs | wc -l`
    done
    ./ratematrix.py --filebase ${filebase0}temp/${i}_${j} --reference 0 $H2 3 $num 4 1 8 $AR --fix 1 $i 2 $j 8 $AR --calculate 0 0 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature  &
  done
done
end=`date +%s%N`
wait

for i in `seq 0 $((4*num))`; do
  for j in `seq 0 $((2*num))`; do
    cat ${filebase0}temp/${i}_${j}out.dat >> ${filebase0}pout.dat
    rm ${filebase0}temp/${i}_${j}out.dat
  done
done


dim=`awk 'BEGIN{n=0}{if(n%3==0){dim+=$1;} n++;}END{print dim}' ${filebase0}pout.dat`
ns=`awk 'BEGIN{n=0}{if(n%3==0){dim=$2;} n++;}END{print dim}' ${filebase0}pout.dat`
cputime=`awk 'BEGIN{n=0}{if(n%3==0){count+=$3}n++;}END{print count}' ${filebase0}pout.dat`
echo "$dim $ns $cputime >> ${filebase0}temp/out.dat
tail -n 2 ${filebase0}pout.dat >> ${filebase0}temp/out.dat

rm ${filebase0}pout.dat

./ratematrix.py --filebase ${filebase0}temp/ --reference 0 $H2 3 $num 4 1 8 $AR --calculate 0 0 --accumulate 1 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature

mv ${filebase0}temp/out.dat $filebase0
mv ${filebase0}temp/multiindices.npy $filebase0
mv ${filebase0}temp/temperatures.npy $filebase0
mv ${filebase0}temp/pressures.npy $filebase0
rm -r ${filebase0}temp/

end=`date +%s%N`
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
