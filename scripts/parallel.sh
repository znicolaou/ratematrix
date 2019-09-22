#!/bin/bash
threads=4

mkdir -p data/h2o2
for num in `seq 3 4`; do

filebase0=data/h2o2/$num
adiabatic=1
temperature=700

mkdir -p $filebase0
rm -r ${filebase0}*
mkdir -p $filebase0

start=`date +%s`
for i in `seq 0 $((4*num))`; do
  for j in `seq 0 $((2*num))`; do
    js=`jobs | wc -l`
    while [ $js -ge $threads ]; do
      sleep 0.01
      js=`jobs | wc -l`
    done
    ./ratematrix.py --filebase ${filebase0}/${i}_${j} --reference 0 $((2*num)) 3 $num  --fix 1 $i 2 $j --calculate 0 0 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature &
  done
done
end=`date +%s`
wait

for i in `seq 0 $((4*num))`; do
  for j in `seq 0 $((2*num))`; do
    cat ${filebase0}/${i}_${j}out.dat >> ${filebase0}pout.dat
    rm ${filebase0}/${i}_${j}out.dat
  done
done


dim=`awk 'BEGIN{n=0}{if(n%3==0){dim+=$2;} n++;}END{print dim}' ${filebase0}pout.dat`
count=`awk 'BEGIN{n=0}{if(n%3==0){count+=$4}n++;}END{print count}' ${filebase0}pout.dat`
level=`awk 'BEGIN{n=0}{if(n%3==0){if($5>max){max=$5}}n++;}END{print max}' ${filebase0}pout.dat`
echo "$((6*num)) $dim $((end-start)) $count $level" >> ${filebase0}out.dat
tail -n 2 ${filebase0}pout.dat >> ${filebase0}out.dat

rm ${filebase0}pout.dat

./ratematrix.py --filebase ${filebase0} --reference 0 $((2*num)) 3 $num  --calculate 0 0 --accumulate 1 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature
rm -r $filebase0

head -n 1 ${filebase0}out.dat >> data/runtimes4out.dat

dim=`head -n 1 ${filebase0}out.dat | awk '{print $2}'`
echo $dim

for start in `seq 0 100 $((dim+100))`; do
echo $start
js=`jobs | wc -l`
while [ $js -ge $threads ]; do
  sleep 1
  js=`jobs | wc -l`
done
./ratematrix.py --filebase ${filebase0} --reference 0 $((2*num)) 3 $num --calculate $start $((start+100)) --accumulate 1 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature &
done
wait

./ratematrix.py --filebase ${filebase0} --reference 0 $((2*num)) 3 $num --calculate 0 0 --accumulate 1 --eigenvalues 100 --adiabatic $adiabatic --temperature $temperature &
done
