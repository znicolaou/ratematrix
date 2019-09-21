#!/bin/bash
threads=4

for num in `seq 3 4`; do
echo $num
filebase0=data/h2o2/$num

dim=`head -n 1 ${filebase0}out.dat | awk '{print $2}'`
start=0

for start in `seq 0 100 $((dim+100))`; do
echo $start
js=`jobs | wc -l`
while [ $js -ge $threads ]; do
  sleep 1
  js=`jobs | wc -l`
done
./ratematrix.py --filebase ${filebase0} --calculate $start $((start+100)) --accumulate 1 --eigenvalues 0 &
done
wait

./add.py --filebase $filebase0 --dimension $dim
rm ${filebase0}*.npz
./ratematrix.py --filebase ${filebase0} --calculate $((dim+1)) $dim --accumulate 1 --eigenvalues 1 &

done
