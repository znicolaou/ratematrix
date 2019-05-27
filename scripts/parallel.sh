#!/bin/bash
threads=16

for num in `seq 3 20`; do

atoms="$num $num $num"
filebase0=data/parallel$num

mkdir -p $filebase0
rm -r ${filebase0}*
mkdir -p $filebase0

start=`date +%s`
for i in `seq 0 $num`; do
  for j in `seq 0 $num`; do
    js=`jobs | wc -l`
    while [ $js -ge $threads ]; do
      sleep 0.01
      js=`jobs | wc -l`
    done
    ./ratematrix.py --filebase ${filebase0}/${i}_${j} --atoms $atoms --fix 1 $i 2 $j --calculate 0 >> ${filebase0}.out &
  done
done
end=`date +%s`
wait

dim=`awk '{dim+=$2}END{print dim}' ${filebase0}.out`
count=`awk '{count+=$4}END{print count}' ${filebase0}.out`
level=`awk '{if($5>max){max=$5}}END{print max}' ${filebase0}.out`

echo "$((3*num)) $dim $((end-start)) $count $level"
done
