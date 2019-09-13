#!/bin/bash
threads=32

for num in `seq 7 10`; do

atoms="$((2*num)) $((4*num)) $num"
filebase0=data/parallel$num

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
    ./ratematrix.py --filebase ${filebase0}/${i}_${j} --atoms $atoms --fix 1 $i 2 $j 8 $num --calculate 0 &
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
echo "$((7*num)) $dim $((end-start)) $count $level" >> ${filebase0}out.dat
tail -n 2 ${filebase0}pout.dat >> ${filebase0}out.dat

#rm ${filebase0}pout.dat

./ratematrix.py --filebase ${filebase0} --atoms $atoms --fix 1 $i 2 $j --calculate 1 --accumulate 1 --plot 0
#rm -r $filebase0

#cat ${filebase0}out.dat >> data/runtimes3out.dat
done
