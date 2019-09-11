#!/bin/bash
threads=100
num=2

atoms="$((4*num)) $((4*num)) $num $((2*num)) $num"
ref="3 $((2*num)) 13 $num 47 $num 48 $num"
filebase0=data/parallelgri$num

mkdir -p $filebase0
rm -r ${filebase0}*
mkdir -p $filebase0

start=`date +%s`
for i in `seq 0 $((4*num))`; do
  for j in `seq 0 $((4*num))`; do
    for k in `seq 0 $((4*num))`; do

      js=`jobs | wc -l`
      while [ $js -ge $threads ]; do
        sleep 0.01
        js=`jobs | wc -l`
      done
      ./ratematrix.py --filebase ${filebase0}/${i}_${j}_${k} --mechanism mechanisms/gri30.cti --atoms $atoms --fix 1 $i 2 $j 4 $k 48 $num --calculate 0 --adiabatic 1 --reference $ref &
    done
  done
done
end=`date +%s`
wait

for i in `seq 0 $((4*num))`; do
    for j in `seq 0 $((4*num))`; do
      for k in `seq 0 $((4*num))`; do
      cat ${filebase0}/${i}_${j}_${k}out.dat >> ${filebase0}pout.dat
      rm ${filebase0}/${i}_${j}_${k}out.dat
    done
  done
done

dim=`awk 'BEGIN{n=0}{if(n%3==0){dim+=$2;} n++;}END{print dim}' ${filebase0}pout.dat`
count=`awk 'BEGIN{n=0}{if(n%3==0){count+=$4}n++;}END{print count}' ${filebase0}pout.dat`
level=`awk 'BEGIN{n=0}{if(n%3==0){if($5>max){max=$5}}n++;}END{print max}' ${filebase0}pout.dat`
echo "$((12*num)) $dim $((end-start)) $count $level" >> ${filebase0}out.dat
tail -n 2 ${filebase0}pout.dat >> ${filebase0}out.dat
# rm ${filebase0}pout.dat

#./ratematrix.py --filebase ${filebase0} --mechanism mechanisms/gri30.cti --atoms $atoms --fix 1 $i 2 $j 4 $k 48 $num --calculate 0 --accumulate 1
# rm -r $filebase0
