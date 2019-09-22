#!/bin/bash
threads=4

mkdir -p data/h2o2
for num in `seq 2 4`; do
echo $num
filebase0=data/h2o2/$num/
adiabatic=0
temperature=1500

mkdir -p ${filebase0}temp

start=`date +%s`
for i in `seq 0 $((4*num))`; do
  for j in `seq 0 $((2*num))`; do
    js=`jobs | wc -l`
    while [ $js -ge $threads ]; do
      sleep 0.01
      js=`jobs | wc -l`
    done
    ./ratematrix.py --filebase ${filebase0}temp/${i}_${j} --reference 0 $((2*num)) 3 $num  --fix 1 $i 2 $j --calculate 0 0 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature &> /dev/null &
  done
done
end=`date +%s`
wait

for i in `seq 0 $((4*num))`; do
  for j in `seq 0 $((2*num))`; do
    cat ${filebase0}temp/${i}_${j}out.dat >> ${filebase0}pout.dat
    rm ${filebase0}temp/${i}_${j}out.dat
  done
done


dim=`awk 'BEGIN{n=0}{if(n%3==0){dim+=$2;} n++;}END{print dim}' ${filebase0}pout.dat`
runtime=`awk 'BEGIN{n=0}{if(n%3==0){count+=$3}n++;}END{print count}' ${filebase0}pout.dat`
count=`awk 'BEGIN{n=0}{if(n%3==0){count+=$4}n++;}END{print count}' ${filebase0}pout.dat`
level=`awk 'BEGIN{n=0}{if(n%3==0){if($5>max){max=$5}}n++;}END{print max}' ${filebase0}pout.dat`
echo "$((6*num)) $dim $runtime $count $level" >> ${filebase0}temp/out.dat
tail -n 2 ${filebase0}pout.dat >> ${filebase0}temp/out.dat
echo "dimension $dim $runtime"

rm ${filebase0}pout.dat

./ratematrix.py --filebase ${filebase0}temp/ --reference 0 $((2*num)) 3 $num  --calculate 0 0 --accumulate 1 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature

mv ${filebase0}temp/out.dat $filebase0
mv ${filebase0}temp/multiindices.npy $filebase0
mv ${filebase0}temp/temperatures.npy $filebase0
mv ${filebase0}temp/pressures.npy $filebase0
rm -r ${filebase0}temp/

head -n 1 ${filebase0}out.dat >> data/runtimes4out.dat

dim=`head -n 1 ${filebase0}out.dat | awk '{print $2}'`

for start in `seq 0 100 $((dim+100))`; do
js=`jobs | wc -l`
while [ $js -ge $threads ]; do
  sleep 1
  js=`jobs | wc -l`
done
./ratematrix.py --filebase ${filebase0} --reference 0 $((2*num)) 3 $num --calculate $start $((start+100)) --accumulate 1 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature >> ${filebase0}_${start}cout.dat &
done
wait

cat ${filebase0}_*cout.dat >> ${filebase0}cout.dat
rm ${filebase0}_*cout.dat

runtime=`awk '{t+=$3}END{print(t)}' ${filebase0}cout.dat`
echo "calculate runtime: $runtime"

evals=`bc <<< "$dim/20"`

./ratematrix.py --filebase ${filebase0} --reference 0 $((2*num)) 3 $num --calculate 0 0 --accumulate 1 --eigenvalues $evals --adiabatic $adiabatic --temperature $temperature

rm -r ${filebase0}rows
rm -r ${filebase0}columns
rm -r ${filebase0}data

done
