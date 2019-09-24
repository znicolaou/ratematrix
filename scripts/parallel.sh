#!/bin/bash
#SBATCH -A p30575
#SBATCH -n 100
#SBATCH -p short
#SBATCH -t 04:00:00
#SBATCH --mem=7500
#SBATCH --output=parallel.out
threads=72

for num in `seq 10 15`; do
mkdir -p data/h2o2

echo $num
filebase0=data/h2o2/$num/
adiabatic=1
temperature=1000

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
    ./ratematrix.py --filebase ${filebase0}temp/${i}_${j} --reference 0 $((2*num)) 3 $num 4 1  --fix 1 $i 2 $j --calculate 0 0 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature &
  done
done
end=`date +%s%N`
wait
sleep 1

for i in `seq 0 $((4*num))`; do
  for j in `seq 0 $((2*num))`; do
    cat ${filebase0}temp/${i}_${j}out.dat >> ${filebase0}pout.dat
    rm ${filebase0}temp/${i}_${j}out.dat
  done
done


dim=`awk 'BEGIN{n=0}{if(n%3==0){dim+=$2;} n++;}END{print dim}' ${filebase0}pout.dat`
cputime=`awk 'BEGIN{n=0}{if(n%3==0){count+=$3}n++;}END{print count}' ${filebase0}pout.dat`
count=`awk 'BEGIN{n=0}{if(n%3==0){count+=$4}n++;}END{print count}' ${filebase0}pout.dat`
level=`awk 'BEGIN{n=0}{if(n%3==0){if($5>max){max=$5}}n++;}END{print max}' ${filebase0}pout.dat`
echo "$((6*num)) $dim $cputime $count $level" >> ${filebase0}temp/out.dat
tail -n 2 ${filebase0}pout.dat >> ${filebase0}temp/out.dat

rm ${filebase0}pout.dat

./ratematrix.py --filebase ${filebase0}temp/ --reference 0 $((2*num)) 3 $num 4 1 --calculate 0 0 --accumulate 1 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature

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
./ratematrix.py --filebase ${filebase0} --reference 0 $((2*num)) 3 $num 4 1 --calculate $start $((start+100)) --accumulate 1 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature &
done
wait

cat ${filebase0}_*cout.dat >> ${filebase0}cout.dat
rm ${filebase0}_*cout.dat

cputime=`awk '{t+=$3}END{print(t)}' ${filebase0}cout.dat`
end=`date +%s%N`
runtime=`bc -l <<< "($end-$starttime)*0.000000001"`
echo "calculate cputime: $cputime"
echo "calculate runtime: $runtime"

evals=`bc <<< "$dim/20"`

./ratematrix.py --filebase ${filebase0} --reference 0 $((2*num)) 3 $num 4 1 --calculate 0 0 --accumulate 1 --eigenvalues 0 --propogate 1 --adiabatic $adiabatic --temperature $temperature
runtime=`awk '{print $1}' ${filebase0}rout.dat`
echo "eigenvalues runtime: $runtime"

rm -r ${filebase0}rows
rm -r ${filebase0}columns
rm -r ${filebase0}data
done
