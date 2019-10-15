#!/bin/bash
#make sure to include the above line (nohup can try to use dash instead of bash sometimes, and do funny things)

threads=72 #number of jobs to run concurrently

mkdir -p data/h2o2
for num in `seq 3 10`; do

#prepare different directories for each job to output to
filebase0=data/h2o2/$num/
mkdir -p ${filebase0}temp

for i in `seq 0 $((4*num))`; do
  for j in `seq 0 $((2*num))`; do
    #only submit $threads jobs at a time, and wait for them to finish
    js=`jobs | wc -l`
    while [ $js -ge $threads ]; do
      sleep 0.01
      js=`jobs | wc -l`
    done
    #submit the script with different parameters in the background (include ampersand)
    #I used argparse package in python to process these arguments easily
    ./ratematrix.py --filebase ${filebase0}temp/${i}_${j} --reference 0 $H2 3 $num 4 1 8 $AR --fix 1 $i 2 $j 8 $AR --calculate 0 0 --eigenvalues 0 --adiabatic $adiabatic --temperature $temperature &
  done
done
wait #include this wait, otherwise the script may terminate before background jobs complete
sleep 1

done
