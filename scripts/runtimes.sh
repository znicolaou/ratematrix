#!/bin/bash
for i in {3..7}; do head -n 1 data/h2o2/$i/out.dat | awk -v num=$i '{print(num, $1, $3)}'; done > statesruntime.dat
for i in {3..7}; do head -n 1 data/h2o2/$i/cout.dat | awk -v num=$i '{print(num, $1, $2)}'; done > calculateruntime.dat
for i in {3..7}; do head -n 1 data/h2o2/$i/eout.dat | awk -v num=$i '{print(num, $1, $2)}'; done > eigenvaluesruntime.dat
for i in {3..7}; do head -n 1 data/h2o2/$i/rout.dat | awk -v num=$i '{print(num, $1, $2)}'; done > propogateruntime.dat
