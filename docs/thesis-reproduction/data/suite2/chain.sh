#!/bin/bash
# $1 = tag, $2 = convection, $3 = extra kw json (without lor)
cd /home/claude/work/var/src
for N in 21 31 41; do
  python3 thesis_suite.py step $N $2 ../../out/suite2/$1_step_$N.json 150 "$(echo $3 | sed 's/}/,"lor":2700}/')" > ../../out/suite2/$1_step_$N.log 2>&1
  python3 thesis_suite.py gaussian $N $2 ../../out/suite2/$1_gaussian_$N.json 150 "$(echo $3 | sed 's/}/,"lor":11500}/')" > ../../out/suite2/$1_gaussian_$N.log 2>&1
done
echo FIN > ../../out/suite2/$1.done
