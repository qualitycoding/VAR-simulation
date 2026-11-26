#!/bin/bash
cd /home/claude/work/var/src
for N in 31 41 61; do
  for s in step gaussian; do
    python3 thesis_suite.py $s $N $1 ../../out/suite/${s}_${N}_$1.json 150 > ../../out/suite/${s}_${N}_$1.log 2>&1
  done
done
