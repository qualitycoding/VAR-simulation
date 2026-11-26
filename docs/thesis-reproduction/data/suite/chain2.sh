#!/bin/bash
cd /home/claude/work/var/src
for N in 21 31 41; do
  for s in step gaussian; do
    python3 thesis_suite.py $s $N flux_upwind ../../out/suite/nodal_${s}_${N}.json 150 '{"profile":"nodal","c0":1e10,"rth":0.1,"rfl":0.5}' > ../../out/suite/nodal_${s}_${N}.log 2>&1
  done
done
