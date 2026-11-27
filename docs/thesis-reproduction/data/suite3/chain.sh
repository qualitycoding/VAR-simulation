#!/bin/bash
cd /home/claude/work/var/src
for N in 61 101 151; do
  for s in step gaussian; do
    if [ $s = step ]; then L=2700; else L=11500; fi
    python3 thesis_suite.py $s $N flux_vanleer ../../out/suite3/${s}_${N}.json 400 "{\"lor\":$L,\"seq\":[200,125,100,85,70],\"save\":\"/home/claude/work/out/suite3/${s}_${N}\"}" > ../../out/suite3/${s}_${N}.log 2>&1
  done
done
echo FIN > ../../out/suite3/all.done
