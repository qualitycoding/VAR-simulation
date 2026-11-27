#!/bin/bash
cd /home/claude/work/phase0
export PYTHONPATH=/home/claude/work/var/src
for i in $(seq 1 40); do
  python3 /home/claude/work/var/src/run_checkpointed.py m=101 n=101 target=12 epochs=3 time_budget=240 >> run.log 2>&1
  tail -n1 run.log | grep -q "STATUS=done" && break
done
echo FINISHED >> run.log
