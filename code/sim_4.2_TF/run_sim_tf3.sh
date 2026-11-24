#!/bin/sh
## run_sim_tf3.sh -- run sim_tf3.R as 4 independent R processes (no forking), 25 replications each
##   sh run_sim_tf3.sh
## Each process loops over smooth -> step -> abrupt for its own replications. Progress for all
## four goes to sim_tf3_progress_<shape>_n1000_J30.log; the final sim_tf3_<shape>_n1000_J30.rds
## is written by whichever process finishes that shape last. Errors: sim_tf3_w1.log ... w4.log.
## A stopped run resumes where it left off: just run this script again.
export OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1   # one thread per process
NPROC=4; B=100; PER=$((B / NPROC))
for p in $(seq 1 $NPROC); do
  from=$(( (p - 1) * PER + 1 )); to=$(( p * PER ))
  nohup caffeinate -i sh -c "for s in smooth step abrupt; do SHAPE=\$s B=$B MC_CORES=1 RHO_MIN=1 REPS=$from:$to Rscript sim_tf3.R; done" > sim_tf3_w$p.log 2>&1 &
done
echo "started $NPROC processes; progress: tail -f sim_tf3_progress_smooth_n1000_J30.log"
