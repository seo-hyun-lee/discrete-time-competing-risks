#!/bin/sh
## run_seer_tf3.sh -- SEER trend-filtering fits as 4 independent R processes
##   sh run_seer_tf3.sh            (then, when all four are done:  Rscript seer_tf3.R summary)
## Each process runs one log-link job (minutes to ~1 h) and then one logit-link job (hours).
## Progress: seer_out3/log_<link>_k?_b?.txt (one line per rho). Errors: seer_tf3_w1.log ... w4.log.
## A stopped run resumes at the next rho: just run this script again.
export OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1
Rscript seer_tf3.R fast || exit 1
for p in 1 2 3 4; do
  nohup caffeinate -i env JOBS=$p,$((p + 4)) Rscript seer_tf3.R fit > seer_tf3_w$p.log 2>&1 &
done
echo "started 4 processes; progress: tail -f seer_out3/log_*.txt"
