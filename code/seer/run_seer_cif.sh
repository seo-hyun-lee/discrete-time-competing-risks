#!/bin/sh
## run_seer_cif.sh -- CIF bootstrap as 4 independent R processes (B replications in total)
##   sh run_seer_cif.sh            (B=100 sh run_seer_cif.sh for a shorter run)
## then, at any time (uses the replications saved so far):  Rscript seer_cif.R figure
## Progress: seer_cif_w1.log ... w4.log (one line per replication). A stopped run resumes:
## just run this script again (finished replications are skipped).
export OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=1
B=${B:-200}; PER=$((B / 4))
for p in 1 2 3 4; do
  from=$(( (p - 1) * PER + 1 )); to=$(( p * PER ))
  nohup caffeinate -i env B=$B BOOT=$from:$to Rscript seer_cif.R boot > seer_cif_w$p.log 2>&1 &
done
echo "started 4 processes, $B replications; progress: tail -n 2 seer_cif_w*.log"
