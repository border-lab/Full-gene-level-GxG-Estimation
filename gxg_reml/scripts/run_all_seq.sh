#!/bin/bash
export OMP_NUM_THREADS=6 OPENBLAS_NUM_THREADS=6 MKL_NUM_THREADS=6 NUMEXPR_NUM_THREADS=6
cd /home/rsb/Dropbox/students/fast-gxg-trace
for s in he_vs_reml check_knobs adversarial effective_markers crlb_floor verify_consistency additive_vs_epi joint_scaling explain_ziyan; do
  echo "######################## $s.py ########################"
  /usr/bin/time -v python3 scripts/$s.py 2>&1 | grep -vE "Elapsed|Maximum resident|^[[:space:]]*(User|System|Percent|Average|Page|Minor|Major|Voluntary|Involuntary|Swaps|File system|Socket|Context|Command being|Exit status|Maximum)" 
  echo "DONE_$s"
done
echo "ALL_SEQ_COMPLETE"
