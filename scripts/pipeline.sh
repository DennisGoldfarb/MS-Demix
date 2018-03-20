#! /bin/bash

# first job - no dependencies
jid1=$(sbatch clean.sh)
jid1="${jid1##* }"

# second job - depends on job 1
jid2=$(sbatch --dependency=afterany:$jid1 demix.sh)
jid2="${jid2##* }"

# third job - depends on job 2
jid3=$(sbatch --dependency=afterany:$jid2 perc.sh)
