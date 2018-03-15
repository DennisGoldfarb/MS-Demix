#! /bin/bash

# first job - no dependencies
jid1=$(sbatch clean.sh)
jid1="${jid1##* }"

# second job - depends on job 1
jid2=$(sbatch --dependency=afterany:$jid1 demix.sh)
