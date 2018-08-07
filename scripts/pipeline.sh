#! /bin/bash

# first job - no dependencies
jid1=$(sbatch clean.sh)
jid1="${jid1##* }"
echo $jid1

# second job - depends on job 1
jid2=$(sbatch --dependency=afterany:$jid1 demix.sh)
jid2="${jid2##* }"
echo $jid2

# third job - depends on job 2
jid3=$(sbatch --dependency=afterany:$jid2 perc.sh)
jid3="${jid3##* }"
echo $jid3
