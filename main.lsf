#!/bin/bash

#BSUB -q q_hpc
#BSUB -n 32
#BSUB -oo %J.out
#BSUB -eo %J.err
#BSUB -R "span[hosts=1]"

module load git

###show current git commit
git log -1
###show current git commit hash in short format
git show -s --format=%h
###show current branch name
git rev-parse --abbrev-ref HEAD

echo "LSB_DJOB_NUMPROC = $LSB_DJOB_NUMPROC"
echo "LSB_DJOB_RANKFILE = $LSB_DJOB_RANKFILE"
echo "LSB_HOSTS = $LSB_HOSTS"
echo "LSB_MCPU_HOSTS = $LSB_MCPU_HOSTS"
echo "LSB_JOBNAME = $LSB_JOBNAME"

cp $LSB_DJOB_RANKFILE .

echo "*** START main.lsf ***"
cat main.lsf
echo "*** END   main.lsf ***"
echo "********************"
echo "***START main.jl***"
cat main.jl
echo "***END main.jl***"

$HOME/bin/julia1.5 --version
$HOME/bin/julia1.5 -t $LSB_DJOB_NUMPROC --project=@. main.jl
