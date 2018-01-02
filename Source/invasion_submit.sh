#!/bin/bash
 
#$ -S /bin/bash 
#$ -N invasion_submit 
#$ -l h_rt=60,h_vmem=1G 
#$ -binding linear:1 

# Output files 
#$ -o invasion.out 
#$ -e invasion.err 
 
 # -o /work/$USER/$JOB_NAME-$JOB_ID
 
# Bash script
for arg in "$@" ; do
  echo $arg
done

# Run script: bash invasion_submit.sh 1 2 3 4