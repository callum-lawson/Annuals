#!/bin/bash

### Inputs
datapath="/gpfs1/data/idiv_brose/lawson/Annuals/"
workpath="/gpfs1/work/lawson/"
# nk=10

### Referencing
cd $datapath

NAME_PREFIX=${1:-ESS_Annuals}
label=$NAME_PREFIX-$(date +%FT%H-%M-%S)
storepath=/work/$USER/$label

mkdir /work/$USER/logs
mkdir $storepath || {
    echo "output dir $storepath already exists" >&2
    exit 1
}

### Input
bash Source/ESS_input.sh \
  $storepath \
  $label 

### Program
#ARRAY_JOB_ID=$(qsub \
#	-terse \
#	-t 1-2 \
#	-wd $datapath \
#	Source/ESS_program.sub \
#	  $storepath \
#	  $label \
#	  $TASK_ID
#)

### Assemble
#qsub \
#    -hold_jid $ARRAY_JOB_ID \
#    -wd $datapath \
#    Source/ESS_assemble.sub \
#  	  $storepath \
#  	  $label
