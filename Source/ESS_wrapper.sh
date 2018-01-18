#!/bin/bash

### Inputs
datapath=/data/idiv_brose/lawson/Annuals
nk=10

### Referencing
cd $datapath

label="nk_"$nk"_"$(date +%F-%H-%M-%S)
storepath=/work/$USER/${1:-ESS}_$label

mkdir $storepath || {
    echo "output dir $storepath already exists" >&2
    exit 1
}
mkdir $storepath/logs

### Input
INPUT=$(qsub \
  -S /bin/bash
  -N ESS_input
  -l h_rt=60,h_vmem=1G
  -binding linear:1
  -o $storepath/logs/$JOB_NAME-$JOB_ID.log
  -j y
  -terse \
  -wd $datapath \
  Source/ESS_input.sub \
    $label 
)

## Program
PROGRAM=$(qsub \
  -S /bin/bash
  -N ESS_program
  -l h_rt=60,h_vmem=1G
  -binding linear:1
  -o $storepath/logs/$JOB_NAME-$JOB_ID-$TASK_ID.log
  -j y
  -hold_jid $INPUT \
  -terse \
  -t 1-2 \
  -wd $datapath \
  Source/ESS_program.sh \
    $storepath \
    $label \
    $TASK_ID
)

## Assemble
qsub \
  -S /bin/bash
  -N ESS_assemble
  -l h_rt=60,h_vmem=1G 
  -binding linear:1
  -o $storepath/logs/$JOB_NAME-$JOB_ID.log
  -j y
  -hold_jid $PROGRAM \
  -wd $datapath \
  Source/ESS_assemble.sh \
    $storepath \
    $label
