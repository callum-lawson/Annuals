#!/bin/bash

### Inputs
datapath=/data/idiv_brose/lawson/Annuals
# nk=10

### Referencing
cd $datapath

NAME_PREFIX=${1:-ESS_Annuals}
label=$NAME_PREFIX-$(date +%F-%H-%M-%S)
storepath=/work/$USER/$label

mkdir $storepath || {
    echo "output dir $storepath already exists" >&2
    exit 1
}
mkdir $storepath/logs

### Input
INPUT=$(qsub \
  -terse \
  -wd $datapath \
  Source/ESS_input.sub \
    $storepath \
    $label 
)

### Program
#PROGRAM=$(qsub \
#  -hold_jid $INPUT \
#  -terse \
#  -t 1-2 \
#  -wd $datapath \
#  Source/ESS_program.sub \
#    $storepath \
#    $label \
#    $TASK_ID
#)

### Assemble
#qsub \
#  -hold_jid $PROGRAM \
#  -wd $datapath \
#  Source/ESS_assemble.sub \
#    $storepath \
#    $label
