#!/bin/bash

### Inputs
datapath=/data/idiv_brose/lawson/Annuals/
nk=${1:-100}

label="nk"$nk"_"$(date +%F-%H-%M-%S)
storepath=/work/$USER/ESS_$label

mkdir $storepath || {
    echo "output dir $storepath already exists" >&2
    exit 1
}
mkdir $storepath/logs

### Input
cd $datapath
module load R
Rscript \
  Source/ESS_input.R \
    $datapath \
    $label \
    $nk

### Program
PROGRAM=$(qsub \
  -terse \
  -t 1-2 \
  -wd $datapath \
  -o $storepath/logs/ESS_program_$label.out \
  -e $storepath/logs/ESS_program_$label.err \
  Source/ESS_program.sub \
    $storepath \
    $label 
)

### Assemble
# qsub \
#   -hold_jid $PROGRAM \
#   -wd $datapath \
#   -o $storepath/logs/ESS_assemble_$label.out \
#   -e $storepath/logs/ESS_assemble_$label.err \
#   Source/ESS_assemble.sub \
#     $storepath \
#     $label
