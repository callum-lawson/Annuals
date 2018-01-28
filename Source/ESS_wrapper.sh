#!/bin/bash

### Inputs
datapath=/data/idiv_brose/lawson/Annuals/
ni=${1:-100}
nr=${2:-100}
nt=${3:-100}
nk=${4:-100}

label="i"$ni"_r"$nr"_t"$nt"_k"$nk"_"$(date +%F-%H-%M-%S)
storepath=/work/$USER/ESS_$label

mkdir $storepath || {
    echo "output dir $storepath already exists" >&2
    exit 1
}
mkdir $storepath/logs
mkdir $datapath/Sims/ESS_$label

### Input
cd $datapath
module load R
Rscript \
  Source/ESS_input.R \
    $datapath \
    $label \
    $ni \
    $nr \
    $nt \
    $nk

### Program
ntasks="$(($ni * 22 * 3))"
ARRAYJOB_ID=$(
  qsub \
    -N $label \
    -terse \
    -t 1-$ntasks \
    -wd $datapath \
    -o $storepath/logs/ESS_program_$label.out \
    -e $storepath/logs/ESS_program_$label.err \
    Source/ESS_program.sub \
      $storepath \
      $label
)

## Assemble
qsub \
  -hold_jid "$label*" \
  -wd $datapath \
  -o $storepath/logs/ESS_assemble_$label.out \
  -e $storepath/logs/ESS_assemble_$label.err \
  Source/ESS_assemble.sub \
    $storepath \
    $label
