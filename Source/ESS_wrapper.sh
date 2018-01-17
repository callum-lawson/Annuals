#!/bin/bash

datapath="/gpfs1/data/idiv_brose/lawson/Annuals/"
workpath="/gpfs1/work/lawson/"
nk=10

cd $datapath

NAME_PREFIX=${1:-ESS_Annuals}
intermediate=/work/$USER/$NAME_PREFIX-$(date +%FT%H-%M-%S)

intermediate/out-1.rds

mkdir /work/$USER/logs
mkdir $intermediate || {
    echo "output dir $intermediate already exists" >&2
    exit 1
}


label=$NAME-$JOB_ID-date +%d%b%Y

bash Source/ESS_input.sub $label $nk

ARRAY_JOB_ID=$(
    qsub \
	-terse \
	-N ESS_program \
	-t 1-2 \
	-wd $datapath \
	Source/ESS_program.sub \
	  $workpath \
	  $label
)

qsub \
    -N ESS_assemble \
    -hold_jid $ARRAY_JOB_ID \
    -wd $datapath \
    Source/ESS_assemble.sub \
      $intermediate \
      $label
