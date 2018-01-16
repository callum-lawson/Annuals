#!/bin/bash

datapath="/gpfs1/data/idiv_brose/lawson/Annuals/"
workpath="/gpfs1/work/lawson/"

label=$JOB_ID-date +%d%b%Y

qsub -N ESS_input ~/Annuals/Source/ESS_input.sub $datapath $label
qsub -N ESS_program -t 1-2 -hold_jid ESS_input ~/Annuals/Source/ESS_program.sub $datapath $workpath $label
qsub -N ESS_assemble -hold_jid ESS_input,ESS_program ~/Annuals/Source/ESS_assemble.sub $datapath $workpath $label