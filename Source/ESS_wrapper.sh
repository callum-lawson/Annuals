#!/bin/bash

curdate=date +%d%b%Y
qsub -N ESS_input ~/Annuals/Source/ESS_input.sub $curdate $JOB_ID
qsub -N ESS_program -t 1-2 -hold_jid ESS_input ~/Annuals/Source/ESS_program.sub $curdate $JOB_ID
qsub -N ESS_assemble -hold_jid ESS_input,ESS_program ~/Annuals/Source/ESS_assemble_submit.sub $curdate $JOB_ID
