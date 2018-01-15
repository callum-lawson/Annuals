#!/bin/bash

trial="trial"
qsub -N ESS_input ~/Annuals/Source/ESS_input.sub $trial 
qsub -N ESS_program -t 1-2 -hold_jid ESS_input ~/Annuals/Source/ESS_program.sub $trial
#qsub ~/Annuals/Source/ESS_assemble_submit.sub -v curdate="trial" -hold_jid ESS_input
