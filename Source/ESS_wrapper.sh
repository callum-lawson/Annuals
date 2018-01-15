curdate="trial"

module load R

Rscript ~/Annuals/Source/ESS_input.R $curdate

qsub Ess_input_submit.sub -t 1-2

Rscript ~/Annuals/Source/ESS_assemble.R $curdate