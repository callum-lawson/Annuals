### Run single instance of ESS calculation

library(optparse)

# Parsing arguments -------------------------------------------------------

parser <- OptionParser(
  usage       = "Rscript %prog sourcepath outpath curdate",
  description = "\ncalculate ESS"
)

cli <- parse_args(parser, positional_arguments = 2)

# Shortcuts ---------------------------------------------------------------

task_id <- cli$args[1]
curdate <- cli$args[2]

# other args: sourcepath, outpath

# Functions ---------------------------------------------------------------

source("~/Annuals/Source/ESS_functions.R")

# Read in data ------------------------------------------------------------

pd <- readRDS(paste0("~/work/lawson/ESS_input_",curdate,".rds"))
pdi <- pd[task_id,]

# Calculate ESS -----------------------------------------------------------

set.seed(task_id)
ess <- evolve(pdi) 

# Save ESS ----------------------------------------------------------------

saveRDS(ess,paste0("~/work/lawson/ESS_output_",task_id,"_",curdate,".rds"))
