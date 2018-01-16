### Run single instance of ESS calculation

# Parsing arguments -------------------------------------------------------

library(optparse)

parser <- OptionParser(
  usage       = "Rscript %prog datapath workpath label task_id",
  description = "\ncalculate ESS"
)

cli <- parse_args(parser, positional_arguments = 4)
# other args: sourcepath, workpath

sourcepath <- paste0(cli$args[1],"Source/")
workpath   <- cli$args[2]
label      <- cli$args[3]
task_id    <- cli$args[4]

# Functions ---------------------------------------------------------------

source(paste0(sourcepath,"ESS_functions.R"))

# Read in data ------------------------------------------------------------

pd <- readRDS(paste0(workpath,"ESS_input_",label,".rds"))
pdi <- pd[task_id,]

# Calculate ESS -----------------------------------------------------------

set.seed(pdi$iteration)
ess <- evolve(pdi) 

# Save ESS ----------------------------------------------------------------

saveRDS(ess,paste0(workpath,"ESS_output_",task_id,"_",label,".rds"))