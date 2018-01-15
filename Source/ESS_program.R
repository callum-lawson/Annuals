### Run single instance of ESS calculation

# Hard inputs -------------------------------------------------------------

sourcepath <- "~/Annuals/Source/"
workpath <- "/gpfs1/work/lawson/"

# Parsing arguments -------------------------------------------------------

library(optparse)

parser <- OptionParser(
  usage       = "Rscript %prog task_id curdate",
  description = "\ncalculate ESS"
)

cli <- parse_args(parser, positional_arguments = 2)
# other args: sourcepath, workpath

task_id <- cli$args[1]
curdate <- cli$args[2]

# Functions ---------------------------------------------------------------

source(paste0(sourcepath,"ESS_functions.R"))

# Read in data ------------------------------------------------------------

pd <- readRDS(paste0(workpath,"ESS_input_",curdate,".rds"))
pdi <- pd[task_id,]

# Calculate ESS -----------------------------------------------------------

set.seed(pdi$iteration)
ess <- evolve(pdi) 

# Save ESS ----------------------------------------------------------------

saveRDS(ess,paste0(workpath,"ESS_output_",task_id,"_",curdate,".rds"))
