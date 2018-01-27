### Run single instance of ESS calculation

# Parsing arguments -------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

storepath <- args[1]
label     <- args[2]
task_id   <- as.numeric(args[3])

# Functions ---------------------------------------------------------------

source("Source/ESS_functions.R")

# Read in data ------------------------------------------------------------

pd <- readRDS(paste0("Sims/ESS_input_",label,".rds"))
pdi <- pd[task_id,]

# Calculate ESS -----------------------------------------------------------

set.seed(pdi$iteration)
ess <- evolve(pdi) 

# Save ESS ----------------------------------------------------------------

saveRDS(ess,paste0(storepath,"/ESS_output_",task_id,"_",label,".rds"))