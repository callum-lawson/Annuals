### Run instance of ESS calculation

# Inputs ------------------------------------------------------------------

# read-in file path (specified in submit script)
# ij <- task number
ij <- 1
cur_date <- format(Sys.Date(),"%d%b%Y")
# Should be passed on from wrapper function
# SAVE FILE PATH (SPECIFIED IN WRAPPER SCRIPT)

# Functions ---------------------------------------------------------------

source("Source/ESS_functions.R")

# Read in data ------------------------------------------------------------

pd <- readRDS(paste0("Sims/ESS_inputs_",cur_date,".rds"))
pdi <- pd[ij,]

# Calculate ESS -----------------------------------------------------------

set.seed(ij)
ess <- evolve(pdi) 

# Save ESS ----------------------------------------------------------------

saveRDS(ess,paste0("Sims/",savefile,"_",cur_date,".rds"))
