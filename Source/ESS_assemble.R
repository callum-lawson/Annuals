### Assemble individual RDS ESS outputs into a single RDS file

# Hard inputs -------------------------------------------------------------

# workpath <- "Sims/"
workpath <- "/gpfs1/work/lawson/"
savepath <- "/gpfs1/data/idiv_brose/lawson/Annuals/Sims/"

# Parsing arguments -------------------------------------------------------

library(optparse)

parser <- OptionParser(
  usage       = "Rscript %prog nstasks curdate",
  description = "\ncombine ESS outputs into single RDS file"
)

cli <- parse_args(parser, positional_arguments = 2)
# other args: sourcepath, workpath

ntasks  <- cli$args[1]
curdate <- cli$args[2]

# Read in files -----------------------------------------------------------

pd <- readRDS(paste0(workpath,"ESS_input_",curdate,".rds"))
finite <- with(pd[1,], nk>0 & nk<Inf)

zw <- with(pd[1,], array(dim=c(nt,2,nj,ni)))
if(finite==FALSE){
  es <- with(pd[1,], array(dim=c(nr,2,nj,ni)))
} 
if(finite==TRUE){
  es <- with(pd[1,], array(dim=c(nr,3,nj,ni)))
} 
  # extra column for extinction in finite populations

for(i in 1:ntasks){
  curlist <- readRDS(paste0(workpath,"ESS_output_",ntasks,"_",curdate,".rds"))
  zw[,,pd$species[i],pd$iteration[i]] <- curlist$zw
  es[,,pd$species[i],pd$iteration[i]] <- as.matrix(curlist$es)
}

# Save output -------------------------------------------------------------

outlist <- list(pd=pd,zw=zw,es=es)
saveRDS(outlist,paste0(workpath,"ESS_output_",curdate,".rds"))
