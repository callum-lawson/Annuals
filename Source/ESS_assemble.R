### Assemble individual RDS ESS outputs into a single RDS file

# Parsing arguments -------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

storepath <- args[1]
label     <- args[2]

# Read in files -----------------------------------------------------------

pd <- readRDS(paste0("Sims/ESS_",label,"/ESS_input_",label,".rds"))
ccat <- interaction(pd[,c("zm","zs")])
climate <- as.numeric(factor(ccat),levels=unique(ccat))
nc <- nlevels(ccat)
ntasks  <- nrow(pd) # = nj * ni * nm
finite <- with(pd[1,], nk>0 & nk<Inf)

zw <- with(pd[1,], array(dim=c(nt,2,nj,ni,nc)))
if(finite==FALSE){
  es <- with(pd[1,], array(dim=c(nr,2,nj,ni,nc)))
} 
if(finite==TRUE){
  es <- with(pd[1,], array(dim=c(nr,3,nj,ni,nc)))
} 
  # extra column for extinction in finite populations

for(i in 1:ntasks){
  filename <- paste0(storepath,"/ESS_output_",i,"_",label,".rds")
  if(file.exists(filename)){
    curlist <- readRDS(filename)
    zw[,,pd$species[i],pd$iteration[i],climate[i]] <- curlist$zw
    es[,,pd$species[i],pd$iteration[i],climate[i]] <- as.matrix(curlist$es)
  }
}

# Save output -------------------------------------------------------------

saveRDS(zw,paste0("Sims/ESS_",label,"/ESS_climate_",label,".rds"))
saveRDS(es,paste0("Sims/ESS_",label,"/ESS_strategies_",label,".rds"))