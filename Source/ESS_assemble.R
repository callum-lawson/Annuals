# Assemble individual RDS ESS outputs into a single file

cur_date <- format(Sys.Date(),"%d%b%Y")

input <- readRDS(paste0("Sims/ESS_input_",cur_date,".rds"))


zw <- with(input[1,], array(nt,2,nj*ni))
zw <- with(input[1,], array(nt,2,nj*ni))

ess <- with(input[1,], array())

# for(i in 1:)