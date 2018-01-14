### Assemble input for ESS analyses

# Hard inputs -------------------------------------------------------------

parpath <- "data/idiv_brose/lawson/Annuals/Models/"
workpath <- "/work/lawson/"

# Parsing arguments -------------------------------------------------------

library(optparse)

parser <- OptionParser(
  usage       = "Rscript %prog curdate",
  description = "\ngenerate inputs for ESS calculations"
)

cli <- parse_args(parser, positional_arguments = 1)
# other args: sourcepath, workpath

curdate <- cli$args[1]

# Input parameters --------------------------------------------------------

ni <- 2
nj <- 22
nk <- Inf
nt <- 20
nb <- 5 # number of "burn-in" timesteps to stabilise resident dynamics
nr <- 5 # number of repeated invasions

cmd <- c(0,0,0)  # amount of sds to add
csd <- c(-1,0,1) # powers
nc <- length(cmd)

# Basic parameters --------------------------------------------------------

plastic <- TRUE
rho <- 0.82
n0 <- 1
ddfun <- "BHS"
Sg <- 1

smut_a <- 5
if(plastic==TRUE)  smut_b <- 5
if(plastic==FALSE) smut_b <- 0

nsmin <- 10^-10
ngmin <- 10^-50

# Model parameters --------------------------------------------------------

pl <- list(
  go = readRDS(paste0(parpath,"go_pars_lnGtnt_BH_25Nov2017.rds")),
  # go = readRDS("Models/go_pars_tdistpois_naspecies_noerr_noGDD_loglik_RICKER_15Oct2017.rds"),
  gs = readRDS(paste0(parpath,"Models/gnzhh_onhh_pars_medians_26Oct2015.rds")),
  # gs = g site level
  # source script: venable_Stan_GO_descriptive_gnzhh_onhh_26Oct2015
  # uses tau_s = 100
  # but tau actually irrelevant because all multiplicative?
  pr = readRDS(paste0(parpath,"Models/pr_pars_yearhet_squared_pc_02Mar2016.rds")),
  rs = readRDS(paste0(parpath,"Models/rs_pars_yearhet_squared_pc_trunc_05Mar2016.rds"))
)
# already permuted

# Starting parameters -----------------------------------------------------
# Creates grid of starting alpha and beta values

am0 <- runif(ni,-5,5)
if(plastic==TRUE)  bm0 <- runif(ni,-5,5)
if(plastic==FALSE) bm0 <- 0

# Climate -----------------------------------------------------------------

pp <- read.csv("Output/prcp_projection_summaries_03Sep2017.csv",header=T)
mpam <- with(pp, median[measure=="mpam" & scenario==60 & yearcat==100])
mpsd <- with(pp, median[measure=="mpsd" & scenario==60 & yearcat==100])
# using projected season precipitation for 
# germination season precipitation change (both very similar)
# year = 2100
# Representative Concentration Pathway 6.0

ncy <- read.csv("Output/ncy_15Jan2016.csv",header=T)
ncy <- subset(ncy,is.na(seasprcp)==F)
# removes first value (missing because no previous winter)

zmo <- mean(log(ncy$seasprcp))
zso <- sd(log(ncy$seasprcp))

wmo <- mean(log(ncy$germprcp))
wso <- sd(log(ncy$germprcp))

zm <- zmo + cmd * zso
zs <- zso ^ csd
wm <- wmo + cmd * wso
ws <- wso ^ csd

# Parameter dataframe -----------------------------------------------------

posd <- expand.grid(ipos=1:ni,jpos=1:nj,mpos=1:nc)

itot <- 10^4 # total iterations in Stan output
ei <- posd$ipos
ec <- with(posd, jpos * itot + ipos)
em <- posd$mpos

pd <- with(pl, data.frame(
  iteration = posd$ipos,
  species = posd$jpos,
  zm = zm[em],
  zs = zs[em],
  wm = wm[em],
  ws = ws[em],
  am0 = am0[ei],
  bm0 = bm0[ei],
  beta_p1 = pr$beta_p[,,1][ec],
  beta_p2 = pr$beta_p[,,2][ec],
  beta_p3 = pr$beta_p[,,3][ec],
  beta_p4 = pr$beta_p[,,4][ec],
  beta_r1 = rs$beta_r[,,1][ec],
  beta_r2 = rs$beta_r[,,2][ec],
  beta_r3 = rs$beta_r[,,3][ec],
  beta_r4 = rs$beta_r[,,4][ec],
  sig_y_p = pr$sig_y_p[ec],
  sig_y_r = rs$sig_y_r[ec],
  sig_s_g = gs$sig_s_g[ec],
  sig_s_p = pr$sig_s_p[ei],
  sig_s_r = rs$sig_s_r[ei],
  sig_o_p = pr$sig_o_p[ei],
  phi_r = rs$phi_r[ei],
  theta_g = gs$theta_g[ec],
  So = exp(-exp(go$alpha_m[ec])),
  m0 = exp(go$alpha_m[ec]),
  m1 = exp(go$beta_m[ec]),
  ni,nj,nk,nt,nb,nr,
  n0,rho,ddfun,Sg,
  smut_a,smut_b,
  nsmin,ngmin,
  plastic
))

# Save RDS file -----------------------------------------------------------

saveRDS(pd,paste0(workpath,"ESS_input_",curdate,".rds"))
