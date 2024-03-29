# Annuals
Analyses of bet-hedging in desert annuals.

Workflow:
1. speciescorrect + plot_sizes + national_climate_process
2. dataprocess
3. reprod + GO_descriptive
4. bootstrap + seedpredict
5. poplevel + prcp_projections
6. popsims
7. traitmeasures

## Data
Master copies of data files used in the analysis.

## Models
Parameters stored from Stan model runs.

## Output
Data produced by R scripts (csv, rds, etc.).
* *aggregate_rainfall*: Daily total precipitation 1946-2014 (10ths of a mm?).
* *cd*: Raw plant reproduction census data.
* *cdpos*: cd for plants that produced at least one seed.
* *census_data_date_errors*: Date errors identified by Ursula (?).
* *csy*: cd aggregated over individuals within plots, then aggregated over plots within years.
* *csyp*: cd aggregated over individuals within plots.
* *csyp_seedests*: csyp with reproduction estimates added for all plots (based on reproduction model predictions).
* *DI_survival*: Table of density-independent seed survival rates as calculated from model.
* *gnzhh_onhh_pars_medians*: Variance components (plot,year,species) for germination and old seed counts (?), based on descriptive Stan model. 
* *medtraits*: Species-averaged vital rates, for use in explaining population simulation results. 
* *msy*: Combined census and seed bank aggregated data, with single count values for each species-year combination. 
* *msy*: msy with estimated new seed densities added (based on reproduction model predictions).
* *ncy*: Total precipitation over germination and growing season periods 1946-2014 (10ths of a mm?).
* *observation_error_byspecies*: observation errors for germinant, old seed, and new seed densities, based on non parametric bootstrap (?). 
* *pr_pars*: Parameters from probability of reproduction model.
* *rs_pars*: Parameters from model of expected number of seeds, given that the plant reproduced. 
* *sb*: Raw seed bank data.
* *Seed_Bank_sharedspecies*: UNCERTAIN.
* *species_list*: List of species abbreviations and full scientific names.
* *ssy*: sb, aggregated across replicates within plots, then plots within years.
* *ssyp*: sb, aggregated across replicates within plots.
* *Tvalues*: Duration of each part of season, estimated from germination and death dates.
* *venablePlots_processed*: Areas for each plant census plot (how have these been processed?).

## Plots
Graphical plots of results.

## Sims
Stored runs of population simulations.

## Source
R code.
* *annual_climate_exploration*: Expolore distribution and time series of rainfall data.
* *dataprocess*: Set up seedling and seed databases.
* *speciescorrect*: Correct names of species in Census and Seed data and remove low-abundance species. 
* *plot_sizes*: Check and edit Venable plot area data.
* *national_climate_process*: Process rainfall data from US to produce single daily rainfall time series
* *prcp_projections*: Running mean of prcp values.
* *outputgraphs_GO*: Create model fit plots for germination and seed survival models.
* *outputgraphs_PrRs*: Create model fit plots for reproduction models.
* *prediction_functions*: Functions for predicting reproduction, germination, and seed survival. 
* *bootstrap_nonparametric*: Non-parametric bootstrap of total counts.
* *seedpredict*: Predictions of new seed density for all plots. 
* *popsims_stochastic*: Simulations of population dynamics for finite area.
* *popsims_stochastic_medpars*: Simulations of population dynamics for finite area, using species-level median values for all parameters except germination ones.
* *popsims_analytical*: Simulations of population dynamics for infinite area.
* *poplevel_binomialG*: Germination and seed survival models in Stan, with binomial germination calculation.
* *poplevel_lnpoisdist*: Germination and seed survival models in Stan, with poisson seed counts.
* *talk_plots*: Population size and germination graphs for presentations.
* *traitmeasures*: Calculate species-level predictive traits (partially obsolete).
* *trait_functions*: Functions for calculating derived trait values from popualtion parameters.
* *simulation_functions_stochastic*: Functions for simulations of population dynamics for finite area.
* *simulation_functions_analytical*: Functions for simulations of population dynamics for infinite area.
* *figure_functions*: Functions to plot figures from Stan model output.
* *reprod*: Estimate per-capita reproduction parameters in Stan.
* *GO_descriptive*: Model to describe G and O values as function of species, year, and site, comparing negative binomial and lognormal-Poisson models (currently saving two versions, not sure which is right)
* *popsims_analytical_old*: Simulations of population dynamics for infinite area, density-independent growth (depricated).
* *mortality_priors*: Simulate gamma distributions to develop reasonable priors for density-independent seed mortality.
* *seasonwindow*: Find best times for open and close of growing season.
* *correlations*: Look at among-species correlations across years and plots.
* *popsize_timeseries*: Plot population sizes of seedlings and seeds through time. Only two versions, debatable whether to keep / update. 
* *queries_reproduction*: Subset and save census data to provide examples of potential problems for Venable to investigate
* *site_germinant_distributions*: How does distributing germinants among discrete sites affect the among-site variance in the number of germinants? 
* *invasion_sims*: Iterated simulations of resident population dynamics and
invasions by new strategies.   
* *invasion_functions*: Calculate ES alpha_G and beta_G by iteratively perturbing parameters and attempting re-invasion. 
* *ESS_explore*: Develop explanations for patterns in ESS germination results.
