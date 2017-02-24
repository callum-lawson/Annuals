# Annuals
Analyses of bet-hedging in desert annuals.

## Data
Master copies of data files used in the analysis.

## Models
Parameters stored from Stan model runs.

## Output
Data produced by R scripts (csv, rds, etc.).
	
## Plots
Graphical plots of results.

## Sims
Stored runs of population simulations.

## Source
R code.
* dataprocess: Set up seedling and seed databases.
* speciescorrect: Correct names of species in Census and Seed data and remove low-abundance species. 
* plot_sizes: Check and edit Venable plot area data.
* national_climate_process: Process rainfall data from US to produce single daily rainfall time series
* annual_climate_exploration: Expolore distribution and time series of rainfall data.
* prcp_projections: Running mean of prcp values.

* outputgraphs_GO: Create model fit plots for germination and seed survival models.
* outputgraphs_PrRs: Create model fit plots for reproduction models.
* prediction_functions: Functions for predicting reproduction, germination, and seed survival. 
* bootstrap_nonparametric: Non-parametric bootstrap of total counts.
* seedpredict: Predictions of new seed density for all plots. 

# Organise the above by workflow, not title