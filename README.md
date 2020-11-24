# Group testing for efficient screening of SARS-CoV-2
This code accompanies our manuscript *Using viral load and epidemic dynamics to optimize pooled testing in resource constrained settings*.

## Getting Started
------------
There are two parts to this project:
  A. Fitting and simulation of viral load trajectories
  B. Assessing pooling designs for individual identification and prevalence estimation using simulated viral loads
  
The first part of this project takes place entirely in R. We used R version 4.0.2. The second part takes place entirely in Python. For part A, continue reading. For part B, please see [this section](code/group_test_simulations/README.md) to get started.


## Part A: Viral kinetics analyses
------------

### 1. Libraries
------------
The project requires the following R packages:
```r
library(tidyverse)
library(data.table)
library(coda)
library(MASS)
library(doParallel)
library(ggpubr)
library(rethinking)
library(extraDistr)
library(patchwork)
library(lazymcmc)
library(deSolve)
library(Matrix)
library(odin)
```
  
There is also a dependency on a custom R package for MCMC sampling, `lazymcmc` available [here](https://github.com/jameshay218/lazymcmc). Note that this package can be installed using `devtools::install_github("jameshay218/lazymcmc")`.

### 2. Git repo structure
------------
Code for this part are in the [code/viral_kinetics/](https://github.com/cleary-lab/covid19-group-tests/tree/master/code/viral_kinetics/) folder. The repo structure is as follows:
  - `chains`: this folder contains the MCMC chains from fitting the viral kinetics model to swab and sputum data. Chains are stored in separate folders for these two scenarios.
  - `data`: folder containing the digitized viral load data (see below) and intermediate files.
  - `fitting_scripts`: scripts used for fitting the viral kinetics model using MCMC and assessing convergence.
  - `functions`: all R files with functions used in the analyses.
  - `pars`: csv files containing data frames that control the MCMC framework.
  - `plots`: folder containing all intermediate plots, paper figures (for the part A related analyses) and *R scripts for figure generation*.
  - `sim_scripts`: scripts used to simulate the SEIR model and viral load trajectories used in part B.
  - `sims`: outputs from running `sim_scripts` scripts. *NOTE* these outputs are not included in the git repo, as they are huge files (tens of gigabytes). They should be reproduced if needed, see below.
  - `test_scripts`: redundant scripts used during the development stage.

### 3. Data
------------
Most data for this part are simulated. However, for model fitting, we used swab and sputum data from [Wölfel, R., Corman, V.M., Guggemos, W. et al. Virological assessment of hospitalized patients with COVID-2019. Nature 581, 465–469 (2020). https://doi-org.ezp-prod1.hul.harvard.edu/10.1038/s41586-020-2196-x](https://www-nature-com.ezp-prod1.hul.harvard.edu/articles/s41586-020-2196-x). Data were extracted using a [web plot digitizer](https://automeris.io/WebPlotDigitizer/). All extracted data are provided in [code/viral_kinetics/data](https://github.com/cleary-lab/covid19-group-tests/tree/master/code/viral_kinetics/data).

### 4. Scripts
------------
R scripts should be run in the following order:
  1. `fitting_scripts/fit_data_multivariate_hinge_swab.R`: reads in the Wölfel et al. swab data and fits the hinge viral kinetics model. Note that the working directories must be changed to your local file system, and there is a flag, `rerun_mcmc` that must be set to TRUE if you wish to rerun the whole thing. Otherwise, the pre-computed chains are used instead.
  2. `fitting_scripts/fit_data_multivariate_hinge_sputum.R`: reads in the Wölfel et al. sputum data and fits the hinge viral kinetics model. Exactly as in the first script, but using different data.
  3. `sim_scripts/simulate_data_hinge_swab_new.R`: uses the MCMC chains generated in script 1 to simulate viral loads from swab samples in a population undergoing a full SEIR epidemic. Again, working directories must be changed. It is also advisable to set `n` at the start of the script to something small, otherwise the script takes a few hours to run. This script will save simulated viral loads automatically to disk as it runs in the `sims` folder.
  4. `sim_scripts/simulate_data_hinge_sputum_new.R`: uses the MCMC chains generated in script 1 to simulate viral loads from sputum in a population undergoing a full SEIR epidemic. Exactly as in script 3, but using the sputum fits.
  5. `sim_scripts/simulate_data_hinge_swab_switch.R`: uses the MCMC chains generated in script 1 to simulate viral loads from swab in a population. Unlike script 3, this instead assumes that R0 changes during the epidemic: it begins at 2.5 drops to 0.8 (representing the implementation of interventions), and then returns to 1.5 later on (representing relaxation of interventions).
  6. `plots`: various figure generating scripts in this folder, named according to the figure they generate.


  


