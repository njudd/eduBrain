# No effect of additional education on long-term brain structure – a preregistered natural experiment in thousands of individuals.
Dr. Nicholas Judd & Prof. Rogier Kievit

[Slides on the results & resources for regression discontinuity](https://www.njudd.com/causalinference/)

See the paper in [Elife](https://elifesciences.org/reviewed-preprints/101526)

See the preregistration [here](https://osf.io/rv38z)


**Recommendation:** I would take a look at [this directory](https://github.com/njudd/EduTelomere) it is a much cleaner implementation of testing the natural experiment (ROSLA) in the UK Biobank. This script uses an outdated package version of RDHonest before covariate correction was added.

![](https://www.njudd.com/assets/proj_imgs/rd/Fig1_4_take2.png)



## 0_functions.R

A collection of functions

1. `vec_to_fence` brings vectors to the boxplot fence
2. `RDjacked` outdated covariate correction for RDHonest as explained [here](https://github.com/kolesarm/RDHonest/issues/7)
3. `MICE_imp` hardcoded does imputation
4. `bayesFIT` fits a stan_glm model used for the map functions
5. `bayes_to_results` grabs results from a stan_glm model object

## 1_ContinuityGlobal.R

All plots and analyses for continuity-based regression discontinuity using the RDHonest package. Calls `RDjacked` to do manual covariate correction, it is recommended to use the latest version of the package. See the directory above for a streamlined script.

## 2_ContinuityROI.R

Does continuity-based regression discontinuity (and plots) for all the neuroimaging regions of interest (e.g., free surfer atlas). Heavily uses [map](https://purrr.tidyverse.org/reference/map.html) to iterate through regions.

## 3_LocalRandGlobal.R

Uses Bayesian local randomization regression discontinuity. Does this by subseting the data to include subjects only 1-months (m1) and 5-months (m6) around the ROSLA cutoff and fitting bayesian models using the function `BayesFIT`. Also includes reviewer initiated additional analyses and plots.














 
