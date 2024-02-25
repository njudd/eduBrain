#### ### ----------------- ### ####
# Nicholas Judd
# Donders Institute - LCD Lab
# 2024-02-14 https://njudd.com 

### Table of Contents
# 3.1 Packages & data from 1_ContinutyGlobal.R
# 3.2 1 month window analysis
# 3.3 1 month window plotting & diagnostics
# 3.4 6 month window analysis
# 3.5 6 month window plotting & diagnostics
# 3.6 misc
#### ### ----------------- ### ####


#### notes

# some covs are meaningless since they have no variance 
# like "summer"; "FLAIR"


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3.1 Loading & cleaning ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

if (!require(pacman)){
  install.packages('pacman')
}

pacman::p_load(tidyverse, lubridate, stringr, RDHonest, fastDummies, mice, ggseg, anytime, kableExtra,
               rstanarm, insight, bayestestR, furrr, bayesplot)

# grab functions
source("~/projects/roslaUKB/main_analysis/0_functions.R") 
# function: vec_to_fence() takes a vector and puts outliers to the fence (i.e., boxplot.stats)
# function: RDjacked() for cov corrected fuzzy RD

# function to fit a simple mod with covs 
# bayesFIT("SA", "ROSLA", b1_covs, b1) # outputs a list of models

# function to get the output I want
# bayes_to_results(list_of_models)

# 1) covariates, 2) covariates to fence & 3)cols with running_var
covs <- c("sex", "visit_day_correct", "visit_day_correct2", "t2_FLAIR" , "summer", "headmotion", #, "t2_FLAIR"* , "summer"
          "imaging_center_11026", "imaging_center_11027", "imaging_center_11028", 
          "dMRI_25922_1", "dMRI_25921_1", "dMRI_25928_1")
covs_fence <- c("dMRI_25922_1", "dMRI_25921_1", "dMRI_25928_1", "headmotion")
cols = c("EduAge16", "running_var", covs)

# loading data
fullset <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240223_fullset.csv")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3.2 1 month window analysis ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
b1 <- fullset[running_var %in% c(-1,0)][
  , ROSLA := running_var >= 0
][
  !is.na(visit_date) #must be an imaging subjects (visit date instance 2)
]

# table(b1$month); table(b1$ROSLA) #checking that it really is the right number and groups

# no variance for the following covs
b1_covs <- covs[!covs %in% c("summer", "t2_FLAIR", "imaging_center_11028")]

# first we need to establish that there are more likely in this narrow window to stay until 16
# First, we established that subjects within this one month window where more likely to stay in school, confirming our instrument.

first_stage <- stan_glm("EduAge16 ~ running_var", data = b1, iter = 40000)
hdi(first_stage)
first_stage_bf <- bayesfactor_parameters(first_stage, null = 0)
first_stage_bf$log_BF[2]

# rather than individually fitting 12 STAN models, we will paralize it with future_map2
# how map2 works
# x <- list(1, 1, 1)
# y <- list(10, 20, 30)
# map2(x, y, \(x, y) x + y)

DVs <- rep(c("SA", "CT", "CSF_norm", "TBV_norm", "WM_hyper", "wFA"), each = 2) #amazing 
IVs <- rep(c("ROSLA", "EduAge"), 6)

plan(multisession, workers = 6) # THIS WORKS
modB1 <- future_map2(DVs, IVs, \(y, x) bayesFIT(y, x, b1_covs, b1),
                   .options = furrr_options(seed = T))

# bayes_to_results(modB1) %>%
#   kbl(caption = "One Month Bandwidth Bayesian Analysis") %>%
#   kable_styling("hover", full_width = F) %>%
#   save_kable("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/1mBayes.html")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3.3 1 month window diagnostics and plotting ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# https://data-se.netlify.app/2021/02/06/plotting-multiple-plots-using-purrr-map-and-ggplot/
# how map2 works
# x <- list(1, 1, 1)
# y <- list(10, 20, 30)
# map2(x, y, \(x, y) x + y)

# first testing the covariates again in the 1mnth window...

b1CovTest <- b1_covs %>%
  future_map(~bayesFIT(., "ROSLA", covs = c(), b1), .options = furrr_options(seed = T)) %>% 
  future_map_dfr(~bayes_to_results(list(.))) %>%
  kbl(caption = "1 Month Bandwidth Covariates Bayesian Test") %>%
  kable_styling("hover", full_width = F) %>%
  save_kable("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/1mBayes_cov.html")


# plotting and saving trace plots
plt_name_trace <- str_c("Trace plot BayesLocal (1 mnth) of", IVs," on ", DVs)
plt_name_trace_s <- str_c("trace_", IVs, "_", DVs)

map2(modB1, plt_name_trace, 
     \(x, y) {mcmc_trace(x) + labs(title = y)}) %>% 
  map2(., plt_name_trace_s, 
       \(q, w) ggsave(paste0("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/plts/diags/", w, ".png"), q, bg = "white"))

# plotting and saving posterior draws
plt_name_draw <- str_c("Posterior draws BayesLocal (1 mnth) of", IVs," on ", DVs)
plt_name_draw_s <- str_c("draw_", IVs, "_", DVs)

map2(modB1, plt_name_draw, 
     \(x, y) {color_scheme_set("pink"); ppc_dens_overlay(y = x$y, yrep = posterior_predict(x, draws = 200)) + labs(title = y)}) %>% 
  map2(., plt_name_draw_s, 
       \(q, w) ggsave(paste0("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/plts/diags/", w, ".png"), q, bg = "white"))


############
# just manually go thru performance::check_model() because changing the base.size is impossible lol
# plt_name_perform_full <- str_c("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/plts/diags/performance_", IVs, "_", DVs, ".png")
# 
# map(modB1, ~performance::check_model(., theme = "ggplot2::theme_minimal(base.size =50)")) %>% 
#   map2(., plt_name_perform_full, 
#        \(q, w) {png(w, width = 3000, height = 1500); plot(q); dev.off()})

# https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html
# ^^^ go thru this!! 



b1SA <- stan_glm(paste0(c("SA ~ ROSLA", b1_covs), collapse=" + "), data = b1, iter = 40000)
# do ploy on scan date & see if anything changes as this orthogonlizes it!
# VERY HIGH VIFs

# plotting the posterior of the IV's with rainclouds...


# you need to raincloud EduYears vs ROSLA...

#### just manually extract them...

sa_post <- data.frame(iv = c(rep(names(get_parameters(modB1[[1]]))[2], dim(get_parameters(modB1[[1]])[2])[1]), rep(names(get_parameters(modB1[[2]]))[2], dim(get_parameters(modB1[[2]])[2])[1])),
                      posterior = c(get_parameters(modB1[[1]])[,2], get_parameters(modB1[[2]])[,2]))

ggplot(sa_post, aes(1, posterior, fill = iv, color = iv)) + 
  ggrain::geom_rain(alpha = .05, rain.side = 'l', boxplot.args.pos = list(
    position = ggpp::position_dodgenudge(x = .1, width = 0.1), width = 0.1)) +
  theme_classic() +
  scale_fill_brewer(palette = 'Dark2') +
  scale_color_brewer(palette = 'Dark2')


# old code that isn't great
# plt_path = "~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/plts/"
# modB1 %>% 
#   map(~data.frame(iv_posterior = get_parameters(.)[,2], 
#                   dv_name = rep(as.character(.$formula[2]), dim(get_parameters(.)[2])[1]),
#                   iv_name = rep(names(get_parameters(.))[2], dim(get_parameters(.)[2])[1]))) %>% 
#   map(~{ggplot(., aes(1, iv_posterior)) + 
#       ggrain::geom_rain(point.args = list(alpha = .05)) + 
#       theme_classic(base_size = 20) +
#       labs(y = str_c("Effect of ", unique(.$iv_name), " on ", unique(.$dv_name)))}) %>% 
#   map(~ggsave(str_c(plt_path, str_replace_all(.$labels$y, " ", "_"), ".png"), ., bg = "white"))
  


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3.4 6 month window analysis ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

b6 <- fullset[running_var %in% c(-6:5)][
  , ROSLA := running_var >= 0
][
  !is.na(visit_date) #must be an imaging subjects (visit date instance 2)
]

covs[!covs %in% b1_covs]

# table(b6$imaging_center_11028) # this is doable...

b6_covs <- c(b1_covs, "imaging_center_11028")

plan(multisession, workers = 6) # THIS WORKS
modB6 <- future_map2(DVs, IVs, \(y, x) bayesFIT(y, x, b6_covs, b6),
                     .options = furrr_options(seed = T))

modB6 %>% 
  future_map_dfr(~bayes_to_results(list(.))) %>%
  kbl(caption = "Six Month Bandwidth Bayesian Analysis") %>%
  kable_styling("hover", full_width = F) %>%
  save_kable("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/6mBayes.html")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3.4 6 month window diagnostics and plotting ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###






### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3.5 Misc ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


### fooling around to see when the age differnces between the groups start having an effect
# essentially the age difference between the groups have an effect because of the running variable

bMAX <- fullset[running_var %in% c(-65:64)][
  , ROSLA := running_var >= 0
]
bMaxSA <- stan_glm(CT ~ ROSLA, data = bMAX, iter = 40000)
bMaxSA_bf <- bayesfactor_parameters(bMaxSA, null = 0)
bMaxSA_bf$log_BF[2]






#  we will use the same covariates used in the local linear analysis yet they will be entered in the model.



# by moving the cutoff up a year (Sept. 1958), up two years (Sept. 1959), down a year (Sept. 1956), and down two years (Sept. 1955). 
# This doesn't make much sense if you don't have an effect...






# Secondly, we will analyze predetermined covariates and placebo outcomes similarly to the local-linear approach

c1 <- stan_glm(sex ~ ROSLA, data = b1, iter = 40000)
c1_bf <- bayesfactor_parameters(c1, null = 0); c1_bf$log_BF[2]

c2.1 <- stan_glm(visit_day_correct ~ ROSLA, data = b1, iter = 40000)
c2.1_bf <- bayesfactor_parameters(c2.1, null = 0); c2.1_bf$log_BF[2]
c2.2 <- stan_glm(visit_day_correct ~ ROSLA, data = b1, iter = 40000)
c2.2_bf <- bayesfactor_parameters(c2.2, null = 0); c2.2_bf$log_BF[2]

c3.1 <- stan_glm(imaging_center_11026 ~ ROSLA, data = b1, iter = 40000)
c3.1_bf <- bayesfactor_parameters(c3.1, null = 0); c3.1_bf$log_BF[2]
c3.2 <- stan_glm(imaging_center_11027 ~ ROSLA, data = b1, iter = 40000)
c3.2_bf <- bayesfactor_parameters(c3.2, null = 0); c3.2_bf$log_BF[2]
# c3.3 <- stan_glm(imaging_center_11028 ~ ROSLA, data = b1, iter = 40000)
# c3.3_bf <- bayesfactor_parameters(c3.3, null = 0); c3.3_bf$log_BF[2]

c4.1 <- stan_glm(dMRI_25922_1 ~ ROSLA, data = b1, iter = 40000)
c4.1_bf <- bayesfactor_parameters(c4.1, null = 0); c4.1_bf$log_BF[2]
c4.2 <- stan_glm(dMRI_25921_1 ~ ROSLA, data = b1, iter = 40000)
c4.2_bf <- bayesfactor_parameters(c4.2, null = 0); c4.2_bf$log_BF[2]
c4.3 <- stan_glm(dMRI_25928_1 ~ ROSLA, data = b1, iter = 40000)
c4.3_bf <- bayesfactor_parameters(c4.3, null = 0); c4.3_bf$log_BF[2]














