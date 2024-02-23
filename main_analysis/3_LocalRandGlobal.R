#### ### ----------------- ### ####
# Nicholas Judd
# Donders Institute - LCD Lab
# 2024-02-14 https://njudd.com 

### Table of Contents
# 3.1 Loading data from 1_ContinutyGlobal
# 3.2 
# 3.3 
# 3.4 
# 3.5 
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
source("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/R_scripts/0_functions.R") 
# function: vec_to_fence() takes a vector and puts outliers to the fence (i.e., boxplot.stats)
# function: RDjacked() for cov corrected fuzzy RD

# function to fit a simple mod with covs 
# bayesFIT("SA", "ROSLA", b1_covs, b1)

# 1) covariates, 2) covariates to fence & 3)cols with running_var
covs <- c("sex", "visit_day_correct", "visit_day_correct2", "t2_FLAIR" , "summer", "headmotion", #, "t2_FLAIR"* , "summer"
          "imaging_center_11026", "imaging_center_11027", "imaging_center_11028", 
          "dMRI_25922_1", "dMRI_25921_1", "dMRI_25928_1")
covs_fence <- c("dMRI_25922_1", "dMRI_25921_1", "dMRI_25928_1", "headmotion")
cols = c("EduAge16", "running_var", covs)

# loading data
fullset <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240223_fullset.csv")

# 1mnth window
b1 <- fullset[running_var %in% c(-1,0)][
  , ROSLA := running_var >= 0
][
  !is.na(visit_date) #must be an imaging subjects (visit date instance 2)
]
table(b1$month)

# no variance for the following covs
b1_covs <- covs[!covs %in% c("summer", "t2_FLAIR", "imaging_center_11028")]

# first we need to establish that there are more likely in this narrow window to stay until 16
# First, we established that subjects within this one month window where more likely to stay in school, confirming our instrument.

first_stage <- stan_glm("EduAge16 ~ running_var", data = b1, iter = 40000)
hdi(first_stage)
first_stage_bf <- bayesfactor_parameters(first_stage, null = 0)
first_stage_bf$log_BF[2]

# you now automate this whole process
# b1SA <- stan_glm(paste0(c("SA ~ ROSLA", b1_covs), collapse=" + "), data = b1, iter = 40000)
# b1CT <- stan_glm(paste0(c("CT ~ ROSLA", b1_covs), collapse=" + "), data = b1, iter = 40000)
# b1CSF <- stan_glm(paste0(c("CSF_norm ~ ROSLA", b1_covs), collapse=" + "), data = b1, iter = 40000)
# b1TBV <- stan_glm(paste0(c("TBV_norm ~ ROSLA", b1_covs), collapse=" + "), data = b1, iter = 40000)
# b1WMh <- stan_glm(paste0(c("WM_hyper ~ ROSLA", b1_covs), collapse=" + "), data = b1, iter = 40000)
# b1wFA <- stan_glm(paste0(c("wFA ~ ROSLA", b1_covs), collapse=" + "), data = b1, iter = 40000)
# b1SA_confound <- stan_glm(paste0(c("SA ~ EduAge", b1_covs), collapse=" + "), data = b1, iter = 40000)
# b1CT_confound <- stan_glm(paste0(c("CT ~ EduAge", b1_covs), collapse=" + "), data = b1, iter = 40000)
# b1CSF_confound <- stan_glm(paste0(c("CSF_norm ~ EduAge", b1_covs), collapse=" + "), data = b1, iter = 40000)
# b1TBV_confound <- stan_glm(paste0(c("TBV_norm ~ EduAge", b1_covs), collapse=" + "), data = b1, iter = 40000)
# b1WMh_confound <- stan_glm(paste0(c("WM_hyper ~ EduAge", b1_covs), collapse=" + "), data = b1, iter = 40000)
# b1wFA_confound <- stan_glm(paste0(c("wFA ~ EduAge", b1_covs), collapse=" + "), data = b1, iter = 40000)

# how map2 works
# x <- list(1, 1, 1)
# y <- list(10, 20, 30)
# map2(x, y, \(x, y) x + y)

IVs <- rep(c("SA", "CT", "CSF_norm", "TBV_norm", "WM_hyper", "wFA"), each = 2) #amazing 
DVs <- rep(c("ROSLA", "EduAge"), 6)

plan(multisession, workers = 6) # THIS WORKS
nkj <- future_map2(IVs, DVs, \(x, y) bayesFIT(x, y, b1_covs, b1),
                   .options = furrr_options(seed = T))

# checking diags

# just save them in a folder 

performance::check_model(nkj[[1]]) # has some issues...
performance::check_model(b1SA)






bayesplot::mcmc_trace(nkj[[1]])

color_scheme_set("brightblue")
mcmc_hist_by_chain(posterior, pars = c("wt", "sigma"))


nkj_results <- bayes_to_results(nkj)

##### THIS WORKS!!!!!!!!^^^^^


# Error: Error : 'data' must be a data frame.
# what a bitch...d




nkj %>% map(~str_split(.$formula, " ")[[1]][1]) 
nkj %>% map(~str_split(.$formula, " ")[[1]][3]) 


nkj%>% map(~hdi(.)$CI_low[2])


nkj %>% map(~mean(get_parameters(.)[,2]))


bayesfactor_parameters(get_parameters(nkj[[1]])[,2], prior = distribution_normal(40000, 0, 1))

bayesfactor_parameters(b1SA, null = 0)
bayesfactor_parameters(get_parameters(b1SA)[,2], prior = distribution_normal(80000, 0, 1), null = 0)

map(str_split(wFA_codes$Field, " "), `[`, -c(1:4))

ct %>% map(`[`) %>% map(~hdi(.)$CI_low[2])







bayesfactor_parameters(ct[[1]])


# b1SA_p <- get_parameters(b1SA)
# b1CT_p <- get_parameters(b1CT)
# b1CSF_p <- get_parameters(b1CSF)
# b1TBV_p <- get_parameters(b1TBV)
# b1WMh_p <- get_parameters(b1WMh)
# ggplot(m1SAp, aes(x = ROSLATRUE)) +
#   geom_density(fill = "red")

bayes_inference <- function(mod, param){
  
  posterior <- eval(parse(text=paste0("get_parameters(", r, ")$",v)))
  h <- hdi(posterior)
  
  bf <- eval(parse(text=paste0("bayesfactor_parameters(", r, ", null =0)")))
  
  df <- data.frame(model = mod, term = param, mean = mean(posterior), 
                   median = median(posterior), 
                   hdi_95 = paste0("[", round(h$CI_low, 2), ", ", round(h$CI_high, 2), "]"),
                   logBF = bf$log_BF[bf$Parameter == param])
  
  return(df)
}

g <- c("b1SA", "b1CT") %>% 
  future_map(~bayes_inference(., "ROSLATRUE")) 

g %>% list_rbind()
  
  


bayes_inference("b1SA", "ROSLATRUE")


b1SA_bf <- bayesfactor_parameters(b1SA, null = 0)
b1CT_bf <- bayesfactor_parameters(b1CT, null = 0)
b1CSF_bf <- bayesfactor_parameters(b1CSF, null = 0)
b1TBV_bf <- bayesfactor_parameters(b1TBV, null = 0)
b1WMh_bf <- bayesfactor_parameters(b1WMh, null = 0)

b1SA_bf$log_BF[2]; b1CT_bf$log_BF[2]; b1CSF_bf$log_BF[2]; b1TBV_bf$log_BF[2]; b1WMh_bf$log_BF[2]

b1SA_bf_CON <- bayesfactor_parameters(b1SA_confound, null = 0)
b1SA_bf_CON$log_BF[2]


# half-year window 6:5

b6 <- fullset[running_var %in% c(-6:5)][
  , ROSLA := running_var >= 0
]

b6SA <- stan_glm(SA ~ ROSLA, data = b6, iter = 40000)
b6CT <- stan_glm(CT ~ ROSLA, data = b6, iter = 40000)
b6CSF <- stan_glm(CSF_norm ~ ROSLA, data = b6, iter = 40000)
b6TBV <- stan_glm(TBV_norm ~ ROSLA, data = b6, iter = 40000)
b6WMh <- stan_glm(WM_hyper ~ ROSLA, data = b6, iter = 40000)

b6SA_bf <- bayesfactor_parameters(b6SA, null = 0)
b6CT_bf <- bayesfactor_parameters(b6CT, null = 0)
b6CSF_bf <- bayesfactor_parameters(b6CSF, null = 0)
b6TBV_bf <- bayesfactor_parameters(b6TBV, null = 0)
b6WMh_bf <- bayesfactor_parameters(b6WMh, null = 0)

b6SA_bf$log_BF[2]; b6CT_bf$log_BF[2]; b6CSF_bf$log_BF[2]; b6TBV_bf$log_BF[2]; b6WMh_bf$log_BF[2]


# so it didn't work with a month bracket, yet does with half a year
b6SA_confound <- stan_glm(SA ~ EduAge, data = b6, iter = 40000)
b6SA_bf_CON <- bayesfactor_parameters(b6SA_confound, null = 0)
b6SA_bf_CON$log_BF[2]






### fooling around to see when the age differnces between the groups start having an effect

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













