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
# map & map2 are used a lot for parallelization; they are pretty much for loops with lists
# map2 takes two lists (of equal lenght & itterates thru them in order)
# https://data-se.netlify.app/2021/02/06/plotting-multiple-plots-using-purrr-map-and-ggplot/
# how map2 works
# x <- list(1, 1, 1)
# y <- list(10, 20, 30)
# map2(x, y, \(x, y) x + y)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3.1 Loading & cleaning ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

if (!require(pacman)){
  install.packages('pacman')
}

pacman::p_load(tidyverse, lubridate, stringr, RDHonest, fastDummies, mice, ggseg, anytime, kableExtra,
               rstanarm, insight, bayestestR, furrr, bayesplot, ggseg, ggsegTracula)

# Enable this universe
# options(repos = c(
#   ggseg = 'https://ggseg.r-universe.dev',
#   CRAN = 'https://cloud.r-project.org'))
# install.packages('ggsegTracula')

plan(multisession, workers = 6) # THIS WORKS
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

# "we will use the same covariates used in the local linear analysis yet they will be entered in the model"
# was not possible for the following covs due to no variance
b1_covs <- covs[!covs %in% c("summer", "t2_FLAIR", "imaging_center_11028")]

# first we need to establish that there are more likely in this narrow window to stay until 16
# First, we established that subjects within this one month window where more likely to stay in school, confirming our instrument.

# no priors on this atm
# first_stage <- stan_glm("EduAge16 ~ running_var", data = b1, iter = 40000)
# hdi(first_stage)
# first_stage_bf <- bayesfactor_parameters(first_stage, null = 0)
# first_stage_bf$log_BF[2]

# rather than individually fitting 12 STAN models, we will paralize it with future_map2
DVs <- rep(c("SA", "CT", "CSF_norm", "TBV_norm", "WM_hyper", "wFA"), each = 2) #amazing 
IVs <- rep(c("ROSLA", "EduAge"), 6)

# https://discourse.mc-stan.org/t/evidence-ratios-when-using-default-priors/4847
# never use the STAN defaults for BayesFactors

# The Bayes Factor package uses a Cauchy(0, .7) as default
# plot(x, dnorm(x, sd = .7),type="l", col = "blue")
# lines(x, dcauchy(x, scale = .7),type="l", col = "purple")
# lines(x, dnorm(x, sd = 1),type="l", col = "red")

# stan_glm can't specify priors to specific terms & I couldn't get scalling (e.g., autoscale) to work in brms

# Hi Markus,
# 
# For the Cauchy, the prior width r equals the interquartile range. So if r=0.707, 
# there is a 50% chance that the true value of effect size lies in the interval from -0.707 to +0.707. 
# I encourage you to read the papers on Bayesian inference (part I and II) here: https://osf.io/m6bi8/
#   
# Finally, people who argue for a lower value of r fail to understand that lowering the width r but keeping the location fixed at 0 has undesirable consequences, 
# as this just makes H1 very similar to H0. If you feel that your effect sizes are around .36, 
# you should have a prior that peaks around .36. Right now the only parameter to adjust is the width 
# (we are working on this). Bottom-line: just use the default.
# 
# Cheers,
# E.J.


# 1 SD scalled priors
modB1_1sd <- future_map2(DVs, IVs, \(y, x) bayesFIT(y, x, b1_covs, b1,BF_prior = normal(location = 0, scale = 1, autoscale = TRUE)), .options = furrr_options(seed = T))
m1Bayes_1SD <- modB1_1sd %>%
  future_map_dfr(~bayes_to_results(list(.)), .options = furrr_options(seed = T))
# .5 SD # prior_summary(modB1_.5sd[[1]])
modB1_.5sd <- future_map2(DVs, IVs, \(y, x) bayesFIT(y, x, b1_covs, b1, BF_prior = normal(location = 0, scale = .5, autoscale = TRUE)), .options = furrr_options(seed = T))
m1Bayes_.5SD <- modB1_.5sd %>%
  future_map_dfr(~bayes_to_results(list(.)), .options = furrr_options(seed = T)) 
# 1.5 SD
modB1_1.5sd <- future_map2(DVs, IVs, \(y, x) bayesFIT(y, x, b1_covs, b1, BF_prior = normal(location = 0, scale = 1.5, autoscale = TRUE)), .options = furrr_options(seed = T))
m1Bayes_1.5SD <- modB1_1.5sd %>%
  future_map_dfr(~bayes_to_results(list(.)), .options = furrr_options(seed = T)) 

m1Bayes_1SD <- m1Bayes_1SD[,c(1,2,3,5,8,9,7,6)]

colnames(m1Bayes_1SD)[c(8,7)] <- c("logBF_1", "BF_1")
# colnames(m1Bayes_1SD)[c(7,8)] <- c("BF_1SD", "logBF_1SD")
m1Bayes_1SD$BF_0.5 <- m1Bayes_.5SD$BF
m1Bayes_1SD$BF_1.5 <- m1Bayes_1.5SD$BF

# credible intervals; getting rid of low, high & puting it in the right spot
m1Bayes_1SD$ci_low <- stringi::stri_c("(", m1Bayes_1SD$ci_low, ", ", m1Bayes_1SD$ci_high, ")")
m1Bayes_1SD <- m1Bayes_1SD[,-6] # getting rid of high
colnames(m1Bayes_1SD)[5] <- "CI" # renaming low

# m1Bayes_1SD %>% 
#   kbl(caption = "One Month Bandwidth Bayesian Analysis") %>%
#   kable_styling(full_width = F) %>% 
#   save_kable("~/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/results/localRand/20240414_1mBayes_correctPriors.html")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3.3 1 month window diagnostics and plotting ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# first testing the covariates again in the 1mnth window...
# "we will analyze predetermined covariates and placebo outcomes similarly to the local-linear approach"

# missing one of the centers for the covariate test, since it is dummy coded
# b1$imaging_center_11025 <- b1$imaging_center == 11025
# b1_covsT <- c(b1_covs, "imaging_center_11025")
# b1CovTest <- b1_covsT %>%
#   future_map(~bayesFIT(., "ROSLA", covs = c(), b1, BF_prior = normal(location = 0, scale = 1, autoscale = TRUE)), .options = furrr_options(seed = T)) %>%
#   future_map_dfr(~bayes_to_results(list(.))) %>%
#   kbl(caption = "1 Month Bandwidth Covariates Bayesian Test") %>%
#   kable_styling("hover", full_width = F) %>%
#   save_kable("~/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/results/localRand/20240414_1mBayes_cov_Prior1SD.html")


# by moving the cutoff up a year (Sept. 1958), up two years (Sept. 1959), down a year (Sept. 1956), and down two years (Sept. 1955). 
# This doesn't make much sense if you don't have an effect... so we will not do it

# plotting and saving trace plots
plt_name_trace <- str_c("Trace plot BayesLocal (1 mnth) of", IVs," on ", DVs)
plt_name_trace_s <- str_c("trace_", IVs, "_", DVs)

# map2(modB1, plt_name_trace, 
#      \(x, y) {mcmc_trace(x) + labs(title = y)}) %>% 
#   map2(., plt_name_trace_s, 
#        \(q, w) ggsave(paste0("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/plts/diags/", w, ".png"), q, bg = "white"))

# plotting and saving posterior draws
plt_name_draw <- str_c("Posterior draws BayesLocal (1 mnth) of", IVs," on ", DVs)
plt_name_draw_s <- str_c("draw_", IVs, "_", DVs)

# map2(modB1, plt_name_draw, 
#      \(x, y) {color_scheme_set("pink"); ppc_dens_overlay(y = x$y, yrep = posterior_predict(x, draws = 200)) + labs(title = y)}) %>% 
#   map2(., plt_name_draw_s, 
#        \(q, w) ggsave(paste0("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/plts/diags/", w, ".png"), q, bg = "white"))


############
# just manually go thru performance::check_model() because changing the base.size is impossible lol
# plt_name_perform_full <- str_c("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/plts/diags/performance_", IVs, "_", DVs, ".png")
# 
# map(modB1, ~performance::check_model(., theme = "ggplot2::theme_minimal(base.size =50)")) %>% 
#   map2(., plt_name_perform_full, 
#        \(q, w) {png(w, width = 3000, height = 1500); plot(q); dev.off()})

# https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html
# ^^^ go thru this!! 

# do ploy on scan date & see if anything changes as this orthogonlizes it!
# VERY HIGH VIFs


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3.4 6 month window analysis ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

b6 <- fullset[running_var %in% c(-5:4)][
  , ROSLA := running_var >= 0
][
  !is.na(visit_date) #must be an imaging subjects (visit date instance 2)
]

# first_stage6 <- stan_glm("EduAge16 ~ running_var", data = b6, iter = 40000)
# this first stage should be the dummy coding (to be similar to first stage 1mnth)
# hdi(first_stage6)
# first_stage_bf6 <- bayesfactor_parameters(first_stage6, null = 0)
# first_stage_bf6$log_BF[2]

# taking as many covs as possible from the cont. approach
covs[!covs %in% b1_covs]
# table(b6$imaging_center_11028) # this is doable...
b6_covs <- c(b1_covs, "imaging_center_11028")

# 1 SD scalled priors
m6_1sd <- future_map2(DVs, IVs, \(y, x) bayesFIT(y, x, b6_covs, b6, BF_prior = normal(location = 0, scale = 1, autoscale = TRUE)), .options = furrr_options(seed = T))
m6Bayes_1SD <- m6_1sd %>%
  future_map_dfr(~bayes_to_results(list(.)), .options = furrr_options(seed = T)) 
# .5 SD # prior_summary(modB1_.5sd[[1]])
m6_.5sd <- future_map2(DVs, IVs, \(y, x) bayesFIT(y, x, b6_covs, b6, BF_prior = normal(location = 0, scale = .5, autoscale = TRUE)), .options = furrr_options(seed = T))
m6Bayes_.5SD <- m6_.5sd %>%
  future_map_dfr(~bayes_to_results(list(.)), .options = furrr_options(seed = T)) 
# 1.5 SD
m6_1.5sd <- future_map2(DVs, IVs, \(y, x) bayesFIT(y, x, b6_covs, b6, BF_prior = normal(location = 0, scale = 1.5, autoscale = TRUE)), .options = furrr_options(seed = T))
m6Bayes_1.5SD <- m6_1.5sd %>%
  future_map_dfr(~bayes_to_results(list(.)), .options = furrr_options(seed = T))

# making the table
m6Bayes_1SD <- m6Bayes_1SD[,c(1,2,3,5,8,9,6,7)]
colnames(m6Bayes_1SD)[c(7,8)] <- c("logBF_1", "BF_1")
m6Bayes_1SD$BF_0.5 <- m6Bayes_.5SD$BF
m6Bayes_1SD$BF_1.5 <- m6Bayes_1.5SD$BF

# credible intervals; getting rid of low, high & puting it in the right spot
m6Bayes_1SD$ci_low <- stringi::stri_c("(", m6Bayes_1SD$ci_low, ", ", m6Bayes_1SD$ci_high, ")")
m6Bayes_1SD <- m6Bayes_1SD[,-6] # getting rid of high
colnames(m6Bayes_1SD)[5] <- "CI" # renaming low
#putting the right order
m6Bayes_1SD <- m6Bayes_1SD[, c(1:6,8,7,9)]

m6Bayes_1SD %>% 
  kbl(caption = "Five Month Bandwidth Bayesian Analysis") %>%
  kable_styling(full_width = F) %>% 
  save_kable("~/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/results/localRand/20240414_5mBayes_correctPriors.html")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3.4 6 month window diagnostics and plotting ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# first testing the covariates again in the 6mnth window...
# "we will analyze predetermined covariates and placebo outcomes similarly to the local-linear approach"

# b6$imaging_center_11025 <- b6$imaging_center == 11025
# b6_covsT <- c(b6_covs, "imaging_center_11025")

# b6_covsT %>%
#   future_map(~bayesFIT(., "ROSLA", covs = c(), b6, BF_prior = normal(location = 0, scale = 1, autoscale = TRUE)), .options = furrr_options(seed = T)) %>%
#   future_map_dfr(~bayes_to_results(list(.)), .options = furrr_options(seed = T)) %>%
#   kbl(caption = "5 Month Bandwidth Covariates Bayesian Test") %>%
#   kable_styling("hover", full_width = F) %>%
#   save_kable("~/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/results/localRand/20240414_5mBayes_cov_Prior1SD.html")

# again no point in moving the window up or down a year since there are no sig results to disprove/test

# plotting posterior Total Surface area & CSF_norm for 5 month Local Rand

# grab posterior ROSLA total SA; # grab posterior EduYear SA;

# sa6_posterior <- data.frame(posterior = c(get_parameters(modB6[[1]])[,2], get_parameters(modB6[[2]])[,2]),
#                             effect = c(rep("Causal", 160000), rep("Observational", 160000)))
# 
# ggplot(sa6_posterior, aes(1, posterior, group = effect)) +
#   ggrain::geom_rain()

#  grab posterior ROSLA CSF_norm; # grab posterior EduYear CSF_norm

# csf6_posterior <- data.frame(posterior = c(get_parameters(modB6[[5]])[,2], get_parameters(modB6[[6]])[,2]),
#                             effect = c(rep("Causal", 160000), rep("Observational", 160000)))
# 
# ggplot(csf6_posterior, aes(posterior, group = effect)) +
#   geom_density()


#### trying a plot to show evidence
# you just want to show assosiational vs causal BFs

# draw arrows of the effect changing 
hold <- m6Bayes_1SD

# rename legend to causal & correlational
m6Bayes_1SD$X[m6Bayes_1SD$X == "ROSLA"] <- "Causal"
m6Bayes_1SD$X[m6Bayes_1SD$X == "EduAge"] <- "Correlational"

# rename and order labels corectly
m6Bayes_1SD$Y[m6Bayes_1SD$Y == "CSF_norm"] <- "CSF"
m6Bayes_1SD$Y[m6Bayes_1SD$Y == "TBV_norm"] <- "TBV"
m6Bayes_1SD$Y[m6Bayes_1SD$Y == "WM_hyper"] <- "WMh"

m6Bayes_1SD$Y <- factor(m6Bayes_1SD$Y, levels = c("CSF", "SA", "wFA", "TBV", "WMh", "CT"))

# think about an ROI plot on the brain. than also do a dot plot connected moving?

# shows how to do subscripts
# axis_labels <- c(expression(paste('100\nBF'["01"])), expression("30\nBF"["01"]), expression("10\nBF"["01"]), expression("1\nBF"["01"]), expression("0\nBF"[" "]), 
#                  expression("1\nBF"["10"]), expression("10\nBF"["10"]), expression("30\nBF"["10"]),expression("100\nBF"["10"]))

plt2 <- ggplot(m6Bayes_1SD, aes(Y, logBF, shape = X))  +
  geom_point(size = 4, color = "NA") +
  ylim(-5.5, 5.5) +
  theme_classic(base_size = 20) +
  annotate("rect", xmin = .5, xmax = 6.5, ymin = 4.6, ymax = 5.5, alpha = .8, fill = heatmaply::RdBu(10)[1]) + # extreme evidence
  annotate("rect", xmin = .5, xmax = 6.5, ymin = 3.4, ymax = 4.6, alpha = .8, fill = heatmaply::RdBu(10)[2]) + # very strong
  annotate("rect", xmin = .5, xmax = 6.5, ymin = 2.3, ymax = 3.4, alpha = .8, fill = heatmaply::RdBu(10)[3]) + # strong
  annotate("rect", xmin = .5, xmax = 6.5, ymin = 1.1, ymax = 2.3, alpha = .8, fill = heatmaply::RdBu(10)[4]) + # substantial
  annotate("rect", xmin = .5, xmax = 6.5, ymin = 1, ymax = 1.1, alpha = .8, fill = heatmaply::RdBu(10)[5]) + # anecdotal
  annotate("rect", xmin = .5, xmax = 6.5, ymin = -1, ymax = 1, alpha = .8, fill = "white") + # no evidence either way
  annotate("rect", xmin = .5, xmax = 6.5, ymin = -1, ymax = -1.1, alpha = .8, fill = heatmaply::RdBu(10)[6]) +
  annotate("rect", xmin = .5, xmax = 6.5, ymin = -1.1, ymax = -2.3, alpha = .8, fill = heatmaply::RdBu(10)[7]) +
  annotate("rect", xmin = .5, xmax = 6.5, ymin = -2.3, ymax = -3.4, alpha = .8, fill = heatmaply::RdBu(10)[8]) +
  annotate("rect", xmin = .5, xmax = 6.5, ymin = -3.4, ymax = -4.6, alpha = .8, fill = heatmaply::RdBu(10)[9]) +
  annotate("rect", xmin = .5, xmax = 6.5, ymin = -4.6, ymax = -5.5, alpha = .8, fill = heatmaply::RdBu(10)[10]) +
  labs(x = "", y = "Bayes Factors", shape = "Parameter") + # 
  geom_point(size = 4) + # a hack to put them on top
  scale_y_continuous(breaks=c(-4.6, -3.4, -2.3,-1, 0, 1, 2.3,3.4,4.6), labels = c("100", "30", "10", "1", "0", "1", "10", "30", "100")) + 
  scale_shape_manual(values = c(16,17)) + # http://www.sthda.com/english/wiki/ggplot2-point-shapes
  theme(axis.line=element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(color = "black"),
        # axis.text.x = element_text(angle= -30, vjust = 1.5, hjust=.5), # angle = 45,
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        legend.position = c(.954, .861),
        legend.box.just = "center",
        legend.justification = c("right", "bottom"),
        legend.margin = margin(0.5,6, 0.5, 1),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', linewidth = .3)) +
  coord_flip() 

plt2

ggsave("~/Google Drive/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/results/plts/Fig2.png", 
       plt2, width = 7, height = 5, bg = "white")
  


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


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 3.5 Figure 3 NOT PreReged whole brain Bayes ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# not taking any inference from this it is just meant to help visualize the null result

b6; b6_covs

saROI <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240222_SAroi.csv")[
  eid %in% b6$eid
][
  b6[, .(eid, headmotion, ROSLA)], on = "eid" # was missing from SA_roi
]

f3_SA_DVs <- colnames(saROI)[as.logical(str_detect(colnames(saROI), "lh_") + str_detect(colnames(saROI), "rh_"))]
f3_SA_IVs <- rep("ROSLA", length(f3_SA_DVs))

f3_SAroi <- future_map2(f3_SA_DVs, f3_SA_IVs, \(y, x) bayesFIT(y, x, b6_covs, saROI, BF_prior = normal(0, 1, autoscale = TRUE)),
                     .options = furrr_options(seed = T))

# was crashing my comp with future_map_dfr & ~bayes_to_results
roi_vec <- c(); bf_vec <- c();rhat_vec <- c() # empty vectors for loop
for(m in 1:length(f3_SAroi)){
  roi_vec <- c(roi_vec, as.character(f3_SAroi[[m]]$formula[2]))
  bf_vec <- c(bf_vec, round(bayesfactor_parameters(f3_SAroi[[m]], null = 0)$log_BF[2], 2))
  rhat_vec <- c(rhat_vec, as.numeric(rhat(f3_SAroi[[m]])[2]))
} # there should be a way to paralize without running out of memory...

f3_SA_BFs <- tibble(
  label = roi_vec,
  bf = bf_vec)
# this took FOREVER!!!
# data.table::fwrite(f3_SA_BFs, "/Volumes/home/lifespan/nicjud/UKB/proc/20240415_SAroiBFs.csv")

f3_SA_BFs <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240415_SAroiBFs.csv")

f3_SA_BFs$bf_log <- f3_SA_BFs$bf
f3_SA_BFs$bf[f3_SA_BFs$bf < 0] <- -1/exp(f3_SA_BFs$bf[f3_SA_BFs$bf < 0])
f3_SA_BFs$bf[f3_SA_BFs$bf > 0] <- exp(f3_SA_BFs$bf[f3_SA_BFs$bf > 0])

# now for the brain plt'n

# this is for a categorical brain plt
# f3_SA_BFs$logBF <- rep(NA, length(f3_SA_BFs$bf_log))
# f3_SA_BFs$logBF[f3_SA_BFs$bf <= -4.6] <- "H0 extreme"
# f3_SA_BFs$logBF[f3_SA_BFs$bf <= -3.4 & f3_SA_BFs$bf >= -4.6] <- "H0 very strong"
# f3_SA_BFs$logBF[f3_SA_BFs$bf <= -2.3 & f3_SA_BFs$bf >= -3.4] <- "H0 strong"
# f3_SA_BFs$logBF[f3_SA_BFs$bf <= -1.1 & f3_SA_BFs$bf >= -2.3] <- "H0 substantial"
# f3_SA_BFs$logBF[f3_SA_BFs$bf <= -1 & f3_SA_BFs$bf >= -1.1] <- "H0 anecdotal"
# f3_SA_BFs$logBF[f3_SA_BFs$bf > -1] <- "No evidence"
# 
# f3_SA_BFs$logBF <- factor(f3_SA_BFs$logBF, levels = c("H1 extreme", "H1 very strong", "H1 strong", "H1 substantial", "H1 anecdotal", "no evidence",
#                                                       "H0 anecdotal", "H0 substantial", "H0 strong", "H0 very strong", "H0 extreme"))

# working on the contineous color scale
#https://stackoverflow.com/questions/58698075/assign-specific-color-to-definite-value-in-bar-plot-using-scale-fill-gradientn/58699545#58699545
# https://stackoverflow.com/questions/67082430/how-to-highlight-a-single-value-on-a-continuous-color-scale-in-ggplot2

plt_f3_SA_cont <- f3_SA_BFs %>%
  ggplot() +
  geom_brain(atlas = dk, color = "black",
             position = position_brain(side ~ hemi),
             aes(fill = bf, color = bf),
             show.legend = TRUE) +
  theme_void() +
  theme(legend.position="bottom",legend.direction = "horizontal", legend.box = "horizontal",
        legend.key.width = unit(2, "cm"),
        text = element_text(family = 'Arial')) +
  scale_fill_gradient2(low = heatmaply::RdBu(10)[9], mid = "white", high = heatmaply::RdBu(10)[2],
                       limits = c(-40,40), breaks=c(-40,-30, -20, -10, 0,  10, 20, 30, 40), labels=c('40','30', '20', '10', '0',  '10', '20', '30', '40'),
                       na.value = "white", name = "Surface Area Bayes Factors",
                       guide = guide_colorbar(title.position = "top", label.position = "bottom", 
                                              title.hjust = 0.5)) # centeres the color bar with text

  
# ggsave("~/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/results/plts/SI_plt3_SA_cont.png",
#        plt_f3_SA_cont, bg = "white", width = 8, height = 5)



##### STOPPING NOT RUNNING FOR OTHER MODALITIES



# now for the brain plt'n
# https://github.com/tidyverse/ggplot2/issues/5728
# plt_f3_SA <- f3_SA_BFs %>%
#   ggplot() +
#   geom_brain(atlas = dk, colour = "black",
#              position = position_brain(side ~ hemi),
#              aes(fill = logBF),
#              show.legend = TRUE) +
#   theme_void() +
#   scale_fill_manual(values = c(heatmaply::RdBu(10)[1], heatmaply::RdBu(10)[2], heatmaply::RdBu(10)[3], heatmaply::RdBu(10)[4], heatmaply::RdBu(10)[5],
#                                "white",
#                                heatmaply::RdBu(10)[6], heatmaply::RdBu(10)[7], heatmaply::RdBu(10)[8], heatmaply::RdBu(10)[9], heatmaply::RdBu(10)[10]),
#                     name="Strength of Evidence \nin Bayes Factors", labels=c("H1 extreme", "H1 very strong", "H1 strong", "H1 substantial", "H1 anecdotal", "No evidence",
#                                                                              "H0 anecdotal", "H0 substantial", "H0 strong", "H0 very strong", "H0 extreme"), 
#                     drop=FALSE,
#                     na.value = "white", na.translate = F)
# 
# ggsave("~/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/results/plts/SI_plt3_SA.png",
#        plt_f3_SA, bg = "white")

# now for CT
rm(list = c("saROI", "f3_SAroi"))

ctROI <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240222_CTroi.csv")[
  eid %in% b6$eid
][
  b6[, .(eid, headmotion, ROSLA)], on = "eid" # was missing from SA_roi
]

f3_CT_DVs <- colnames(ctROI)[as.logical(str_detect(colnames(ctROI), "lh_") + str_detect(colnames(ctROI), "rh_"))]
f3_CT_IVs <- rep("ROSLA", length(f3_CT_DVs))

f3_CTroi <- future_map2(f3_CT_DVs, f3_CT_IVs, \(y, x) bayesFIT(y, x, b6_covs, ctROI, BF_prior = normal(0, 1, autoscale = TRUE)),
                        .options = furrr_options(seed = T))

# was crashing my comp with future_map_dfr & ~bayes_to_results
roi_vec <- c(); bf_vec <- c(); rhat_vec <- c() # empty vectors for loop
for(m in 1:length(f3_CTroi)){
  roi_vec <- c(roi_vec, as.character(f3_CTroi[[m]]$formula[2]))
  bf_vec <- c(bf_vec, round(bayesfactor_parameters(f3_CTroi[[m]], null = 0)$log_BF[2], 2))
  rhat_vec <- c(rhat_vec, as.numeric(rhat(f3_CTroi[[m]])[2]))
} # there should be a way to paralize without running out of memory...

f3_CT_BFs <- tibble(
  label = roi_vec,
  bf = bf_vec)
# this took FOREVER!!!
# data.table::fwrite(f3_CT_BFs, "/Volumes/home/lifespan/nicjud/UKB/proc/20240417_CTroiBFs.csv")
# f3_CT_BFs <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240417_CTroiBFs.csv")


f3_CT_BFs$bf_log <- f3_CT_BFs$bf
f3_CT_BFs$bf[f3_CT_BFs$bf < 0] <- -1/exp(f3_CT_BFs$bf[f3_CT_BFs$bf < 0])
f3_CT_BFs$bf[f3_CT_BFs$bf > 0] <- exp(f3_CT_BFs$bf[f3_CT_BFs$bf > 0])

# making a contineous brain plt for CT 
plt_f3_CT <- f3_CT_BFs %>%
  ggplot() +
  geom_brain(atlas = dk, color = "black",
             position = position_brain(side ~ hemi),
             aes(fill = bf, color = bf),
             show.legend = TRUE) +
  theme_void() +
  theme(legend.position="bottom",legend.direction = "horizontal", legend.box = "horizontal",
        legend.key.width = unit(2, "cm"),
        text = element_text(family = 'Arial')) +
  scale_fill_gradient2(low = heatmaply::RdBu(10)[9], mid = "white", high = heatmaply::RdBu(10)[2],
                       limits = c(-40,40), breaks=c(-40,-30, -20, -10, 0,  10, 20, 30, 40), labels=c('40','30', '20', '10', '0',  '10', '20', '30', '40'),
                       na.value = "white", name = "Mean Cortical Thickness Bayes Factors",
                       guide = guide_colorbar(title.position = "top", label.position = "bottom", 
                                              title.hjust = 0.5)) # centeres the color bar with text

# ggsave("~/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/results/plts/SI_plt2_CT.png",
#        plt_f3_CT, bg = "white")



# categorical plotting...
# https://github.com/tidyverse/ggplot2/issues/5728
# f3_CT_BFs$logBF <- rep(NA, length(f3_CT_BFs$bf))
# f3_CT_BFs$logBF[f3_CT_BFs$bf <= -4.6] <- "H0 extreme"
# f3_CT_BFs$logBF[f3_CT_BFs$bf <= -3.4 & f3_CT_BFs$bf >= -4.6] <- "H0 very strong"
# f3_CT_BFs$logBF[f3_CT_BFs$bf <= -2.3 & f3_CT_BFs$bf >= -3.4] <- "H0 strong"
# f3_CT_BFs$logBF[f3_CT_BFs$bf <= -1.1 & f3_CT_BFs$bf >= -2.3] <- "H0 substantial"
# f3_CT_BFs$logBF[f3_CT_BFs$bf <= -1 & f3_CT_BFs$bf >= -1.1] <- "H0 anecdotal"
# f3_CT_BFs$logBF[f3_CT_BFs$bf > -1] <- "No evidence"
# 
# f3_CT_BFs$logBF <- factor(f3_CT_BFs$logBF, levels = c("H1 extreme", "H1 very strong", "H1 strong", "H1 substantial", "H1 anecdotal", "no evidence",
#                                                       "H0 anecdotal", "H0 substantial", "H0 strong", "H0 very strong", "H0 extreme"))
# 
# # now for the brain plt'n
# plt_f3_CT <- f3_CT_BFs %>%
#   ggplot() +
#   geom_brain(atlas = dk, colour = "black",
#              position = position_brain(side ~ hemi),
#              aes(fill = logBF),
#              show.legend = TRUE) +
#   theme_void() +
#   scale_fill_manual(values = c(heatmaply::RdBu(10)[1], heatmaply::RdBu(10)[2], heatmaply::RdBu(10)[3], heatmaply::RdBu(10)[4], heatmaply::RdBu(10)[5],
#                                "white",
#                                heatmaply::RdBu(10)[6], heatmaply::RdBu(10)[7], heatmaply::RdBu(10)[8], heatmaply::RdBu(10)[9], heatmaply::RdBu(10)[10]),
#                     name="Strength of Evidence \nin Bayes Factors", labels=c("H1 extreme", "H1 very strong", "H1 strong", "H1 substantial", "H1 anecdotal", "No evidence",
#                                             "H0 anecdotal", "H0 substantial", "H0 strong", "H0 very strong", "H0 extreme"), 
#                     drop=FALSE,
#                     na.value = "white", na.translate = F)
# 
# ggsave("~/My Drive/Assembled Chaos/10 Projects/10.02 ROSLA UK BioBank/results/plts/Plt3_CT.png",
#        plt_f3_CT, bg = "white")

# now for FA
rm(list = c("ctROI", "f3_CTroi"))

faROI <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240222_wFAGLOBAL.csv")[
  eid %in% b6$eid
][
  b6[, .(eid, headmotion, ROSLA)], on = "eid" # was missing from SA_roi
]

f3_FA_DVs <- colnames(faROI)[34:60]
f3_FA_IVs <- rep("ROSLA", length(f3_FA_DVs))

f3_FAroi <- future_map2(f3_FA_DVs, f3_FA_IVs, \(y, x) bayesFIT(y, x, b6_covs, faROI),
                        .options = furrr_options(seed = T))

# was crashing my comp with future_map_dfr & ~bayes_to_results
roi_vec <- c(); bf_vec <- c(); rhat_vec <- c() # empty vectors for loop
for(m in 1:length(f3_FAroi)){
  roi_vec <- c(roi_vec, as.character(f3_FAroi[[m]]$formula[2]))
  bf_vec <- c(bf_vec, round(bayesfactor_parameters(f3_FAroi[[m]], null = 0)$log_BF[2], 2))
  rhat_vec <- c(rhat_vec, as.numeric(rhat(f3_FAroi[[m]])[2]))
} # there should be a way to paralize without running out of memory...

f3_FA_BFs <- tibble(
  label = roi_vec,
  bf = bf_vec)
# this took FOREVER!!!
data.table::fwrite(f3_FA_BFs, "/Volumes/home/lifespan/nicjud/UKB/proc/20240405_FAroiBFs.csv")


# it is not this atlas so just ignore for now...
# as.data.frame(tracula$data)[,1:4]
# it is named differently....also dim(tracula$data) is 38...








# f3_FA_BFs %>% 
#   ggplot() +
#   geom_brain(atlas = tracula, 
#              #position = position_brain(hemi ~ side),
#              aes(fill = logBF),
#              show.legend = TRUE)
# 
# 
# plot(tracula) +
#   theme(legend.position = "bottom", 
#         legend.text = element_text(size = 9)) +
#   guides(fill = guide_legend(ncol = 3))



# now for subcortical
rm(list = c("faROI", "f3_FAroi"))

sROI <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240222_sROI.csv")[
  eid %in% b6$eid
][
  b6[, .(eid, headmotion, ROSLA)], on = "eid" # was missing from SA_roi
]

f3_sub_DVs <- colnames(sROI)[as.logical(str_detect(colnames(sROI), "lh_") + str_detect(colnames(sROI), "rh_"))]
f3_sub_IVs <- rep("ROSLA", length(f3_sub_DVs))

f3_SUBroi <- future_map2(f3_sub_DVs, f3_sub_IVs, \(y, x) bayesFIT(y, x, b6_covs, sROI),
                        .options = furrr_options(seed = T))

# was crashing my comp with future_map_dfr & ~bayes_to_results
roi_vec <- c(); bf_vec <- c(); rhat_vec <- c() # empty vectors for loop
for(m in 1:length(f3_SUBroi)){
  roi_vec <- c(roi_vec, as.character(f3_SUBroi[[m]]$formula[2]))
  bf_vec <- c(bf_vec, round(bayesfactor_parameters(f3_SUBroi[[m]], null = 0)$log_BF[2], 2))
  rhat_vec <- c(rhat_vec, as.numeric(rhat(f3_SUBroi[[m]])[2]))
} # there should be a way to paralize without running out of memory...

f3_SUB_BFs <- tibble(
  label = roi_vec,
  bf = bf_vec)
# this took FOREVER!!!
data.table::fwrite(f3_SUB_BFs, "/Volumes/home/lifespan/nicjud/UKB/proc/20240405_SUBroiBFs.csv")

# making a df of the atlas to edit
aseg_edit <- aseg$data
# subsetting relevant regions
aseg_edit <- aseg_edit[aseg_edit$hemi %in% c("left", "right"),]
aseg_edit <- aseg_edit[!is.na(aseg_edit$label),]

# trying to make a col similar to f3_SUB_BFs to merge by
aseg_edit$hemi[aseg_edit$hemi == "left"] <- rep("lh_", length(aseg_edit$hemi[aseg_edit$hemi == "left"]))
aseg_edit$hemi[aseg_edit$hemi == "right"] <- rep("rh_", length(aseg_edit$hemi[aseg_edit$hemi == "right"]))


aseg_edit$hemi <- str_c(aseg_edit$hemi, str_remove(aseg_edit$region, " "))

data.frame(aseg_edit$hemi, tolower(f3_SUB_BFs$label[-c(18,9)]))

aseg_edit$hemi[!aseg_edit$hemi %in% tolower(f3_SUB_BFs$label[-c(18,9)])]


f3_SUB_BFs$label <- tolower(f3_SUB_BFs$label)

# reworking
f3_SUB_BFs <- f3_SUB_BFs[,-3]
colnames(f3_SUB_BFs)[1] <- "hemi"

SUBplt_data <- left_join(data.frame(hemi = aseg_edit$hemi, label = aseg_edit$label), f3_SUB_BFs)

colnames(SUBplt_data)[1] <- "NOThemi"

SUBplt_data %>%
  ggplot() +
  geom_brain(atlas = aseg, 
             #position = position_brain(hemi ~ side),
             aes(fill = bf))

ggplot() +
  geom_brain(atlas = aseg)


SUBplt_data %>%
  ggplot() +
  geom_brain(atlas = aseg, 
             # position = position_brain(hemi ~ side),
             aes(color = bf))



