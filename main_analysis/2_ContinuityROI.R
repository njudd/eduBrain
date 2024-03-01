#### ### ----------------- ### ####
# Nicholas Judd
# Donders Institute - LCD Lab
# 2024-02-14 https://njudd.com 

### Table of Contents
# 2.1 Loading
# 2.2 SA regional fuzzy RD
# 2.3 CT regional fuzzy RD
# 2.4 FA roi fuzzy RD
# 2.5 Subcortical fuzzy RD
# 2.6 misc
#### ### ----------------- ### ####


### worknotes
# check *'s


# https://stackoverflow.com/questions/10323817/r-unexpected-results-from-p-adjust-fdr

### ***FLAIR is removed due to to little obs


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 2.1 Loading & cleaning ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

if (!require(pacman)){
  install.packages('pacman')
}

pacman::p_load(tidyverse, lubridate, stringr, RDHonest, fastDummies, mice, ggseg, kableExtra, rddensity)

# grab functions
source("~/projects/roslaUKB/main_analysis/0_functions.R") 
# function: vec_to_fence() takes a vector and puts outliers to the fence (i.e., boxplot.stats)
# function: RDjacked() for cov corrected fuzzy RD

# 1) covariates, 2) covariates to fence & 3)cols with running_var
covs <- c("sex", "summer", "visit_day_correct", "visit_day_correct2", "headmotion", #, "t2_FLAIR"*
          "imaging_center_11026", "imaging_center_11027", "imaging_center_11028", 
          "dMRI_25922_1", "dMRI_25921_1", "dMRI_25928_1")
covs_fence <- c("dMRI_25922_1", "dMRI_25921_1", "dMRI_25928_1", "headmotion")
cols = c("EduAge16", "running_var", covs)

# loading data
ukb <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/raw/current/31_brain_IDPs.csv")
fullset <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/proc/20240223_fullset.csv")
showcase <- data.table::fread("/Volumes/home/lifespan/nicjud/UKB/raw/Data_Dictionary_Showcase.csv")

# remove all the cols for the imaging followup (instance 3; -3)
followupCols <- colnames(ukb)[str_detect(colnames(ukb), "-3")]
ukb[, (followupCols):=NULL]


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 2.2 SA regional fuzzy RD ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# making sub steps so I don't mess up

# x.1 grab regional Surface Area
sa_codes <- showcase[Category==192,][str_detect(Field, "Area"),]

# x.2 grabbing the 3rd element in the string; "Area of X..."
sa_codes$ggSeg <- map_chr(str_split(sa_codes$Field, " "), 3)
# detecting left or right and adding it in a similar manor to ggseg # dk$data$label
sa_codes$ggSeg[str_detect(sa_codes$Field, "left")] <- str_c("lh_", sa_codes$ggSeg[str_detect(sa_codes$Field, "left")])
sa_codes$ggSeg[str_detect(sa_codes$Field, "right")] <- str_c("rh_", sa_codes$ggSeg[str_detect(sa_codes$Field, "right")])

# x.3 Check the string & seeing we don't need the left & right total SA
sa_codes <- sa_codes[!str_detect(Field, "Total"),]

# seems like now it matches with DK atlas
# sa_codes$ggSeg %in% dk$data$label

# x.4 selecting the cols from the UKB brain IDs that match our codes
SAroi <- ukb %>% 
  select(matches(as.character(sa_codes$FieldID))) # selecting based off the key; orders as well it seems

# making 100% sure the key is in the same order as the colnames
# UPDATE: automatically matches above in order...
# colnames(SAroi) <- str_remove_all(colnames(SAroi), "-2.0") # removing the instance code
# sa_codes <- sa_codes[match(colnames(SAroi), sa_codes$FieldID),]

# x.5 renmaing the cols to the key; this looks sketch up matches above orders it already
colnames(SAroi) <- sa_codes$ggSeg 

# x.6 giving IDs back
SAroi$eid <- ukb$eid 

# x.7 joining fullset for the running var & covs
SAroi <- fullset[SAroi, on = 'eid']

# x.8 removing subjects with missing Y's, running and IV
SAout <- c(sa_codes$ggSeg, "running_var", "EduAge16")
SAroi <- SAroi[complete.cases(SAroi[, ..SAout]),]

# x.9 fencing all relevant vars; regional Y's and covs to fence
SAroi <- SAroi[,c(sa_codes$ggSeg, covs_fence) := lapply(.SD, vec_to_fence), .SDcols=c(sa_codes$ggSeg, covs_fence)] 

# x.10 removing missing covs; for now...
SAroi <- SAroi[complete.cases(SAroi[, ..covs]),]

# str(SAroi[, c(covs), with = FALSE])

# **only two observations without T2 flair
# SAroi <- SAroi[t2_FLAIR ==1,]

# x.10 saving output
# data.table::fwrite(SAroi, "/Volumes/home/lifespan/nicjud/UKB/proc/20240222_SAroi.csv")

# x.11 looping thru each ROI and fitting RDjacked, while FDR correcting
# old loop, now using map
# saROI_results <- c()
# for (i in sa_codes$ggSeg) {
#   if (!is.data.frame(saROI_results)){
#     saROI_results <- RDjacked(i, "running_var", fuzzy = 'EduAge16', df = SAroi, covs = covs)
#   } else
#     saROI_results <- rbind(saROI_results,
#                            RDjacked(i, "running_var", fuzzy = 'EduAge16', df = SAroi, covs = covs))
# }
# saROI_results$pfdr <- p.adjust(saROI_results$p.value, method = "fdr")

# more efficient code

sa_results <- sa_codes$ggSeg %>%
  map_df(~RDjacked(., "running_var", fuzzy = 'EduAge16', df = SAroi, covs = covs)) %>% # looping thru the regions applying function and getting a DF
  mutate(pfdr = p.adjust(p.value, method = "fdr")) # adjusting the p-vals

sum(sa_results$pfdr < .05) # anything sig?

# x12 save results into html table

# sa_results %>% 
#   kbl(caption = "Surface Area ROI results") %>%
#   kable_styling("hover", full_width = F) %>% 
#   save_kable("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/RDcont/sa_results.html")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 2.3 CT regional fuzzy RD ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# x.1 grab regional Cortical Thickness
ct_codes <- showcase[Category==192,][str_detect(Field, "thickness"),]

# x.2 grabbing the 4rd element in the string; "Mean thickness of X..."
ct_codes$ggSeg <- map_chr(str_split(ct_codes$Field, " "), 4)
# detecting left or right and adding it in a similar manor to ggseg # dk$data$label
ct_codes$ggSeg[str_detect(ct_codes$Field, "left")] <- str_c("lh_", ct_codes$ggSeg[str_detect(ct_codes$Field, "left")])
ct_codes$ggSeg[str_detect(ct_codes$Field, "right")] <- str_c("rh_", ct_codes$ggSeg[str_detect(ct_codes$Field, "right")])

# x.3 we don't need the left & right "global mean thickness"
ct_codes <- ct_codes[!str_detect(Field, "Global"),]

# x.4 selecting the cols from the UKB brain IDs that match our codes
CTroi <- ukb %>% 
  select(matches(as.character(ct_codes$FieldID))) #selecting based off the key
# x.5 renmaing the cols to the key
colnames(CTroi) <- ct_codes$ggSeg
# x.6 giving IDs back
CTroi$eid <- ukb$eid 
# x.7 joining fullset for the running var & covs
CTroi <- fullset[CTroi, on = 'eid']
# x.8 removing subjects with missing Y's, running and IV
CTout <- c(ct_codes$ggSeg, "running_var", "EduAge16")
CTroi <- CTroi[complete.cases(CTroi[, ..CTout]),]
# x.9 fencing all relevant vars
CTroi <- CTroi[,c(ct_codes$ggSeg, covs_fence) := lapply(.SD, vec_to_fence), .SDcols=c(ct_codes$ggSeg, covs_fence)] 
# x.10 removing missing covs; for now...
CTroi <- CTroi[complete.cases(CTroi[, ..covs]),]

# x.10 saving output
# data.table::fwrite(CTroi, "/Volumes/home/lifespan/nicjud/UKB/proc/20240222_CTroi.csv")

# x.11 looping thru each ROI and fitting RDjacked, while FDR correcting
ct_results <- ct_codes$ggSeg %>%
  map_df(~RDjacked(., "running_var", fuzzy = 'EduAge16', df = CTroi, covs = covs)) %>% # looping thru the regions applying function and getting a DF
  mutate(pfdr = p.adjust(p.value, method = "fdr")) # adjusting the p-vals

sum(ct_results$pfdr < .05)



# x12 save results into html table
# ct_results %>% 
#   kbl(caption = "Cortical thickness ROI results") %>%
#   kable_styling("hover", full_width = F) %>% 
#   save_kable("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/RDcont/ct_results.html")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 2.4 weighted FA regional fuzzy RD ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# x.1 grabbing regional weighted FA
wFA_codes <- showcase[showcase$Category==135,][str_detect(Field, "Weighted-mean FA"),]

# x.2 grabbing the element in the string; more complicated cause it varies in length
wFA_codes$ggSeg <- map(str_split(wFA_codes$Field, " "), `[`, -c(1:4)) %>% #split the char on spaces & get rid of the first 4 elements
  map(str_flatten, "_") %>% #merge it back together, with a _ inbetween words
  flatten_chr() # flatten the list into a char vector
# there are some without left & right; so I won't code for this

# x.3 I think I need everything there, so I wont remove anything
#yet, the naming is a problem RDjacked thinks it is a function
wFA_codes$ggSeg <- str_remove(wFA_codes$ggSeg, "\\(")
wFA_codes$ggSeg <- str_remove(wFA_codes$ggSeg, "\\)")
wFA_codes$ggSeg <- str_remove(wFA_codes$ggSeg, "-")

# x.4 selecting the cols from the UKB brain IDs that match our codes
wFAroi <- ukb %>% 
  select(matches(as.character(wFA_codes$FieldID))) #selecting based off the key and matches order
# x.5 renmaing the cols to the key
colnames(wFAroi) <- wFA_codes$ggSeg
# x.6 giving IDs back
wFAroi$eid <- ukb$eid 
# x.7 joining fullset for the running var & covs
wFAroi <- fullset[wFAroi, on = 'eid']
# x.8 removing subjects with missing Y's, running and IV
wFAout <- c(wFA_codes$ggSeg, "running_var", "EduAge16")
wFAroi <- wFAroi[complete.cases(wFAroi[, ..wFAout]),]
# x.9 fencing all relevant vars; regional Y's and covs to fence
wFAroi <- wFAroi[,c(wFA_codes$ggSeg, covs_fence) := lapply(.SD, vec_to_fence), .SDcols=c(wFA_codes$ggSeg, covs_fence)] 
# x.10 removing missing covs; for now...
wFAroi <- wFAroi[complete.cases(wFAroi[, ..covs]),]

# x.10 saving output; also for global analysis
# y <- wFA_codes$ggSeg
# wFAroi$wFA <- rowMeans(wFAroi[, ..y])
# data.table::fwrite(wFAroi, "/Volumes/home/lifespan/nicjud/UKB/proc/20240222_wFAGLOBAL.csv")

# x.11 looping thru each ROI and fitting RDjacked, while FDR correcting
wFA_results <- wFA_codes$ggSeg %>%
  map_df(~RDjacked(., "running_var", fuzzy = 'EduAge16', df = wFAroi, covs = covs)) %>% # looping thru the regions applying function and getting a DF
  mutate(pfdr = p.adjust(p.value, method = "fdr")) # adjusting the p-vals

sum(wFA_results$pfdr < .05) # anything sig?

# x12 save results into html table
# wFA_results %>% 
#   kbl(caption = "weighted FA ROI results") %>%
#   kable_styling("hover", full_width = F) %>% 
#   save_kable("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/RDcont/wfa_results.html")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
#### 2.5 subcortical ROIs fuzzy RD ####
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# x.1 grabbing subcortical ROIs; we care about volumne ofc
sub_codes <- showcase[showcase$Category==190,][str_detect(Field, "Volume"),]

# x.2 grabbing the element in the string; more complicated cause there a bunch of stuff
# nothing whole brain, no ventricles, no vessels
sub_codes <- sub_codes[str_detect(sub_codes$Field, "whole brain|Ventricle|vent|Vent|vessel", negate = TRUE),] # nothing whole brain
# grabbing 3rd element "Volume of X..."
sub_codes$ggSeg <- map_chr(str_split(sub_codes$Field, " "), 3)
# adding left & right tags
sub_codes$ggSeg[str_detect(sub_codes$Field, "left")] <- str_c("lh_", sub_codes$ggSeg[str_detect(sub_codes$Field, "left")])
sub_codes$ggSeg[str_detect(sub_codes$Field, "right")] <- str_c("rh_", sub_codes$ggSeg[str_detect(sub_codes$Field, "right")])

# x.3 removed elements in x.2 but there is still unwanted stuff
sub_codes <- sub_codes[str_detect(sub_codes$Field, " Cortex|Cerebellum-White-Matter|CerebralWhiteMatter", negate = TRUE),]
sub_codes$ggSeg <- str_remove(sub_codes$ggSeg, "-") # doesn't play nice with RDHonest i.e., as.formula()

# x.4 selecting the cols from the UKB brain IDs that match our codes
sROI <- ukb %>%
  select(matches(as.character(sub_codes$FieldID))) #selecting based off the key and matches order
# x.5 renmaing the cols to the key
colnames(sROI) <- sub_codes$ggSeg
# x.6 giving IDs back
sROI$eid <- ukb$eid 
# x.7 joining fullset for the running var & covs
sROI <- fullset[sROI, on = 'eid']
# x.8 removing subjects with missing Y's, running and IV
Sout <- c(sub_codes$ggSeg, "running_var", "EduAge16")
sROI <- sROI[complete.cases(sROI[, ..Sout]),]
# x.9 fencing all relevant vars; regional Y's and covs to fence
sROI <- sROI[,c(sub_codes$ggSeg, covs_fence) := lapply(.SD, vec_to_fence), .SDcols=c(sub_codes$ggSeg, covs_fence)] 
# x.10 removing missing covs; for now...
sROI <- sROI[complete.cases(sROI[, ..covs]),]
# x.10 saving output
# data.table::fwrite(sROI, "/Volumes/home/lifespan/nicjud/UKB/proc/20240222_sROI.csv")

# x.11 looping thru each ROI and fitting RDjacked, while FDR correcting
sROI_results<- sub_codes$ggSeg %>%
  map_df(~RDjacked(., "running_var", fuzzy = 'EduAge16', df = sROI, covs = covs)) %>% # looping thru the regions applying function and getting a DF
  mutate(pfdr = p.adjust(p.value, method = "fdr")) # adjusting the p-vals

sum(sROI_results$pfdr < .05) # anything sig?

# x12 save results into html table
# sROI_results %>% 
#   kbl(caption = "Subcortical Vol results") %>%
#   kable_styling("hover", full_width = F) %>% 
#   save_kable("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/RDcont/sub_results.html")



# plot of number of observations

sa_results$modality <- rep("SA", dim(sa_results)[1])
ct_results$modality <- rep("CT", dim(ct_results)[1])
wFA_results$modality <- rep("wFA", dim(wFA_results)[1])
sROI_results$modality <- rep("Subcortical", dim(sROI_results)[1])

largedf <- rbind(sa_results, ct_results, wFA_results, sROI_results)

ggplot(largedf, aes(modality, eff.obs, fill = modality)) +
  ggrain::geom_rain() +
  theme_minimal()

ggplot(largedf, aes(modality, Y, fill = modality)) +
  ggrain::geom_rain() +
  theme_minimal()


SI_3.1 <- ggplot(largedf, aes(modality, eff.obs, fill = modality)) +
  ggrain::geom_rain() +
  labs(y = "Effective n", x = "") +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = 'Dark2')
# ggsave("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/plts/SI_3.1_regionN.png", bg = "white", width = 10, height = 6)
# 
# 
SI_3.2 <- ggplot(largedf, aes(modality, p.value, fill = modality)) +
  ggrain::geom_rain() +
  labs(y = "p value (uncorrected)", x = "") +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none") +
  scale_fill_brewer(palette = 'Dark2')

SI_3 <- SI_3.1 / SI_3.2 + plot_annotation(tag_levels = 'a')

# ggsave("~/My_Drive/life/10 Projects/10.02 ROSLA UK BioBank/results/plts/SI_3.png", SI_3, bg = "white", width = 14, height = 14)











