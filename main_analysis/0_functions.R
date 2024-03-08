# Nicholas Judd
# LCD Lab - Donders
# 2024-02-12
# njudd.com

# Script for a RDHonest fuzzy RD with cov treatment and missing covs

# add nessesary packages

if (!require("pacman")) install.packages("pacman")
pacman::p_load(mice)
pacman::p_load_gh("kolesarm/RDHonest")

# adapted from github
# https://github.com/kolesarm/RDHonest/issues/7
# NOTE! Issues in his lm subsetting bounds; only subseting one side of the bounds



# bring a single vector to the fence using the Tukey Method
vec_to_fence <- function(vec){
  stats <- boxplot.stats(vec)$stats
  vec[vec < stats[1]] <- stats[1]
  vec[vec > stats[5]] <- stats[5]
  return(vec)
}



# fuzzy RD with covs corrected
RDjacked <- function(y, running, fuzzy, df, covs = NULL, report_efficiancy = FALSE) {
  
  print("Only handles fuzzy RD atm; no missingness")
  
  df = as.data.frame(df) # make sure df is a data.frame
  
  # subsetting the DF, so the complete.cases args later works out
  if (!is.null(covs)){
    df = df[, c(y, running, fuzzy, covs)] #all the data we need
  } else {
    df = df[, c(y, running, fuzzy)]
  }
  
  equation = paste0(y, "~", fuzzy, "|", running)
  rdpre <- RDHonest(equation, data = df) # no covs, estimate M + bandwidth
  
  rd1 <-  RDHonest(equation, data = df,
                     T0 = rdpre$coefficients$estimate) # as seen in the RDHonest tutorial
  
  if (is.null(covs)){
    print("Covariates not supplied")
    
    rd1 <- rd1$coefficients
    rd1$term <- paste(y)
    rd1$CI <- paste0("(", as.character(round(rd1$conf.low, 2)), ", ", as.character(round(rd1$conf.high, 2)), ")")
    
    names(rd1)[2:3] <- c("uY", "se")
    
    # selecting params
    rd1 <- rd1[c("term", "eff.obs", "bandwidth", "uY", "CI","p.value", "first.stage")]
    rd1[c("eff.obs", "bandwidth", "uY", "p.value", "first.stage")] <- round(rd1[c("eff.obs", "bandwidth", "uY", "p.value", "first.stage")], 3)
    
    return(rd1)
    
    
  } else {
    
    print("Warning: covs need to be dummy coded see fastDummies::dummy_cols")
    
    if(sum(colnames(df) %in% "y_adj") >= 1){
      print("WARNING y_adj present in dataframe; overwritting")
    }
    
    if(dim(df)[1] == sum(complete.cases(df))){
      print("No missing data; proceeding with Adjustment")
      
      # only on the earlier bandwidth
      sSet <- eval(parse(text = paste0(rlang::expr(df), "$", running))) <= rd1$coefficients$bandwidth & eval(parse(text = paste0(rlang::expr(df), "$", running))) >= -rd1$coefficients$bandwidth
      df2 = df[sSet,]
      
      # Step 1 Fuzzy
      f1 <- paste(c(paste0(fuzzy, "~", running, "*I(", running, ">=0)"), covs), collapse=" + ")
      r1 <- lm(f1, data = df2) 
      df$fuzzy_adj <- drop(df[,fuzzy]-as.matrix(df[, covs]) %*% r1$coefficients[covs])
      
      # Step 2 
      f2 <- paste(c(paste0(y, "~", running, "*I(", running, ">=0)"), covs), collapse=" + ")
      r2 <- lm(f2, data = df2)
      df$y_adj <- drop(df[,y]-as.matrix(df[, covs]) %*% r2$coefficients[covs])
      
      # he estimates new M's, lets do the same! https://github.com/kolesarm/RDHonest/issues/7
      
      rd2 <- RDHonest(paste0("y_adj ~ fuzzy_adj |", running), data = df, 
                      T0 = rdpre$coefficients$estimate, M = c(rd1$coefficients$M.rf, rd1$coefficients$M.fs))
      
      if (report_efficiancy == TRUE) {
        ##### now I am just making a nice output & showing the difference from cov correction
        # rd2 <- rd2$coefficients[c('term', 'estimate', 'conf.low', 'conf.high', 'bandwidth', 'eff.obs', 'p.value')]
        # names(rd2)[2] <- "Y"
        # 
        # Y_un <- rd1$coefficients$estimate # estimate uncorrected
        # 
        # dY <- abs(rd2$Y - Y_un)
        # dY_pi <- Y_un/abs(rd2$Y)*100
        # 
        # dCI_un <- abs(rd1$coefficients$conf.low - rd1$coefficients$conf.high)
        # 
        # dCI <- abs(rd2$conf.low - rd2$conf.high) - dCI_un
        # dCI_pi <- dCI/abs(rd2$conf.low - rd2$conf.high)*100
        # 
        # rd2 <- cbind(rd2, Y_un, dY, dY_pi, dCI, dCI_pi)
        # # Y uncorrected, abs change of the Y's, perfect change of Y's, raw CI change, percent CI change
        # 
        # rd2$CI <- paste0("(", as.character(round(rd2$conf.low, 2)), ", ", as.character(round(rd2$conf.high, 2)), ")")
        # 
        # rd2 <- rd2[c("term", "eff.obs", "bandwidth", "Y", "Y_un", "dY", "dY_pi","p.value", "CI", "dCI", "dCI_pi")]
        # 
        # rd2[c("eff.obs", "bandwidth", "Y", "Y_un", "dY", "dY_pi", "dCI", "dCI_pi")] <- round(rd2[c("eff.obs", "bandwidth", "Y", "Y_un", "dY", "dY_pi", "dCI", "dCI_pi")], 3)
        # 
        # rd2$term <- paste(y)
        # 
        # return(rd2)
      } else {
        
        # above compares corrected vs uncorrected, here we just get the corrected results
        rd2 <- rd2$coefficients
        rd2$term <- paste(y)
        rd2$CI <- paste0("(", as.character(round(rd2$conf.low, 2)), ", ", as.character(round(rd2$conf.high, 2)), ")")
        
        names(rd2)[2:3] <- c("Y", "se")
        
        # selecting params
        rd2 <- rd2[c("term", "eff.obs", "bandwidth", "Y", "CI","p.value", "first.stage")]
        
        rd2[c("eff.obs", "bandwidth", "Y", "p.value", "first.stage")] <- round(rd2[c("eff.obs", "bandwidth", "Y", "p.value", "first.stage")], 3)
        
        return(rd2)
      }
      

    } else if(dim(df)[1] != sum(complete.cases(df))){
      
      print("MISSING data issue, you have no value to mulitple the beta with to subtract from Y")
      print("ERROR: not continuiting atm")
      
      
      # # Step 1 formula
      # f1 <- paste(c(paste0(y, "~", running, "*I(", running, ">=0)"), covs), collapse=" + ")
      # # subsetting to the nessesary bandwidth
      # sSet <- eval(parse(text = paste0(rlang::expr(df), "$", running))) <= rd1$coefficients$bandwidth & eval(parse(text = paste0(rlang::expr(df), "$", running))) >= -rd1$coefficients$bandwidth
      # df2 = df[sSet,]
      # 
      # # the df is already subsetted by the nessesary cols
      # ini <- mice(df2, maxit=0, print=F)
      # pred <- ini$pred
      # 
      # pred[ ,fuzzy] <- 0
      # pred[ ,running] <- 0
      # imp <- mice(df2, pred=pred, print=F, m = 10)
      # 
      # fit = with(imp, f1)
      # 
      # r1 <- lm(f1, data = df2) 
      # # r1 <- lm(f1, data = data, subset = sSet)
      # # OMFG what a nightmare; gonna just subset a new df... this doesn't work in very interesting ways.
      # 
      # ### now this is tricky as you need to assign it
      # 
      # df$y_adj <- drop(df[,y]-as.matrix(df[, covs]) %*% r1$coefficients[covs])
      # 
      # rd2 <- RDHonest(paste0("y_adj ~", fuzzy, "|", running), data = df, 
      #                 T0 = rdpre$coefficients$estimate, # as in the vinjette for fuzzy RD's
      #                 M = c(rd1$coefficients$M.rf, rd1$coefficients$M.fs))
      # 
      # 
      # ##### now I am just making a nice output & showing the difference from cov correction
      # rd2 <- rd2$coefficients[c('term', 'estimate', 'conf.low', 'conf.high', 'bandwidth', 'eff.obs', 'p.value')]
      # names(rd2)[2] <- "Y"
      # 
      # Y_un <- rd1$coefficients$estimate # estimate uncorrected
      # 
      # dY <- abs(rd2$Y - Y_un)
      # dY_pi <- Y_un/abs(rd2$Y)*100
      # 
      # dCI_un <- abs(rd1$coefficients$conf.low - rd1$coefficients$conf.high)
      # 
      # dCI <- abs(rd2$conf.low - rd2$conf.high) - dCI_un
      # dCI_pi <- dCI/abs(rd2$conf.low - rd2$conf.high)*100
      # 
      # rd2 <- cbind(rd2, Y_un, dY, dY_pi, dCI, dCI_pi)
      # # Y uncorrected, abs change of the Y's, perfect change of Y's, raw CI change, percent CI change
      # 
      # rd2$CI <- paste0("(", as.character(round(rd2$conf.low, 2)), ", ", as.character(round(rd2$conf.high, 2)), ")")
      # 
      # rd2 <- rd2[c("term", "eff.obs", "bandwidth", "Y", "Y_un", "dY", "dY_pi","p.value", "CI", "dCI", "dCI_pi")]
      # 
      # rd2[c("eff.obs", "bandwidth", "Y", "Y_un", "dY", "dY_pi", "dCI", "dCI_pi")] <- round(rd2[c("eff.obs", "bandwidth", "Y", "Y_un", "dY", "dY_pi", "dCI", "dCI_pi")], 3)
      # 
      # rd2$term <- paste(y, rd2$term)
      # 
      # return(rd2)
      # 
      
      
      
      
      }
    }
}

# function to fit a simple Bayes model
# function to fit a simple mod with covs 
# bayesFIT("SA", "ROSLA", b1_covs, b1)

# allows parallel
# IVs <- rep(c("SA", "CT", "CSF_norm", "TBV_norm", "WM_hyper", "wFA"), each = 2) #amazing 
# DVs <- rep(c("ROSLA", "EduAge"), 6)
# plan(multisession, workers = 6) # THIS WORKS
# nkj <- future_map2(IVs, DVs, \(x, y) bayesFIT(x, y, b1_covs, b1), .options = furrr_options(seed = T))

bayesFIT <- function(dv, iv, covs, df){
  # fit the model
  
  f <- eval(parse(text = paste0(c(paste(dv, "~", iv), covs), collapse=" + ")))
  #^^^ this is needed for bayesfactor_parameters in bayes_to_results()
  m <- stan_glm(f, data = df, iter = 80000, refresh=0) #40000
    # assign it
  # assign(paste0("m", iv,"_", dv), m) #, envir = globalenv()) #  envir = parent.frame()
  
  # output the models
  # print(paste0("The following model is finished:", "m", iv,"_", dv))

  # logBF <- bayesfactor_parameters(m, null = 0)
  # logBF <- bayesfactor_parameters(m, null = 0)$log_BF[2]
  return(m) #list(m, logBF)
}

# a function that outputs results
# having issues with calling bayesfactor_parameters on the list object...

bayes_to_results <- function(listofMODs){
  
  print("Warning this function assumes the first IV in your model is the IV of interest!!!")
  
  effobs <- listofMODs %>% map_int(~dim(.$model)[1])
  binsize <- listofMODs %>% map_int(~length(unique(.$data$running_var))/2)
  
  # as.character(.formula[2])
  # Y <- listofMODs %>% map_chr(~str_split(.$formula, " ")[[1]][1])
  Y <- listofMODs %>% map_chr(~as.character(.$formula[2]))
  
  # as.character(.$formula[3])
  X <- listofMODs %>% map_chr(~str_split(as.character(.$formula[3]), " ")[[1]][1])
  median <- listofMODs %>% map_dbl(~round(mean(get_parameters(.)[,2]), 2))
  # # mean <- listofMODs %>% map(~mean(get_parameters(.)[,2]))
  ci_low <- listofMODs %>% map_dbl(~round(hdi(.)$CI_low[2],2))
  ci_high <- listofMODs %>% map_dbl(~round(hdi(.)$CI_high[2],2))
  
  # ci95 <- stringr::str_c("[", ci_low, ", ", ci_high, "]")
  
  logBF <- listofMODs %>% map_dbl(~round(bayesfactor_parameters(., null = 0)$log_BF[2], 2))
  
  df <- data.frame(Y, X, effobs, binsize, median, logBF, ci_low, ci_high) #, median, ci95, logBF)
  return(df)
}







