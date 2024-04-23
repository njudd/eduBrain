





############## Missingness Trouble Shooting

rd0 <- RDHonest(CT ~ EduAge16 | running_var, data = ct)
rd1 <- RDHonest(CT ~ EduAge16 | running_var, data = ct,
                T0 = rd0$coefficients$estimate) # no covs, estimate M + bandwidth

## Step 0 allow missingness in the covs

ini <- mice(ct, maxit=0, print=F)
pred <- ini$pred

pred[ ,"EduAge16"] <- 0
pred[ ,"running_var"] <- 0
imp <- mice(ct, pred=pred, print=F, m = 10)

## Step 1 Estimate the effects of the covs, use the same bandwidth as above
f1 <- paste(c("SA ~ running_var*I(running_var>=0)", covs), collapse=" + ")


data = rlang::expr(sa)
run ="running_var"

ct <- eval(parse(text = paste0(data, "$", run))) <= rd1$coefficients$bandwidth & eval(parse(text = paste0(data, "$", run))) >= -rd1$coefficients$bandwidth

r1 <- lm(f1, data = sa, subset = eval(parse(text = paste0(data, "$", run))) <= rd1$coefficients$bandwidth & eval(parse(text = paste0(data, "$", run))) >= -rd1$coefficients$bandwidth)
r1 <- lm(f1, data = sa, subset = ct)



fit = with(imp, lm(as.formula(f1), subset = (running_var <= rd1$coefficients$bandwidth)))
pool.fit <- pool(fit)
summary(pool.fit)

ct <- as.data.frame(ct)


ct$CT_adj2 <-drop(ct$CT-as.matrix(ct[, covs]) %*% pool.fit$pooled$estimate[pool.fit$pooled$term %in% covs])

## Step 2 Create covariate-adjusted outcome

ct$CT_adj <- drop(ct$CT-as.matrix(ct[, covs]) %*% r1$coefficients[covs])


cor(sa$SA_adj, sa$SA_adj2)


# Step 3 run RDHonest as normal with Y adjusted, keep same M as above
rd2 <- RDHonest(CT_adj~running_var, data = ct, M = rd1$coefficients$M)

#### some random stuff
ggplot(sa, aes(SA, fill = as.character(sex))) +
  geom_density(alpha = .5)


summary(lm(SA ~ EduAge16, data = sa))


temp = RDjackedFuzzy(y = 'SA', running = 'running_var', fuzzy = 'EduAge16', data = sa, covs = covs)






############## missingness playspace
############## missingness playspace






# poly(SAroi$visit_day_correct, 2)
# do everything with orthogonal polys?!??!







running = "running_var"
fuzzy = "EduAge16"
y = "SA"

df = as.data.frame(SAroi)

df = df[, c(y, running, fuzzy, covs)] 

rd1 <- RDHonest("SA ~ EduAge16 | running_var", df)

temp <- df$visit_day_correct


df$visit_day_correct <- poly(temp, 2)[,1]
df$visit_day_correct2 <- poly(temp, 2)[,2]




# Step 1 formula
f1 <- paste(c(paste0(y, "~", running, "*I(", running, ">=0)"), covs), collapse=" + ")
# subsetting to the nessesary bandwidth
sSet <- eval(parse(text = paste0(rlang::expr(df), "$", running))) <= rd1$coefficients$bandwidth & eval(parse(text = paste0(rlang::expr(df), "$", running))) >= -rd1$coefficients$bandwidth
df2 = df[sSet,]

# the df is already subsetted by the nessesary cols
ini <- mice(df2, maxit=0, print=F)
pred <- ini$pred

pred[ ,fuzzy] <- 0
pred[ ,running] <- 0
imp <- mice(df2, pred=pred, print=F, m = 10, method = "cart")
# https://stackoverflow.com/questions/48355250/do-imputation-in-r-when-mice-returns-error-that-system-is-computationally-singu
fit = with(imp, f1)


# https://github.com/amices/mice/issues/82
# predict than combine

r1 <- lm(f1, data = df2) 
# r1 <- lm(f1, data = data, subset = sSet)
# OMFG what a nightmare; gonna just subset a new df... this doesn't work in very interesting ways.

### now this is tricky as you need to assign it

df$y_adj <- drop(df[,y]-as.matrix(df[, covs]) %*% r1$coefficients[covs])

rd2 <- RDHonest(paste0("y_adj ~", fuzzy, "|", running), data = df, 
                T0 = rdpre$coefficients$estimate, # as in the vinjette for fuzzy RD's
                M = c(rd1$coefficients$M.rf, rd1$coefficients$M.fs))


##### now I am just making a nice output & showing the difference from cov correction
rd2 <- rd2$coefficients[c('term', 'estimate', 'conf.low', 'conf.high', 'bandwidth', 'eff.obs', 'p.value')]
names(rd2)[2] <- "Y"

Y_un <- rd1$coefficients$estimate # estimate uncorrected

dY <- abs(rd2$Y - Y_un)
dY_pi <- Y_un/abs(rd2$Y)*100

dCI_un <- abs(rd1$coefficients$conf.low - rd1$coefficients$conf.high)

dCI <- abs(rd2$conf.low - rd2$conf.high) - dCI_un
dCI_pi <- dCI/abs(rd2$conf.low - rd2$conf.high)*100

rd2 <- cbind(rd2, Y_un, dY, dY_pi, dCI, dCI_pi)
# Y uncorrected, abs change of the Y's, perfect change of Y's, raw CI change, percent CI change

rd2$CI <- paste0("(", as.character(round(rd2$conf.low, 2)), ", ", as.character(round(rd2$conf.high, 2)), ")")

rd2 <- rd2[c("term", "eff.obs", "bandwidth", "Y", "Y_un", "dY", "dY_pi","p.value", "CI", "dCI", "dCI_pi")]

rd2[c("eff.obs", "bandwidth", "Y", "Y_un", "dY", "dY_pi", "dCI", "dCI_pi")] <- round(rd2[c("eff.obs", "bandwidth", "Y", "Y_un", "dY", "dY_pi", "dCI", "dCI_pi")], 3)

rd2$term <- paste(y, rd2$term)

return(rd2)



############## missingness playspace
############## missingness playspace







