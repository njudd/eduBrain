





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





