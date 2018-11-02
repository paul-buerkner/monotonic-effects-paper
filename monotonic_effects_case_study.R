library(tidyverse)
library(brms)

# load the data
data("ICFCoreSetCWP", package = "ordPens")
cwp <- ICFCoreSetCWP %>%
  select(-starts_with("e"))

# specify priors
prior1 <- prior(normal(0, 10), class = "b") +
  prior(dirichlet(1, 1, 1), class = "simo", coef = "mod4501") +
  prior(dirichlet(1, 1, 1, 1), class = "simo", coef = "mod4551")

# fit and summarize the first model
fit1 <- brm(phcs ~ mo(d450) + mo(d455), data = cwp, prior = prior1)
summary(fit1)
marginal_effects(fit1)

# fit and summarize the second model
prior2 <- prior(normal(0, 2.5), class = "b")
fit2 <- brm(phcs ~ d450 + d455, data = cwp, prior = prior2)
summary(fit2)
marginal_effects(fit2)

# compare models using stacking
model_weights(fit1, fit2, weights = "loo2")

# prepare data for the fully models
cwp_NA <- rbind(cwp, c(rep(4, 51), phcs = NA))

# prior for the full models
prior3 <- prior(horseshoe(par_ratio = 0.1), class = "b")

# fit the third model
# takes some time
pred <- setdiff(names(cwp_NA), "phcs")
formula3 <- paste0("phcs | mi() ~ ", paste0("mo(", pred, ")", collapse = "+"))
formula3 <- as.formula(formula3)
fit3 <- brm(formula3, data = cwp_NA, prior = prior3) 

# plot size parameters
ps_bsp <- posterior_samples(fit3, pars = "^bsp_")
names(ps_bsp) <- sub("^bsp_mo", "", names(ps_bsp))
medians_bsp <- apply(ps_bsp, 2, median)
ps_bsp <- ps_bsp[, order(medians_bsp)]
bayesplot::mcmc_intervals(ps_bsp, prob_outer = 0.95)  

# fit the full model of linear effects
fit4 <- brm(phcs | mi() ~ ., data = cwp_NA, prior = prior3)
summary(fit4)

# compare models
model_weights(fit3, fit4, weights = "loo2", newdata = cwp) 
