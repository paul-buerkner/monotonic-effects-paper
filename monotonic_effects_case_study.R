library(tidyverse)
library(brms)

# sequential difference coding of factors
sdif_coding <- function(x) {
  x <- as.factor(x)
  contrasts(x) <- MASS::contr.sdif(levels(x))
  x
}

# load the data
data("ICFCoreSetCWP", package = "ordPens")
cwp <- ICFCoreSetCWP %>%
  select(-starts_with("e")) %>%
  mutate(
    # sequential coding for use in categorical models
    d450c = sdif_coding(d450),
    d455c = sdif_coding(d455)
  )

# we will include the following predictors:
# d450: impairments in walking
# d455: impairments in moving around

# investigate correlation between the two predictors
cor(cwp$d450, cwp$d455, method = "kendall")

# we start by analysing the effect of d455 only

# specify priors
prior_b <- prior(normal(0, 2.5), class = "b")
prior_s1 <- prior(dirichlet(1, 1, 1, 1), class = "simo", coef = "mod4551")

# fit and summarize the monotonic model
fit_mo <- brm(phcs ~ mo(d455), data = cwp, prior = prior_b + prior_s1)
summary(fit_mo)
marginal_effects(fit_mo)

# fit and summarize the linear model
fit_lin <- brm(phcs ~ d455, data = cwp, prior = prior_b)
summary(fit_lin)
marginal_effects(fit_lin)

# fit and summarize the categorical model
fit_cat <- brm(phcs ~ d455c, data = cwp, prior = prior_b)
summary(fit_cat)
marginal_effects(fit_cat)

# compare models using approximate LOO-CV
loo_compare(loo(fit_mo), loo(fit_lin), loo(fit_cat))
model_weights(fit_mo, fit_lin, fit_cat, weights = "loo")


# we will now include both d455 and d450 into the model

# add a prior for the simplex parameter of d450
prior_s2 <- 
  prior(dirichlet(1, 1, 1), class = "simo", coef = "mod4501") +
  prior(dirichlet(1, 1, 1, 1), class = "simo", coef = "mod4551")

# fit and summarize models including two predictors
fit_mo2 <- brm(phcs ~ mo(d455) + mo(d450), data = cwp, 
               prior = prior_b + prior_s2)
summary(fit_mo2)
plot(marginal_effects(fit_mo2), ask = FALSE)

fit_lin2 <- brm(phcs ~ d455 + d450, data = cwp, prior = prior_b)
summary(fit_lin2)
plot(marginal_effects(fit_lin2), ask = FALSE)

fit_cat2 <- brm(phcs ~ d455c + d450c, data = cwp, prior = prior_b)
summary(fit_cat2)
plot(marginal_effects(fit_cat2), ask = FALSE)

# compare models using approximate LOO-CV
loo_compare(loo(fit_mo2), loo(fit_lin2), loo(fit_cat2))
model_weights(fit_mo2, fit_lin2, fit_cat2, weights = "loo")

