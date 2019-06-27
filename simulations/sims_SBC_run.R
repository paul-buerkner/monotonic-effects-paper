mo_trans <- function(x, zeta) {
  .mo_trans <- function(x, zeta) {
    if (x == 0) {
      out <- 0
    } else {
      out <- length(zeta) * sum(zeta[1:x])
    }
    return(out)
  }
  sapply(x, .mo_trans, zeta = zeta)
}

get_dummy_model <- function(cond) {
  require(brms)
  if (cond$pred == "main") {
    x <- sample(0:cond$D, size = cond$nobs, TRUE)
    df <- data.frame(x, y = 1)
    bform <- bf(y ~ mo(x))
    bprior <- prior(normal(0, 1), class = "Intercept") +
      prior(normal(0, 1), class = "b")
    # simo parameters have flat default priors
  } else if (cond$pred == "interaction") {
    x <- sample(0:cond$D, cond$nobs, TRUE)
    z <- sample(0:cond$D, cond$nobs, TRUE)
    df <- data.frame(x, z, y = 1)
    bform <- bf(y ~ mo(x)*mo(z))
    bprior <- prior(normal(0, 1), class = "Intercept") +
      prior(normal(0, 1), class = "b")
    # simo parameters have flat default priors
  }
  if (cond$likelihood == "gaussian") {
    bfamily <- brmsfamily("gaussian")
    bprior <- bprior + 
      # half-normal prior
      prior(normal(0, 1), class = "sigma")
  } else if (cond$likelihood == "poisson") {
    bfamily <- brmsfamily("poisson")
  }
  out <- brm(
    bform, data = df, family = bfamily, 
    prior = bprior, chains = 0
  )
  out
}

run_trial <- function(cond, dummy_model, seed = NULL) {
  require(brms)
  if (!is.null(seed)) {
    # message("Using seed ", seed)
    set.seed(seed)
  }
  out <- list(truth = list())
  alpha <- rep(1, cond$D)
  if (cond$pred == "main") {
    out$truth$b_Intercept <- rnorm(1, 0, 1)
    out$truth$bsp_mox <- rnorm(1, 0, 1) / cond$D
    out$truth$simo_mox1 <- as.vector(rdirichlet(1, alpha))
  } else if (cond$pred == "interaction") {
    out$truth$b_Intercept <- rnorm(1, 0, 1)
    out$truth$bsp_mox <- rnorm(1, 0, 1) / cond$D
    out$truth$bsp_moz <- rnorm(1, 0, 1) / cond$D
    out$truth$`bsp_mox:moz` <- rnorm(1, 0, 1) / cond$D^2
    out$truth$simo_mox1 <- as.vector(rdirichlet(1, alpha))
    out$truth$simo_moz1 <- as.vector(rdirichlet(1, alpha))
    out$truth$`simo_mox:moz1` <- as.vector(rdirichlet(1, alpha))
    out$truth$`simo_mox:moz2` <- as.vector(rdirichlet(1, alpha))
  }
  if (cond$likelihood == "gaussian") {
    # half normal prior
    out$truth$sigma <- abs(rnorm(1, 0, 1))
  }
  if (cond$pred == "main") {
    x <- sample(0:cond$D, cond$nobs, replace = TRUE)
    eta <- out$truth$b_Intercept + 
      out$truth$bsp_mox * mo_trans(x, out$truth$simo_mox1)
    out$data <- data.frame(x, eta)
  } else if (cond$pred == "interaction") {
    x <- sample(0:cond$D, size = cond$nobs, replace = TRUE)
    z <- sample(0:cond$D, size = cond$nobs, replace = TRUE)
    eta <- out$truth$b_Intercept + 
      out$truth$bsp_mox * mo_trans(x, out$truth$simo_mox1) +
      out$truth$bsp_moz * mo_trans(z, out$truth$simo_moz1) +
      out$truth$`bsp_mox:moz` * 
      mo_trans(x, out$truth$`simo_mox:moz1`) *
      mo_trans(z, out$truth$`simo_mox:moz2`)
    out$data <- data.frame(x, z, eta)
  }
  if (cond$likelihood == "gaussian") {
    y <- rnorm(cond$nobs, eta, out$truth$sigma)
  } else if (cond$likelihood == "poisson") {
    # creates unrealisticly large and variable data 
    # because of the vague and uncorrelated priors
    # does not make much sense to analyse this via SBC
    y <- rpois(cond$nobs, exp(eta))
  }
  out$data$y <- y
  fit <- suppressMessages(update(
    dummy_model, newdata = out$data, recompile = FALSE, 
    warmup = cond$ndraws, iter = 2 * cond$ndraws, chains = 1,
    refresh = 0
  ))
  out$ps <- posterior_samples(fit)
  out
}


file <- "simulations/SBC_results.rds"
if (!grepl("2018_monotonic_effects$", getwd())) {
  file <- paste0("2018_monotonic_effects/", file)
}

if (file.exists(file)) {
  SBC_results <- readRDS(file)
} else {
  SBC_results <- expand.grid(
    ntrials = 1000,
    D = c(4, 10, 50),
    ndraws = 500,
    nobs = c(50, 200, 1000),
    likelihood = "gaussian",
    pred = c("main", "interaction")
  )
  SBC_results$cond <- seq_len(nrow(SBC_results))
  SBC_results$res <- list(list())
}

I <- seq_len(nrow(SBC_results))
# only run conditions for which no results exist yet
I <- intersect(I, which(!lengths(SBC_results$res)))


library(foreach)
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
for (i in I) {
  message("Simulating condition ", i)
  cond <- SBC_results[i, ]
  dummy_model <- get_dummy_model(cond)
  J <- seq_len(cond$ntrials)
  SBC_results$res[[i]] <- 
    foreach(j = J, .packages = "brms") %dopar% 
    run_trial(cond, dummy_model, seed = j)
  print(warnings())
  saveRDS(SBC_results, file = file)
}
stopCluster(cl)

saveRDS(SBC_results, file = file)

