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

# invert the sign of half of the elements of x
invert_half <- function(x) {
  len <- length(x)
  signs <- c(rep(1, ceiling(len / 2)), rep(-1, floor(len / 2)))
  x * signs[sample(seq_len(len))]
}

# sequential difference coding of factors
sdif_coding <- function(x) {
  x <- as.factor(x)
  contrasts(x) <- MASS::contr.sdif(levels(x))
  x
}

# compile a dummy brms model to avoid recompliation in each trial
get_dummy_model <- function(cond, effect) {
  require(brms)
  effect <- match.arg(effect, c("lin", "mo", "cat"))
  if (cond$pred == "main") {
    x <- sample(0:cond$D, size = cond$nobs, TRUE)
    xf <- sdif_coding(x)
    df <- data.frame(x, xf, y = 1)
    if (effect == "lin") {
      bform <- bf(y ~ x)
    } else if (effect == "mo") {
      bform <- bf(y ~ mo(x))
    } else if (effect == "cat") {
      bform <- bf(y ~ xf) 
    }
  } else if (cond$pred == "interaction") {
    x <- sample(0:cond$D, cond$nobs, TRUE)
    z <- sample(0:cond$D, cond$nobs, TRUE)
    xf <- sdif_coding(x)
    zf <- sdif_coding(z)
    df <- data.frame(x, xf, z, zf, y = 1)
    if (effect == "lin") {
      bform <- bf(y ~ x * z)
    } else if (effect == "mo") {
      bform <- bf(y ~ mo(x) * mo(z))
    } else if (effect == "cat") {
      bform <- bf(y ~ xf * zf) 
    }
  }
  if (cond$likelihood == "gaussian") {
    bfamily <- brmsfamily("gaussian")
  }
  # very wide priors considering the scale of the data
  bprior <- prior(normal(0, 10), class = "b")
  # use the default priors of brms for now
  out <- brm(bform, data = df, family = bfamily, 
             prior = bprior, chains = 0)
  out
}

run_trial <- function(cond, dummy_models, seed = NULL) {
  require(brms)
  require(bigsplines)
  if (!is.null(seed)) {
    # message("Using seed ", seed)
    set.seed(seed)
  }
  out <- list(truth = list())
  alpha <- rep(1, cond$D)
  if (cond$pred == "main") {
    out$truth$b_Intercept <- rnorm(1, 0, 1) 
    out$truth$bsp_mox <- rnorm(1, 0, 1) / cond$D
    if (cond$effect == "lin") {
      out$truth$simo_mox1 <- alpha / sum(alpha)
    } else if (cond$effect == "mo") {
      out$truth$simo_mox1 <- as.vector(rdirichlet(1, alpha)) 
    } else if (cond$effect == "cat") {
      out$truth$simo_mox1 <- invert_half(as.vector(rdirichlet(1, alpha)))
    }
  } else if (cond$pred == "interaction") {
    out$truth$b_Intercept <- rnorm(1, 0, 1)
    out$truth$bsp_mox <- rnorm(1, 0, 1) / cond$D
    out$truth$bsp_moz <- rnorm(1, 0, 1) / cond$D
    out$truth$`bsp_mox:moz` <- rnorm(1, 0, 1) / cond$D^2
    if (cond$effect == "lin") {
      out$truth$simo_mox1 <- alpha / sum(alpha)
      out$truth$simo_moz1 <- alpha / sum(alpha)
      out$truth$`simo_mox:moz1` <- alpha / sum(alpha)
      out$truth$`simo_mox:moz2` <- alpha / sum(alpha)
    } else if (cond$effect == "mo") {
      out$truth$simo_mox1 <- as.vector(rdirichlet(1, alpha))
      out$truth$simo_moz1 <- as.vector(rdirichlet(1, alpha))
      out$truth$`simo_mox:moz1` <- as.vector(rdirichlet(1, alpha))
      out$truth$`simo_mox:moz2` <- as.vector(rdirichlet(1, alpha))
    } else if (cond$effect == "cat") {
      out$truth$simo_mox1 <- invert_half(as.vector(rdirichlet(1, alpha)))
      out$truth$simo_moz1 <- invert_half(as.vector(rdirichlet(1, alpha)))
      out$truth$`simo_mox:moz1` <- invert_half(as.vector(rdirichlet(1, alpha)))
      out$truth$`simo_mox:moz2` <- invert_half(as.vector(rdirichlet(1, alpha)))
    }
  }
  if (cond$likelihood == "gaussian") {
    out$truth$sigma <- abs(rnorm(1, 0, 1))
  }
  if (cond$pred == "main") {
    x <- sample(0:cond$D, cond$nobs, TRUE)
    xf <- sdif_coding(x)
    eta <- out$truth$b_Intercept + 
      out$truth$bsp_mox * mo_trans(x, out$truth$simo_mox1)
    out$data <- data.frame(x, xf, eta)
  } else if (cond$pred == "interaction") {
    x <- sample(0:cond$D, size = cond$nobs, TRUE)
    z <- sample(0:cond$D, size = cond$nobs, TRUE)
    xf <- sdif_coding(x)
    zf <- sdif_coding(z)
    eta <- out$truth$b_Intercept + 
      out$truth$bsp_mox * mo_trans(x, out$truth$simo_mox1) +
      out$truth$bsp_moz * mo_trans(z, out$truth$simo_moz1) +
      out$truth$`bsp_mox:moz` * 
      mo_trans(x, out$truth$`simo_mox:moz1`) *
      mo_trans(z, out$truth$`simo_mox:moz2`)
    out$data <- data.frame(x, xf, z, zf, eta)
  }
  if (cond$likelihood == "gaussian") {
    y <- rnorm(cond$nobs, eta, out$truth$sigma)
  }
  out$data$y <- y
  
  # otherwise computing time explodes or errors occur for some trials
  skip.iter <- TRUE
  
  if (cond$pred == "main") {
    fit_lin <- lm(y ~ x, data = out$data)
    eta_lin <- fitted(fit_lin)
    
    fit_cat <- lm(y ~ xf, data = out$data)
    eta_cat <- fitted(fit_cat)
    
    # set knots for spline models
    nxlev <- length(unique(out$data$x))
    knotid <- binsamp(out$data$x, nmbin = nxlev)
    
    fit_os <- bigssp(
      y ~ x, data = out$data, type = c(x = "ord"),
      nknots = knotid, skip.iter = skip.iter, rseed = NULL
    )
    eta_os <- as.vector(predict(fit_os))
    
    fit_ls <- bigssp(
      y ~ x, data = out$data, type = c(x = "lin"),
      nknots = knotid, skip.iter = skip.iter, rseed = NULL
    )
    eta_ls <- as.vector(predict(fit_ls))
    
    fit_cs <- bigssp(
      y ~ x, data = out$data, type = c(x = "cub"),
      nknots = knotid, skip.iter = skip.iter, rseed = NULL
    )
    eta_cs <- as.vector(predict(fit_cs))
    
    # isoreg always assumes the relationship to be monotonically increasing
    x <- out$data$x
    y <- out$data$y
    is_neg <- mean(y[x == min(x)]) > mean(y[x == max(x)])
    sign <- if (is_neg) -1 else 1 
    
    fit_iso <- isoreg(x = x, y = sign * y)
    eta_iso <- sign * fitted(fit_iso)
    
    fit_osmo <- ordspline(
      x = x, y = sign * y, monotone = TRUE, 
      knots = unique(x)
    )
    eta_osmo <- sign * as.vector(predict(fit_osmo))
    
    out$data <- cbind(out$data, 
      Estimate_lin = eta_lin,
      Estimate_cat = eta_cat,
      Estimate_os = eta_os,
      Estimate_ls = eta_ls,
      Estimate_cs = eta_cs,
      Estimate_iso = eta_iso,
      Estimate_osmo = eta_osmo
    )
  } else if (cond$pred == "interaction") {
    fit_lin <- lm(y ~ x * z, data = out$data)
    eta_lin <- fitted(fit_lin)
    
    fit_cat <- lm(y ~ xf * zf, data = out$data)
    eta_cat <- fitted(fit_cat)
    
    # set knots for ordinal spline models
    nxlev <- length(unique(out$data$x))
    nzlev <- length(unique(out$data$z))
    pred_mat <- cbind(out$data$x, out$data$z)
    knotid <- binsamp(pred_mat, nmbin = c(nxlev, nzlev))
    
    fit_os <- bigssp(
      y ~ x * z, data = out$data, type = c(x = "ord", z = "ord"),
      nknots = knotid, skip.iter = skip.iter, rseed = NULL
    )
    eta_os <- as.vector(predict(fit_os))
    
    fit_ls <- bigssp(
      y ~ x * z, data = out$data, type = c(x = "lin", z = "lin"),
      nknots = knotid, skip.iter = skip.iter, rseed = NULL
    )
    eta_ls <- as.vector(predict(fit_ls))
    
    fit_cs <- bigssp(
      y ~ x * z, data = out$data, type = c(x = "cub", z = "cub"),
      nknots = knotid, skip.iter = skip.iter, rseed = NULL
    )
    eta_cs <- as.vector(predict(fit_cs))
    
    out$data <- cbind(out$data,
      Estimate_lin = eta_lin,
      Estimate_cat = eta_cat,
      Estimate_os = eta_os,
      Estimate_ls = eta_ls,
      Estimate_cs = eta_cs,
      # no version available for interactions
      Estimate_iso = NA,
      Estimate_osmo = NA
    )
  }
  
  # fit Bayesian models via brms
  # fit_lin <- suppressMessages(update(
  #   dummy_models[["lin"]], newdata = out$data, recompile = FALSE, 
  #   warmup = cond$ndraws, iter = 2 * cond$ndraws, chains = 1,
  #   refresh = 0
  # ))
  # eta_lin <- as.data.frame(fitted(fit_lin))
  # names(eta_lin) <- paste0(names(eta_lin), "_lin")
  
  fit_mo <- suppressMessages(update(
    dummy_models[["mo"]], newdata = out$data, recompile = FALSE, 
    warmup = cond$ndraws, iter = 2 * cond$ndraws, chains = 1,
    refresh = 0
  ))
  eta_mo <- as.data.frame(fitted(fit_mo))
  names(eta_mo) <- paste0(names(eta_mo), "_mo")
  out$data <- cbind(out$data, eta_mo) 
  
  # if (cond$D == 50 && cond$pred == "interaction" && cond$nobs <= 200) {
  #   # too many parameters for too few observations
  #   # which the unpenalized categorical model cannot handle
  #   eta_cat <- NA
  # } else {
  #   fit_cat <- suppressMessages(update(
  #     dummy_models[["cat"]], newdata = out$data, recompile = FALSE, 
  #     warmup = cond$ndraws, iter = 2 * cond$ndraws, chains = 1,
  #     refresh = 0
  #   ))
  #   eta_cat <- as.data.frame(fitted(fit_cat))
  #   names(eta_cat) <- paste0(names(eta_cat), "_cat") 
  # }
  # out$data <- cbind(out$data, eta_lin, eta_mo, eta_cat) 
  
  out
}

file <- "simulations/comp_results.rds"
if (!grepl("2018_monotonic_effects$", getwd())) {
  file <- paste0("2018_monotonic_effects/", file)
}

if (file.exists(file)) {
  comp_results <- readRDS(file)
} else {
  comp_results <- expand.grid(
    ntrials = 1000,
    D = c(4, 10, 50),
    ndraws = 500,
    nobs = c(50, 200, 1000),
    likelihood = "gaussian",
    effect = c("lin", "mo", "cat"),
    pred = c("main", "interaction"),
    stringsAsFactors = FALSE
  )
  comp_results$cond <- seq_len(nrow(comp_results))
  comp_results$res <- list(list())
}

I <- seq_len(nrow(comp_results))
# only run conditions for which no results exist yet
I <- intersect(I, which(!lengths(comp_results$res)))

library(foreach)
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
for (i in I) {
  message("Simulating condition ", i)
  cond <- comp_results[i, ]
  dummy_models <- list(lin = NULL, mo = NULL, cat = NULL)
  # for (eff in names(dummy_models)) {
  #   dummy_models[[eff]] <- get_dummy_model(cond, eff) 
  # }
  dummy_models[["mo"]] <- get_dummy_model(cond, "mo") 
  J <- seq_len(cond$ntrials)
  comp_results$res[[i]] <- 
    foreach(j = J, .packages = "brms") %dopar% 
    run_trial(cond, dummy_models, seed = j)
  print(warnings())
  saveRDS(comp_results, file = file)
}
stopCluster(cl)

saveRDS(comp_results, file = file)

