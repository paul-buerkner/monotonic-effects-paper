rank_stat <- function(samples, truth, L = 10, B = L + 1) {
  samples <- as.matrix(samples)
  out <- rep(NA, NCOL(samples))
  for (i in seq_along(out)) {
    out[i] <- sum(samples[, i] < truth[i]) 
  }
  if (B < L + 1) {
    stopifnot((L + 1) %% B == 0)
    div <- (L + 1) %/% B
    out <- out %/% div
  }
  out
}

get_pars <- function(ps) {
  pars <- colnames(ps)
  pars <- setdiff(pars, "lp__")
  pars <- sub("\\[.+", "", pars)
  unique(pars)
}

SBC_analysis <- function(cond, L = 10, B = L + 1) {
  require(brms)
  out <- cond$res[[1]]
  if (!length(out)) {
    return(out)
  }
  pars <- get_pars(out[[1]]$ps)
  nsamples <- nrow(out[[1]]$ps)
  ntrials <- cond$ntrials
  for (i in seq_len(cond$ntrials)) {
    out[[i]]$ranks <- list(cond = cond$cond)
    subset <- seq.int(1, nsamples, length.out = L)
    for (par in pars) {
      ps_par <- brms:::get_samples(out[[i]]$ps, paste0("^", par, "($|\\[)"))
      ps_par <- brms:::p(ps_par, subset)
      true_par <- brms:::p(out[[i]]$truth[[par]])
      out[[i]]$ranks[[par]] <- rank_stat(ps_par, true_par, L = L, B = B)
    }
  }
  out
}

library(dplyr)
library(tidyr)
library(purrr)

path <- "simulations/"
if (!grepl("2018_monotonic_effects$", getwd())) {
  path <- paste0("2018_monotonic_effects/", path)
}
SBC_results <- readRDS(paste0(path, "SBC_results.rds"))
I <- seq_len(nrow(SBC_results))
L <- 219
B <- 11

for (i in I) {
  SBC_results$res[[i]] <- SBC_analysis(SBC_results[i, ], L = L, B = B)
}

SBC_ranks <- SBC_results$res %>%
  map(
    ~ map(., ~ unlist(.$ranks)) %>%
      do.call(bind_rows, args = .)
  ) %>%
  map(~ if (nrow(.)) gather(., "par", "value", -cond)
  ) %>%
  bind_rows() %>%
  inner_join(SBC_results %>% select(-res), by = "cond")

saveRDS(SBC_ranks, file = paste0(path, "SBC_ranks.rds"))
