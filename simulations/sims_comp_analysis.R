rmse <- function(x, y) {
  sqrt(mean((x - y)^2))
}

comp_analysis <- function(cond) {
  require(brms)
  out <- cond$res[[1]]
  if (!length(out)) {
    return(out)
  }
  varnames <- names(out[[1]]$data)
  estimates <- varnames[grepl("^Estimate_", varnames)]
  models <- sub("^Estimate_", "", estimates)
  ntrials <- cond$ntrials
  for (i in seq_len(cond$ntrials)) {
    out[[i]]$pred <- data.frame(
      cond = cond$cond, 
      trial = i,
      model = models,
      rmse = NA
    )
    for (j in seq_len(nrow(out[[i]]$pred))) {
      out[[i]]$pred$rmse[j] <- rmse(
        out[[i]]$data[[estimates[j]]], 
        out[[i]]$data$eta
      )
    }
  }
  out
}

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)

path <- "simulations/"
if (!grepl("2018_monotonic_effects$", getwd())) {
  path <- paste0("2018_monotonic_effects/", path)
}
comp_results <- readRDS(paste0(path, "comp_results.rds"))
I <- seq_len(nrow(comp_results))

for (i in I) {
  comp_results$res[[i]] <- comp_analysis(comp_results[i, ])
}

comp_preds <- comp_results$res %>%
  map(
    ~ map(., ~ .$pred) %>%
      do.call(bind_rows, args = .) %>%
      identity()
  ) %>%
  bind_rows() %>%
  inner_join(comp_results %>% select(-res), by = "cond")

# TODO: triton
# saveRDS(comp_preds, file = paste0(path, "comp_preds.rds"))


comp_preds %>%
  filter(cond == 21) %>%
  filter(!model %in% c("lin")) %>%
  ggplot(aes(model, rmse)) +
  geom_boxplot() +
  facet_wrap("cond", scales = "free_y")
  
