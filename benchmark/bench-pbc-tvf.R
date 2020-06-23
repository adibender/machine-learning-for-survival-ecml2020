library(pammtools)
library(dplyr)
library(tidyr)
library(purrr)
library(pammtools)
library(pem.xgb)
library(survival)
library(pec)
library(prodlim)
library(batchtools)

## Here we have to redo the pbc analysis because initially the TVF setting was
# not comparable with the time-constant covariates setting
# (different number of observations, different CV indices, hyperparams, etc.)

## functions to generate problem instances and wrappers for different algorithms
source("problems.R")
source("algorithms.R")

## data sets + parameters
benchmark_data <- readRDS("benchmark_data.Rds")
load("pec_params.RData")
data = list(
  benchmark_data = benchmark_data,
  pec_params = pec_params
)

# custom adoptation to make sure everything is equal at estimation time
# for both modeling options
pam_xgb_wrapper_pbc <- function(
  data,
  job,
  instance,
  budget = 20,
  param_space = list(
    max_depth        = c(1L, 20L),
    gamma            = c(0, 5),
    lambda           = c(1, 3),
    colsample_bytree = c(.5, 1),
    min_child_weight = c(5L, 50L),
    subsample        = c(.5, 1)),
  nrounds    = 5000L,
  eta        = .05,
  early_stopping_rounds = 50L,
  nthread = 1L) {

  ped_params = list(formula = Surv(time, status) ~ .)
  ped_params_tvf = list(
    formula = Surv(time, status)~. +
      concurrent(bili, albumin, ast, protime, lag = 0, tz_var = "day"),
    id = "id")

  ## for data sets with huge p, reset colsample_bytree
  if(class(instance$data$train) == "list") {
    p <- ncol(instance$data$train[[1]])
  } else {
    p <- ncol(instance$data$train) - 2
  }

  param_df <- purrr::map_dfr(
    seq_len(budget),
    ~ pem.xgb::random_params(param_space))
  param_df$eta <- eta

  tune_res_x_t1 <- xgb.tune_pam(
    param_df              = param_df,
    data                  = instance$data$train[[1]],
    nrounds               = nrounds,
    cv_indices            = instance$cv_indices,
    split_frac            = split_frac,
    ped_params            = ped_params,
    verbose               = TRUE,
    print_every_n         = 500L,
    nthread               = nthread,
    early_stopping_rounds = early_stopping_rounds
  )
  tune_res_x_tvf <- xgb.tune_pam(
    param_df              = param_df,
    data                  = instance$data$train,
    nrounds               = nrounds,
    cv_indices            = instance$cv_indices,
    split_frac            = split_frac,
    ped_params            = ped_params_tvf,
    verbose               = TRUE,
    print_every_n         = 500L,
    nthread               = nthread,
    early_stopping_rounds = early_stopping_rounds
  )

  best_set_x_t1 <- tune_res_x_t1 %>%
    filter(IBS == min(IBS))
  best_set_x_tvf <- tune_res_x_tvf %>%
    filter(IBS == min(IBS))

  ped_params$data <- instance$data$train[[1]]
  ped_params_tvf$data <- instance$data$train
  # ped_params$data <- instance$data$train[split_idx$ind_train, , drop = FALSE]
  pam_xgb_x_t1 <- xgb.train.ped(
    params  = as.list(best_set_x_t1$param_set[[1]]),
    data    = do.call(as_ped, ped_params),
    nrounds = best_set_x_t1$nrounds,
    vebose  = TRUE,
    nthread = nthread)
  pam_xgb_x_tvf <- xgb.train.ped(
    params  = as.list(best_set_x_tvf$param_set[[1]]),
    data    = do.call(as_ped, ped_params_tvf),
    nrounds = best_set_x_tvf$nrounds,
    vebose  = TRUE,
    nthread = nthread)
  ibs <- get_ibs(
    object         = list(pam_xgb = pam_xgb_x_t1),
    data           = instance$data$test[[1]],
    pec_params     = data$pec_params,
    keep_reference = TRUE
  )
  ibs_tvf <- get_ibs(
    object         = list(pam_xgb = pam_xgb_x_tvf),
    data           = instance$data$test,
    pec_params     = data$pec_params,
    keep_reference = TRUE
  ) %>%
  mutate(method = ifelse(method == "pam_xgb", "pam_xgb_tvf", method))

  attr(ibs, "best_set") <- best_set_x_t1
  attr(ibs, "pam_xgb") <- pam_xgb_x_t1

  attr(ibs_tvf, "best_set") <- best_set_x_tvf
  attr(ibs_tvf, "pam_xgb") <- pam_xgb_x_tvf

  list(ibs, ibs_tvf)

}

# pbc_tvf <- benchmark_data[["pbc_tvf"]]
set.seed(251052)
instances <- map(
  1:25,
  ~{
    problem_wrapper(data, .x, "pbc_tvf")
  })
# saveRDS(instances_pbc, "instances/pbc.Rds")

library(furrr)
plan("multicore")
options(mc.cores = 25L)

res_pbc <- future_map(
  .x = seq_along(instances),
  ~{
    inst <- instances[[1]]
    res <- pam_xgb_wrapper_pbc(data, .x, inst, nthread = 1L)
    saveRDS(res, paste0("instances/results/xgb/xgb_pbc_tvf_", .x, ".Rds"))
    res
  }
)
xgb_ibs <- map_dfr(
  res_pbc,
  ~{
    .x[[2]] <- .x[[2]] %>% filter(method != "Reference")
    do.call(rbind, .x)
  }
)
# saveRDS(res_xgb, "instances/results/pred_pbc_tvf_xgb.Rds")

## orsf
res_orsf <- future_map(
  .x = seq_along(instances),
  ~{
    inst <- instances[[1]]
    inst$data$train <- inst$data$train[[1]]
    inst$data$test <- inst$data$test[[1]]
    res <- orsf_wrapper(data, .x, inst)
    saveRDS(res, paste0("instances/results/xgb/orsf_pbc", .x, ".Rds"))
    res
  }
)

orsf_ibs <- do.call(rbind, res_orsf)

## coxph
res_cph <- future_map(
  .x = seq_along(instances),
  ~{
    inst <- instances[[1]]
    inst$data$train <- inst$data$train[[1]]
    inst$data$test <- inst$data$test[[1]]
    res <- cph_wrapper(data, .x, inst)
    saveRDS(res, paste0("instances/results/xgb/cph_pbc", .x, ".Rds"))
    res
  }
)

cph_ibs <- do.call(rbind, res_cph)


ibs <- rbind(xgb_ibs, orsf_ibs, cph_ibs)

saveRDS(ibs, "../paper/res-ibs-pbc.Rds")
