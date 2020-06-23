library(pammtools)
library(dplyr)
library(tidyr)
library(purrr)
library(pammtools)
library(pem.xgb)
library(survival)
library(pec)
library(prodlim)
devtools::load_all("../")

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


## seed
init_seed <- 1012019
repls <- 25
instances <- map(
  .x = seq_len(repls),
  ~{
    set.seed(init_seed + .x)
    problem_wrapper(data, .x, "mgus2")
  })

saveRDS(instances, "instances/mgus2.Rds")
library(furrr)
plan("multicore")

res_xgb <- future_map(
  .x = seq_along(instances),
  ~ {
    print(.x)
    inst <- instances[[.x]]
    res <- pam_xgb_wrapper(data, .x, inst, nthread = 1, budget = 20)
    saveRDS(res, paste0("instances/results/xgb/xgb_mgus2_", .x, ".Rds"))
    pam_xgb <- attr(res, "pam_xgb")
    eval_times_1 <- pem.xgb:::times_list(inst$data$test, status_value = 1)$eval_times
    eval_times_2 <- pem.xgb:::times_list(inst$data$test, status_value = 2)$eval_times
    cause1 <- predictEventProb(pam_xgb, inst$data$test, eval_times_1, cause = 1)
    cause2 <- predictEventProb(pam_xgb, inst$data$test, eval_times_2, cause = 2)
    list(cause1 = cause1, cause2 = cause2)
  }
)
saveRDS(res_xgb, "instances/results/pred_mgus2.Rds")

# res_temp <- map(1:4, ~{
#   inst <- instances[[.x]]
#   name <- paste0("instances/results/xgb/xgb_mgus2_", .x, ".Rds")
#   res <- readRDS(name)
#   pam_xgb <- attr(res, "pam_xgb")
#   eval_times_1 <- pem.xgb:::times_list(inst$data$test, status_value = 1)$eval_times
#   eval_times_2 <- pem.xgb:::times_list(inst$data$test, status_value = 2)$eval_times
#   cause1 <- predictEventProb(pam_xgb, inst$data$test, eval_times_1, cause = 1)
#   cause2 <- predictEventProb(pam_xgb, inst$data$test, eval_times_2, cause = 2)
#   list(cause1 = cause1, cause2 = cause2)
# })
# saveRDS(res_temp, "instances/results/pred_mgus2_tmp.Rds")

## coxph
res_csc <- purrr::map(
  seq_len(repls),
  ~{
    print(.x)
    inst <- instances[[.x]]
    res_i <- riskRegression::CSC(Hist(time, status)~ age + sex + dxyr +
      hgb + creat + mspike, data= inst$data$train)
    eval_times_1 <- pem.xgb:::times_list(inst$data$test, status_value = 1)$eval_times
    eval_times_2 <- pem.xgb:::times_list(inst$data$test, status_value = 2)$eval_times
    cause1 <- predictEventProb(res_i, inst$data$test, eval_times_1, cause = 1)
    cause2 <- predictEventProb(res_i, inst$data$test, eval_times_2, cause = 2)
    list(cause1 = cause1, cause2 = cause2)

})

saveRDS(res_csc, "instances/results/pred_mgus2_coxph.Rds")
