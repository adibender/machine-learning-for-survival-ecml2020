library(pammtools)
library(dplyr)
library(tidyr)
library(purrr)
library(pammtools)
library(pem.xgb)
library(survival)
library(pec)
library(prodlim)

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
    problem_wrapper(data, .x, "metabric")
  })

saveRDS("instances/metabric.Rds")

res_xgb <- map(
  .x = seq_along(instances),
  ~ {
    inst <- instances[[.x]]
    res <- pam_xgb_wrapper(data, .x, inst, nthread = 30)
    saveRDS(res, paste0("instances/results/xgb/xgb_metabric_", .x, ".Rds"))
    pam_xgb <- attr(res, "pam_xgb")
    eval_times <- pem.xgb:::times_list(inst$data$test)$eval_times
    predictSurvProb(pam_xgb, inst$data$test, eval_times)
  }
)
saveRDS(res_xgb, "instances/results/pred_metabric_xgb.Rds")

## coxph
res_csc <- purrr::map(
  seq_len(repls),
  ~{
    print(.x)
    inst <- instances[[.x]]
    # Note: need to exclude some variables, otherwise NA are produced in some
    # coefficients and prediction doesn't work
    res_i <- coxph(Surv(time, status)~. - grade_2 - histological_11 - ER_IHC_status_1 -
      ER_Expr_1 - PR_Expz_1 - HER2_IHC_status_3 - HER2_SNP6_state_2 - Her2_Expr_1 -
      Treatment_7 - inf_men_status_1 - group_3 - group_4 - cellularity_2 - histological_8 -
      Pam50_Subtype_5 - int_clust_memb_9 - site_4 - Genefu_3,
    data= inst$data$train,
    x = TRUE)
    eval_times <- pem.xgb:::times_list(inst$data$test)$eval_times
    predictSurvProb(res_i, inst$data$test, eval_times)

})

saveRDS(res_csc, "instances/results/pred_metabric_coxph.Rds")
