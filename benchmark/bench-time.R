library(dplyr)
library(tidyr)
library(purrr)
library(pammtools)
library(survival)
library(pec)
library(prodlim)
library(pem.xgb)

## functions to generate problem instances and wrappers for different algorithms
source("../simulation/sim-funs.R")
load("../benchmark/pec_params.RData")
pec_params <- list(
  formula = Surv(time, status) ~ 1,
  exact = FALSE,
  start =0.01)


# wrapper that performs estimation based on different choices w.r.t. cut-point
# selection (using all event times vs. subsample)
pam_xgb_wrapper_time <- function(
  instance,
  budget = 20,
  param_space = list(
    max_depth        = c(1L, 20L),
    gamma            = c(0, 5),
    lambda           = c(1, 3),
    colsample_bytree = c(.5, 1),
    min_child_weight = c(5L, 50L),
    subsample        = c(.5, 1)),
  ped_params = list(formula = Surv(time, status) ~ .),
  pec_params = list(formula = Surv(time, status) ~ 1, exact = FALSE, start =0.01),
  nrounds    = 3200L,
  eta        = .05,
  early_stopping_rounds = 50L,
  return_model = FALSE,
  nthread = 1L) {

  param_df <- purrr::map_dfr(
    seq_len(budget),
    ~ pem.xgb::random_params(param_space))
  param_df$eta <- eta

  st1 <- system.time({

  tune_res <- xgb.tune_pam(
    param_df              = param_df,
    data                  = instance$data$train,
    nrounds               = nrounds,
    cv_indices            = instance$cv_indices,
    split_frac            = split_frac,
    ped_params            = ped_params,
    verbose               = TRUE,
    print_every_n         = 400L,
    nthread               = nthread,
    early_stopping_rounds = early_stopping_rounds
  )

  best_set <- tune_res %>%
    filter(IBS == min(IBS))

  # split_idx <- get_split_indices(
  #   nrow(instance$data$train),
  #   split_frac = split_frac)
  ped_params$data <- instance$data$train
  # ped_params$data <- instance$data$train[split_idx$ind_train, , drop = FALSE]
  pam_xgb <- xgb.train.ped(
    params  = as.list(best_set$param_set[[1]]),
    data    = do.call(as_ped, ped_params),
    nrounds = best_set$nrounds,
    vebose  = TRUE,
    nthread = nthread)

  })

  ibs <- get_ibs(
    object         = list(pam_xgb = pam_xgb),
    data           = instance$data$test,
    pec_params     = pec_params,
    keep_reference = TRUE
  ) %>%
  mutate(
    time = st1["elapsed"],
    strategy = "full")

  ## with cut points based on subsampling
  st2 <- system.time({

  samp_id <- sample(seq_len(nrow(instance$data$train)), 200, replace = FALSE)
  cut <- pammtools:::get_cut.default(instance$data$train[samp_id, ], Surv(time, status)~.)
  ped_params$cut <- cut

  tune_res <- xgb.tune_pam(
    param_df              = param_df,
    data                  = instance$data$train,
    nrounds               = nrounds,
    cv_indices            = instance$cv_indices,
    split_frac            = split_frac,
    ped_params            = ped_params,
    verbose               = TRUE,
    print_every_n         = 400L,
    nthread               = nthread,
    early_stopping_rounds = early_stopping_rounds
  )

  best_set <- tune_res %>%
    filter(IBS == min(IBS))

  ped_params$data <- instance$data$train
  pam_xgb <- xgb.train.ped(
    params  = as.list(best_set$param_set[[1]]),
    data    = do.call(as_ped, ped_params),
    nrounds = best_set$nrounds,
    vebose  = TRUE,
    nthread = nthread)

  })


  ibs2 <- get_ibs(
    object         = list(pam_xgb = pam_xgb),
    data           = instance$data$test,
    pec_params     = pec_params,
    keep_reference = TRUE
  ) %>%
  mutate(
    time = st2["elapsed"],
    strategy = "subsample")

  rbind(ibs, ibs2)

}


## seed
init_seed <- 261019
repls <- 10
instances400 <- map(
  .x = seq_len(repls),
  ~{
    set.seed(init_seed + .x)
    sim_wrapper(n = 400)
  })
instances800 <- map(
  .x = seq_len(repls),
  ~{
    set.seed(init_seed + .x)
    sim_wrapper(n = 800)
  })

instances1600 <- map(
  .x = seq_len(repls),
  ~{
    set.seed(init_seed + .x)
    sim_wrapper(n = 1600)
  })

instances3200 <- map(
  .x = seq_len(repls),
  ~{
    set.seed(init_seed + .x)
    sim_wrapper(n = 3200)
  })

instances <- list(
  n400 = instances400,
  n800 = instances800,
  n1600 = instances1600,
  n3200 = instances3200)

saveRDS(instances, "instances/instances-sim-time.Rds")

library(furrr)
plan("multicore")
options(mc.cores = 10)
res_400 <- future_map(instances400,
    ~pam_xgb_wrapper_time(.x, budget = 10, nthread = 1L,
      pec_params = pec_params))

saveRDS(res_400, "bench-time-400.Rds")

res_800 <- future_map(instances800,
    ~pam_xgb_wrapper_time(.x, budget = 10, nthread = 1L,
      pec_params = pec_params))

saveRDS(res_800, "bench-time-800.Rds")

res_1600 <- future_map(instances1600,
    ~pam_xgb_wrapper_time(.x, budget = 10, nthread = 1L,
      pec_params = pec_params))
saveRDS(res_1600, "bench-time-1600.Rds")

set.seed(311916)
res_3200 <- future_map(instances3200,
    ~pam_xgb_wrapper_time(.x, budget = 10, nthread = 1L,
      pec_params = pec_params))
saveRDS(res_3200, "bench-time-3200.Rds")

## eval
res_400 <- readRDS("bench-time-400.Rds") %>% bind_rows() %>%
  mutate(n = 400)
res_800 <- readRDS("bench-time-800.Rds") %>% bind_rows() %>%
  mutate(n = 800)
res_1600 <- readRDS("bench-time-1600.Rds") %>% bind_rows() %>%
  mutate(n = 1600)
res_3200 <- readRDS("bench-time-3200.Rds") %>% bind_rows() %>%
  mutate(n = 3200)


res <- bind_rows(res_400, res_800, res_1600, res_3200) %>%
  filter(method == "pam_xgb")

saveRDS(res, "../paper/results-scaling-benchmark.Rds")
