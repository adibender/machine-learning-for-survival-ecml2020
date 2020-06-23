library(batchtools)

reg <- loadRegistry("output/benchmark-registry")
ids <- findExperiments(prob.name = "benchmark")
pars <- unwrap(getJobPars(ids)) %>% as_tibble()
ids_df_unique <- pars %>%
  filter(data_name != "icu", data_name != "synthetic") %>%
  filter(algorithm == "pam_xgb")
ids_data_list <- split(ids_df_unique, f = ids_df_unique$data_name)
names(ids_data_list) <- purrr::map(ids_data_list, ~.x$data_name[1])

purrr::walk(
  ids_data_list,
  ~{
    data_name <- .x$data_name[1]
    out <- purrr::map(.x$job.id, function(id) {
      instance_i <- makeJob(id)$instance
    })
    saveRDS(out, paste0("instances/", data_name, ".Rds"))
    invisible()
  })


reg <- loadRegistry("../simulation/output/tve-registry")
ids <- findExperiments(prob.name = "syntheticTVE")
pars <- unwrap(getJobPars(ids)) %>% as_tibble()
ids_df_unique <- pars %>%
  filter(algorithm == "pam_xgb")

out <- purrr::map(
  ids_df_unique$job.id,
  ~ makeJob(.x)$instance
)
saveRDS(out, "instances/syntheticTVE.Rds")


##################### Competing risks ##########################################
## time-varying effects for cause 1, linear effects for cause 2
reg <- loadRegistry("../simulation/output/tve-registry")
ids <- findExperiments(prob.name = "syntheticTVE_CR")
pars <- unwrap(getJobPars(ids)) %>% as_tibble()
ids_df_unique <- pars %>%
  filter(algorithm == "pam_xgb")

out <- purrr::map(
  ids_df_unique$job.id,
  ~ makeJob(.x)$instance
)
saveRDS(out, "instances/syntheticTVE_CR.Rds")

## save predictEventProb results for comparison with DeepHit
library(batchtools)
source("algorithms.R")
benchmark_data <- readRDS("benchmark_data.Rds")
pec_params <- list(formula = Hist(time, status) ~ 1, exact = FALSE, start = .01)
data <- list(
  benchmark_data = benchmark_data,
  pec_params     = pec_params)
reg <- loadRegistry("../simulation/output/tve-registry")
pars <- unwrap(getJobPars()) %>% as_tibble()
pars <- filter(pars, algorithm == "pam_xgb", problem == "syntheticTVE_CR")
ids <- pars$job.id
res <- purrr::map(
  seq_along(ids),
  ~{
    res_i <- loadResult(ids[.x])
    best_set <- attributes(res_i)$best_set$param_set
    nrounds <- attributes(res_i)$best_set$nrounds
    instance_i <- makeJob(ids[.x])$instance
    ped_params = list(formula = Surv(time, status) ~ .)
    ped_params$data <- instance_i$data$train
    pam_xgb <- xgb.train.ped(
      params  = as.list(best_set[[1]]),
      data    = do.call(as_ped, ped_params),
      nrounds = nrounds,
      vebose  = TRUE)
    eval_times_1 <- pem.xgb:::times_list(instance_i$data$test, status_value = 1)$eval_times
    eval_times_2 <- pem.xgb:::times_list(instance_i$data$test, status_value = 2)$eval_times
    cause1 <- predictEventProb(pam_xgb, instance_i$data$test, eval_times_1, cause = 1)
    cause2 <- predictEventProb(pam_xgb, instance_i$data$test, eval_times_2, cause = 2)
    list(cause1 = cause1, cause2 = cause2)
})
saveRDS(res, "instances/results/pred_syntheticTVE_CR_xgb.Rds")

synth_tve_cr_csc <- filter(pars, algorithm == "csc", problem == "syntheticTVE_CR")
ids_csc <- synth_tve_cr_csc$job.id
res_csc <- purrr::map(
  seq_along(ids_csc),
  ~{
    print(.x)
    res_i <- loadResult(ids_csc[.x])
    instance_i <- makeJob(ids[.x])$instance
    form <- paste0("Hist(time, status)~", paste0("x", 1:13, collapse = "+"))
    res_i <- riskRegression::CSC(as.formula(form), data= instance_i$data$train)
    eval_times_1 <- pem.xgb:::times_list(instance_i$data$test, status_value = 1)$eval_times
    eval_times_2 <- pem.xgb:::times_list(instance_i$data$test, status_value = 2)$eval_times
    cause1 <- predictEventProb(res_i, instance_i$data$test, eval_times_1, cause = 1)
    cause2 <- predictEventProb(res_i, instance_i$data$test, eval_times_2, cause = 2)
    list(cause1 = cause1, cause2 = cause2)
})
saveRDS(res_csc, "instances/results/pred_syntheticTVE_CR_csc.Rds")

############################### Metabric #######################################
## save predictEventProb results for comparison with DeepHit
library(batchtools)
source("algorithms.R")
benchmark_data <- readRDS("benchmark_data.Rds")
pec_params <- list(formula = Hist(time, status) ~ 1, exact = FALSE, start = .01)
data <- list(
  benchmark_data = benchmark_data,
  pec_params     = pec_params)
reg <- loadRegistry("output/benchmark-registry")
pars <- getJobPars()
pars$algo.pars <- NULL
### XGB (PEM)
pars <- unwrap(pars) %>% as_tibble()
pars <- filter(pars, algorithm == "pam_xgb", data_name == "metabric")
ids <- pars$job.id[1:25] # only first 25 estimated with deephit
instances <- purrr::map(ids, ~ makeJob(.x)$instance)
res <- purrr::map(
  seq_along(ids),
  ~{
    print(.x)
    res_i      <- loadResult(ids[.x])
    best_set   <- attributes(res_i)$best_set$param_set
    nrounds    <- attributes(res_i)$best_set$nrounds
    instance_i <- makeJob(ids[.x])$instance
    ped_params <- list(formula = Surv(time, status) ~ .)
    ped_params$data <- instance_i$data$train
    pam_xgb <- xgb.train.ped(
      params  = as.list(best_set[[1]]),
      data    = do.call(as_ped, ped_params),
      nrounds = nrounds,
      nthreads = 5,
      verbose  = TRUE)
    eval_times <- pem.xgb:::times_list(instance_i$data$test, status_value = 1)$eval_times
    predictSurvProb(pam_xgb, instance_i$data$test, eval_times_1)
})
saveRDS(res, "instances/results/pred_metabric_xgb.Rds")

## coxph
pars <- getJobPars()
pars$algo.pars <- NULL
pars <- unwrap(pars) %>% as_tibble()
pars <- filter(pars, algorithm == "cph", problem == "metabric")
ids <- pars$job.id
res <- purrr::map(
  seq_along(ids),
  ~{
    print(.x)
    res_i <- loadResult(ids[.x])
    instance_i <- makeJob(ids[.x])$instance
    form <- paste0("Hist(time, status)~", paste0("x", 1:13, collapse = "+"))
    res_i <- riskRegression::CSC(as.formula(form), data= instance_i$data$train)
    eval_times_1 <- pem.xgb:::times_list(instance_i$data$test, status_value = 1)$eval_times
    eval_times_2 <- pem.xgb:::times_list(instance_i$data$test, status_value = 2)$eval_times
    cause1 <- predictEventProb(res_i, instance_i$data$test, eval_times_1, cause = 1)
    cause2 <- predictEventProb(res_i, instance_i$data$test, eval_times_2, cause = 2)
    list(cause1 = cause1, cause2 = cause2)
})

saveRDS(res_csc, "instances/results/pred_metabric_coxph.Rds")
