library(ggplot2)
theme_set(theme_bw())
library(batchtools)
source("../benchmark/algorithms.R")
load("../benchmark/pec_params.RData")
source("sim-funs.R")

if(!checkmate::test_directory_exists("output/tve-registry")) {

  reg <- makeExperimentRegistry(
    "output/tve-registry",
    packages = c("dplyr", "tidyr", "pammtools", "pem.xgb", "obliqueRSF",
      "survival", "pec", "prodlim"),
    seed     = 22112019)
  reg <- loadRegistry("output/tve-registry", writeable = TRUE)
  reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = 10)
  addProblem(
    name = "syntheticTVE",
    data = list(pec_params = pec_params),
    fun  = sim_wrapper,
    seed = 23112019)

  addAlgorithm(name = "pam_xgb", fun = pam_xgb_wrapper)
  addAlgorithm(name = "orsf", fun = orsf_wrapper)
  addAlgorithm(name = "coxph", fun = cph_wrapper)

  sim_df  <- data.frame(n = c(1000))
  algo_df <- data.frame(budget = 20)
  addExperiments(
    prob.design  = list(syntheticTVE = sim_df),
    algo.designs = list(pam_xgb = algo_df, orsf = algo_df, coxph = data.frame(x = TRUE)),
    repls        = 25L)
  ids_tve <- findExperiments(prob.name = "syntheticTVE")

  # system.time({
    # j1 <- makeJob(ids_tve$job.id[1])
    # inst1 <- j1$instance
    # prob <- testJob(ids_tve$job.id[1])
    # test1 <- testJob(id = ids_tve$job.id[51])
  # })
  # test1
  # plot(test1)
  # crps(test1, times = c(1, 4, 9))

  submitJobs(findNotDone(ids_tve))
  waitForJobs()

  ############################### competing risks
  addProblem(
    name = "syntheticTVE_CR",
    data = list(pec_params = pec_params),
    fun  = sim_wrapper_cr,
    seed = 23112019)

  addAlgorithm(name = "pam_xgb", fun = pam_xgb_wrapper)
  addAlgorithm(name = "csc", fun = csc_wrapper)

  sim_df <- data.frame(n = c(500))

  addExperiments(
    prob.designs = list(syntheticTVE_CR = sim_df),
    algo.designs = list(
      pam_xgb = data.frame(budget = 20),
      csc = data.frame(budget = 10)
    ),
    repls = 25L)

  ids_tve_cr <- findExperiments(prob.name = "syntheticTVE_CR")
  # system.time({
  #   test1 <- testJob(id = ids_tve_cr[1])
  #   test26 <- testJob(id = ids_tve_cr[11])
  # })
  # test1
  # plot(test1)
  # crps(test1, times = c(1, 4, 9))

  submitJobs(ids = findNotDone(ids_tve_cr))
  waitForJobs()

}

library(batchtools)
library(ggplot2)
theme_set(theme_bw())
reg     <- loadRegistry("output/tve-registry", writeable = TRUE)
reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = 10L)
# reg$cluster.functions = makeClusterFunctionsInteractive()
ids     <- findExperiments(prob.name="syntheticTVE")
# censoring percentage
mean_cens <- purrr::map_dbl(ids[[1]][1:25], ~{
  1 - mean(makeJob(.x)$instance$data$train$status)
})
mean(mean_cens)

## IBS evaluation
pars    <- unwrap(getJobPars()) %>% as_tibble()
res_tve <- reduceResultsDataTable(ids=findDone(ids)) %>%
  as_tibble() %>%
  tidyr::unnest() %>%
  left_join(pars) %>%
  mutate(method = factor(method, levels = c("Reference", "coxph", "orsf", "pam_xgb")))
saveRDS(res_tve, "res-synthetic-tve.Rds")

## competing risks setting
ids     <- findExperiments(prob.name="syntheticTVE_CR")
# censoring percentage
mean_cens <- purrr::map_dbl(ids[[1]][1:25], ~{
  mean(makeJob(.x)$instance$data$train$status == 0)
})
mean(mean_cens)

## IBS evaluation (but comparison with deephit was based on metrics defined there)
pars    <- unwrap(getJobPars()) %>% as_tibble()
res_tve <- reduceResultsDataTable(ids=findDone(ids)) %>%
  as_tibble() %>%
  tidyr::unnest() %>%
  left_join(pars) %>%
  mutate(method = factor(method, levels = c("Reference", "csc", "pam_xgb")))
# res_tve <- filter(res_tve, IBS < 1)
saveRDS(res_tve, "res-synthetic-tve-cr.Rds")
