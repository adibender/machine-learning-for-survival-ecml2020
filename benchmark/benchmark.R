library(batchtools)
library(ggplot2)
theme_set(theme_bw())
# devtools::load_all("../")

benchmark_data <- readRDS("benchmark_data.Rds")
pec_params <- list(formula = Surv(time, status) ~ 1, exact = FALSE, start = .01)
save(pec_params, file = "pec_params.RData")

source("problems.R")
source("algorithms.R")

if(!checkmate::test_directory_exists("output/benchmark-registry")) {

  reg <- makeExperimentRegistry(
    "output/benchmark-registry",
    packages = c("dplyr", "tidyr", "pammtools", "pem.xgb", "obliqueRSF",
      "survival", "pec", "prodlim"),
    load = "pec_params.RData",
    seed = 30112019)
  reg <- loadRegistry("output/benchmark-registry", writeable = TRUE)
  reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = 5)
  addProblem(
    name = "benchmark",
    data = list(
      benchmark_data = benchmark_data,
      pec_params     = pec_params),
    fun  = problem_wrapper,
    seed = 1012020)

  prob_df <- tibble::tibble(data_name = names(benchmark_data)) %>%
    dplyr::filter(!(data_name %in%
      c("pbc_tvf", "mgus2", "metabric", "churn"))) # these are estimated in separate benchmarks
  algo_df <- data.frame(budget = 20)

  addAlgorithm(name = "pam_xgb", fun = pam_xgb_wrapper)
  addAlgorithm(name = "orsf", fun = orsf_wrapper)
  addAlgorithm(name = "coxph", fun = cph_wrapper)

  addExperiments(
    prob.design = list(benchmark = prob_df),
    algo.designs = list(
      pam_xgb = algo_df,
      orsf    = algo_df,
      coxph   = algo_df),
    repls = 25L)

  # system.time({
    test1 <- testJob(id = 31)
  # })
  ids_bench <- findExperiments(prob.name = "benchmark")
  submitJobs(findNotDone(ids_bench))
  waitForJobs()

}

library(batchtools)
reg  <- loadRegistry("output/benchmark-registry", writeable = TRUE)
reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = 15)
# # reg$cluster.functions = makeClusterFunctionsInteractive()
# submitJobs(findNotDone())
# waitForJobs()
ids  <- findExperiments(prob.name = "benchmark")
pars <- unwrap(getJobPars()) %>% as_tibble()

## collect results
res <- reduceResultsDataTable(ids=findDone()) %>%
  as_tibble() %>%
  tidyr::unnest() %>%
  left_join(pars) %>%
  mutate(method = factor(method, levels = c("Reference", "coxph", "orsf", "pam_xgb")))

## save results
saveRDS(res, "res_benchmark.Rds")

## plot results
p_benchmark <- ggplot(res, aes(x = factor(quantile), y = IBS, fill = method)) +
  geom_boxplot() +
  facet_wrap(~data_name)

ggsave("benchmark.png", p_benchmark, width = 12, height = 9)
