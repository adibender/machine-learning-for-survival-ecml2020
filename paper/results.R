library(dplyr)
library(purrr)
library(batchtools)
library(ggplot2)
theme_set(theme_bw())

######################### comparison with ORSF #################################
# synthetic data
reg_synth <- loadRegistry("../simulation/output/tve-registry")
ids_tve <- findExperiments(prob.name = "syntheticTVE")
res_tve <- reduceResultsDataTable(ids_tve) %>% as_tibble() %>% unnest()
smry_tve <- res_tve %>%
  mutate(quantile = paste0("Q", quantile * 100)) %>%
  mutate(data_name = "synthetic (TVE)") %>%
  group_by(data_name, method, quantile) %>%
  summarise(value = mean(IBS, na.rm = TRUE)) %>%
  mutate(value = round(value, 3) * 100) %>%
  rename(data = data_name, algorithm = method) %>%
  ungroup() %>%
  mutate(
    algorithm = factor(algorithm,
      levels = c("Reference", "coxph", "orsf", "pam_xgb"),
      labels = c("Kaplan-Meier", "Cox-PH", "ORSF", "XGB (PEM)"))
  )

## benchmark
res_bench <- readRDS("../benchmark/res_benchmark.Rds")
res_smry <- res_bench %>%
  filter(data_name != "pbc") %>% # pbc evaluated separately (b/c of TVF)
  mutate(quantile = paste0("Q", quantile * 100)) %>%
  group_by(data_name, method, quantile) %>%
  summarise(value = mean(IBS, na.rm = TRUE)) %>%
  mutate(value = round(value, 3) * 100) %>%
  rename(data = data_name, algorithm = method) %>%
  ungroup() %>%
  mutate(algorithm = factor(algorithm,
    levels = c("Reference", "coxph", "orsf", "pam_xgb"),
    labels = c("Kaplan-Meier", "Cox-PH", "ORSF", "XGB (PEM)")))

## pbc tvf
res_pbc_tvf <- readRDS("res-ibs-pbc.Rds")
res_pbc_smry <- res_pbc_tvf %>%
  mutate(data_name = "pbc") %>%
  mutate(quantile = paste0("Q", quantile * 100)) %>%
  group_by(data_name, method, quantile) %>%
  summarise(value = mean(IBS, na.rm = TRUE)) %>%
  mutate(value = round(value, 3) * 100) %>%
  rename(data = data_name, algorithm = method) %>%
  ungroup() %>%
  mutate(algorithm = factor(algorithm,
    levels = c("Reference", "coxph", "orsf", "pam_xgb", "pam_xgb_tvf"),
    labels = c("Kaplan-Meier", "Cox-PH", "ORSF", "XGB (PEM)", "XGB (PEM, TVF)")))



library(tables)
tab_bench <- tabular(
    RowFactor(data) * Heading() * Factor(quantile) ~
    algorithm *  Heading() * value * Heading() * identity,
  data = rbind(res_smry, res_pbc_smry, smry_tve))
latex(tab_bench)


###################### comparison with deephit #################################
## metabric
res_metabric <- readRDS("results_metabric.RDS")
res_smry_meta <- res_metabric %>%
  group_by(algorithm, index, quantile, cause) %>%
  summarise(value = round(mean(value), 3) * 100) %>%
  mutate(data = "metabric")

## mgus 2
res_mgus <- readRDS("results_mgus2.RDS")
res_smry_mgus <- res_mgus %>%
  group_by(algorithm, index, quantile, cause) %>%
  summarise(value = round(mean(value), 3) * 100) %>%
  mutate(data = "MGUS 2")



## synthetic TVE CR
synth_tve_cr <- readRDS("../benchmark/results_syntheticTVE_CR.RDS")

smry  <- synth_tve_cr %>%
  filter(value < 1) %>% # couple of values for cph in iteration 12
  group_by(algorithm, index, quantile, cause) %>%
  summarise(value= round(mean(value), 3) * 100) %>%
  mutate(data = "synthetic (TVE, CR)")


## combine
res_smry <- do.call(rbind, list(smry, res_smry_meta, res_smry_mgus))
library(tables)
tab <- tabular(
  RowFactor(data) * index * algorithm ~
    Factor(cause) * Heading() * quantile * Heading() * value * Heading() * identity,
  data = res_smry)
latex(tab)


################################# scaling benchmark ############################
res_time <- readRDS("results-scaling-benchmark.Rds")

res_time_smry <- res_time %>%
  filter(quantile == .5) %>%
  group_by(n, strategy) %>%
  summarise(time = mean(time), IBS = mean(IBS)) %>%
  mutate(IBS = 100 * round(IBS, 3)) %>%
  mutate(time = round(time/3600, 2)) %>%
  tidyr::pivot_longer(cols=one_of(c("time", "IBS")), names_to="measure")

tab_time <- tabular(
  RowFactor(measure, "")*RowFactor(strategy) ~
    Factor(n)*Heading()*identity*Heading()*value, data = res_time_smry)
latex(tab_time)
