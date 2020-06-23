library(magrittr)
library(pammtools)

## prepare data sets

# pbc tvf
data("pbc", package = 'survival')
pbc_tvf <- pbc %>% dplyr::filter(.data$id <= 312) %>%
  dplyr::mutate_if(is.numeric, ~ifelse(is.na(.), mean(., na.rm = TRUE), .)) %>%
  dplyr::select(-ascites, -hepato, -spiders) %>%
  dplyr::filter(!is.na(stage)) %>%
  dplyr::mutate(status = (.data$status == 2) * 1L)
pbcseq <- pbcseq %>%
  dplyr::select(id, day, bili, albumin, ast, protime)
pbc_tvf <- list(pbc_tvf, pbcseq)



data("Breast", package = 'biospear')
names(Breast)[4:ncol(Breast)]=paste('var',4:ncol(Breast),sep = '_')


data('GBSG2', package = 'pec')
GBSG2 <- GBSG2 %>%
  dplyr::rename(status = cens) %>%
  dplyr::mutate(tgrade = factor(tgrade, ordered = FALSE))


data("tumor", package = "pammtools")
tumor <- tumor %>%
  dplyr::rename(
    time = days,
    status = status,
  )

## customer churn
churn <- readr::read_csv("../datasets/WA_Fn-UseC_-Telco-Customer-Churn.csv")
churn <- churn %>%
  dplyr::mutate_if(is.character, as.factor) %>%
  select(-customerID, -TotalCharges, -MonthlyCharges, -Contract) %>%
  rename(time = tenure, status = Churn) %>%
  mutate(status = 1L * (status == "Yes")) %>%
  filter(time != 0)
churn <- as.data.frame(churn)



## DeepHit data sets
# Metabric (single event), available from DeepHit Repo
mb <- readr::read_csv("../datasets/cleaned_features_final.csv")
mb_outcome <- readr::read_csv("../datasets/label.csv") %>%
  dplyr::rename(time = event_time, status = label)
mb <- cbind(mb_outcome, mb)

# mgus2
data("mgus", package = "survival")
mgus2$time <- with(mgus2, ifelse(pstat == 0, futime, ptime))
mgus2$status <- with(mgus2, ifelse(pstat == 0, 2 * death, 1))
mgus2 <- mgus2 %>%
  select(-id, -ptime, -futime, -death, -pstat) %>%
   dplyr::mutate_if(is.numeric, ~ifelse(is.na(.), mean(., na.rm = TRUE), .))




## combine and save
data_list <- list(
  pbc_tvf   = pbc_tvf,
  breast    = Breast,
  gbsg2     = GBSG2,
  tumor     = tumor,
  metabric  = mb,
  mgus2     = mgus2
)


saveRDS(data_list, "benchmark_data.Rds")
