# setup
library(pammtools) # requires github version
library(dplyr)
library(mgcv)
library(pem.xgb)

library(ggplot2)
theme_set(theme_bw())

## daten simulation
n <- 500
# create data set with variables which will affect the hazard rate.
df <- cbind.data.frame(x1 = runif (n, -3, 3), x2 = runif (n, 0, 6)) %>%
 as_tibble()
# the formula which specifies how covariates affet the hazard rate
f0 <- function(t) {
 dgamma(t, 8, 2) * 6
}
form <- ~ - 3.5 + f0(t) - 0.5 * x1 + sqrt(x2)

sim_df <- sim_pexp(form, df, seq(0, 10, by = .2)) %>%
  mutate(time = round(time, 5))
head(sim_df)


####################### Step 1. Data Trafo #####################################
head(sim_df, 3)
# as_ped transforms data into PED format
ped <- sim_df %>% as_ped(formula = Surv(time, status) ~ .)

ped %>%
  group_by(id) %>%
  slice(1:2, n())


################## Step 2. Model Fit ###########################################

# fit PAM
pam <- pamm(ped_status ~ s(tend, k = 20) + x1 + s(x2), data = ped)

# fit xgboost model
xgb_params <- list(
  max_depth        = 3,
  eta              = .01,
  colsample_bytree = .7,
  min_child_weight = 10,
  subsample        = .7)

xgb_pam <- xgboost.ped(
  ped,
  params        = xgb_params,
  nrounds       = 500,
  print_every_n = 500)



################## Step 3. Postprocessing
# 3.1 model evaluation
get_ibs(
  list(pam, xgb_pam),
  data = sim_df[1:100,],
  exact = FALSE)

# 3.1. visualization
## add hazard and surv_prob for xgb model
pred_df  <- ped %>% make_newdata(tend = unique(tend), x1 = c(-.5, 1))
pred_df  <- pred_df %>%
  add_hazard2(xgb_pam) %>%
  add_hazard(pam)
## add hazard and surv_prob for pam
pred_df <- pred_df %>%
  group_by(x1) %>%
  add_surv_prob2(xgb_pam) %>%
  add_surv_prob(pam)
# add true hazard and survival probability
pred_df <- pred_df %>%
  mutate(
    truth_hazard = exp(-3.5 + f0(tend) -0.5 * x1 + sqrt(x2)),
    truth        = exp(-cumsum(truth_hazard * intlen)))

# Hazard
ggplot(pred_df, aes(x = tend)) +
  geom_stephazard(aes(y = hazard_pam_xgb), col = "green") +
  geom_stephazard(aes(y = hazard), col = "black") +
  geom_line(aes(y = truth_hazard), col = "blue", lty = 3) +
  facet_wrap(~x1)

# Survival probability
ggplot(pred_df, aes(x = tend)) +
  geom_line(aes(y = surv_prob_pam_xgb), col = "green") +
  geom_line(aes(y = surv_prob), col = "black") +
  geom_line(aes(y = truth), col = "blue", lty = 3) +
  scale_y_continuous(limits = c(0, 1)) +
  facet_wrap(~x1)
