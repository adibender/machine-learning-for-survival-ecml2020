sim_wrapper <- function(
  data,
  job,
  n          = 500,
  time_grid  = seq(0, 10, by = 0.1),
  split_frac = .7,
  nfold      = 4) {

  # create data set with covariates
  df <- tibble::tibble(
    x0 = sample(c(-1,1), n, .3),
    x1 = runif(n, -3, 3),
    x2 = runif(n, -3, 3),
    x3 = runif(n, -3, 3))
  # baseline hazard
  f0 <- function(t, x0) { x0 * dgamma(t, 8, 2) * 6 }
  fx1 <- function(x1, t) { -0.1 * x1}
  fx2 <- function(x2, t) { dnorm(x2) * 10 * (1/(t + 1)) }
  fx3 <- function(x3, t) {2 * x3 * cos(t / max(time_grid) * pi)}

  instance <- pammtools::sim_pexp(
    formula = ~ -3.5 + f0(t, x0)  + fx1(x1, t) + fx2(x2, t) + fx3(x3, t),
    data    = df,
    cut     = time_grid)

  # noise variables
  df_noise <- purrr::map_dfc(1:20, ~runif(nrow(instance)), -3, 3)
  instance <- cbind(instance, df_noise)
  # add censoring
  cens_times <- runif(nrow(instance), 0, 20)
  status <- instance$time <= cens_times
  time <- pmin(instance$time, cens_times)
  status[time == 10] <- 0
  instance$time <- time
  instance$status <- status
  instance$id <- NULL
  instance$x0 <- as.factor(instance$x0)

  train_test_idx <- pem.xgb::get_split_indices(nrow(instance), split_frac = split_frac)
  train <- instance[train_test_idx[["ind_train"]], , drop = FALSE]
  test <- instance[train_test_idx[["ind_eval"]], , drop = FALSE]

  list(
    data = list(train = train, test = test),
    cv_indices = pem.xgb::get_cv_indices(train$status, nfold = nfold))

}

sim_wrapper_cr <- function(
  data,
  job,
  n          = 1000,
  time_grid  = seq(0, 10, by = 0.1),
  split_frac = .7,
  nfold      = 4) {

  # create data set with covariates
  df <- tibble::tibble(
    x0 = sample(c(-1,1), n, .3),
    x1 = runif(n, -3, 3),
    x2 = runif(n, -3, 3),
    x3 = runif(n, -3, 3))
  df2 <- mvtnorm::rmvnorm(n = nrow(df), mean=rep(0, 10))
  colnames(df2) <- paste0("x", 4:(ncol(df2)+3))
  eff <- sample(c(-1, 1), ncol(df2), replace = TRUE)
  df <- cbind(df, df2)
  # baseline hazard


  instance <- sim_pexp_cr(
    formula = ~ -3.5 + ( x0 * dgamma(t, 8, 2) * 6)  -0.1 * x1 +
      dnorm(x2) * 10 * (1/(t + 1)) + 2 * x3 * cos(t / max(seq(0, 10, by = 0.1)) * pi) |
      -3.5 + ( x0 * dgamma(t, 8, 2) * 6) + 2 * x4 -.1 * x5,
    data    = df,
    cut     = time_grid)
  instance$status <- NULL
  instance <- instance %>%
    rename(status = type)
  # add censoring
  cens_times <- runif(nrow(instance), 0, 20)
  cens <- instance$time > cens_times
  instance$time <- pmin(instance$time, cens_times)
  instance$status[cens] <- 0
  instance$status[instance$time == 10] <- 0
  instance <- as.data.frame(instance)
  instance$id <- instance$hazard1 <- instance$hazard2 <- NULL
  instance$x0 <- as.factor(instance$x0)

  train_test_idx <- pem.xgb::get_split_indices(nrow(instance), split_frac = split_frac)
  train <- instance[train_test_idx[["ind_train"]], , drop = FALSE]
  test <- instance[train_test_idx[["ind_eval"]], , drop = FALSE]

  list(
    data = list(train = train, test = test),
    cv_indices = pem.xgb::get_cv_indices(train$status, nfold = nfold))

}
