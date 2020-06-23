cph_wrapper <- function(
  data,
  job,
  instance,
  formula = Surv(time, status) ~.,
  x = TRUE,
  budget = 20) {

  train <- instance$data$train
  if (nrow(train) > ncol(train)) {

    mod <- coxph(formula, data = instance$data$train, x = x)
    ibs <- get_ibs(
      object         = list(coxph = mod),
      data           = instance$data$test,
      pec_params     = data$pec_params,
      keep_reference = FALSE)

  } else {

    ibs <- data.frame(
      method   = rep("coxph", 3),
      time     = quantile(train$time, c(.25, .5, .75)),
      quantile = c(.25, .5, .75),
      IBS      = rep(NA, 3))

  }

  ibs

}

pam_xgb_wrapper <- function(
  data,
  job,
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
  cut_subsamp = NULL,
  nrounds    = 5000L,
  eta        = .05,
  early_stopping_rounds = 50L,
  return_model = FALSE,
  nthread = 1L) {

  ## for data sets with huge p, reset colsample_bytree
  if(class(instance$data$train) == "list") {
    p <- ncol(instance$data$train[[1]])
  } else {
    p <- ncol(instance$data$train) - 2
  }

  if (p > 100) {
    param_space$colsample_bytree[2] <- sqrt(p) / p
    param_space$colsample_bytree[1] <- (sqrt(p) / p) / 10
  }

  param_df <- purrr::map_dfr(
    seq_len(budget),
    ~ pem.xgb::random_params(param_space))
  param_df$eta <- eta

  if(!is.null(cut_subsamp)) {
    samp_id <- sample(seq_len(nrow(instance$data$train)), cut_subsamp, replace = F)
    cut <- pammtools:::get_cut.default(instance$data$train[samp_id, ], Surv(time, status)~.)
    ped_params$cut <- cut
  }
  tune_res <- xgb.tune_pam(
    param_df              = param_df,
    data                  = instance$data$train,
    nrounds               = nrounds,
    cv_indices            = instance$cv_indices,
    split_frac            = split_frac,
    ped_params            = ped_params,
    verbose               = TRUE,
    print_every_n         = 500L,
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

  if(return_model) {
    return(pam_xgb)
  }

  ibs <- get_ibs(
    object         = list(pam_xgb = pam_xgb),
    data           = instance$data$test,
    pec_params     = data$pec_params,
    keep_reference = TRUE
  )
  attr(ibs, "best_set") <- best_set
  attr(ibs, "pam_xgb") <- pam_xgb

  ibs

}



orsf_wrapper <- function(
  data,
  job,
  instance,
  budget = 20,
  param_set = list(
    alpha                    = c(0, 1),
    gamma                    = c(.25, .75),
    min_events_to_split_node = c(5L, 20L),
    min_obs_to_split_node    = c(10L, 40L)
  )) {

  train <- instance$data$train

  t_med <- quantile(train$time[train$status == 1], .5)

  p <- ncol(train) - 2
  p2 <- ceiling(sqrt(p))
  param_set$mtry <- c(min(2, p2), min(2*p2, p))

  param_list <- purrr::map(seq_len(budget), function(i) random_params(param_set))
  tune_res <- purrr::map_dfr(
    param_list,
    ~{
      cv_res <- purrr::map(
        instance$cv_indices,
        function(cv_ind) {
          orsf_params <- .x
          orsf_params$data <- train[cv_ind$ind_train, , drop = FALSE]
          orsf_params$verbose <- FALSE
          do.call(obliqueRSF::ORSF, orsf_params)
        }
      )
      pec_list <- purrr::map2(
        cv_res,
        instance$cv_indices,
        function(res, cv_ind) {
          pec::pec(
            object = list(orsf = res),
            formula = Surv(time, status) ~ 1,
            times   = seq(.01, t_med, length.out = 500),
            data    = instance$data$train[cv_ind$ind_test, , drop = FALSE],
            start   = .01,
            exact   = FALSE
          )
        }
      )
      ibs_df <- purrr::map_dfr(
        .x = pec_list,
        .f = ~ as.data.frame(pec::crps(.x, times = t_med)) %>%
        dplyr::filter(method != "Reference") %>%
        dplyr::mutate(time = "Q50"),
        .id = "fold"
      )
    },
    .id = "id"
  )

  tune_res <- tune_res %>%
    dplyr::group_by(id) %>%
    dplyr::summarize(IBS = mean(IBS)) %>%
    dplyr::mutate(param_set = param_list)

  best_set <- tune_res %>%
    dplyr::filter(IBS == min(IBS))

  # final training
  params_orsf <- best_set$param_set[[1]]
  params_orsf$data <- instance$data$train
  params_orsf$verbose <- FALSE

  orsf <- do.call(obliqueRSF::ORSF, params_orsf)

  suppressMessages({
    ibs <- get_ibs(
      object         = list(orsf = orsf),
      data           = instance$data$test,
      pec_params     = data$pec_params,
      keep_reference = FALSE
    )
  })

  attr(ibs, "best_set") <- best_set
  attr(ibs, "orsf") <- orsf

  ibs

}


csc_wrapper <- function(data, job, instance, ...) {

  library(riskRegression)
  cnames <- colnames(instance$data$train)
  cnames <- cnames[!(cnames %in% c("time", "status"))]
  form <- as.formula(paste0("Hist(time, status)~", paste0(cnames, collapse = "+")))

  train <- instance$data$train
  data$pec_params$formula <- Hist(time, status) ~1
  if (nrow(train) > ncol(train)) {

    mod1 <- CSC(form, data = train, cause = 1)
    mod2 <- CSC(form, data = train, cause = 2)
    ibs <- get_ibs(
      object         = list(csc = mod1),
      data           = instance$data$test,
      pec_params     = data$pec_params,
      keep_reference = FALSE)

  } else {

    ibs <- data.frame(
      method   = rep("coxph", 3),
      time     = quantile(train$time, c(.25, .5, .75)),
      quantile = c(.25, .5, .75),
      IBS      = rep(NA, 3))

  }

  ibs

}
