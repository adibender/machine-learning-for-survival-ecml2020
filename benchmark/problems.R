problem_wrapper <- function(data, job, data_name, split_frac = .7, nfold = 4,
  tvf_form = NULL) {

  instance <- data$benchmark_data[[data_name]]
  if (class(instance) == "list") {
    dat <- instance
    instance <- instance[[1]]
  }

  train_test_idx <- get_split_indices(nrow(instance), split_frac = split_frac)
  train <- instance[train_test_idx[["ind_train"]], , drop = FALSE]
  test <- instance[train_test_idx[["ind_eval"]], , drop = FALSE]

  if (exists("dat")) {
    out <- list(
      data = list(train = list(train, dat[[2]]), test = list(test, dat[[2]])),
      cv_indices = get_cv_indices(train$status, nfold = nfold)
    )
  } else {
    out <- list(
      data = list(train = train, test = test),
      cv_indices = get_cv_indices(train$status, nfold = nfold)
    )
  }

  out

}
