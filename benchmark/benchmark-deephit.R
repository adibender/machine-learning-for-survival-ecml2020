############## Run Benchmark with DeepHit #############

source("calling_deephit.R")
library(pem.xgb)
library(purrr)


###### other functions for data prep #######
problem_wrapper <- function(data, job, data_name, split_frac = .7, nfold = 4) {

  instance <- data$benchmark_data[[data_name]]

  train_test_idx <- get_split_indices(nrow(instance), split_frac = split_frac)
  train <- instance[train_test_idx[["ind_train"]], , drop = FALSE]
  test <- instance[train_test_idx[["ind_eval"]], , drop = FALSE]

  list(
    data = list(train = train, test = test),
    cv_indices = get_cv_indices(train$status, nfold = nfold))

}

prep_data_for_deephit <- function(data_from_instance){

  this_data <- model.matrix(~ 0 + ., data =
                              data_from_instance[,setdiff(colnames(data_from_instance),
                                                          c("time", "status"))])

  make_data_for_deephit(matrix(as.integer(data_from_instance$time), ncol=1),
                        matrix(as.integer(data_from_instance$status), ncol=1),
                        this_data)

}


###################################################################

datasets <- list.files("instances/", full.names = T)

for(dataname in datasets[5]){

  instance_list <- readRDS(dataname)
  dtn <- gsub(".*\\/\\/(.*?)\\.Rds","\\1",dataname)
  if(dtn %in% c("breast", "syntheticTVE")){

    instance_list <- lapply(instance_list, function(il){

      # round times of breast and syntheticTVE
      il$data$train$time <- round(il$data$train$time, 2)*100
      il$data$test$time <- round(il$data$test$time, 2)*100

      return(il)

    })

  }

  cat("==========================================================\n")
  cat("==========================================================\n")
  cat("==========================================================\n")
  cat("================ RUNNING DATASET ", dtn, "=================\n")
  cat("==========================================================\n")
  cat("==========================================================\n")
  cat("==========================================================\n")

  for(simiter in 1:length(instance_list)){

    cat("==========================================================\n")
    cat("================= RUNNING ITERATION", simiter, "===================\n")
    cat("==========================================================\n")

    instance <- instance_list[[simiter]]

    train <- prep_data_for_deephit(instance$data$train)
    test <- prep_data_for_deephit(instance$data$test)
    split_inds <- instance$cv_indices

    q_last = 0.8
    q_eval = c(0.25, 0.5, 0.75)
    times <- times_list(instance$data$test, q_last = q_last, q_eval = q_eval)
    times_to_predict <- times$eval_times

    ret <- DeepHit_RS_wrapper(RS_budget = 20,
                              train = train,
                              test = test,
                              split_inds = split_inds,
                              times_to_predict = times_to_predict,
                              facCheck = ifelse(dtn=="syntheticCR",1000,100))

    saveRDS(ret, file=paste0("results/result_", dtn,
                             "_iter_", simiter, ".RDS"))

  }

}
