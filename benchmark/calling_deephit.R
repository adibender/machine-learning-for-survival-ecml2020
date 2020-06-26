############## Random Search DeepHit Wrapper ############

# load reticulate
library("reticulate")

# create virtual environment in this folder using
#
# virtualenv deephit
# source deephit/bin/activate
# pip3 install numpy pandas tensorflow sklearn lifelines
#
# using the terminal. Then install the necessary packages

# # Set the path to the Python executable file
#
use_virtualenv("deephit")
# use_python("./deephit/bin/python", required = T)

# install packages
# -> actually works better using pip from terminal
# virtualenv_install("deephit",
#                    packages = c("numpy", "pandas", "tensorflow", "sklearn", "lifelines"),
#                    ignore_installed = T)

DH <- import("class_DeepHit")
ID <- import("import_data")
UE <- import("utils_eval")

source("deep_c_index.R")

# # data <- ID$import_dataset_SYNTHETIC()
# data <- ID$import_dataset_METABRIC()
# 
# # this is now a list of
# # - ncol(covariates)
# # - list with
# # - - covariates
# # - - time
# # - - label
# # - list of masks
# # - - mask1 (observation x number_events x number_categories)
# # - - mask2
# str(data,2)
# # however the model needs a reordering, namely
# data[[2]] <- data[[2]][c(1,3,2)] # covariates, label, time
# 
# actual_mask <- data[[3]]
# actual_data <- data[[2]]

library("tensorflow")
library("dplyr")
library("keras")

##### CREATE DEEPHIT NETWORK

make_data_for_deephit <- function(time, label, data){

  num_Category    = as.integer(round(max(time) * 1.2))
  num_Event       = as.integer(round(length(unique(label)) - 1))

  x_dim           = dim(data)[2]

  mask1           = ID$f_get_fc_mask2(time, label, num_Event, num_Category)
  mask2           = ID$f_get_fc_mask3(time, -1, num_Category)

  list(x_dim, list(data, label, time), list(mask1, mask2))

}

DeepHit_RS_wrapper <- function(
  lr = 1e-4,
  hyperparams = list(h_dim_shared = c(50L, 100L, 200L, 300L),
                     num_layers_shared = c(1L,2L,3L,5L),
                     h_dim_CS = c(50L, 100L, 200L, 300L),
                     num_layers_CS = c(1L,2L,3L,5L),
                     active_fn = c("relu", "elu", "tanh"),
                     mb_size = c(32L, 64L, 128L),
                     iteration = c(50000),
                     keep_prob = c(0.6),
                     lr_train = c(lr),
                     alpha =  c(1),
                     beta =  c(0.1, 0.5, 1, 3, 5),
                     gamma = c(0)),
  RS_budget = 50,
  train, # preprate_data_for_deephit
  test, # preprate_data_for_deephit
  eval_time = NULL,
  this_seed = 1337,
  split_inds,
  nr_splits = length(split_inds),
  times_to_predict,
  facCheck = 100
)
{

  set.seed(this_seed)

  # use first in case each parameter is a list of values
  # hyperparams <- lapply(hyperparams, function(x) x[1])
  hyperparams_all_combs <- expand.grid(hyperparams, stringsAsFactors = FALSE)
  hyperparams_RS_index <- sample(1:nrow(hyperparams_all_combs), RS_budget)
  best_comb <- hyperparams_all_combs[1,]
  hyperparams_all_combs <- hyperparams_all_combs[hyperparams_RS_index,]
  
  data <- train[[2]]
  mask <- train[[3]]
  
  global_max_valid = -99
  global_best_iter = 0

  x_dim = NCOL(data[[1]])
  num_Event = dim(mask[[1]])[2]
  num_Category  = dim(mask[[1]])[3]
  
  for(j in 1:RS_budget){
  
    hyperparams <- hyperparams_all_combs[j,]
    parameters = hyperparams[c("alpha","beta","gamma")]
    
    cat( "PARAMETER SETTING ", j, ": ","\n")
    print(hyperparams)
    cat(" ==================================\n")
    
    this_setting_max_valid <- 0
    this_global_best_iter <- 0
    
    for(k in 1:nr_splits){
      
      idx_train <- split_inds[[k]]$ind_train
      idx_val <- split_inds[[k]]$ind_test
      
      # reset TF graph
      tf$compat$v1$reset_default_graph()
      config = tf$compat$v1$ConfigProto()
      sess = tf$compat$v1$Session(config=config)
      
      model = DH$Model_DeepHit(sess, "DeepHit",
                               input_dims = list('x_dim' = x_dim,
                                                 'num_Event' = num_Event,
                                                 'num_Category'  = num_Category),
                               network_settings =
                                 c(hyperparams[c('h_dim_shared',
                                                 'num_layers_shared',
                                                 'h_dim_CS',
                                                 'mb_size',
                                                 'num_layers_CS')],
                                   list('active_fn' = switch(hyperparams$active_fn,
                                                             relu = tf$nn$relu,
                                                             elu = tf$nn$elu,
                                                             tanh = tf$nn$tanh)),
                                   list('initial_W' = 
                                          tf$contrib$layers$xavier_initializer())))
      
      saver = tf$train$Saver()
      
      sess$run(tf$global_variables_initializer())
      
      this_val_data <- lapply(data, function(d) d[idx_val, , drop=F])
      this_val_mask <- lapply(mask, function(m){ 
        if(length(dim(m))==3) m[idx_val,,,drop=F] else
        m[idx_val,]})
      
      this_train_data <- lapply(data, function(d) d[idx_train, , drop=F])
      this_train_mask <- lapply(mask, function(m){ 
        if(length(dim(m))==3) m[idx_train,,,drop=F] else
        m[idx_train,]})
      
      if(is.null(eval_time)) eval_time <- quantile(this_train_data[[3]], c(0.25, 0.5, 0.75))
      
      max_valid = 0
      best_iter = 0
      stop_flag = 0
      
      cat( "MAIN TRAINING ...")
      cat("EVALUATION TIMES: ", paste(eval_time, collapse = ", "),"\n")
      
      avg_loss = 0
      for(i in 1:hyperparams$iteration){
        
        if(stop_flag > 5){ #for faster early stopping
          break
        }else{
          
          idx <- sample(1:nrow(this_train_data[[1]]), hyperparams$mb_size)
          this_data <- lapply(this_train_data, function(d) d[idx, , drop=F])
          this_mask <- lapply(this_train_mask, function(m)
            if(length(dim(m))==3) m[idx,,,drop=F] else m[idx,])
          
          loss_curr =  model$train(this_data, this_mask, unlist(parameters), 
                                   hyperparams[['keep_prob']],
                                   hyperparams[['lr_train']])
          avg_loss = avg_loss +  loss_curr[[2]]/facCheck
          
          if(i%%facCheck==0){
            
            cat('|| ITR: ', i, ' | Loss: ', avg_loss,"\n")
            avg_loss = 0
            
            pred = model$predict(this_val_data[[1]])
            
            tmp_valid <- weighted_c_index(
              pred,
              this_train_data[[3]],
              this_train_data[[2]][,1],
              this_val_data[[3]], #va_time
              this_val_data[[2]][,1],
              eval_time = eval_time)
            
            if(tmp_valid >  max_valid){
              stop_flag = 0
              max_valid = tmp_valid
              best_iter = i
            
            }else{
              stop_flag <- stop_flag + 1
            }
            
            
          }
          
        }
        
      }
      
      this_setting_max_valid <- this_setting_max_valid + max_valid
      this_global_best_iter <- this_global_best_iter + best_iter
      
    }
    
    cat( 'updated.... average c-index = ', tmp_valid,"\n")
    if(this_setting_max_valid/nr_splits > global_max_valid){
      global_best_iter <- this_global_best_iter/nr_splits
      global_max_valid <- this_setting_max_valid/nr_splits
      best_comb <- hyperparams
      cat( 'updated.... best c-index = ', tmp_valid,"\n")
    }
    
    
  }
  
  tf$compat$v1$reset_default_graph()
  config = tf$compat$v1$ConfigProto()
  sess = tf$compat$v1$Session(config=config)
  
  cat(" ================================================\n")
  cat(" REFITTING THE BEST MODEL FOUND IN RANDOM SEARCH...\n")
  
  data <- train[[2]]
  mask <- train[[3]]
  
  x_dim = NCOL(data[[1]])
  num_Event = dim(mask[[1]])[2]
  num_Category  = dim(mask[[1]])[3]
  
  hyperparams <- best_comb
  hyperparams$iterations <- global_best_iter
  
  parameters = hyperparams[c("alpha","beta","gamma")]
  
  model = DH$Model_DeepHit(sess, "DeepHit",
                           input_dims = list('x_dim' = x_dim,
                                             'num_Event' = num_Event,
                                             'num_Category'  = num_Category),
                           network_settings =
                             c(hyperparams[c('h_dim_shared',
                                             'num_layers_shared',
                                             'h_dim_CS',
                                             'mb_size',
                                             'num_layers_CS')],
                               list('active_fn' = switch(hyperparams$active_fn,
                                                         relu = tf$nn$relu,
                                                         elu = tf$nn$elu,
                                                         tanh = tf$nn$tanh)),
                               list('initial_W' = tf$contrib$layers$xavier_initializer())))
  
  saver = tf$train$Saver()
  sess$run(tf$global_variables_initializer())
  
  train_loss <- c()
  
  for(i in 1:global_best_iter){
    train_loss =  model$train(data, mask, unlist(parameters), hyperparams[['keep_prob']],
                                 hyperparams[['lr_train']])
  }
  
  pred <- model$predict(test[[2]][[1]])
  
  res <- list(
    pred = pred,
    T_train = data[[3]],
    Y_train = data[[2]],
    T_test = test[[2]][[3]], #va_time
    Y_test = test[[2]][[2]],
    eval_time = times_to_predict, 
    also_brier = T,
    apply_fun = function(x)x)
  
  return(res)
  
}

