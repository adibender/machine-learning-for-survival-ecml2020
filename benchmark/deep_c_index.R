########### R version of the C-index implemented in DeepHit #############

library("reticulate")
use_virtualenv("deephit")
UE <- import("utils_eval")

weighted_c_index <- function(pred, # matrix [#test_obs. x events x times]
                             T_train, # matrix [train_time x 1]
                             Y_train, # matrix [train_ind x 1] with events
                             T_test, # matrix [test_time x 1]
                             Y_test, # matrix [test_ind x 1] with events
                             eval_time, # times for evaluation
                             num_Event = dim(pred)[2], # number of modelled time
                             num_Category = dim(pred)[3],
                             also_brier = FALSE,
                             apply_fun = mean
                             ){


  va_result1 = matrix(0, nrow = num_Event, ncol = length(eval_time))
  va_result2 = matrix(0, nrow = num_Event, ncol = length(eval_time))

  for(t in 1:length(eval_time)){

    eval_horizon = round(eval_time[t])

    if(eval_horizon > num_Category){
      stop('ERROR: evaluation horizon is out of range\n')
    }else{
      if(dim(pred)[2]==1)
        risk = rowSums(pred[,1,1:(eval_horizon),drop=F]) else
          risk = rowSums(pred[,,1:(eval_horizon),drop=F], dims=2) #risk score until eval_time
        for(k in 1:num_Event){

          if(num_Event==1) this_risk <- matrix(risk, ncol=1) else this_risk <- risk[,k,drop=F]

          va_result1[k, t] =
            UE$weighted_c_index(T_train , #time
                                matrix((Y_train == k)*1,ncol=1), #label
                                this_risk,
                                T_test,
                                matrix((Y_test == k)*1,ncol=1),
                                eval_horizon)
          
          va_result2[k,t] = 
            UE$weighted_brier_score(T_train , #time
                                    matrix((Y_train == k)*1,ncol=1), #label
                                    this_risk,
                                    T_test,
                                    matrix((Y_test == k)*1,ncol=1),
                                    eval_horizon)

        }
    }
  }

  if(!also_brier) return(apply_fun(va_result1)) else return(list(apply_fun(va_result1),
                                                                 apply_fun(va_result2)))

}
