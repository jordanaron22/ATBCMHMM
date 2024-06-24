################## Intro ################## 

library(Rcpp)
library(RcppArmadillo)
library(matrixStats)
library(MASS)
library(tidyverse)

RE_type <- "norm"

real_data <- T
set_seed <- T

sim_num <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
if (is.na(sim_num)){sim_num <- 9}
if (set_seed){set.seed(sim_num)}

RE_num <- as.numeric(commandArgs(TRUE)[1])
# sim_size <- as.numeric(commandArgs(TRUE)[2])
# RE_type <- as.character(commandArgs(TRUE)[3])
# print(paste("Sim Seed:",sim_num,"Size",sim_size,"RE type",RE_type,"Clust Num:",RE_num))
print(paste("Sim Seed:",sim_num,"HMM Num:",RE_num))

if(is.na(RE_num)){RE_num <- 5}
# if(is.na(sim_size)){sim_size <- 0}
# if(is.na(RE_type)){RE_type <- "norm"}

################## Functions ################## 

readCpp <- function(path) {
  tryCatch(
    {
      sourceCpp(file = path)
    },
    error = function(cond) {
      message("Wrong environment")
      # Choose a return value in case of error
      NA
    },
    warning = function(cond) {
      message("Wrong environment")
      # Choose a return value in case of warning
      NULL
    },
    finally = {
      message("Done")
    }
  )
}


ChooseCovar <- function(covar_vec){
  return(which(covar_vec == 1))
}

SimulateMC <- function(day_length,init,tran,covar_ind){
  hidden_states <- numeric(day_length)
  for (i in 1:day_length){
    
    if (i == 1) {
      hidden_states[1] <- rbinom(1,1,init[covar_ind,2])
    } else {
      hidden_states[i] <- rbinom(1,1,tran[hidden_states[i-1] + 1,2])
    }
  }
  
  
  return(hidden_states)
}

SimulateHMM <- function(num_of_people,init,tran_array,emit,pi_l,missing_perc){
  
  day_length <- 25
  mixture_num <- dim(emit)[3]
  covar_mat <- t(rmultinom(num_of_people,1,pi_l))
  covar_vec <- apply(covar_mat,1,ChooseCovar)
  
  depression_array <- array(NA,dim = c(day_length,num_of_people,7))
  hidden_states_matrix <- matrix(NA,ncol = num_of_people,nrow = day_length)
  
  for (ind in 1:num_of_people){
    depression <- matrix(NA,25,7)
    covar_ind <- covar_vec[ind]
    hidden_states <- SimulateMC(day_length,init,tran_array[,,covar_ind],covar_ind)
    
    for (i in 1:day_length){
      
      probd <- emit[hidden_states[i] + 1,2,] 
      depression[i,] <- rbinom(7,1,probd)
    }
    
    dep_missing_ind <- rbinom(day_length,1,missing_perc)
    depression[dep_missing_ind==1,] <- NA
    
    depression_array[,ind,] <- depression
    hidden_states_matrix[,ind] <- hidden_states
    

  }
  
  return(list(hidden_states_matrix,depression_array,covar_vec))
}

CondMarginalize <- function(alpha,beta,pi_l){
  alpha_beta <- simplify2array(alpha) + simplify2array(beta)
  
  
  for (ind in 1:dim(alpha_beta)[4]){
    for (re_ind in 1:dim(alpha_beta)[3]){
      alpha_beta[,,re_ind,ind] <- alpha_beta[,,re_ind,ind] + log(pi_l[re_ind])
    }
  }
  
  ind_like_mat <- apply(alpha_beta,c(1,4),logSumExp)
  
  weight_array_wake <- array(0, dim = c(dim(alpha_beta)[1],dim(alpha_beta)[4],dim(alpha_beta)[3]))
  weight_array_sleep <- array(0, dim = c(dim(alpha_beta)[1],dim(alpha_beta)[4],dim(alpha_beta)[3]))
  for (ind in 1:dim(alpha_beta)[4]){
    for (t in 1:dim(alpha_beta)[1]){
      weight_array_wake[t,ind,] <- alpha_beta[t,1,,ind] - ind_like_mat[t,ind]
      weight_array_sleep[t,ind,] <- alpha_beta[t,2,,ind] - ind_like_mat[t,ind]
    }
  }
  
  return(list(weight_array_wake,weight_array_sleep))
}

CalcInit <- function(alpha, beta,pi_l,log_sweights_vec){
  
  num_obs <- dim(alpha[[1]][,,1])[1]
  time <- 1
  init_0_vec <- matrix(0,length(alpha),length(pi_l))
  init_1_vec <- matrix(0,length(alpha),length(pi_l))
  init_mat <- matrix(0,length(pi_l),2)
  
  ind_like_vec <- CalcLikelihoodIndVec(alpha,pi_l)
  
  
  for(ind in 1:length(alpha)){ 
    ind_like <- ind_like_vec[ind]
    
    init_0_vec[ind,] <- alpha[[ind]][time,1,] + beta[[ind]][time,1,] + log(pi_l) - ind_like + log_sweights_vec[ind]
    init_1_vec[ind,] <- alpha[[ind]][time,2,] + beta[[ind]][time,2,] + log(pi_l) - ind_like + log_sweights_vec[ind]
    
  }
  
  for (re_ind in 1:length(pi_l)){
    init_0 <- logSumExp(init_0_vec[,re_ind])
    init_1 <- logSumExp(init_1_vec[,re_ind])
    init_vec <- exp(c(init_0,init_1) - logSumExp(c(init_0,init_1)))
    init_mat[re_ind,] <- init_vec
  }
  
  return(init_mat)
  
}

CalcProbRE <- function(alpha,pi_l){
  
  len <- dim(alpha[[1]])[1]
  re_len <- dim(alpha[[1]])[3]
  re_weight_vec <- numeric(re_len)
  re_weights <- matrix(0,nrow = length(alpha),ncol = re_len)
  
  
  for (ind in 1:length(alpha)){
    for (re_ind in 1:re_len){
      re_weights[ind,re_ind] <- logSumExp(alpha[[ind]][len,,re_ind]) + log(pi_l[re_ind])
    }
    re_weights[ind,] <- exp(re_weights[ind,] - logSumExp(c(re_weights[ind,])))
    
  }
  
  return(re_weights)
  
}

CalcTranC <- function(alpha,beta,dep,emit,pi_l){
  
  len <- dim(dep)[1]
  re_num <- dim(alpha[[1]])[3]
  tran_array_working <- array(0,c(2,2,re_num))
  
  ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))
  
  tran_vals_re_00 <- CalcTranHelper(init_state = 0,new_state = 0,depression = dep,
                                     tran_array = tran_array,emit_nd_vec = emit[1,1,],
                                     ind_like_vec = ind_like_vec,
                                     alpha = alpha,beta = beta, pi_l = pi_l)  
  
  tran_vals_re_01 <- CalcTranHelper(init_state = 0,new_state = 1,depression = dep,
                                     tran_array = tran_array,emit_nd_vec = emit[2,1,],
                                     ind_like_vec = ind_like_vec,
                                     alpha = alpha,beta = beta, pi_l = pi_l)
  
  tran_vals_re_10 <- CalcTranHelper(init_state = 1,new_state = 0,depression = dep,
                                     tran_array = tran_array,emit_nd_vec = emit[1,1,],
                                     ind_like_vec = ind_like_vec,
                                     alpha = alpha,beta = beta, pi_l = pi_l)  
  
  tran_vals_re_11 <- CalcTranHelper(init_state = 1,new_state = 1,depression = dep,
                                     tran_array = tran_array,emit_nd_vec = emit[2,1,],
                                     ind_like_vec = ind_like_vec,
                                     alpha = alpha,beta = beta, pi_l = pi_l)  
  
  for (init_state in 1:2){
    for (new_state in 1:2){
      
      
      
      if (init_state == 1 & new_state == 1){tran_vals <- tran_vals_re_00}
      if (init_state == 1 & new_state == 2){tran_vals <- tran_vals_re_01}
      if (init_state == 2 & new_state == 1){tran_vals <- tran_vals_re_10}
      if (init_state == 2 & new_state == 2){tran_vals <- tran_vals_re_11}
      
      for(re_ind in 1:re_num){
        tran_array_working[init_state,new_state,re_ind] <- sum(tran_vals[,,re_ind])

      }
    }
  }
  
  
  for (re_ind in 1:re_num){
    for (i in 1:2){
      row_sum <- sum(tran_array_working[i,,re_ind])
      for (j in 1:2){
        tran_array_working[i,j,re_ind] <- tran_array_working[i,j,re_ind] / row_sum
      }
    }
  }
  
  
  
  return(tran_array_working)
}


CalcEmit <- function(alpha,beta,dep,emit,pi_l){
  
  n <- dim(dep)[2]
  len <- dim(dep)[1]
  re_num <- dim(alpha[[1]])[3]
  emit_working <- array(0,c(2,2,7))
  
  ind_like_vec <- unlist(lapply(c(1:length(alpha)),IndLike,alpha = alpha, pi_l = pi_l, len = len))
  
  for (ind in 1:n){
    alpha_beta <- alpha[[ind]] + beta[[ind]] - ind_like_vec[ind]

    for (data_ind in 1:dim(emit_working)[3]){
      
      dep_vec <- dep[,ind,data_ind]
      zinds <- which(dep_vec == 0)
      oinds <- which(dep_vec == 1)
      
      for (mc_state in 1:2){
        for (re_ind in 1:re_num){
              
          emit_working[mc_state,1,data_ind] <- logSumExp(c(emit_working[mc_state,1,data_ind],alpha_beta[zinds,mc_state,re_ind] + log(pi_l[re_ind])))
          emit_working[mc_state,2,data_ind] <- logSumExp(c(emit_working[mc_state,2,data_ind],alpha_beta[oinds,mc_state,re_ind] + log(pi_l[re_ind])))
          
        }
      }
    }
  }
  
  emit_working <- exp(emit_working)
  
  for (data_ind in 1:dim(emit_working)[3]){
    for (i in 1:2){
      row_sum <- sum(emit_working[i,,data_ind])
      for (j in 1:2){
        emit_working[i,j,data_ind] <- emit_working[i,j,data_ind] / row_sum
      }
    }
  }
  
  return(emit_working)
  
}
  
CalcLikelihood <- function(alpha,pi_l){
  return(sum(CalcLikelihoodIndVec(alpha,pi_l)))
}

CalcLikelihoodIndVec <- function(alpha,pi_l){
  num_obs <- dim(alpha[[1]][,,1])[1]
  like_vec <- c()
  #i is number of people
  for (i in 1:length(alpha)){
    ind_like <- logSumExp(c(SumOverREIndTime(alpha,pi_l,i,num_obs)))
    like_vec <- c(like_vec,ind_like)
  }
  return(like_vec)
}

SumOverREIndTime <- function(fb,pi_l,ind,time, add_re = T){
  
  fb_ind <- fb[[ind]]
  
  fb_sum <- numeric(2)
  if (add_re){
    fb_sum[1] <- logSumExp(c(fb_ind[time,1,] + log(pi_l)))
    fb_sum[2] <- logSumExp(c(fb_ind[time,2,] + log(pi_l)))
  } else {
    fb_sum[1] <- logSumExp(c(fb_ind[time,1,]))
    fb_sum[2] <- logSumExp(c(fb_ind[time,2,]))
  }
  
  return(fb_sum)
}

IndLike <- function(alpha,pi_l,ind,len){
  likelihood <- logSumExp(SumOverREIndTime(alpha,pi_l,ind,len))
  return(likelihood)
}

UniqueIDHelp <- function(num){
  return(c(1:num))
}

CalcBIC <- function(new_likelihood,RE_num,dep){
  
  #init tran emit pi
  num_of_param <- RE_num + RE_num*2 + 7*2 + (RE_num-1)

  bic <- num_of_param * log(sum(!is.na(dep))) - (2 * new_likelihood)
  return(bic)
}

Forward <- function(dep,init,tran_array,emit_nd0_vec,emit_nd1_vec){
  alpha <- list()
  len <- dim(dep)[1]
  n <- dim(dep)[2]
  num_re <- dim(tran_array)[3]
  
  for (ind in 1:n){
    alpha_array <- array(NA,dim = c(len,2,num_re))
    dep_ind <- dep[,ind,]
    
    for (clust_i in 1:num_re){
      tran <- tran_array[,,clust_i]
      init_vec <- init[clust_i,]
      
      alpha_array[,,clust_i] <- ForwardIndC(dep_ind,init_vec,tran,emit_nd0_vec,emit_nd1_vec)
    }
    
    alpha[[ind]] <- alpha_array
  }
  
  return(alpha)
}

Backward <- function(dep,tran_array,emit_nd0_vec,emit_nd1_vec){
  beta <- list()
  len <- dim(dep)[1]
  n <- dim(dep)[2]
  num_re <- dim(tran_array)[3]
  
  for (ind in 1:n){
    beta_array <- array(NA,dim = c(len,2,num_re))
    dep_ind <- dep[,ind,]
    
    for (clust_i in 1:num_re){
      tran <- tran_array[,,clust_i]
      
      beta_array[,,clust_i] <- BackwardIndC(dep_ind,tran,emit_nd0_vec,emit_nd1_vec)
    }
    
    beta[[ind]] <- beta_array
  }
  
  return(beta)
}

CalcTranHelper <- function(init_state, new_state,
                           depression, tran_array, 
                           emit_nd_vec, ind_like_vec,
                           alpha, beta, pi_l){
  
  init_state <- init_state + 1
  new_state <- new_state + 1
  
  num_of_people <- dim(depression)[2]
  len <- dim(depression)[1]
  num_re <- dim(tran_array)[3]
  
  tran_val_re_array <- array(NA, dim = c(len-1, num_of_people,num_re))
  
  for (clust_i in 1:num_re){
    tran_vals_re_mat <- matrix(0,len-1,num_of_people)
    tran_val <- tran_array[init_state,new_state,clust_i]
    
    for (ind in 1:num_of_people){
      alpha_ind <- alpha[[ind]]
      beta_ind <- beta[[ind]]
      likelihood <- ind_like_vec[ind]
      
      depression_ind_m1 <- depression[2:len,ind,]
      class_vec <- logClassificationC(depression_ind_m1,emit_nd_vec)
      
      alpha_ind_slice <- alpha_ind[1:(len-1),init_state,clust_i]
      beta_ind_slice <- beta_ind[2:len,new_state,clust_i]
      
      tran_vals_re_ind <- exp(alpha_ind_slice + beta_ind_slice + log(tran_val) + log(pi_l[clust_i]) + class_vec - likelihood)
      
      tran_vals_re_mat[,ind] <- tran_vals_re_ind
    }
    
    tran_val_re_array[,,clust_i] <- tran_vals_re_mat;
  }
  
  return(tran_val_re_array)
}

ViterbiInd <- function(dep_ind,est_params,assignment){
  
  init_ind <- est_params[[1]][assignment,]
  tran_ind <- est_params[[2]][,,assignment]
  emit <- est_params[[3]]
  
  len <- dim(dep_ind)[1]
  
  lc_nd <- logClassificationC(dep_ind,emit[1,1,])
  lc_d <- logClassificationC(dep_ind,emit[2,1,])
  
  viterbi_mat <- matrix(NA,2,len)
  viterbi_mat[1,1] <- log(init_ind[1]) + lc_nd[1]
  viterbi_mat[2,1] <- log(init_ind[2]) + lc_d[1]
  
  traceback <- matrix(NA,2,len)
  
  for (time in 2:len){
    
    
    viterbi_mat[1,time] <- lc_nd[time]+ max(viterbi_mat[1,time-1] + log(tran_ind[1,1]),
                                            viterbi_mat[2,time-1] + log(tran_ind[2,1]))
    
    
    viterbi_mat[2,time] <- lc_d[time] + max(viterbi_mat[1,time-1] + log(tran_ind[1,2]),
                                            viterbi_mat[2,time-1] + log(tran_ind[2,2]))
    
    
    traceback[1,time] <-  which.max(c(viterbi_mat[1,time-1] + log(tran_ind[1,1]),
                                      viterbi_mat[2,time-1] + log(tran_ind[2,1])))
    
    
    traceback[2,time] <- which.max(c(viterbi_mat[1,time-1] + log(tran_ind[1,2]),
                                     viterbi_mat[2,time-1] + log(tran_ind[2,2])))
    
    
  }
  
  
  best_path <- numeric(len)
  best_path[len] <- which.max(viterbi_mat[, len])
  for (time in (len - 1):1) {
    best_path[time] <- traceback[best_path[time + 1], time + 1]
  }
  
  
  return(best_path-1)
} 


Viterbi <- function(dep,est_params,assign_vec){
  
  len <- dim(dep)[1]
  n <- dim(dep)[2]
  decoded_dep = matrix(NA,len,n)
  
  for (ind in 1:n){
    decoded_dep[,ind] <- ViterbiInd(dep[,ind,],est_params,assign_vec[ind])
  }
  
  return(decoded_dep)
}

################## EM Setup ################## 

readCpp( "cFunctions.cpp" )
readCpp( "../Rcode/cFunctions.cpp" )

print("loaded cfunctions")
if (!real_data){
  
  init_true <- matrix(NA,ncol = 2,nrow = RE_num)
  init_true[,1] <- seq(.95,.1,length.out=RE_num)
  init_true[,2] <- 1 - init_true[,1]

  emit_true <- array(NA, c(2,2,7))
  emit_true[1,1,] <- seq(.9,.6,length.out = 7)
  emit_true[1,2,] <- 1 - emit_true[1,1,]
  emit_true[2,2,] <- seq(.6,.9,length.out = 7)
  emit_true[2,1,] <- 1 - emit_true[2,2,]
  
  tran_array_true <- array(NA, c(2,2,3))
  tran_array_true[,,1] <- matrix(c(.85,.15,.3,.7),2,2,byrow = T)
  tran_array_true[,,2] <- matrix(c(.6,.4,.45,.55),2,2,byrow = T)
  tran_array_true[,,3] <- matrix(c(.4,.6,.65,.35),2,2,byrow = T)
  
  pi_l_true <- c(.28,.32,.4)
  
  
  simulated_hmm <- SimulateHMM(2000,init_true,tran_array_true,emit_true,pi_l_true,.2)
  
  mc <- simulated_hmm[[1]]
  dep <- simulated_hmm[[2]]
  covar_vec <- simulated_hmm[[3]]
  
  
} else {
  
  atbc_data <- read.csv("nATBCAnalyticData_HMM_All_Updated.csv")
  id_reps <- c(table(atbc_data$ID))
  unique_id_vec <- unlist(lapply(id_reps,UniqueIDHelp))
  atbc_data <- atbc_data %>% mutate(uniqueID = unique_id_vec)
  
  atbc_wide <- atbc_data %>% pivot_wider(names_from = c(uniqueID), values_from = c(Anxiety,Depression,Memory,Concentrate,Fatigue,Appetite,Insomnia))
  
  id_vec <- atbc_wide[,1]
  atbc_wide <- atbc_wide[,-1]
  
  dep <- array(NA, dim = c(25,29133,7))
  
  dep[,,1] <- t(atbc_wide[,1:25])
  dep[,,2] <- t(atbc_wide[,26:50])
  dep[,,3] <- t(atbc_wide[,51:75])
  dep[,,4] <- t(atbc_wide[,76:100])
  dep[,,5] <- t(atbc_wide[,101:125])
  dep[,,6] <- t(atbc_wide[,126:150])
  dep[,,7] <- t(atbc_wide[,151:175])
}


###### Initial Settings ###### 
init <- matrix(NA,ncol = 2,nrow = RE_num)
init[,1] <- seq(.75,.25,length.out=RE_num)
init[,2] <- 1 - init[,1]

pi_l <- rep(1/RE_num,RE_num)


emit <- array(NA, c(2,2,7))
emit[1,1,] <- .95
emit[1,2,] <- .05
emit[2,2,] <- .95
emit[2,1,] <- .05

tran_array <- array(NA, c(2,2,RE_num))
tran_array[1,1,] <- seq(.95,.65,length.out=RE_num)
tran_array[1,2,] <- 1 - tran_array[1,1,]
tran_array[2,2,] <- seq(.95,.65,length.out=RE_num)
tran_array[2,1,] <- 1 - tran_array[2,2,]
dim(tran_array)[3] <- RE_num

log_sweights_vec <- numeric(dim(dep)[2])


# init <- init_true
# emit <- emit_true
# tran_array <- tran_array_true
# pi_l <- pi_l_true

################## EM ################## 
time_vec <- c()

print("pre_alpha")
alpha <- Forward(dep,init,tran_array,emit[1,1,],emit[2,1,])
print("post alpha")
beta <- Backward(dep,tran_array,emit[1,1,],emit[2,1,])


new_likelihood <- CalcLikelihood(alpha,pi_l)
likelihood_vec <- c(new_likelihood)
likelihood <- -Inf
like_diff <- new_likelihood - likelihood

# apply(alpha[[1]][,,1]+beta[[1]][,,1],1,logSumExp)

while(like_diff > 1e-3){
  
  start_time <- Sys.time()
  likelihood <- new_likelihood
  
  ##### MC Param  #####
  # init <- init_old
  # tran_array <- tran_array_old
  # emit <- emit_old
  # pi_l <- pi_l_old
  # alpha <- alpha_old
  # beta <- beta_old
  
  if (exists("init_old")){
    init_old_old <- init_old
    tran_array_old_old  <- tran_array_old
    emit_old_old  <- emit_old
    pi_l_old_old  <- pi_l_old
    alpha_old_old  <- alpha_old
    beta_old_old  <- beta_old
  }
  
  init_old <- init
  tran_array_old <- tran_array
  emit_old <- emit
  pi_l_old <- pi_l
  alpha_old <- alpha
  beta_old <- beta
  
  re_prob <- CalcProbRE(alpha,pi_l)
  pi_l <- colSums(re_prob)/dim(dep)[2]
  
  init <- CalcInit(alpha,beta,pi_l,log_sweights_vec)
  tran_array <- CalcTranC(alpha,beta,dep,emit,pi_l)
  emit <- CalcEmit(alpha,beta,dep,emit,pi_l)
  
  ##### Reorder #####
  #Reorder to avoid label switching
  #Cluster means go from small to large by activity
  if (RE_num > 1){
    # reord_inds <- order(tran_array[1,2,])
    # init <- init[reord_inds,]
    # tran_array <- tran_array[,,reord_inds]
    # pi_l <- pi_l[reord_inds]
  }
  
  ##### #####
  
  alpha <- Forward(dep,init,tran_array,emit[1,1,],emit[2,1,])
  beta <- Backward(dep,tran_array,emit[1,1,],emit[2,1,])
  
  new_likelihood <- CalcLikelihood(alpha,pi_l)
  like_diff <- new_likelihood - likelihood
  print(paste("RE num:",RE_num,"Like:",round(like_diff,6)))
  
  if (like_diff < 0){
    alpha <- Forward(dep,init_old,tran_array_old,emit_old[1,1,],emit_old[2,1,])
    beta <- Backward(dep,tran_array_old,emit_old[1,1,],emit_old[2,1,])
    new_likelihood <- CalcLikelihood(alpha,pi_l)
    like_diff <- new_likelihood - likelihood
    print(paste("RE num:",RE_num,"Like:",round(like_diff,6)))
    like_diff <- .0001
  }
  
  likelihood_vec <- c(likelihood_vec,new_likelihood)
  
  end_time <- Sys.time()
  time_vec <- c(time_vec,as.numeric(difftime(end_time, start_time, units = "secs")))
  
}

bic <- CalcBIC(new_likelihood,RE_num,dep)
est_params <- list(init,tran_array,emit,pi_l,re_prob)
assign_vec <- apply(re_prob,1,which.max)
decoded_dep <- Viterbi(dep,est_params,assign_vec)

to_save <- list(est_params,bic,decoded_dep)

save(to_save,file = paste0("atbcMHMM",RE_num,".rda"))

