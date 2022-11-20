rm(list=ls())

current_time <- Sys.time()
library(mvtnorm)
library(dplyr)
library(truncnorm)

set.seed(0202)

setwd("C:/Users/leeu9c/Downloads/R_simul")
dataset <- read.csv("actigraph_csv.csv")

subject1 <- subset(dataset, Index==1)
data <- subject1$axis1
plot(data,type="h")
# plot(data,type="l")

len <- length(data)/10
newdata <- numeric(len)

for(i in 1:len){
  newdata[i] <- mean(data[(10*i-9):(10*i)])
}

y <- round(newdata)
t <- 1:length(y)
L <- t[length(t)]+0.001
plot(y,type="h")

########################## functions ##########################
########################################################################################################
########################################################################################################
########################################################################################################

log_zip_lkh <- function(y_input, pi_input, lambda_input){
  y_0 <- y_input[y_input==0]
  y_not_0 <- y_input[y_input!=0]
  
  log_lkh_0 <- length(y_0) * log(pi_input + (1-pi_input)* exp(-lambda_input))
  log_lkh_not_0 <- length(y_not_0) * (log(1-pi_input)) + sum(dpois(y_not_0,lambda_input,log=T))
  
  log_lkh <- log_lkh_0 + log_lkh_not_0
  return(log_lkh)
}
# test
# log_zip_lkh(y[1:10],300,0.8)

kernel_ftn <- function(t_input, l_input, g_input){
  K <- exp(-(outer(t_input,t_input,FUN = "-")^2)/l_input)+g_input*diag(length(t_input))
  return(K)
} # covariance matrix = sigma^2 * kernel_ftn(t,l,g)


log_mvn_lkh <- function(y_input, t_input, mu_input, sigma2_input, l_input, g_input){
  d <- length(y_input)
  k_t_t_prime <- kernel_ftn(t_input, l_input, g_input)
  log_lkh <- dmvnorm(y_input, mean = rep(mu_input, d), sigma = sigma2_input * k_t_t_prime ,log=T )
  return(log_lkh)
}
# test
# log_mvn_lkh(t[1:10],y[1:10],400,200,1,1)


# case 1: a - r
# case 2: a - a
# case 3: r - a
# case 4: r - r

sum_log_mvn_lkh <- function(y_input, t_input, mu_input, sigma2_input, l_input, g_input, active_input, rest_input, n_active_input, case_input){
  
  sum <- 0
  for(I in 1:n_active_input){
    r <- case_when(
      case_input==1|case_input==2 ~ rest_input[I],
      case_input==3|case_input==4 ~ rest_input[I+1]
    )
    t_set <- t_input[ t_input > active_input[I] & t_input < r]
    y_set <- y_input[ t_input > active_input[I] & t_input < r ]
    
    if(length(y_set)!=0){
      sum <- sum + log_mvn_lkh(y_set, t_set, mu_input, sigma2_input, l_input, g_input)
    }
    
  }
  
  return(sum)
  
}
# sum_log_mvn_lkh(y,t,mu,sigma2,l,g,active,rest,n_active,case)


acc.birth.ar <- function(y_input, t_input, mvn_inputs, zip_inputs, k_input, start_input, end_input, L_input, p_birth_input, p_death_input, n_active_input, n_rest_input){
  
  # mvn_inputs : (mu, sigma^2, l, g)
  mu <- mvn_inputs[1]
  sigma2 <- mvn_inputs[2]
  l <- mvn_inputs[3]
  g <- mvn_inputs[4]
  
  # zip_inputs (pi, lambda)
  pi <- zip_inputs[1]
  lambda <- zip_inputs[2]
  
  log_lkh_ratio <- log_mvn_lkh(y_input, t_input, mu, sigma2, l, g) - log_zip_lkh(y_input, pi, lambda)
  log_step_ratio <- log(2*k_input+2)+log(2*k_input+1)-2*log(L_input)
  log_proposal_ratio <- log(p_death_input) + log(1/(n_active_input+1)) - log(p_birth_input) - log(1/n_rest_input) - log(2/(end_input-start_input)^2)
  
  log_acc <- log_lkh_ratio + log_step_ratio + log_proposal_ratio
  accept <- ifelse(log(runif(1))<log_acc,1,0)
  
  print(paste0("birth ar lkh is ", log_lkh_ratio))
  print(paste0("birth ar pp is ", log_proposal_ratio))
  
  return(accept)
}

acc.birth.ra <- function(y_input, t_input, mvn_inputs, zip_inputs, k_input, start_input, end_input, L_input, p_birth_input, p_death_input, n_active_input, n_rest_input){
  
  # mvn_inputs : (mu, sigma^2, l, g)
  mu <- mvn_inputs[1]
  sigma2 <- mvn_inputs[2]
  l <- mvn_inputs[3]
  g <- mvn_inputs[4]
  
  # zip_inputs (pi, lambda)
  pi <- zip_inputs[1]
  lambda <- zip_inputs[2]
  
  log_lkh_ratio <- log_zip_lkh(y_input, pi, lambda)-log_mvn_lkh(y_input, t_input, mu, sigma2, l, g) 
  log_step_ratio <- log(2*k_input+2)+log(2*k_input+1)-2*log(L_input)
  log_proposal_ratio <- log(p_death_input) + log(1/(n_rest_input+1)) - log(p_birth_input) - log(1/n_active_input) - log(2/(end_input-start_input)^2) 
  
  log_acc <- log_lkh_ratio + log_step_ratio + log_proposal_ratio
  accept <- ifelse(log(runif(1))<log_acc,1,0)
  
  print(paste0("birth ra lkh is ", log_lkh_ratio))
  print(paste0("birth ra pp is ", log_proposal_ratio))
  
  return(accept)
}


acc.death.ar <- function(y_input, t_input, mvn_inputs, zip_inputs, k_input, start_input, end_input, L_input, p_birth_input, p_death_input, n_active_input, n_rest_input){
  
  # mvn_inputs : (mu, sigma^2, l, g)
  mu <- mvn_inputs[1]
  sigma2 <- mvn_inputs[2]
  l <- mvn_inputs[3]
  g <- mvn_inputs[4]
  
  # zip_inputs (pi, lambda)
  pi <- zip_inputs[1]
  lambda <- zip_inputs[2]
  
  log_lkh_ratio <- log_zip_lkh(y_input, pi, lambda) - log_mvn_lkh(y_input, t_input, mu, sigma2, l, g) 
  log_step_ratio <- 2*log(L_input)-log(2*k_input+2)-log(2*k_input+1)
  log_proposal_ratio <- log(p_birth_input) + log(1/n_rest_input) + log(2/(end_input-start_input)^2) - log(1/n_active_input) - log(p_birth_input)
  
  log_acc <- log_lkh_ratio + log_step_ratio + log_proposal_ratio
  accept <- ifelse(log(runif(1))<log_acc,1,0) # 1:acc, 0: reject
  
  print(paste0("death ar lkh is ", log_lkh_ratio))
  print(paste0("death ar pp is ", log_proposal_ratio))
  return(accept)
}

acc.death.ra <- function(y_input, t_input, mvn_inputs, zip_inputs, k_input, start_input, end_input, L_input, p_birth_input, p_death_input, n_active_input, n_rest_input){
  
  # mvn_inputs : (mu, sigma^2, l, g)
  mu <- mvn_inputs[1]
  sigma2 <- mvn_inputs[2]
  l <- mvn_inputs[3]
  g <- mvn_inputs[4]
  
  # zip_inputs (pi, lambda)
  pi <- zip_inputs[1]
  lambda <- zip_inputs[2]
  
  log_lkh_ratio <- log_mvn_lkh(y_input, t_input, mu, sigma2, l, g) - log_zip_lkh(y_input, pi, lambda) 
  log_step_ratio <- 2*log(L_input)-log(2*k_input+2)-log(2*k_input+1)
  log_proposal_ratio <- log(p_birth_input) + log(1/n_active_input) + log(2/(end_input-start_input)^2) - log(1/n_rest_input) - log(p_birth_input)
  
  log_acc <- log_lkh_ratio + log_step_ratio + log_proposal_ratio
  accept <- ifelse(log(runif(1))<log_acc,1,0) # 1:acc, 0: reject
  print(paste0("death ra lkh is ", log_lkh_ratio))
  print(paste0("death ra pp is ", log_proposal_ratio))
  
  return(accept)
}

########################################################################################################
################################ initial values ##############################################
######################################################################### ###############################

initial_length <- 20 # initial_length/2 = # of active = # of rest for now
initial_cuts <- seq(1-0.01,length(t)+0.01,length=initial_length)

# initial case = 1 -> a start r end
case <- 1
active <- initial_cuts[(1:length(initial_cuts))%%2==1]
rest <- initial_cuts[(1:length(initial_cuts))%%2==0]
print(active)
print(rest)

abline(v=active,col="grey")
abline(v=rest,col="grey")

# active[1] < rest[1]  -> active,rest / else -> rest,active
position <- c(rbind(active,rest))
group_indicator <- as.numeric(cut(t,position))

#  active[1] < rest[1] -> odds : active / else -> odds:rest
y_rest <- y[group_indicator%%2==0]; y_active <- y[group_indicator%%2==1]
t_rest <- t[group_indicator%%2==0]; t_active <- t[group_indicator%%2==1]

# length of active & rest location
a_len <- length(active)
r_len <- length(rest)

# number of active set 
n_active <- a_len # or length(active)-1
n_rest <- r_len - 1

# MVN paramter (mu, sigma^2, l, g)
mu <- mean(y[y!=0])
sigma2 <- var(y[y!=0])
l <- 1
g <- 1
# ZIP paramter (pi, lambda)
pi <- 0.7
lambda <- mean(y[y!=0 & y<100])

# iteration
iter <- 50000
# store
mu_store <- sigma2_store <- l_store <- g_store <- pi_store <- lambda_store <- case_store <- numeric(iter)
active_list <- rest_list <- list()

mu_store[1] <- mu
sigma2_store[1] <- sigma2
l_store[1] <- l
g_store[1] <- g
pi_store[1] <- pi
lambda_store[1] <- lambda

# initial case
# case 1: a - r
# case 2: a - a
# case 3: r - a
# case 4: r - r
# if(max(active,rest)==tail(rest,1)){ # end with r
#   case <- ifelse(min(active,rest)==active[1],1,4)
# }else{
#   case <- ifelse(min(active,rest)==active[1],2,3)
# }
y_active_list <- t_active_list <- y_rest_list <- t_rest_list <- list()
acc_b_ar <- acc_b_ra <- acc_d_ar <- acc_d_ra <- numeric(iter)
########################################################################################################
################################ MCMC run ##############################################
########################################################################################################
for(i in 2:iter){
  
  #print("start")
  ############### MVN paramters update ##################
  # 1. mu
  mu_q <- rnorm(1, mu, 1.5)
  log_acc <- sum_log_mvn_lkh(y, t, mu_q, sigma2, l, g, active, rest, n_active, case) - sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active, case)
  
  mu <- ifelse(log(runif(1)) < log_acc, mu_q, mu)
  #print("mu")
  # 2. sigma2 >0 
  log.sigma2_q <- rnorm(1, log(sigma2), 0.1)
  sigma2_q <- exp(log.sigma2_q)
  
  log_acc <- sum_log_mvn_lkh(y, t, mu, sigma2_q, l, g, active, rest, n_active, case) + log.sigma2_q - sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active, case)-log(sigma2)
  sigma2 <- ifelse(log(runif(1)) < log_acc, sigma2_q, sigma2)
  #print("sigma")
  
  # 3. l > 0
  log.l_q <- rnorm(1, log(l), 0.1)
  l_q <- exp(log.l_q)
  
  log_acc <- sum_log_mvn_lkh(y, t, mu, sigma2, l_q, g, active, rest, n_active, case) + log.l_q - sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active, case)-log(l)
  l <- ifelse(log(runif(1)) < log_acc, l_q, l)
  #print("l")
  
  # 4. g
  log.g_q <- rnorm(1, log(g), 0.1)
  g_q <- exp(log.g_q)
  
  log_acc <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g_q, active, rest, n_active, case) + log.g_q - sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active, case)-log(g)
  g <- ifelse(log(runif(1)) < log_acc, g_q, g)
  
  
  ############### ZIP parameters update #####################
  # 1. pi (0 < pi <1)
  pi_q <- rtruncnorm(1,0,1,pi,0.5)
  log_acc <- log_zip_lkh(y_rest,pi_q,lambda) + log(dtruncnorm(pi,0,1,pi_q,1)) - log_zip_lkh(y_rest,pi,lambda) - log(dtruncnorm(pi_q,0,1,pi,1)) 
  pi <- ifelse(log(runif(1)) < log_acc, pi_q, pi)
  #print("pi")
  # 2. lambda
  log.lambda_q <- rnorm(1, log(lambda), 0.1)
  lambda_q <- exp(log.lambda_q)
  log_acc <- log_zip_lkh(y_rest,pi,lambda_q) + log.lambda_q - log_zip_lkh(y_rest,pi,lambda) - log(lambda)
  
  lambda <- ifelse(log(runif(1)) < log_acc, lambda_q, lambda)
  #print("lambda")
  
  ############### position change/ fix endpoints ########################
  
  if(length(active)!=1 & length(rest)!=1){
    # update position
    
    #print(position)
    k_len <- length(position)
    for(k in 2:(k_len-1)){
      
      current <- position[k]
      trunc_var <- 3
      proposal <- rtruncnorm(1, position[k-1], position[k+1], current, trunc_var)
      
      position_q <- position
      position_q[k] <- proposal
      
      
      c_indi <- as.numeric(cut(t,position))
      p_indi <- as.numeric(cut(t,position_q))
      
      subcase <- ifelse(case==1|case==2,1,2)
      
      switch(subcase,
             {y_rest <- y[c_indi%%2==0]
             y_rest_q <- y[p_indi%%2==0]},
             {y_rest <- y[c_indi%%2==1]
             y_rest_q <- y[p_indi%%2==1]})
      
      if(k%%2==0){
        
        if(subcase==1){ # rest position update
          
          rest_q <- rest
          rest_q[k/2] <- proposal
          
          new_lkh <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest_q, n_active,case) + log_zip_lkh(y_rest_q, pi, lambda)
          old_lkh <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active,case) + log_zip_lkh(y_rest, pi, lambda)
          
          log_acc <- new_lkh +  log(dtruncnorm(current, position[k-1], position[k+1], proposal, trunc_var)) - old_lkh - log(dtruncnorm(proposal, position[k-1], position[k+1], current, trunc_var))
          if(log(runif(1)) < log_acc){
            rest <- rest_q
          }
        }else{ # active position update
          
          active_q <- active
          active_q[k/2] <- proposal
          
          new_lkh <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active_q, rest, n_active,case) + log_zip_lkh(y_rest_q, pi, lambda)
          old_lkh <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active,case) + log_zip_lkh(y_rest, pi, lambda)
          
          log_acc <- new_lkh + log(dtruncnorm(current, position[k-1], position[k+1], proposal, 1)) - old_lkh - log(dtruncnorm(proposal, position[k-1], position[k+1], current, 1) )
          if(log(runif(1)) < log_acc){
            active <- active_q
          }
        }
      }else{
        
        if(subcase==1){ # active
          
          active_q <- active
          active_q[(k+1)/2] <- proposal
          
          new_lkh <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active_q, rest, n_active,case) + log_zip_lkh(y_rest_q, pi, lambda)
          old_lkh <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active,case) + log_zip_lkh(y_rest, pi, lambda)
          
          log_acc <- new_lkh + log(dtruncnorm(current, position[k-1], position[k+1], proposal, 1)) - old_lkh - log(dtruncnorm(proposal, position[k-1], position[k+1], current, 1) )
          if(log(runif(1)) < log_acc){
            active <- active_q
          }
        }else{ # rest
          
          rest_q <- rest
          rest_q[(k+1)/2] <- proposal
          
          new_lkh <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest_q, n_active,case) + log_zip_lkh(y_rest_q, pi, lambda)
          old_lkh <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active,case) + log_zip_lkh(y_rest, pi, lambda)
          
          log_acc <- new_lkh +  log(dtruncnorm(current, position[k-1], position[k+1], proposal, trunc_var)) - old_lkh - log(dtruncnorm(proposal, position[k-1], position[k+1], current, trunc_var))
          if(log(runif(1)) < log_acc){
            rest <- rest_q
          }
        }
        
      }
      
      position <- switch (case,
                          c(rbind(active,rest)),
                          c(rbind(active[1:(a_len-1)],rest),active[a_len]),
                          c(rbind(rest,active)),
                          c(rbind(rest[1:(r_len-1)],active),rest[r_len]),
      )
      #print(position)
      
    }# position update
    #print("position")
  }
  
  
  
  ############### birth of death ########################
  
  if(length(active)==1 & length(rest)==1){ # only one period -> do birth
    
    if(active<rest){ # r,a birth
      
      j <- sample(1:n_active, 1) 
      start <- active[j]
      end <- case_when(
        case==1|case==2 ~ rest[j],
        case==3|case==4 ~ rest[j+1]
      ) 
      
      # temp <- runif(4,start,end) %>% sort
      # round_t <- round(temp) 
      # 
      # r_star <- temp[which.min(y[round_t])]
      # a_star <- runif(1,r_star,end)

      
      temp <- runif(2, start, end)
      r_star <- min(temp); a_star <- max(temp) # proposed (r*, a*)
      
      t_star <- t[t>r_star & t<a_star]
      y_star <- y[t>r_star & t<a_star]
      
      k <- length(c(active,rest))
      mvn_parameters <- c(mu, sigma2, l, g)
      zip_parameters <- c(pi,lambda)
      if(length(y_star)==0){
        birth_indicator <-0 # no birth
      }else{
        birth_indicator <- acc.birth.ra(y_star, t_star, mvn_parameters, zip_parameters, k, start, end, L, p_birth, p_death, n_active, n_rest)
      }
      # print(birth_indicator) 
      if(birth_indicator==1){ # index
        
        active <- c(active,a_star)
        rest <- c(r_star,rest)
        
      }# index end

    }else{ # a,r birth
      j <- sample(1:n_rest, 1) 
      start <- rest[j]
      end <- case_when(
        case==1|case==2 ~ active[j+1],
        case==3|case==4 ~ active[j]
      ) # selected rest interval = (start, end) = r[j],a[j] or r[j],a[j+1]
      
      temp <- runif(2, start, end)
      a_star <- min(temp); r_star <- max(temp) # proposed (a*, r*)
      
      t_star <- t[t>a_star & t<r_star]
      y_star <- y[t>a_star & t<r_star]
      
      k <- length(c(active,rest))
      mvn_parameters <- c(mu, sigma2, l, g)
      zip_parameters <- c(pi,lambda)
      if(length(y_star)==0){
        birth_indicator <-0 # no birth
      }else{
        birth_indicator <- acc.birth.ar(y_star, t_star, mvn_parameters, zip_parameters, k, start, end, L, p_birth, p_death, n_active, n_rest)
      }
      # print(birth_indicator) # 1: accept, 2:reject
      
      if(birth_indicator==1){ # if birth, new index
        
        active <- c(a_star,active)
        rest <- c(rest,r_star)
        
      }# index end
      
    }# (a,r) birth end
    
  }else{ # birth(a,r)/ birth(r,a)/ death(a,r)/ death(r,a)
    
    p_birth <- 0.5;p_death <- 1-p_birth # prob
    
    #(a,r) birth
    j <- sample(1:n_rest, 1) 
    start <- rest[j]
    end <- case_when(
      case==1|case==2 ~ active[j+1],
      case==3|case==4 ~ active[j]
    ) # selected rest interval = (start, end) = r[j],a[j] or r[j],a[j+1]
    
    temp <- runif(2, start, end)
    a_star <- min(temp); r_star <- max(temp) # proposed (a*, r*)
    
    t_star <- t[t>a_star & t<r_star]
    y_star <- y[t>a_star & t<r_star]
    
    k <- length(c(active,rest))
    mvn_parameters <- c(mu, sigma2, l, g)
    zip_parameters <- c(pi,lambda)
    if(length(y_star)==0){
      birth_indicator <-0 # no birth
    }else{
      birth_indicator <- acc.birth.ar(y_star, t_star, mvn_parameters, zip_parameters, k, start, end, L, p_birth, p_death, n_active, n_rest)
    }
    # print(birth_indicator) # 1: accept, 2:reject
    
    if(birth_indicator==1){ # if birth, new index
      
      acc_b_ar[i] <- 1
      
      if(j==1){
        
        rest <- c(rest[1:j],r_star,rest[(j+1):r_len])
        active <- switch(case,
                         c(active[1:j],a_star,active[(j+1):a_len]),
                         c(active[1:j],a_star,active[(j+1):a_len]),
                         c(a_star,active),
                         c(a_star,active)
        )
      }else if(j==n_rest){
        
        rest <- switch(case,
                       c(rest[1:j],r_star,rest[(j+1):r_len]),
                       c(rest,r_star),
                       c(rest,r_star),
                       c(rest[1:j],r_star,rest[(j+1):r_len])
        )
        active <- switch(case,
                         c(active[1:j],a_star,active[(j+1):a_len]),
                         c(active[1:j],a_star,active[(j+1):a_len]),
                         c(active[1:(j-1)],a_star,active[j:a_len]),
                         c(active[1:(j-1)],a_star,active[j:a_len])
        )
      }else{
        
        rest <- c(rest[1:j],r_star,rest[(j+1):r_len])
        active <- switch(case,
                         c(active[1:j],a_star,active[(j+1):a_len]),
                         c(active[1:j],a_star,active[(j+1):a_len]),
                         c(active[1:(j-1)],a_star,active[j:a_len]),
                         c(active[1:(j-1)],a_star,active[j:a_len])
        )
        
      }
      
    }# (a,r) birth end
    
    # new case
    if(max(active,rest)==tail(rest,1)){ # end with r
      case <- ifelse(min(active,rest)==active[1],1,4)
    }else{
      case <- ifelse(min(active,rest)==active[1],2,3)
    }
    #print("here1")
    
    a_len <- length(active); r_len <- length(rest)
    
    position <- switch (case,
                        c(rbind(active,rest)),
                        c(rbind(active[1:(a_len-1)],rest),active[a_len]),
                        c(rbind(rest,active)),
                        c(rbind(rest[1:(r_len-1)],active),rest[r_len]),
    )
    group_indicator <- as.numeric(cut(t,position))
    
    if(case==1|case==2){
      y_rest <- y[group_indicator%%2==0]; y_active <- y[group_indicator%%2==1]
      t_rest <- t[group_indicator%%2==0]; t_active <- t[group_indicator%%2==1]
      
    }else{
      y_rest <- y[group_indicator%%2==1]; y_active <- y[group_indicator%%2==0]
      t_rest <- t[group_indicator%%2==1]; t_active <- t[group_indicator%%2==0]
    }
    #print("here2")
    # n_active, n_rest, a_len, r_len
    # case 1: a - r
    # case 2: a - a
    # case 3: r - a
    # case 4: r - r
    
    
    if(case==2|case==3){
      n_active <- a_len -1
      n_rest <- r_len
    }else{
      n_active <- a_len
      n_rest <- r_len-1
    }
    
    #(r,a) birth
    j <- sample(1:n_active, 1) 
    start <- active[j]
    end <- case_when(
      case==1|case==2 ~ rest[j],
      case==3|case==4 ~ rest[j+1]
    ) 
    
    # temp <- runif(4,start,end) %>% sort
    # round_t <- round(temp) 
    # 
    # r_star <- temp[which.min(y[round_t])]
    # a_star <- runif(1,r_star,end)
    
    temp <- runif(2, start, end)
    r_star <- min(temp); a_star <- max(temp) # proposed (r*, a*)
    # 
    t_star <- t[t>r_star & t<a_star]
    y_star <- y[t>r_star & t<a_star]
    
    k <- length(c(active,rest))
    mvn_parameters <- c(mu, sigma2, l, g)
    zip_parameters <- c(pi,lambda)
    if(length(y_star)==0){
      birth_indicator <-0 # no birth
    }else{
      birth_indicator <- acc.birth.ra(y_star, t_star, mvn_parameters, zip_parameters, k, start, end, L, p_birth, p_death, n_active, n_rest)
    }
    # print(birth_indicator) 
    if(birth_indicator==1){ # index
      
      acc_b_ra[i] <- 1
      # case 1: a - r
      # case 2: a - a
      # case 3: r - a
      # case 4: r - r
      
      if(j==1){
        
        rest <- switch(case,
                       c(r_star,rest),
                       c(r_star,rest),
                       c(rest[1:j],r_star,rest[(j+1):r_len]),
                       c(rest[1:j],r_star,rest[(j+1):r_len]))
        
        active <- c(active[1:j],a_star,active[(j+1):a_len])
        
      }else if(j==n_active){
        rest <- switch(case,
                       c(rest[1:(j-1)],r_star,rest[j:r_len]),
                       c(rest[1:(j-1)],r_star,rest[j:r_len]),
                       c(rest[1:j],r_star,rest[(j+1):r_len]),
                       c(rest[1:j],r_star,rest[(j+1):r_len]))
        
        active <- switch(case,
                         c(active,a_star),
                         c(active[1:j],a_star,active[(j+1):a_len]),
                         c(active[1:j],a_star,active[(j+1):a_len]),
                         c(active,a_star))
        
      }else{
        rest <- switch(case,
                       c(rest[1:(j-1)],r_star,rest[j:r_len]),
                       c(rest[1:(j-1)],r_star,rest[j:r_len]),
                       c(rest[1:j],r_star,rest[(j+1):r_len]),
                       c(rest[1:j],r_star,rest[(j+1):r_len]))
        
        active <- c(active[1:j],a_star,active[(j+1):a_len])
      }
    }# r.a index end
    # new case
    if(max(active,rest)==tail(rest,1)){ # end with r
      case <- ifelse(min(active,rest)==active[1],1,4)
    }else{
      case <- ifelse(min(active,rest)==active[1],2,3)
    }
    #print("here1")
    
    a_len <- length(active); r_len <- length(rest)
    
    position <- switch (case,
                        c(rbind(active,rest)),
                        c(rbind(active[1:(a_len-1)],rest),active[a_len]),
                        c(rbind(rest,active)),
                        c(rbind(rest[1:(r_len-1)],active),rest[r_len]),
    )
    group_indicator <- as.numeric(cut(t,position))
    
    if(case==1|case==2){
      y_rest <- y[group_indicator%%2==0]; y_active <- y[group_indicator%%2==1]
      t_rest <- t[group_indicator%%2==0]; t_active <- t[group_indicator%%2==1]
      
    }else{
      y_rest <- y[group_indicator%%2==1]; y_active <- y[group_indicator%%2==0]
      t_rest <- t[group_indicator%%2==1]; t_active <- t[group_indicator%%2==0]
    }
    #print("here2")
    # n_active, n_rest, a_len, r_len
    # case 1: a - r
    # case 2: a - a
    # case 3: r - a
    # case 4: r - r
    
    
    if(case==2|case==3){
      n_active <- a_len -1
      n_rest <- r_len
    }else{
      n_active <- a_len
      n_rest <- r_len-1
    }
    
    #(a,r) death
    if(length(rest)!=1|length(active)!=1){
      j <- sample(1:n_active, 1)
      start <- active[j]
      end <- case_when(
        case==1|case==2 ~ rest[j],
        case==3|case==4 ~ rest[j+1]
      ) # selected active interval = (start, end) = a[j],r[j] or a[j],r[j+1]
      
      a_star <- start; r_star <- end
      
      t_star <- t[t>a_star & t<r_star]
      y_star <- y[t>a_star & t<r_star]
      
      k <- length(c(active,rest))
      mvn_parameters <- c(mu, sigma2, l, g)
      zip_parameters <- c(pi,lambda)
      
      j_prime <- sample(1:n_rest, 1)
      start <- rest[j_prime]
      end <- case_when(
        case==1|case==2 ~ active[j_prime+1],
        case==3|case==4 ~ active[j_prime]
      ) 
      
      if(length(y_star)==0){
        death_indicator <- 0
      }else{
        death_indicator <- acc.death.ar(y_star, t_star, mvn_parameters, zip_parameters, k, start, end, L, p_birth, p_death, n_active, n_rest)
      }
      
      if(death_indicator==1){ # accept -> indexing
        
        acc_d_ar[i] <- 1
        active <- active[-j]
        if(case==1|case==2){
          rest <- rest[-j]
        }else{
          rest <- rest[-(j+1)]
        }
        # print("here")
        
        if(j==1){
          if(case==1|case==2){
            rest <- c(0.99,rest)
          }
        }else if(j==n_active){
          if(case==1|case==4){
            active <- c(active,L)
          }
        }
        
      } # index end
    }# (a,r) death end
    
    # new case
    if(max(active,rest)==tail(rest,1)){ # end with r
      case <- ifelse(min(active,rest)==active[1],1,4)
    }else{
      case <- ifelse(min(active,rest)==active[1],2,3)
    }

    a_len <- length(active); r_len <- length(rest)
    
    position <- switch (case,
                        c(rbind(active,rest)),
                        c(rbind(active[1:(a_len-1)],rest),active[a_len]),
                        c(rbind(rest,active)),
                        c(rbind(rest[1:(r_len-1)],active),rest[r_len]),
    )
    group_indicator <- as.numeric(cut(t,position))
    
    if(case==1|case==2){
      y_rest <- y[group_indicator%%2==0]; y_active <- y[group_indicator%%2==1]
      t_rest <- t[group_indicator%%2==0]; t_active <- t[group_indicator%%2==1]
      
    }else{
      y_rest <- y[group_indicator%%2==1]; y_active <- y[group_indicator%%2==0]
      t_rest <- t[group_indicator%%2==1]; t_active <- t[group_indicator%%2==0]
    }

    
    
    if(case==2|case==3){
      n_active <- a_len -1
      n_rest <- r_len
    }else{
      n_active <- a_len
      n_rest <- r_len-1
    }
    # (r,a) death
    if(length(rest)!=1|length(active)!=1){
      
      j <- sample(1:n_rest, 1)
      start <- rest[j]
      end <- case_when(
        case==1|case==2 ~ active[j+1],
        case==3|case==4 ~ active[j]
      ) 
      
      r_star <- start; a_star <- end
      
      t_star <- t[t>r_star & t<a_star]
      y_star <- y[t>r_star & t<a_star]
      
      k <- length(c(active,rest))
      mvn_parameters <- c(mu, sigma2, l, g)
      zip_parameters <- c(pi,lambda)
      
      j_prime <- sample(1:n_rest, 1)
      start <- rest[j_prime]
      end <- case_when(
        case==1|case==2 ~ active[j_prime+1],
        case==3|case==4 ~ active[j_prime]
      ) 
      
      if(length(y_star)==0){
        death_indicator <- 0
      }else{
        death_indicator <- acc.death.ra(y_star, t_star, mvn_parameters, zip_parameters, k, start, end, L, p_birth, p_death, n_active, n_rest)
      }
      
      if(death_indicator==1){ # index
        acc_d_ra[i] <- 1
        rest <- rest[-j]
        if(case==1|case==2){
          active <- active[-(j+1)]
        }else{
          active <- active[-j]
        }
        
        if(j==1){
          if(case==3|case==4){
            active <- c(0.99, active)
          }
        }else if(j==n_rest){
          if(case==2|case==3){
            rest <- c(rest,L)
          }
        }
        
        
      }# index end
    }# (r,a) death end
    # new case
    if(max(active,rest)==tail(rest,1)){ # end with r
      case <- ifelse(min(active,rest)==active[1],1,4)
    }else{
      case <- ifelse(min(active,rest)==active[1],2,3)
    }
    #print("here1")
    
    a_len <- length(active); r_len <- length(rest)
    
    position <- switch (case,
                        c(rbind(active,rest)),
                        c(rbind(active[1:(a_len-1)],rest),active[a_len]),
                        c(rbind(rest,active)),
                        c(rbind(rest[1:(r_len-1)],active),rest[r_len]),
    )
    group_indicator <- as.numeric(cut(t,position))
    
    if(case==1|case==2){
      y_rest <- y[group_indicator%%2==0]; y_active <- y[group_indicator%%2==1]
      t_rest <- t[group_indicator%%2==0]; t_active <- t[group_indicator%%2==1]
      
    }else{
      y_rest <- y[group_indicator%%2==1]; y_active <- y[group_indicator%%2==0]
      t_rest <- t[group_indicator%%2==1]; t_active <- t[group_indicator%%2==0]
    }
    #print("here2")
    # n_active, n_rest, a_len, r_len
    # case 1: a - r
    # case 2: a - a
    # case 3: r - a
    # case 4: r - r
    
    
    if(case==2|case==3){
      n_active <- a_len -1
      n_rest <- r_len
    }else{
      n_active <- a_len
      n_rest <- r_len-1
    }
    


      

    
    
  }# birth or death end
  
  
  
  
  
  ################## new case, update y_rest, y_active, t,rest, t_active ##################
  
  # new case
  if(max(active,rest)==tail(rest,1)){ # end with r
    case <- ifelse(min(active,rest)==active[1],1,4)
  }else{
    case <- ifelse(min(active,rest)==active[1],2,3)
  }
  #print("here1")
  
  a_len <- length(active); r_len <- length(rest)
  
  position <- switch (case,
                      c(rbind(active,rest)),
                      c(rbind(active[1:(a_len-1)],rest),active[a_len]),
                      c(rbind(rest,active)),
                      c(rbind(rest[1:(r_len-1)],active),rest[r_len]),
  )
  group_indicator <- as.numeric(cut(t,position))
  
  if(case==1|case==2){
    y_rest <- y[group_indicator%%2==0]; y_active <- y[group_indicator%%2==1]
    t_rest <- t[group_indicator%%2==0]; t_active <- t[group_indicator%%2==1]
    
  }else{
    y_rest <- y[group_indicator%%2==1]; y_active <- y[group_indicator%%2==0]
    t_rest <- t[group_indicator%%2==1]; t_active <- t[group_indicator%%2==0]
  }
  #print("here2")
  # n_active, n_rest, a_len, r_len
  # case 1: a - r
  # case 2: a - a
  # case 3: r - a
  # case 4: r - r
  
  
  if(case==2|case==3){
    n_active <- a_len -1
    n_rest <- r_len
  }else{
    n_active <- a_len
    n_rest <- r_len-1
  }
  
  #print("here3")
  ############### save #####################
  active_list[[i]] <- active
  rest_list[[i]] <- rest
  mu_store[i] <- mu
  sigma2_store[i] <- sigma2
  l_store[i] <- l
  g_store[i] <- g
  pi_store[i] <- pi
  lambda_store[i] <- lambda
  
  y_active_list[[i]] <- y_active
  y_rest_list[[i]] <- y_rest
  t_active_list[[i]] <- t_active
  t_rest_list[[i]] <- t_rest
  
  print(active)
  print(rest)
  
  print(i)
  
  
  plot(t_active,y_active,type="h",col='red', xlim=c(0,length(t)),ylim=c(0,max(y_active)),xlab="",ylab="")
  par(new=T)
  plot(t_rest,y_rest,type="h",col='blue', xlim=c(0,length(t)),lwd=3, ylim=c(0,max(y_active)),xlab="",ylab="")
  
  # 
  if(i%%500==0){
    setwd("C:/Users/leeu9c/Downloads/R_simul/result")
    png(file=paste0("Fig",i,".png"),width=500,height=400,pointsize=20)
    plot(t_active,y_active,type="h",col='red', xlim=c(0,length(t)),ylim=c(0,max(y_active)),xlab="",ylab="")
    par(new=T)
    plot(t_rest,y_rest,type="h",col='blue', xlim=c(0,length(t)),lwd=3, ylim=c(0,max(y_active)),xlab="",ylab="")
    dev.off()

    #save to rda file
    save(t_active_list, t_rest_list, file = paste0("locations_",i,".rda"))
  }

  
  
}## MCMC

sum(acc_b_ar)/i
sum(acc_b_ra)/i
sum(acc_d_ar)/i
sum(acc_d_ra)/i



## t
t_post <- numeric(length(t))
for(i in 1:length(t)){
  if(sum(unlist(t_active_list)==i)>=sum(unlist(t_rest_list)==i))t_post[i] <- 1
  
}
table(t_post)
t_active_p <- t[t_post==1]
t_rest_p <- t[t_post==0]
y_active_p <- y[t_post==1]
y_rest_p <- y[t_post==0]
par(mfrow=c(1,1))
plot(t_active_p,y_active_p,type="h",col='red', xlim=c(0,length(t)),ylim=c(0,max(y_active_p)),xlab="",ylab="")
par(new=T)
plot(t_rest_p,y_rest_p,type="h",col='blue', lwd=2,xlim=c(0,length(t)),ylim=c(0,max(y_active_p)),xlab="",ylab="")
abline(v=active,col="green",lwd=2)
abline(v=rest,col="blue")



end_time <-Sys.time()
end_time-current_time

par(mfrow=c(2,3))
plot(mu_store[1:23950],type="l")
plot(sigma2_store[1:23950],type="l")
plot(l_store[1:23950],type="l")
plot(g_store[1:23950],type="l")
plot(pi_store[1:23950],type="l")
plot(lambda_store[1:23950],type="l")


par(mfrow=c(3,1))
plot(l_store[1:4000],type="l")
plot(g_store[1:4000],type="l")
plot(sigma2_store[1:4000],type="l")

