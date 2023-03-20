## 03/15/23 : DP priors version
rm(list=ls())

current_time <- Sys.time()
library(mvtnorm)
library(dplyr)
library(truncnorm)

set.seed(0202)
dataset <- read.csv("/Users/wonny/Downloads/actigraph_csv.csv")
source('/Users/wonny/Documents/Reversible_Jump_MCMC/rjmcmc_ftn_0315.R')

subject1 <- subset(dataset, Index==2)
data <- (subject1$axis1 + subject1$axis2 + subject1$axis3)/3

len <- round(length(data)/10)
newdata <- numeric(len)

for(i in 1:len){
  newdata[i] <- mean(data[(10*i-9):(10*i)])
}

y <- round(newdata)
t <- 1:length(y)
L <- t[length(t)]+0.5


################################ initial values ##############################################

Init <- initial_position(y, quantile(y)[3]+0.1)
active <- Init$active
rest <- Init$rest

# active[1] < rest[1]  -> active,rest / else -> rest,active
position <- c(rbind(active,rest))
group_indicator <- as.numeric(cut(t,position))

#  active[1] < rest[1] -> odds : active / else :rest
y_rest <- y[group_indicator%%2==0]; y_active <- y[group_indicator%%2==1]
t_rest <- t[group_indicator%%2==0]; t_active <- t[group_indicator%%2==1]


# length of active & rest vectors
a_len <- length(active)
r_len <- length(rest)

# number of active & rest intervals 
n_active <- a_len 
n_rest <- r_len - 1

# MVN parameters (mu, sigma^2, l, g)
mu <- mean(y[y!=0])
sigma2 <- var(y[y!=0])
l <- 1
g <- 10^-4 # fix 

Lambda <- seq(1,251,10)
N <- length(Lambda) ## maximum possible number of components
indi <- numeric(length(y_rest)) # indicator
Pi <- rep(1/N, N)  # initial values for Pi


for(j in 1:length(y_rest)){
  
  temp_prob <- numeric(N)
 for(k in 1:N){
   temp_prob[k] <- log(Pi[k]) + dpois(y_rest[j], Lambda[k], log=T)
 }

  indi[j] <- which.max(temp_prob) # initial indicators
}

# parameter of the base distribution G0

alpha0 <- quantile(y)[3]/2 # shape
beta0 <- 2 # scale

dp_alpha <- 1.5

iter <- 10^4

# store
mu_store <- sigma2_store <- l_store <- n_active_store <- n_rest_store <- n_components_store <-  numeric(iter)
active_list <- rest_list <- list()

lambda_store <- pi_store <- matrix(0, iter, N)

mu_store[1] <- mu
sigma2_store[1] <- sigma2
l_store[1] <- l
lambda_store[1,] <- Lambda

y_active_list <- t_active_list <- y_rest_list <- t_rest_list <- indi_list <- list()
acc_b_ar <- acc_b_ra <- acc_d_ar <- acc_d_ra <- death_j <- numeric(iter)


### MCMC ####
for(i in 2:iter){
  
  u <- runif(1)
  
  p_posit <- 0.3
  p_birth <- 0.3
  p_death <- 1 - p_posit - p_birth
  
  # parameter updates
  
  ############### MVN parameters update ##################
  # 1. mu
  mu_q <- rnorm(1, mu, 1.5)
  log_acc <- sum_log_mvn_lkh(y, t, mu_q, sigma2, l, g, active, rest, n_active) - sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active)
  
  mu <- ifelse(log(runif(1)) < log_acc, mu_q, mu)
  # 2. sigma2 >0 
  log.sigma2_q <- rnorm(1, log(sigma2), 0.1)
  sigma2_q <- exp(log.sigma2_q)
  
  log_acc <- sum_log_mvn_lkh(y, t, mu, sigma2_q, l, g, active, rest, n_active) + log.sigma2_q - sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active)-log(sigma2)
  sigma2 <- ifelse(log(runif(1)) < log_acc, sigma2_q, sigma2)
  
  # 3. l > 0
  log.l_q <- rnorm(1, log(l), 0.1)
  l_q <- exp(log.l_q)
  
  log_acc <- sum_log_mvn_lkh(y, t, mu, sigma2, l_q, g, active, rest, n_active) + log.l_q - sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active)-log(l)
  l <- ifelse(log(runif(1)) < log_acc, l_q, l)
  
  # 4. g
  g <- 10^-4
  
  mvn_parameters <- c(mu, sigma2, l, g)
  
  ###################################################################################
  # lambda update
  K <- unique(indi)
  K <- K[order(K)]
  
  for(j in K){
    
    Lambda[j] <- rgamma(1, shape = alpha0 + sum(y_rest[indi==j]), scale = beta0/(length(y_rest[indi==j])*beta0 + 1))
    
  }

  for(j in setdiff(c(1:N),K)){
    Lambda[j] <- rgamma(1, shape = alpha0, scale = beta0)
  }
  print(Lambda)
  #####################################################################################
  # indicator update
  
  for(j in 1:length(y_rest)){
    
    temp_prob <- numeric(N)
    for(k in 1:N){
      temp_prob[k] <- log(Pi[k]) + dpois(y_rest[j], Lambda[k], log=T)
    }
    
    indi[j] <- which.max(temp_prob) 
  }
  
  
  #####################################################################################
  # pi update
  M <- U <- numeric(N)
  for(j in 1:N){
    M[j] <- sum(indi==j)
  }
  
  for(j in 1:(N-1)){
    U[j] <- rbeta(1, 1+M[j], dp_alpha + sum(M[(j+1):N]))
  }
  U[N] <- 1
  
  Pi[1] <- U[1]
  for(j in 2:(N-1)){
    Pi[j] <- U[j] * prod(1-U[1:(j-1)])
  }
  Pi[N] <- 1- sum(Pi[1:(N-1)])
  
  #####################################################################################
  ## position change / birth / death 
  
  if(u < p_posit){ # position updates 
    print("posit")
    ############### position change/ fix endpoints ########################
    
    if(length(active)!=1 & length(rest)!=1){
      # update position
      
      k_len <- length(position)
      for(k in 2:(k_len-1)){
        
        current <- position[k]
        candidates <- seq(position[k-1]+1, position[k+1]-1, by=1)
        if(length(candidates)==1){
          proposal <- candidates
        }else{
          proposal <- sample(candidates, 1)
          
        }
        
        if(current!=proposal){ # if == -> no need to update
          
          if(k%%2==0){
            
            rest_q <- rest
            rest_q[k/2] <- proposal
            
            if(proposal < current){
              
              # A -> R
              
              y_star <-  y[t> proposal & t<current]
              t_star <- t[t>proposal & t<current]
              
              if(length(y_star)!=0){
                log_num <- log_den_y(y_star, Pi, Lambda)
                log_den <- log_mvn_lkh(y_star, t_star, mu, sigma2, l, g)
                
                if(log(runif(1)) < (log_num - log_den)){
                  rest <- rest_q
                }
              }
                
            }else{
              
              # R -> A
              y_star <-  y[t> current & t< proposal]
              t_star <- t[t> current & t< proposal]
              
              if(length(y_star)!=0){
                log_num <- log_mvn_lkh(y_star, t_star, mu, sigma2, l, g)
                log_den <- log_den_y(y_star, Pi, Lambda)
                
                if(log(runif(1)) < ( log_num - log_den )){
                  rest <- rest_q
                }
              }

              
            } # rest case end
            
          }else{
            
            active_q <- active
            active_q[(k+1)/2] <- proposal
            
            if(proposal < current){
              
              # R -> A
              
              y_star <-  y[t> proposal & t<current]
              t_star <- t[t>proposal & t<current]
              
              if(length(y_star)!=0){
                log_num <- log_mvn_lkh(y_star, t_star, mu, sigma2, l, g)
                log_den <- log_den_y(y_star, Pi, Lambda)
                
                if(log(runif(1)) < (log_num - log_den)){
                  active <- active_q
                }
              }
          
            }else{
              
              # A -> R
              y_star <-  y[t> current & t< proposal]
              t_star <- t[t> current & t< proposal]
              
              if(length(y_star)!=0){
                log_num <- log_den_y(y_star, Pi, Lambda)
                log_den <- log_mvn_lkh(y_star, t_star, mu, sigma2, l, g)
                
                if(log(runif(1)) < (log_num - log_den)){
                  active <- active_q
                }
              }

            }
            
          } # active case end
          
          position <- c(rbind(active,rest))
          
        }
        
      }
    }# position update
    
  }else if(u < p_posit + p_birth){ # birth updates
    print("birth")
    ##################### (a,r) birth #######################################
    
    if(n_rest >= 1){
      temp <- acc.birth.ar(y, t, mvn_parameters, Pi, Lambda, p_birth, p_death, active, rest, n_active, n_rest)
      
      active <- temp$active
      rest <- temp$rest
      n_active <- temp$n_active
      n_rest <- temp$n_rest
      acc_b_ar[i] <- temp$accept
      
    }
    
    
    ##################### (r,a) birth #######################################
    
    if(n_active >= 1){
      temp <- acc.birth.ra(y, t, mvn_parameters, Pi, Lambda, p_birth, p_death, active, rest, n_active, n_rest)
      
      active <- temp$active
      rest <- temp$rest
      n_active <- temp$n_active
      n_rest <- temp$n_rest
      acc_b_ra[i] <- temp$accept
      
      #print("r,a birth")
      #print(lambda1)
      #print(lambda2)
    }
    
    
  }else{ # death updates
    print("death")
    ##################### (a,r) death #######################################
    #print("death")
    if(n_active>2){
      temp <- acc.death.ar(y, t, mvn_parameters, Pi, Lambda, p_birth, p_death, active, rest, n_active, n_rest)
      
      active <- temp$active
      rest <- temp$rest
      n_active <- temp$n_active
      n_rest <- temp$n_rest
      acc_d_ar[i] <- temp$accept
      death_j[i] <- temp$j
    }
    
    ##################### (r,a) death #######################################
    
    if(n_rest>1){
      temp <- acc.death.ra(y, t, mvn_parameters, Pi, Lambda, p_birth, p_death, active, rest, n_active, n_rest)
      
      active <- temp$active
      rest <- temp$rest
      n_active <- temp$n_active
      n_rest <- temp$n_rest
      acc_d_ra[i] <- temp$accept
      
      
    }
    
  } # death update end
  
  position <- c(rbind(active,rest))
  group_indicator <- as.numeric(cut(t,position))
  
  y_rest <- y[group_indicator%%2==0]; y_active <- y[group_indicator%%2==1]
  t_rest <- t[group_indicator%%2==0]; t_active <- t[group_indicator%%2==1]
  
  indi <- numeric(length(y_rest))
  for(j in 1:length(y_rest)){
    
    temp_prob <- numeric(N)
    for(k in 1:N){
      temp_prob[k] <- log(Pi[k]) + dpois(y_rest[j], Lambda[k], log=T)
    }
    
    indi[j] <- which.max(temp_prob) 
  }
  
  ############### save #####################
  active_list[[i]] <- active
  rest_list[[i]] <- rest
  
  mu_store[i] <- mu
  sigma2_store[i] <- sigma2
  l_store[i] <- l
  lambda_store[i,] <- Lambda
  pi_store[i,] <- Pi
  
  n_active_store[i] <- n_active
  n_rest_store[i] <- n_rest
  n_components_store[i] <- length(unique(indi))
  
  y_active_list[[i]] <- y_active
  y_rest_list[[i]] <- y_rest
  t_active_list[[i]] <- t_active
  t_rest_list[[i]] <- t_rest
  
  #print(Lambda)
  print(i)

  
  plot(t_active,y_active,type="h",col='red', xlim=c(0,length(t)),ylim=c(0,max(y_active)),xlab="",ylab="")
  par(new=T)
  plot(t_rest,y_rest,type="h",col='blue', xlim=c(0,length(t)), ylim=c(0,max(y_active)),xlab="",ylab="")
  
  
  
  
  
  
} # MCMC end 


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
plot(t_rest_p,y_rest_p,type="h",col='blue',xlim=c(0,length(t)),ylim=c(0,max(y_active_p)),xlab="",ylab="")


# plot(mu_store,type = "l")
# plot(sigma2_store,type = "l")
# plot(l_store,type = "l")
# 
# ## correlation matrix
# 
# # > kernel_ftn(t_rest, mean(l_store[2000:iter]), 10^-4)[1:4,1:4]
# # [,1]        [,2]        [,3]         [,4]
# # [1,] 1.0001000000 0.290483686 0.084380772 0.0006008008
# # [2,] 0.2904836861 1.000100000 0.290483686 0.0020682772
# # [3,] 0.0843807719 0.290483686 1.000100000 0.0071201147
# # [4,] 0.0006008008 0.002068277 0.007120115 1.0001000000
# 
# 
# 
# plot(t_active_p,y_active_p,type="h",col='tomato', xlim=c(0,length(t)),ylim=c(0,max(y_active_p)),xlab="",ylab="")
# par(new=T)
# plot(t_rest_p,y_rest_p,type="h",col='royalblue',xlim=c(0,length(t)),ylim=c(0,max(y_active_p)),xlab="",ylab="")
# 
# 
# mean(lambda1_store[100:iter]);mean(lambda2_store[100:iter]);mean(lambda3_store[100:iter])
# # [1] 0.4092578
# # [1] 6.938948
# # [1] 30.42742
# par(mfrow=c(1,1))
# plot(lambda1_store[1:1000],type="l")
# 
# ## posterior active, rest
# temp_list <- list()
# 
# for(i in 2:(length(t_post)-1)){
#   
#   if(abs(t_post[i+1]-t_post[i])==1){
#     temp_list[[i]] <- i+1
#   }
#   
#   
# }
# 
# temp_vec <- unlist(temp_list) 
# length(temp_vec) # 92
# active_post <- temp_vec[seq(2, length(temp_vec), 2)] - 0.5
# rest_post <- temp_vec[seq(1, length(temp_vec)-1, 2)] - 0.5
# 
# active_post <- c(0.5, active_post)
# rest_post <- c(rest_post, L + 0.5)
# 
# # length(active_post);length(rest_post)
# 
# active_span <- numeric(length(active_post))
# for(i in 1:length(active_post)){
#   active_span[i] <- rest_post[i] - active_post[i]
# }
# 
