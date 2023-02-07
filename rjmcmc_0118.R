rm(list=ls())

current_time <- Sys.time()
library(mvtnorm)
library(dplyr)
library(truncnorm)

set.seed(0202)
dataset <- read.csv("/Users/wonny/Downloads/actigraph_csv.csv")
source('/Users/wonny/Documents/Reversible_Jump_MCMC/rjmcmc_ftn_0118.R')

subject1 <- subset(dataset, Index==1)
data <- (subject1$axis1 + subject1$axis2 + subject1$axis3)/3

len <- length(data)/10
newdata <- numeric(len)

for(i in 1:len){
  newdata[i] <- mean(data[(10*i-9):(10*i)])
}

y <- round(newdata)
t <- 1:length(y)
L <- t[length(t)]+0.5
# plot(y,type="h")


## temporarily use a portion of data
# y <- y[1:200]
# t <- t[1:200]
# L <- t[length(t)]+0.5
# plot(y,type="h")




################################ initial values ##############################################

Init <- initial_position(y, 50)
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

# ZIP parameters (pi, lambda)
pi <- 0.7
lambda <- mean(y[y!=0 & y<100])

# iteration
iter <- 2 * 10^4

# store
mu_store <- sigma2_store <- l_store  <- pi_store <- lambda_store <- n_active_store <- n_rest_store <- numeric(iter)
active_list <- rest_list <- list()

mu_store[1] <- mu
sigma2_store[1] <- sigma2
l_store[1] <- l
pi_store[1] <- pi
lambda_store[1] <- lambda

y_active_list <- t_active_list <- y_rest_list <- t_rest_list <- list()
acc_b_ar <- acc_b_ra <- acc_d_ar <- acc_d_ra <- numeric(iter)


################################ MCMC run ##############################################

for(i in 2:iter){
  
  u <- runif(1)

  p_posit <- 1/3
  p_birth <- 1/3
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
    
    
    ############### ZIP parameters update #####################
    # 1. pi (0 < pi <1)
    pi_q <- rtruncnorm(1,0,1,pi,0.5)
    log_acc <- log_zip_lkh(y_rest,pi_q,lambda) + log(dtruncnorm(pi,0,1,pi_q,1)) - log_zip_lkh(y_rest,pi,lambda) - log(dtruncnorm(pi_q,0,1,pi,1)) 
    pi <- ifelse(log(runif(1)) < log_acc, pi_q, pi)
    # 2. lambda
    log.lambda_q <- rnorm(1, log(lambda), 0.1)
    lambda_q <- exp(log.lambda_q)
    log_acc <- log_zip_lkh(y_rest,pi,lambda_q) + log.lambda_q - log_zip_lkh(y_rest,pi,lambda) - log(lambda)
    
    lambda <- ifelse(log(runif(1)) < log_acc, lambda_q, lambda)

    mvn_parameters <- c(mu, sigma2, l, g)
    zip_parameters <- c(pi,lambda)
    
  if(u < p_posit){ # position updates 
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
        position_q <- position
        position_q[k] <- proposal
        
        c_indi <- as.numeric(cut(t,position))
        p_indi <- as.numeric(cut(t,position_q))
        
        y_rest <- y[c_indi%%2==0]
        y_rest_q <- y[p_indi%%2==0]
        
        if(k%%2==0){
          
          rest_q <- rest
          rest_q[k/2] <- proposal
          
          new_lkh <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest_q, n_active) + log_zip_lkh(y_rest_q, pi, lambda)
          old_lkh <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active) + log_zip_lkh(y_rest, pi, lambda)
          
          log_acc <- new_lkh  - old_lkh
          if(log(runif(1)) < log_acc){
            rest <- rest_q
          }
          
        }else{
          
          active_q <- active
          active_q[(k+1)/2] <- proposal
          
          new_lkh <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active_q, rest, n_active) + log_zip_lkh(y_rest_q, pi, lambda)
          old_lkh <- sum_log_mvn_lkh(y, t, mu, sigma2, l, g, active, rest, n_active) + log_zip_lkh(y_rest, pi, lambda)
          
          log_acc <- new_lkh - old_lkh
          if(log(runif(1)) < log_acc){
            active <- active_q
          }
        }
        
        position <- c(rbind(active,rest))
        
        
      }# position update
    }
    
  }else if(u < p_posit + p_birth){ # birth updates
    ##################### (a,r) birth #######################################
    
    if(n_rest >= 1){
      temp <- acc.birth.ar(y, t, mvn_parameters, zip_parameters, p_birth, p_death, active, rest, n_active, n_rest)
      
      active <- temp$active
      rest <- temp$rest
      n_active <- temp$n_active
      n_rest <- temp$n_rest
      acc_b_ar[i] <- temp$accept
    }


    ##################### (r,a) birth #######################################
    
    if(n_active >= 1){
      temp <- acc.birth.ra(y, t, mvn_parameters, zip_parameters, p_birth, p_death, active, rest, n_active, n_rest)
      
      active <- temp$active
      rest <- temp$rest
      n_active <- temp$n_active
      n_rest <- temp$n_rest
      acc_b_ra[i] <- temp$accept
      
    }


  }else{ # death updates
    
    ##################### (a,r) death #######################################
    
    if(n_active>2){
      temp <- acc.death.ar(y, t, mvn_parameters, zip_parameters, p_birth, p_death, active, rest, n_active, n_rest)
      
      active <- temp$active
      rest <- temp$rest
      n_active <- temp$n_active
      n_rest <- temp$n_rest
      acc_d_ar[i] <- temp$accept
      
    }
    
    ##################### (r,a) death #######################################
    
    if(n_rest>1){
      temp <- acc.death.ra(y, t, mvn_parameters, zip_parameters, p_birth, p_death, active, rest, n_active, n_rest)
      
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

  
  ############### save #####################
  active_list[[i]] <- active
  rest_list[[i]] <- rest
  mu_store[i] <- mu
  sigma2_store[i] <- sigma2
  l_store[i] <- l
  pi_store[i] <- pi
  lambda_store[i] <- lambda
  n_active_store[i] <- n_active
  n_rest_store[i] <- n_rest
  
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
abline(v=active,col="green",lwd=1)
abline(v=rest,col="blue")



end_time <-Sys.time()
end_time-current_time

par(mfrow=c(3,2))
plot(mu_store, type="l")
plot(sigma2_store, type="l")
plot(l_store, type="l")
plot(pi_store, type="l")
plot(lambda_store, type="l")
plot(n_active_store, type="l")

# save(t_post, mu_store, sigma2_store, l_store, pi_store, lambda_store, n_active_store, n_rest_store, file="result0107.rda")

