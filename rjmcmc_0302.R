rm(list=ls())

current_time <- Sys.time()
library(mvtnorm)
library(dplyr)
library(truncnorm)

set.seed(0202)
dataset <- read.csv("/Users/wonny/Downloads/actigraph_csv.csv")
source('/Users/wonny/Documents/Reversible_Jump_MCMC/rjmcmc_ftn_0302.R')

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

## rest group
y_rest_group <- list()
id <- which(t_rest[-1] - t_rest[-length(t_rest)]!=1)

for(i in 1:(length(id)+1)){
  
  if(i==1){
    y_rest_group[[i]] <- y[t_rest[1]:t_rest[id[1]]]
  }else if(i==(length(id)+1)){
    y_rest_group[[i]] <- y[t_rest[id[(i-1)]+1]:t_rest[length(t_rest)]]
  }else{
    y_rest_group[[i]] <- y[t_rest[id[(i-1)]+1]:t_rest[id[i]]]
  }
  
}

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

# HP parameters (lambda_1,...lambda_r, alpha, beta)
lambda1 <- rep(3, n_rest)
lambda2 <- rep(mean(y_rest[y_rest>quantile(y_rest)[4]]), n_rest)

if(lambda2[1] > mu){
  lambda2 <- rep(mu - sqrt(lambda2[1])-5, n_rest)
}

## indicator for y_rest
indi <- list()
for(i in 1:n_rest){
  
  temp_y <- y_rest_group[[i]]
  temp_indi <- numeric(length(temp_y))
  for(j in 1:length(temp_y)){
    p1 <- dpois(temp_y[j], lambda1[i])
    p2 <- dpois(temp_y[j],lambda2[i])
    prob1 <- p1/(p1+p2)
    prob2 <- 1-prob1
    temp_indi[j] <- ifelse(prob1<prob2, 2, 1)
  }
  indi[[i]] <- temp_indi
  
  
}

indi_v <- unlist(indi)

# Hyperparameters of alpha1, beta1, alpha2, beta2
EY1 <- mean(y_rest[y_rest<=10])
VarY1 <- var(y_rest[y_rest<=10])
EY2 <- mean(y_rest[y_rest>10])
VarY2 <- var(y_rest[y_rest>10])

temp_ftn <- function(ey, vary){
  
  hyperprior_var <- 10
  app_b1 <- vary/ey
  app_a1 <- ey/app_b1
  app_D1 <- hyperprior_var/app_a1
  app_C1 <- app_a1/app_D1
  app_D2 <- hyperprior_var/app_b1
  app_C2 <- app_b1/app_D2
  
  return(list(C1=app_C1, D1=app_D1, C2=app_C2, D2=app_D2))
  
}
C11 <- temp_ftn(EY1, VarY1)$C1
D11 <- temp_ftn(EY1, VarY1)$D1
C12 <- temp_ftn(EY1, VarY1)$C2
D12 <- temp_ftn(EY1, VarY1)$D2

C21 <- temp_ftn(EY2, VarY2)$C1
D21 <- temp_ftn(EY2, VarY2)$D1
C22 <- temp_ftn(EY2, VarY2)$C2
D22 <- temp_ftn(EY2, VarY2)$D2


alpha1 <- 2
beta1 <- 1
alpha2 <- 40
beta2 <- 1
# alpha1 <- rgamma(1, 2, scale=1)
# beta1 <- rgamma(1, 40, scale=1)
# alpha2 <- rgamma(1, C21, scale=D21)
# beta2 <- rgamma(1, C22, scale=D22)
# iteration
iter <- 10^4

# store
mu_store <- sigma2_store <- l_store  <- beta1_store <- alpha1_store <- beta2_store <- alpha2_store <- n_active_store <- n_rest_store <- numeric(iter)
active_list <- rest_list <- lambda1_list <- lambda2_list <- list()

mu_store[1] <- mu
sigma2_store[1] <- sigma2
l_store[1] <- l
alpha1_store[1] <- alpha1
beta1_store[1] <- beta1
alpha2_store[1] <- alpha2
beta2_store[1] <- beta2
lambda1_list[[1]] <- lambda1
lambda2_list[[1]] <- lambda2

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
  
  #####################################################################################
  
  
  ## 1st group update ## 
  
  # lambda1
  lambda1 <- sapply(1:length(lambda1), function(x){rgamma(1, shape = alpha1 + sum(y_rest_group[[x]][indi[[x]]==1]), scale = beta1 / (length(y_rest_group[[x]][indi[[x]]==1])*beta1 +1))})
  
  # alpha1
  # alpha1_q <- rtruncnorm(1, a = 0, b= Inf, mean= alpha1, sd=0.5)
  # 
  # num <- sum(dgamma(lambda1, alpha1_q, scale = beta1, log=T)) + dgamma(alpha1_q, C11, scale= D11, log=T ) + log(dtruncnorm(alpha1, 0, Inf, alpha1_q, 0.5))
  # den <- sum(dgamma(lambda1, alpha1, scale = beta1, log=T)) + dgamma(alpha1, C11, scale= D11, log=T ) + log(dtruncnorm(alpha1_q, 0, Inf, alpha1, 0.5))
  # 
  # log_acc <- num - den
  # 
  # if(log(runif(1)) < log_acc) alpha1 <- alpha1_q
  # 
  # # beta1
  # beta1_q <- rtruncnorm(1, 0, Inf, beta1, sd=0.5)
  # 
  # num <- sum(dgamma(lambda1, alpha1, scale = beta1_q, log=T)) + dgamma(beta1_q, C12, scale= D12, log=T ) + log(dtruncnorm(beta1, 0, Inf, beta1_q, 0.5))
  # den <- sum(dgamma(lambda1, alpha1, scale = beta1, log=T)) + dgamma(beta1, C12, scale= D12, log=T ) + log(dtruncnorm(beta1_q, 0, Inf, beta1, 0.5))
  # 
  # 
  # log_acc <- num - den
  # if(log(runif(1)) < log_acc) beta1 <- beta1_q
  # 
  # 
  
  #####################################################################################
  
  
  ## lambda2 update ##
  
  ## 2nd group update ## 
  
  # lambda2
  origin_lambda2 <- lambda2
  lambda2 <- sapply(1:length(lambda2), function(x){rgamma(1, shape = alpha2 + sum(y_rest_group[[x]][indi[[x]]==2]), scale = beta2 / (length(y_rest_group[[x]][indi[[x]]==2])*beta2 +1))})
  
  if(length(which(lambda2 < lambda1 | lambda2 > (mu-sqrt(lambda2))))!=0){
    
    lambda2[which(lambda2 < lambda1 | lambda2 > (mu-sqrt(lambda2)))] <- origin_lambda2[which(lambda2 < lambda1 | lambda2 > (mu-sqrt(lambda2)))]
    
  }
  
  
  # alpha1
  # alpha2_q <- rtruncnorm(1, a = 0, b= Inf, mean= alpha2, sd=0.5)
  # 
  # num <- sum(dgamma(lambda2, alpha2_q, scale = beta2, log=T)) + dgamma(alpha2_q, C21, scale= D21, log=T ) + log(dtruncnorm(alpha2, 0, Inf, alpha2_q, 0.5))
  # den <- sum(dgamma(lambda2, alpha2, scale = beta2, log=T)) + dgamma(alpha2, C21, scale= D21, log=T ) + log(dtruncnorm(alpha2_q, 0, Inf, alpha2, 0.5))
  # 
  # log_acc <- num - den
  # 
  # if(log(runif(1)) < log_acc) alpha2 <- alpha2_q
  # 
  # # beta1
  # beta2_q <- rtruncnorm(1, 0, Inf, beta2, sd=0.5)
  # 
  # num <- sum(dgamma(lambda2, alpha2, scale = beta2_q, log=T)) + dgamma(beta2_q, C22, scale= D22, log=T ) + log(dtruncnorm(beta2, 0, Inf, beta2_q, 0.5))
  # den <- sum(dgamma(lambda2, alpha2, scale = beta2, log=T)) + dgamma(beta2, C22, scale= D22, log=T ) + log(dtruncnorm(beta2_q, 0, Inf, beta2, 0.5))
  # 
  # 
  # log_acc <- num - den
  # if(log(runif(1)) < log_acc) beta2 <- beta2_q
  # 
  # 
  
  #####################################################################################
  hyperparameters <- c(alpha1, beta1, alpha2, beta2)
  ## position change / birth / death 
  
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
        
        if(current!=proposal){ # if == -> no need to update
          
          if(k%%2==0){
            
            rest_q <- rest
            rest_q[k/2] <- proposal
            
            if(proposal < current){
              
              # A -> R
              
              y_star <-  y[t> proposal & t<current]
              t_star <- t[t>proposal & t<current]
              
              # define the indicators
              temp_indi <- Indi_ftn(y_star, lambda1[k/2], lambda2[k/2])
              
              indi1 <- ifelse(sum(temp_indi==1)==0, 0, 1)
              indi2 <- ifelse(sum(temp_indi==2)==0, 0, 1)
              
              log_num <- sum(indi1 * dpois(y_star[temp_indi==1], lambda1[k/2], log=T)) + sum(indi2 * dpois(y_star[temp_indi==2], lambda2[k/2], log=T))
              
              log_den <- log_mvn_lkh(y_star, t_star, mu, sigma2, l, g)
              
              if(log(runif(1)) < (log_num - log_den)){
                rest <- rest_q
              }
              
            }else{
              
              # R -> A
              y_star <-  y[t> current & t< proposal]
              t_star <- t[t> current & t< proposal]
              
              
              temp_indi <- Indi_ftn(y_star, lambda1[k/2], lambda2[k/2])
              
              indi1 <- ifelse(sum(temp_indi==1)==0, 0, 1)
              indi2 <- ifelse(sum(temp_indi==2)==0, 0, 1)
              
              log_num <- log_mvn_lkh(y_star, t_star, mu, sigma2, l, g)
              
              log_den <- sum(indi1*dpois(y_star[temp_indi==1], lambda1[k/2], log=T))+sum(indi2*dpois(y_star[temp_indi==2], lambda2[k/2], log=T))
              
              if(log(runif(1)) < ( log_num - log_den )){
                rest <- rest_q
              }
            } # rest case end
            
          }else{
            
            active_q <- active
            active_q[(k+1)/2] <- proposal
            
            if(proposal < current){
              
              # R -> A
              
              y_star <-  y[t> proposal & t<current]
              t_star <- t[t>proposal & t<current]
              
              temp_indi <- Indi_ftn(y_star, lambda1[(k-1)/2], lambda2[(k-1)/2])
              
              indi1 <- ifelse(sum(temp_indi==1)==0, 0, 1)
              indi2 <- ifelse(sum(temp_indi==2)==0, 0, 1)
              
              log_num <- log_mvn_lkh(y_star, t_star, mu, sigma2, l, g)
              
              log_den <- sum(indi1*dpois(y_star[temp_indi==1], lambda1[(k-1)/2], log=T)) + sum(indi2*dpois(y_star[temp_indi==2], lambda2[(k-1)/2], log=T))
              
              if(log(runif(1)) < (log_num - log_den)){
                active <- active_q
              }
            }else{
              
              # A -> R
              y_star <-  y[t> current & t< proposal]
              t_star <- t[t> current & t< proposal]
              
              temp_indi <- Indi_ftn(y_star, lambda1[(k-1)/2], lambda2[(k-1)/2])
              
              indi1 <- ifelse(sum(temp_indi==1)==0, 0, 1)
              indi2 <- ifelse(sum(temp_indi==2)==0, 0, 1)
              
              log_num <- sum(indi1*dpois(y_star[temp_indi==1], lambda1[(k-1)/2], log=T)) + sum(indi2*dpois(y_star[temp_indi==2], lambda2[(k-1)/2], log=T))
              
              log_den <- log_mvn_lkh(y_star, t_star, mu, sigma2, l, g)
              
              if(log(runif(1)) < (log_num - log_den)){
                active <- active_q
              }
            }
            
          } # active case end
          
          position <- c(rbind(active,rest))
          
        }
        
      }
    }# position update
    
  }else if(u < p_posit + p_birth){ # birth updates
    ##################### (a,r) birth #######################################
    
    if(n_rest >= 1){
      temp <- acc.birth.ar(y, t, mvn_parameters, lambda1, lambda2, hyperparameters, p_birth, p_death, active, rest, n_active, n_rest)
      
      active <- temp$active
      rest <- temp$rest
      n_active <- temp$n_active
      n_rest <- temp$n_rest
      acc_b_ar[i] <- temp$accept
      lambda1 <- temp$lambda1_new
      lambda2 <- temp$lambda2_new
      #print("a,r birth")
      #print(lambda1)
      #print(lambda2)
    }
    
    
    ##################### (r,a) birth #######################################
    
    if(n_active >= 1){
      temp <- acc.birth.ra(y, t, mvn_parameters, lambda1, lambda2, hyperparameters, p_birth, p_death, active, rest, n_active, n_rest)
      
      active <- temp$active
      rest <- temp$rest
      n_active <- temp$n_active
      n_rest <- temp$n_rest
      acc_b_ra[i] <- temp$accept
      lambda1 <- temp$lambda1_new
      lambda2 <- temp$lambda2_new
      #print("r,a birth")
      #print(lambda1)
      #print(lambda2)
    }
    
    
  }else{ # death updates
    
    ##################### (a,r) death #######################################
    #print("death")
    if(n_active>2){
      temp <- acc.death.ar(y, t, mvn_parameters, lambda1, lambda2, hyperparameters, p_birth, p_death, active, rest, n_active, n_rest)
      
      active <- temp$active
      rest <- temp$rest
      n_active <- temp$n_active
      n_rest <- temp$n_rest
      acc_d_ar[i] <- temp$accept
      lambda1 <- temp$lambda1_new
      lambda2 <- temp$lambda2_new
      death_j[i] <- temp$j
    }
    
    ##################### (r,a) death #######################################
    
    if(n_rest>1){
      temp <- acc.death.ra(y, t, mvn_parameters, lambda1, lambda2, hyperparameters, p_birth, p_death, active, rest, n_active, n_rest)
      
      active <- temp$active
      rest <- temp$rest
      n_active <- temp$n_active
      n_rest <- temp$n_rest
      acc_d_ra[i] <- temp$accept
      lambda1 <- temp$lambda1_new
      lambda2 <- temp$lambda2_new
      
    }
    
  } # death update end
  
  position <- c(rbind(active,rest))
  group_indicator <- as.numeric(cut(t,position))
  
  y_rest <- y[group_indicator%%2==0]; y_active <- y[group_indicator%%2==1]
  t_rest <- t[group_indicator%%2==0]; t_active <- t[group_indicator%%2==1]
  
  # new y_rest_group
  ## rest group
  y_rest_group <- list()
  id <- which(t_rest[-1] - t_rest[-length(t_rest)]!=1)
  
  for(ii in 1:(length(id)+1)){
    
    if(ii==1){
      y_rest_group[[ii]] <- y[t_rest[1]:t_rest[id[1]]]
    }else if(ii==(length(id)+1)){
      y_rest_group[[ii]] <- y[t_rest[id[(ii-1)]+1]:t_rest[length(t_rest)]]
    }else{
      y_rest_group[[ii]] <- y[t_rest[id[(ii-1)]+1]:t_rest[id[ii]]]
    }
    
  }
  
  indi <- list()
  for(k in 1:n_rest){
    
    temp_y <- y_rest_group[[k]]
    temp_indi <- numeric(length(temp_y))
    for(j in 1:length(temp_y)){
      p1 <- dpois(temp_y[j], lambda1[k],log=T)
      p2 <- dpois(temp_y[j],lambda2[k],log=T)
      prob1 <- 1/(1+exp(p2-p1))
      prob2 <- 1-prob1
      temp_indi[j] <- ifelse(prob1<prob2, 2, 1)
    }
    indi[[k]] <- temp_indi
    
    
  }
  
  ############### save #####################
  active_list[[i]] <- active
  rest_list[[i]] <- rest
  mu_store[i] <- mu
  sigma2_store[i] <- sigma2
  l_store[i] <- l
  lambda1_list[[i]] <- lambda1
  alpha1_store[i] <- alpha1
  beta1_store[i] <- beta1
  lambda2_list[[i]] <- lambda2
  alpha2_store[i] <- alpha2
  beta2_store[i] <- beta2
  
  
  n_active_store[i] <- n_active
  n_rest_store[i] <- n_rest
  
  y_active_list[[i]] <- y_active
  y_rest_list[[i]] <- y_rest
  t_active_list[[i]] <- t_active
  t_rest_list[[i]] <- t_rest
  indi_list[[i]] <- unlist(indi)
  
  # print(active)
  # print(rest)
  
  print(i)
  #print(lambda1)
  #print(lambda2)
  
  plot(t_active,y_active,type="h",col='red', xlim=c(0,length(t)),ylim=c(0,max(y_active)),xlab="",ylab="")
  par(new=T)
  plot(t_rest,y_rest,type="h",col='blue', xlim=c(0,length(t)),lwd=2, ylim=c(0,max(y_active)),xlab="",ylab="")
  
  
  
  
  
  
} # MCMC end 
