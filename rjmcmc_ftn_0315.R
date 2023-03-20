initial_position <- function(y_input, cutoff){
  
  Signs <- sign(y_input - cutoff)
  Num  <- length(y_input)
  changePoints <- Indi <- numeric(Num)
  
  for(i in 2:Num){
    if(Signs[i]!=Signs[i-1]){
      changePoints[i] <- i
      
      Indi[i] <- ifelse(Signs[i] > Signs[i-1], "A" , "R")
      
    }
  }
  
  CPs <- changePoints[changePoints!=0]  
  Indi <- Indi[Indi!="0"]
  LthCP <- length(CPs)
  
  if(LthCP%%2 == 0){
    
    if(Indi[1]=="A"){
      active <- CPs[seq(1, LthCP-1, by =2)]
      rest <- CPs[seq(2, LthCP, by=2)]
      
      active <- active[-1]
      if(Indi[length(Indi)]=="R"){
        rest <- rest[-length(rest)]
      }
    }else{
      rest <- CPs[seq(1, LthCP-1, by =2)]
      active <- CPs[seq(2, LthCP, by=2)]
      if(Indi[length(Indi)]=="R"){
        rest <- rest[-length(rest)]
      }
    }
    
  }else{
    
    if(Indi[1]=="A"){
      active <- CPs[seq(1, LthCP, by =2)]
      rest <- CPs[seq(2, LthCP-1, by=2)]
      
      active <- active[-1]
      if(Indi[length(Indi)]=="R"){
        rest <- rest[-length(rest)]
      }
      
    }else{
      rest <- CPs[seq(1, LthCP, by =2)]
      active <- CPs[seq(2, LthCP-1, by=2)]
      if(Indi[length(Indi)]=="R"){
        rest <- rest[-length(rest)]
      }
    }
    
  }
  
  active <- c(0.5, active-0.5)
  rest <- c(rest-0.5, L)
  
  return(list(active=active, rest=rest))
  
}

log_den_y <- function(y_input, Pi_input, Lambda_input){
  
  
  if(length(y_input)==1){
    value <- log(sum(Pi_input * dpois(y_input, Lambda_input)))
  }else{
    
    value <- 0
    for(t in 1:length(y_input)){
      value <- value + log(sum(Pi_input * dpois(y_input[t], Lambda_input)))
    }
    
  }
  
  return(value)
  }


kernel_ftn <- function(t_input, l_input, g_input){
  K <- exp(-abs(outer(t_input,t_input,FUN = "-"))/l_input)+g_input*diag(length(t_input))
  return(K)
} # covariance matrix = sigma^2 * kernel_ftn(t,l,g)


log_mvn_lkh <- function(y_input, t_input, mu_input, sigma2_input, l_input, g_input){
  d <- length(y_input)
  k_t_t_prime <- kernel_ftn(t_input, l_input, g_input)
  log_lkh <- dmvnorm(y_input, mean = rep(mu_input, d), sigma = sigma2_input * k_t_t_prime ,log=T )
  return(log_lkh)
}

sum_log_mvn_lkh <- function(y_input, t_input, mu_input, sigma2_input, l_input, g_input, active_input, rest_input, n_active_input){
  
  sum <- 0
  for(I in 1:n_active_input){
    
    r <- rest_input[I]
    t_set <- t_input[ t_input > active_input[I] & t_input < r]
    y_set <- y_input[ t_input > active_input[I] & t_input < r ]
    
    if(length(y_set)!=0){
      sum <- sum + log_mvn_lkh(y_set, t_set, mu_input, sigma2_input, l_input, g_input)
    }
    
  }
  
  return(sum)
  
}

##########################################################################################

acc.birth.ar <- function(y_input, t_input, mvn_inputs, Pi_input, Lambda_input, p_birth_input, p_death_input, active_input, rest_input, n_active_input, n_rest_input){
  
  
  L <- t_input[length(t_input)]+0.5
  accept <- 0
  a_len <- length(active_input); r_len <- length(rest_input)
  
  # choose one (r,a) where (a*,r*) is proposed based on the prob_rest
  ar_prob <- prob_return_birth(active_input, rest_input, n_active_input, n_rest_input)$rest_prob
  j <- sample(1:n_rest_input, 1, prob=ar_prob) 
  start <- rest_input[j]
  end <- active_input[j+1]
  
  # selected rest interval = (start, end) = r[j],a[j+1]
  
  if((end-start)>2){
    
    temp <- sample( seq(start+1, end-1, by=1), 2 )
    a_star <- min(temp); r_star <- max(temp) # proposed (a*, r*)
    
    t_star <- t_input[t_input>a_star & t_input<r_star]
    y_star <- y_input[t_input>a_star & t_input<r_star]
    
    if(length(y_star)!=0){
      k <- length(c(active_input,rest_input))
      
      # mvn_inputs : (mu, sigma^2, l)
      mu <- mvn_inputs[1]
      sigma2 <- mvn_inputs[2]
      l <- mvn_inputs[3]
      g <- mvn_inputs[4]
      
      
      # for the proposal ratio
      den <- (ar_prob[j]/sum(ar_prob)) * 2 * (1/(end-start-1))^2  * p_birth_input
      num <- (r_star - a_star)/(sum(prob_return_death(active_input, rest_input, n_active_input, n_rest_input)$active_prob) + (r_star - a_star)) * p_death_input
      
      # indicator
      log_lkh_ratio <- log_mvn_lkh(y_star, t_star, mu, sigma2, l, g) - log_den_y(y_star, Pi_input, Lambda_input)
      log_step_ratio <- log(k+1) + k*(log(k)-log(k+2))-log(k+2)
      log_proposal_ratio <- log(num)-log(den)
      
      log_acc <- log_lkh_ratio + log_step_ratio + log_proposal_ratio
      accept <- ifelse(log(runif(1))<log_acc,1,0) # 1:acc, 0: reject
      
    }
    
    if(accept==1){
      if(j==1){
        
        rest <- c(rest_input[1:j],r_star,rest_input[(j+1):r_len])
        active <- c(active_input[1:j],a_star,active_input[(j+1):a_len])
        
        
      }else if(j==n_rest_input){
        
        rest <- c(rest_input[1:j],r_star,rest_input[(j+1):r_len])
        
        active <- c(active_input[1:j],a_star,active_input[(j+1):a_len])
        
      }else{
        
        rest <- c(rest_input[1:j],r_star,rest_input[(j+1):r_len])
        active <- c(active_input[1:j],a_star,active_input[(j+1):a_len])
        
      }
      a_len <- length(active); r_len <- length(rest)
      n_active <- a_len
      n_rest <- r_len-1
      
      
    }
    
  }
  return(list(active=active, rest=rest, n_active=n_active, n_rest=n_rest, accept=accept))
}

##########################################################################################

acc.death.ar <- function(y_input, t_input, mvn_inputs, Pi_input, Lambda_input, p_birth_input, p_death_input, active_input, rest_input, n_active_input, n_rest_input){
  
  
  
  accept <- 0
  a_len <- length(active_input); r_len <- length(rest_input)
  L <- t_input[length(t_input)]+0.5
  ar_prob <- prob_return_death(active_input, rest_input, n_active_input, n_rest_input)$active_prob
  
  if(n_active_input==3){
    j <- 2
  }else{
    j <- sample(2:(n_active_input-1), 1, prob=ar_prob)
  }
  
  start <- active_input[j]
  end <- rest_input[j]
  # selected active interval = (start, end) = a[j],r[j]
  
  a_star <- start; r_star <- end
  
  t_star <- t_input[t_input>a_star & t_input<r_star]
  y_star <- y_input[t_input>a_star & t_input<r_star]
  
  if(length(y_star)!=0){
    
    k <- length(c(active_input,rest_input))
    
    mu <- mvn_inputs[1]
    sigma2 <- mvn_inputs[2]
    l <- mvn_inputs[3]
    g <- mvn_inputs[4]
    
    # for the proposal ratio
    
    new_sum <- (sum(prob_return_birth(active_input, rest_input, n_active_input, n_rest_input)$rest_prob)
                - (active_input[j] - rest_input[j-1])
                - (active_input[j+1] - rest_input[j]) 
                + (active_input[j+1]-rest_input[j-1]))
    
    den <- ar_prob[j-1] / sum(ar_prob) * p_death_input
    num <- (active_input[j+1] - rest_input[j-1])/new_sum * 2 * (1/(active_input[j+1] - rest_input[j-1]-1))^2 * p_birth_input
    
    # indicator

    log_lkh_ratio <-  log_den_y(y_star, Pi_input, Lambda_input) - log_mvn_lkh(y_star, t_star, mu, sigma2, l, g) 
    log_step_ratio <- -log(k)-log(k+1)+k*(log(k)-log(k-2))+2*log(k-2)
    log_proposal_ratio <- log(num)- log(den)
    
    log_acc <- log_lkh_ratio + log_step_ratio + log_proposal_ratio
    accept <- ifelse(log(runif(1))<log_acc,1,0) 
    
    if(accept==1){
      active <- active_input[-j]
      rest <- rest_input[-j]
      
      
      if(j==1){
        rest <- c(0.5,rest_input)
      }else if(j==n_active_input){
        active <- c(active_input,L)
      }
      
      a_len <- length(active); r_len <- length(rest)
      n_active <- a_len
      n_rest <- r_len-1
      
    }
    
  }
  
  
  return(list(active=active, rest=rest, n_active=n_active, n_rest=n_rest, accept=accept, j=j))
}

##########################################################################################

acc.birth.ra <- function(y_input, t_input, mvn_inputs, Pi_input, Lambda_input, p_birth_input, p_death_input, active_input, rest_input, n_active_input, n_rest_input){
  
  accept <- 0
  L <- t_input[length(t_input)]+0.5
  a_len <- length(active_input); r_len <- length(rest_input)
  
  ra_prob <- prob_return_birth(active_input, rest_input, n_active_input, n_rest_input)$active_prob
  j <- sample(1:n_active_input, 1, prob=ra_prob) 
  start <- active_input[j]
  end <- rest_input[j]
  if((end-start)>2){
    
    temp <- sample( seq(start+1, end-1, by=1), 2 )
    r_star <- min(temp); a_star <- max(temp) # proposed (r*, a*)
    
    t_star <- t_input[t_input>r_star & t_input<a_star]
    y_star <- y_input[t_input>r_star & t_input<a_star]
    
    if(length(y_star)!=0){
      k <- length(c(active_input,rest_input))
      
      mu <- mvn_inputs[1]
      sigma2 <- mvn_inputs[2]
      l <- mvn_inputs[3]
      g <- mvn_inputs[4]
      
      
      # for the proposal ratio
      den <- (ra_prob[j]/sum(ra_prob)) * 2 * (1/(end-start-1))^2 * p_birth_input
      num <- (a_star - r_star) / (sum(prob_return_death(active_input, rest_input, n_active_input, n_rest_input)$rest_prob) + (a_star-r_star)) * p_death_input
      

      log_lkh_ratio <- log_den_y(y_star, Pi_input, Lambda_input) - log_mvn_lkh(y_star, t_star, mu, sigma2, l, g) 
      log_step_ratio <- log(k+1) + k*(log(k)-log(k+2))-log(k+2)
      log_proposal_ratio <- log(num) - log(den)
      
      log_acc <- log_lkh_ratio + log_step_ratio + log_proposal_ratio
      accept <- ifelse(log(runif(1))<log_acc,1,0)
      
    }
    
    if(accept==1){
      
      
      if(j==1){
        
        rest <- c(r_star,rest_input)
        
        active <- c(active_input[1:j],a_star,active_input[(j+1):a_len])
        
      }else if(j==n_active_input){
        rest <- c(rest_input[1:(j-1)],r_star,rest_input[j:r_len])
        
        active <- c(active_input,a_star)
        
      }else{
        rest <- c(rest_input[1:(j-1)],r_star,rest_input[j:r_len])
        
        active <- c(active_input[1:j],a_star,active_input[(j+1):a_len])
      }
      
      a_len <- length(active); r_len <- length(rest)
      n_active <- a_len
      n_rest <- r_len-1
      
      
    }
    
  }
  
  return(list(active=active, rest=rest, n_active=n_active, n_rest=n_rest, accept=accept))
}

##########################################################################################

acc.death.ra <- function(y_input, t_input, mvn_inputs, Pi_input, Lambda_input , p_birth_input, p_death_input, active_input, rest_input, n_active_input, n_rest_input){
  
  
  L <- t_input[length(t_input)]+0.5
  accept <- 0
  ra_prob <- prob_return_death(active_input, rest_input, n_active_input, n_rest_input)$rest_prob
  j <- sample(1:n_rest_input, 1, prob=ra_prob)
  
  start <- rest_input[j]
  end <- active_input[j+1]
  
  r_star <- start; a_star <- end
  
  t_star <- t_input[t_input>r_star & t_input<a_star]
  y_star <- y_input[t_input>r_star & t_input<a_star]
  
  if(length(y_star)!=0){
    
    k <- length(c(active_input,rest_input))
    
    mu <- mvn_inputs[1]
    sigma2 <- mvn_inputs[2]
    l <- mvn_inputs[3]
    g <- mvn_inputs[4]
    
    
    # for the proposal ratio
    new_sum <- (sum(prob_return_birth(active_input, rest_input, n_active_input, n_rest_input)$active_prob)
                -(rest_input[j] - active_input[j])
                -(rest_input[j+1] - active_input[j+1])
                + rest_input[j+1]- active_input[j])
    
    den <- ra_prob[j]/sum(ra_prob) * p_birth_input
    num <- (rest_input[j+1]- active_input[j])/new_sum * 2 * (1/(rest_input[j+1]-active_input[j]))^2 * p_death_input
    
    
    log_lkh_ratio <- log_mvn_lkh(y_star, t_star, mu, sigma2, l, g) -log_den_y(y_star, Pi_input, Lambda_input)
    log_step_ratio <- -log(k)-log(k+1)+k*(log(k)-log(k-2))+2*log(k-2)
    log_proposal_ratio <- log(num) - log(den)
    
    log_acc <- log_lkh_ratio + log_step_ratio + log_proposal_ratio
    accept <- ifelse(log(runif(1))<log_acc,1,0)
    
    if(accept==1){
      
      print(log_lkh_ratio)
      rest <- rest_input[-j]
      active <- active_input[-(j+1)]
      
      a_len <- length(active); r_len <- length(rest)
      n_active <- a_len
      n_rest <- r_len-1
      
    }
    
  }
  return(list(active=active, rest=rest, n_active=n_active, n_rest=n_rest, accept=accept))
}


prob_return_birth<- function(active_input, rest_input, n_active_input, n_rest_input){
  
  
  active_length <- numeric(n_active_input)
  rest_length <- numeric(n_rest_input)
  
  for(I in 1:n_active_input){
    active_length[I] <- rest_input[I]-active_input[I]
  }
  
  for(I in 1:n_rest_input){
    rest_length[I] <- active_input[I+1]-rest_input[I]
  }
  
  active_prob <- active_length/ sum(active_length)
  rest_prob <- rest_length/ sum(rest_length)
  
  return(list(active_prob=active_prob, rest_prob=rest_prob))
  
  
}



prob_return_death<- function(active_input, rest_input, n_active_input, n_rest_input){
  
  
  active_length <- numeric(n_active_input-2)
  rest_length <- numeric(n_rest_input)
  
  for(I in 2:(n_active_input-1)){
    active_length[I-1] <- rest_input[I]-active_input[I]
  }
  
  for(I in 1:n_rest_input){
    rest_length[I] <- active_input[I+1]-rest_input[I]
  }
  
  ## same but inverse
  active_length <- 1/active_length
  rest_length <- 1/rest_length
  
  active_prob <- active_length/ sum(active_length)
  rest_prob <- rest_length/ sum(rest_length)
  
  return(list(active_prob=active_prob, rest_prob=rest_prob))
  
}