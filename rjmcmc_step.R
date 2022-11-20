rm(list=ls())

##### data generate #####
n <- 50000
L <- 1
num <- c(10,12,13,10)
total_n <- sum(num)
tau <- c(0, 0.3, 0.5, 0.85, 1)
mu <- c(0.1,0.4 ,0.8, 1.2)
sigma <- 0.1
t <- list(); y<- list()
for(i in 1:4){
  t[[i]] <- runif(num[i],tau[i],tau[i+1])
  y[[i]] <- rnorm(num[i],mu[i],sigma)
}
t <- sort(unlist(t)); y <- unlist(y)
########################


##### hyper parameters #####
alpha <- 1 
beta <- 1.5
lambda <- 3
constant <- 3.6/7
########################



##### functions #####
lkh <- function(l,u,height){
  temp_y <- y[l<=t & t<u]
  re <- sum(log(dnorm(temp_y,height,sigma)))
  return(re)
}
birth <- function(k){
  return(constant*min(1,lambda/(k+1))) # k -> k+1
}
death<-function(k){
  prob <- ifelse(k==0,0,constant*min(1,k/lambda)) # k -> k-1
  return(prob)
}
########################


## height change ##
acpt.height <- function(h,h_prime,s,j){

  log_like_ratio <- lkh(s[j],s[j+1],h_prime)-lkh(s[j],s[j+1],h[j])
  log_prior_ratio <- alpha*(log(h_prime)-log(h[j]))-beta*(h_prime-h[j])
  log_pi_ratio <- log_like_ratio + log_prior_ratio
  return(exp(log_pi_ratio))
  
}
########################


## position change ##
acpt.position <- function(s,s_prime,j,h){
  log_like_ratio <- lkh(s[j-1],s_prime,h[j-1])+lkh(s_prime,s[j+1],h[j])-lkh(s[j-1],s[j],h[j-1])-lkh(s[j],s[j+1],h[j])
  log_prior_ratio <- log(s[j+1]-s_prime)+log(s_prime-s[j-1])-log(s[j+1]-s[j])-log(s[j]-s[j-1])
  log_pi_ratio <- log_like_ratio + log_prior_ratio
  return(exp(log_pi_ratio))
}
########################


## birth of a new step ##
acpt.birth <- function(s,h,s_star,h1_prime,h2_prime,j,k){
  m1_prime <- length(t[s[j]<t & t<s_star])
  m2_prime <- length(t[s_star<t & t<s[j+1]])
  mj <- m1_prime + m2_prime
  log_like_ratio <- lkh(s[j],s_star,h1_prime) + lkh(s_star,s[j+1],h2_prime)  - lkh(s[j],s[j+1],h[j])
  log_prior_ratio <- log(lambda)-log(k+1)+log(2*(k+1)*(2*k+3))-2*log(L)+log(s_star-s[j])+log(s[j+1]-s_star)-log(s[j+1]-s[j])
                      +alpha*log(beta)-log(factorial(alpha))+(alpha-1)*(log(h1_prime)+log(h2_prime)-log(h[j]))-beta*(h1_prime+h2_prime-h[j])
  log_proposal_ratio <- log(death(k+1))-log(birth(k))+log(L)-log(k+1)
  log_jacobian <- 2*log(h1_prime+h2_prime)-log(h[j])
  log_pi_ratio <- log_like_ratio + log_prior_ratio + log_proposal_ratio + log_jacobian
  return(exp(log_pi_ratio))
}
########################


## death of a step ##
acpt.death <- function(s,h_prime,h1,h2,j,k){
  m1_prime <- length(t[s[j-1]<t & t<s[j]])
  m2_prime <- length(t[s[j]<t & t<s[j+1]])
  mj <- m1_prime + m2_prime
  log_like_ratio <- lkh(s[j-1],s[j+1],h_prime)-  lkh(s[j-1],s[j],h1) - lkh(s[j],s[j+1],h2)
  log_prior_ratio <- log(k+1)-log(lambda)-log(2*(k+1)*(2*k+3))+2*log(L)+log(s[j+1]-s[j-1])-log(s[j+1]-s[j])-log(s[j]-s[j-1])-alpha*log(beta)+log(factorial(alpha))+(alpha-1)*(log(h_prime)-log(h1)-log(h2))+beta*(h1+h2-h_prime)
  log_proposal_ratio <- log(birth(k-1))-log(death(k))-log(L)+log(k)
  log_jacobian <- log(h_prime)-2*log(h1+h2)
  log_pi_ratio <- log_like_ratio + log_prior_ratio + log_proposal_ratio + log_jacobian
  return(exp(log_pi_ratio))
}
########################

h0 <- mean(y)
s0 <- c(0,1)
k0 <- 0

  H <- list()
  S <- list()
  K <- vector(length=n)
  
  K[1] <- k0
  H[[1]] <- h0
  S[[1]] <- s0
  
  for(i in 2:n){
    k <- K[i-1]
    h <- H[[i-1]]
    s <- S[[i-1]]
    
    position_prob <- ifelse(k<=0, 0, 0.5*(1-birth(k)-death(k)))
    height_prob <- 1-birth(k)-death(k)-position_prob
    type <- runif(1)
    # height change
    if(type>1-height_prob){
      j <- sample(1:(k+1),size=1)
      u <- runif(1,-1,1)
      h_tilde <- h[j]*exp(u)
      U <- runif(1)
      h[j] <- ifelse(U<=acpt.height(h,h_tilde,s,j),h_tilde,h[j])
    }

    # position change
    if(type<=1-height_prob && type>1-height_prob-position_prob){
      j <- sample(1:k,size=1)+1
      s_tilde <- runif(1,s[j-1],s[j+1])
      U <- runif(1)
      s[j] <- ifelse(U<=acpt.position(s,s_tilde,j,h),s_tilde,s[j])
    }
    
    # death of a step
    if(type>=birth(k) && type<=birth(k)+death(k)){
      j <- sample(1:k,size=1)+1
      r <- (s[j+1]-s[j]) / (s[j+1]-s[j-1])
      hj_prime <- h[j-1]^(1-r) * h[j]^r
      U <- runif(1)
      if(U<=acpt.death(s,hj_prime,h[j-1],h[j],j,k)){
        k <- k-1
        h[j-1] <- hj_prime
        h <- h[-j]
        s <- s[-j]
      }
      
    }
    
    # birth of a new step
    if(type<=birth(k)){
      j <- sample(1:(k+1),size=1)
      s_star <- runif(1,s[j],s[j+1])
      u <- runif(1)
      r <- exp(log(s[j+1]-s_star)-log(s[j+1]-s[j]))
      h1_prime <- h[j]*exp(r*(log(u)-log(1-u)))
      h2_prime <- h[j]*exp((1-r)* (log(1-u)-log(u)))
      U <- runif(1)
      if(U<= acpt.birth(s,h,s_star,h1_prime,h2_prime,j,k)){
        s <- c(s[1:j],s_star,s[(j+1):(k+2)])
        if(j==1){
          left <- c()
        }
        else{
          left <- h[1:(j-1)]
        }
        if((j-1)==k){
          right <- c()
        }
        else{
          right <- h[(j+1):(k+1)]
        }
        h <- c(left,h1_prime,h2_prime,right)
        k <- k+1
      }
    }
    
    K[i] <- k
    H[[i]] <- h
    S[[i]] <- s
    
    print(i)
    
  }

table(K) # 3

# step: 3
temp1 <- unlist(lapply(S,function(x){length(x)}))
temp2 <- S[which(temp1==5)]
temp3 <- unlist(temp2)
temp4 <- matrix(temp3,length(which(temp1==5)),5,byrow=T)
s_mat <- temp4[,-c(1,5)]
s_est <- apply(temp4,2,mean)


# height
temp1 <- unlist(lapply(H,function(x){length(x)}))
temp2 <- H[which(temp1==4)]
temp3 <- unlist(temp2)
temp4 <- matrix(temp3,length(which(temp1==4)),4,byrow=T)
h_mat <- temp4
h_est <- apply(temp4,2,mean)
h_est<- c(h_est,h_est[4])

plot(t,y, xlim=c(0,0.99),ylim=c(0,1.4),xlab="position",ylab="height")
par(new=T)
plot(s_est,h_est,type="s",col="red", xlim=c(0,0.99),ylim=c(0,1.4),xlab="position",ylab="height",main="u ~ U[-1,1]")


  
  


