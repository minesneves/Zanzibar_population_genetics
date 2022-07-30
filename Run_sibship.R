
###################
## packages
####################

## library for using the Sterling number function
library(gmp)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)
library(spatstat)
library(dplyr)
library(data.table)
library(ggpmisc)
library(plyr)


#####################
## global parameters
####################
Nmax <- 250  #Cheever 1977 
nsim <- 1000

###################
## functions
####################

# falling factorial
fallingfac <- function(N,n){
  if (N<=170) {
    exp( lfactorial(N) - lfactorial(N-n))
  } else { # Stirling approximation: its an equation
    #that allows you to calculate the natural logarithm for large factorials;
    exp ( (N*log(N)-N) - ((N-n)*log(N-n)-(N-n)) )
  }
}

fallingfac <- Vectorize(fallingfac, "N")


# unique items distribution (UID, from Mendelsen et al 2016)
#N(k)/N^A*Stirling approximation    # gives the probability of finding n unique items from a sample of m miracidia obtained from
# an original number N of female worms inside a host.
prob <- function(n, N, m) {
  fallingfac(N,n)/N^m*as.numeric(Stirling2(m, n))
}

# expectation of UID
expect <- function(N, m) {
  (N^m - (N-1)^m)/(N^(m-1))
}

## variance of UID
vari <- function(N, m) {
  N*(N-1)*(1-2/N)^m + N*(1-1/N)^m - N^2*(1-1/N)^(2*m)
}

## posterior approximation by sampling importance resampling
f <- function(par)
{
  n <- par[1]
  m <- par[2]
  ## mean of prior NBD
  mean <- par[3]
  ## overdispersion of prior NBD
  overdisp <- par[4]
  prior <- par[5]
  
  ## SIR algorthim 
  # 1. sample a bunch of values (nsim) from a uniform dist with min = n worms and max= 350 worms
  # Nstar - needs to be an approximate of the target distribution 
  
  Nstar <- round(runif(nsim, n, Nmax)) #gives a random sample of our initial "N" worms - this is the correct number of female worms in a host
  
  # Latin-hypercube sampling
  # Sampling from a multivariate normal distribution
  # Normalising function needs to be different
  
  
  # 2. calculate importance weights/importance ratios for each Xi
  # these are used as resampling weights to select the sample
  
  if (prior==1) {
    w <- prob(n, Nstar, m)*dnbinom(Nstar, mu=mean, size=overdisp)/   
      dunif(Nstar, n, Nmax)  # Negative binomial prior: Aqui multiplicamos duas probabilidades; 
    #a de encontrar n unique items numa amostra de m sampled do Nstar,
    # vezes a density probability function da negative binomial distribution.
  }
  else {  #Uniform prior        # Assumindo que temos uma range (Nstar) que vai de n (num identificado com sibship) até 350 , 
    #se tirarmos m samples dessa range, qual é a probabilidade de encontrarmos esse mesmo n?
    
    w <- prob(n, Nstar, m)/     # UID - gives the probability of finding n unique genotypes from a sample of m miracidia from Nstar
      dunif(Nstar, n, Nmax)     # (which is our correct number of female worms inside the host)
  }                             # Where is the assumption of egg laying by female worms???
  
  # 3. resample
  samp <- (sample(Nstar, size=nsim, prob=w, replace=T))
  
  
  df <- data.frame(expectN=mean(samp), bias=n-mean(samp),
                   varN = var(samp), percentbias = (n-mean(samp))/mean(samp),
                   mn = m/n,
                   lwr = quantile(samp, probs=c(0.025)), 
                   upr = quantile(samp, probs=c(0.975)))
  df
  
}



########################################################
## 2. explore posterior of N for different n, m and priors
########################################################


# Use uniform prior ----------------------------------------------------------------

df2 <- data.frame(n=S_h$n, m=S_h$m, mean=0, overdisp=0, prior=0)


tmp <- vector("list", nrow(df2))
for (i in 1:nrow(df2)) {
  if (df2[i,1]>df2[i,2]) {
    tmp[[i]] <- NA
  } else {
    tmp[[i]] <- f(as.numeric(df2[i,]))  #runs SIR algorithm
  }
}

df2 <- cbind(df2, do.call(rbind, tmp))


# Combining with Sh data using all data ------------------------------------------

combeggs<- cbind(S_h, df2[,-c(3:5)])


# add + 1 to eggs --------------------------------------

combeggs$eggs<- combeggs$eggs+1

