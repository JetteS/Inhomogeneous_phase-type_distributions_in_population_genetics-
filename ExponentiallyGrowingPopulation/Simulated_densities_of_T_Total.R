## Creating plots of the simulated density of T_Total
## under an exponentially growing population model

##################################################
## The inverse of the integrated intensity under 
## a model with an exponentially growing population
##################################################
## Name : LambdaInvEGP
## 
## Author: Jette Steinbach
## 
## Purpose : Compute the value of the inverse of
##           the integrated intensity 
##           Lambda(t) = int_0^t lambda(u) du.
##
## Input:
## - t : a vector holding the time points measured
##       in units of N generations at which 
##       the function should be evaluated.
## - gamma : a number between 0 and N (the population 
##           size in the current generation). This 
##           number is the important part of the 
##           factor 1/(1-gamma/N) by which the 
##           population size increases. Defaults to 0.
## Output: 
## - The value(s) of the the inverse of the integrated
##   intensity evaluated at the point(s) (in) t.

LambdaInvEGP <- function(t, gamma=0){
  ## Checking whether t is a positive real vector
  if(sum(t <= 0)>0) stop("All time points in t must be positive real numbers.")
  ## Checking whether gamma is a positive real number
  if(gamma < 0) stop("gamma must be positive.")
  
  ## Computing the inverse of the integrated intensity
  ## for all values of t
  if(gamma==0){
    
    res <- t
  }else{
    
    res <- 1/gamma*log(gamma*t +1)
  }
  return(res)
}

##################################################
## Plotting the simulated density of the total
## branch length under a model with an
## exponentially growing population size
##################################################
## Name : plot_simfTotal
## 
## Author: Jette Steinbach
## 
## Purpose : Plotting the simulated probability
##           density function of the total branch
##           length under a model with an
##           exponentially growing population size.
##
## Input:
## - n : a natural number representing the 
##       sample size. Must be greater 
##       than one. 
## - gamma : a number between 0 and N (the population 
##           size in the current generation). This 
##           number is the important part of the 
##           factor 1/(1-gamma/N) by which the 
##           population size increases. Defaults to 0.
## Output: 
## - A ggplot of the simulated probability density
##   function of the total branch length.
##
## Remark: Requires the package gglplot2

plot_simfTotal <- function(n, gamma=0){
  
  ## Simulating the waiting time T_Total once
  Dx <- data.frame(TTotal = simT(n=n, out = 'wtimes', LambdaInv = LambdaInvEGP, gamma=gamma)[2])
  ## and 10000 times more
  for(i in 1:10000){
    
    Dx <- rbind(Dx, simT(n=n, out = 'wtimes', LambdaInv = LambdaInvEGP, gamma=gamma)[2])
  }
  
  ## Creating a simple ggplot
  p <- ggplot(Dx, aes(x=TTotal)) + 
    ggtitle(expression(paste("Probability density function for ", tau["Total"])), 
            subtitle=TeX(paste("n =",n,"$\\gamma = ",gamma,"$")))+ 
    xlab("x") + 
    ylab(expression(paste(f[tau["Total"]], "(x)"))) +
    theme_minimal() + 
    geom_density(color = "red") +
    geom_histogram(aes(y=..density..), bins = 60, alpha=0.5) +
    ylim(0,1) +
    xlim(0,15)
  
  print(p)
}

## Plotting the simulated density for n=2,5,10,15 and gamma=0
plot_simfTotal(n=2, gamma = 0)
plot_simfTotal(n=5, gamma = 0)
plot_simfTotal(n=10, gamma = 0)
plot_simfTotal(n=15, gamma = 0)

## plotting the simulated density for n=2,10 and  gamma=0.1,0.5,1,2.
plot_simfTotal(n=2, gamma = 0.1)
plot_simfTotal(n=10, gamma = 0.1)
plot_simfTotal(n=2, gamma = 0.5)
plot_simfTotal(n=10, gamma = 0.5)
plot_simfTotal(n=2, gamma = 1)
plot_simfTotal(n=10, gamma = 1)
plot_simfTotal(n=2, gamma = 2)
plot_simfTotal(n=10, gamma = 2)
