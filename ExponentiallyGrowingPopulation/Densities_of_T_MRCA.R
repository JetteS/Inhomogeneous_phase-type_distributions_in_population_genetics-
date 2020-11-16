## Creating plots of the density of T_MRCA
## under an exponentially growing population model

##################################################
## Density of the time to the most recent common
## ancestor under a model with an exponentially
## growing population size
##################################################
## Name : fTMRCA
## 
## Author: Jette Steinbach
## 
## Purpose : Compute the value of the probability
##          density function of the time to the 
##          most recent common ancestor under a 
##          model with an exponentially growing 
##          population size.
##
## Input:
## - x : a vector of quantiles at 
##       which the probability density function 
##       should be evaluated. All quantiles must be 
##       positive real numbers.
## - n : a natural number representing the sample
##       size. Must be greater than one. 
##       Defaults to 10.
## - gamma : a number between 0 and N (the population 
##           size in the current generation). This 
##           number is the important part of the 
##           factor 1/(1-gamma/N) by which the 
##           population size increases. Defaults to 0.
##
## Output: 
## - The value(s) of the probability density function 
##   evaluated at the point(s) (in) x.
##
## Remark: Requires the package expm.

fTMRCA <- function(x, n=10, gamma=0){
  ## Checking whether n is a positive natural number:
  if(n <=1 | abs(n - round(n)) > .Machine$double.eps^0.5) stop("The sample size must be a natural number greater than 1.")
  ## Checking whether gamma is a positive real number
  if(gamma < 0) stop("gamma must be positive.")
  ## Checking whether x is a positive quantile
  if(sum(x <= 0)>0) stop("All quantiles in x must be positive real numbers.")
  
  ## Defining the vector holding the densities
  dens <- replicate(length(x),NA)
  
  if(n==2){
    ## For a sample size of n=2, there can only be one 
    ## single coalescent event. Hence, all vectors and matrices 
    ## are just numbers. The initial distribution is equal to 
    ## one and the sub-intensity rate matrix is equal to minus
    ## one.
    ## Compute the density at all points in x
    
    if(gamma==0){
      
      dens <- exp(-x)
      
    }else{
      
      for (i in 1:length(x)) {
        
        dens[i] <- exp(gamma*x[i]-1/gamma*(exp(gamma*x[i])-1))
      }
    }
    
  }else{
    ## Defining the initial distribution
    initDist <- c(1,replicate(n=n-2,0))
    ## Defining the sub-intensity rate matrix T
    Tmat <- diag(choose(n:3,2),nrow = n-2, ncol = n-2)
    Tmat <- cbind(0,Tmat)
    Tmat <- rbind(Tmat,0)
    diag(Tmat) <- -choose(n:2,2)
    
    ## Computing the value of the probability density 
    ## function for all values in x
    if(gamma==0){
      
      for (i in 1:length(x)) {
        
        dens[i] <- -initDist%*%expm(x[i]*Tmat)%*%rowSums(Tmat)
      }
      
    }else{
      
      for (i in 1:length(x)) {
        
        dens[i] <- -exp(gamma*x[i])*initDist%*%expm(1/gamma*(exp(gamma*x[i])-1)*Tmat)%*%rowSums(Tmat)
      }
    }
  }
  return(dens)
}

##################################################
## Plotting the density of the time to the most 
## recent common ancestor under a model with an 
## exponentially growing population size
##################################################
## Name : plot_fTMRCA
## 
## Author: Jette Steinbach
## 
## Purpose : Plot the probability
##          density function of the time to the 
##          most recent common ancestor under a 
##          model with an exponentially growing
##          population size.
##
## Input:
## - x : a vector of quantiles at 
##       which the probability density function 
##       should be evaluated. All quantiles must be 
##       positive real numbers.
## - n : a natural number representing the sample
##       size. Must be greater than one.
## - gamma : a number between 0 and N (the population 
##           size in the current generation). This 
##           number is the important part of the 
##           factor 1/(1-gamma/N) by which the 
##           population size increases. Defaults to 0.
##
## Output: 
## - A ggplot of the probability density function of 
##   the time to the most recent common ancestor.
##
## Remark: Requires the package gglplot2, ggpubr,
## RColorBrewer, tidyr and latex2exp

plot_fTMRCA <- function(x, n, gamma=0){
  
  for(j in gamma){
    
    ## Creating the data frame that will hold the 
    ## results
    Dx <- data.frame(x)
    
    ## Computing the probability density function for all
    ## values of x and all sample sizes
    Dx <- cbind(Dx, sapply(n, fTMRCA, x=x, gamma=j))
    ## Changing the format of the data frame to long
    Dx <- pivot_longer(Dx,cols = 2:(length(n)+1))
    
    ## Creating a simple ggplot
    p <- ggplot(Dx, aes(x=x, y = value, colour = name)) + 
      ggtitle(expression(paste("Probability density function for ", tau["MRCA"])), 
              subtitle=TeX(paste("$\\gamma = ",j,"$"))) + 
      xlab("x") + 
      ylab(expression(paste(f[tau["MRCA"]], "(x)"))) +
      theme_minimal() + 
      geom_line() +
      scale_colour_brewer(name = "Sample size", labels = as.character(n), palette = "Dark2") +
      ylim(0,2)
    
    ## Printing the ggplot
    print(p)
  }
}

## Creating the vector of quantiles
x <- seq(0.001,5, by = 0.001)
## and plotting the densities for gamma= 0, 0.1, 0.5, 1 and 2
plot_fTMRCA(x,n=c(2,5,10,15), gamma = c(0,0.1,0.5,1,2))
