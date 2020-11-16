## Creating plots of the density of T_MRCA under a 
## suddenly expanding population growth model

##################################################
## Density of the time to the most recent common
## ancestor under a model with a suddenly
## expanding population size
##################################################
## Name : fTMRCA
## 
## Author: Jette Steinbach
## 
## Purpose : Compute the value of the probability
##          density function of the time to the 
##          most recent common ancestor under a 
##          model with a suddenly expanding 
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
## - a : a number between 0 and 1 that represents
##       the factor by which the population size
##       decreases. Defaults to 1.
## - b : The generation scaled in units of N 
##       generations where the expansion occurred.
##       Defaults to 0.05.
## Output: 
## - The value(s) of the probability density function 
##   evaluated at the point(s) (in) x.
##
## Remark: Requires the package expm.

fTMRCA <- function(x, n=10, a=1, b=1){
  ## Checking whether n is a positive natural number:
  if(n <=1 | abs(n - round(n)) > .Machine$double.eps^0.5) stop("The sample size must be a natural number greater than 1.")
  ## Checking whether a is a real number between 0 and 1
  if(a <= 0) stop("The factor a must be strictly positive.")
  if(a > 1) stop("The factor a must be less than or equal to one.")
  ## Checking whether b is a positive real number
  if(b <= 0) stop("b must be a positive real number.")
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
    for (i in 1:length(x)) {
      
      if(x[i] <= b){
        
        dens[i] <- exp(-x[i])
      }else{
        dens[i] <- (1/a)*exp(-(a*b-b+x[i])/a)
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
    for (i in 1:length(x)) {
      
      if(x[i] <= b){
        
        dens[i] <- initDist%*%expm(Tmat*x[i])%*%(-rowSums(Tmat))
      }else{
        
        dens[i] <- (1/a)*initDist%*%expm(Tmat*((a*b-b+x[i])/a))%*%(-rowSums(Tmat))
      }
    }
  }
  return(dens)
}

##################################################
## Plotting the density of the time to the most 
## recent common ancestor under a model with a 
## suddenly expanding population size
##################################################
## Name : plot_fTMRCA
## 
## Author: Jette Steinbach
## 
## Purpose : Plot the probability density
##           function of the time to the 
##           most recent common ancestor under a 
##           model with a suddenly expanding 
##           population size.

## Input:
## - x : a vector of quantiles at 
##       which the probability density function 
##       should be evaluated. All quantiles must be 
##       positive real numbers.
## - n : a vector of natural number representing the 
##       sample sizes. All entries must be greater 
##       than one. 
## - a : a number between 0 and 1 that represents
##       the factor by which the population size
##       decreases. Defaults to 1.
## - b : A vector of real numbers that represents the
##       generation scaled in units of N 
##       generations where the expansions occurred.
## Output: 
## - A ggplot of the probability density function of 
##   the time to the most recent common ancestor.
##
## Remark: Requires the package gglplot2, ggpubr,
## RColorBrewer and tidyr

plot_fTMRCA <- function(x, n, a=1, b){
  
  for(j in b){
    
    ## Creating the data frame that will hold the 
    ## results
    Dx <- data.frame(x)
    
    ## Computing the probability density function for all
    ## values of x and all sample sizes
    Dx <- cbind(Dx, sapply(n, fTMRCA, x=x, a=a,b=j))
    ## Changing the format of the data frame to long
    Dx <- pivot_longer(Dx,cols = 2:(length(n)+1))
    
    ## Creating a simple ggplot
    p <- ggplot(Dx, aes(x=x, y = value, colour = name)) + 
      ggtitle(expression(paste("Probability density function for ", tau["MRCA"])), 
              subtitle=paste("b =",j, "and a =", a)) + 
      xlab("x") + 
      ylab(expression(paste(f[tau["MRCA"]], "(x)"))) +
      theme_minimal() + 
      geom_line() +
      scale_colour_brewer(name = "Sample size", labels = as.character(n), palette = "Dark2") +
      ylim(0,1)
    
    ## Printing the ggplot 
    print(p)
  }
}

## Creating the vector of quantiles
x <- seq(0.001,5, by = 0.001)
## and plotting the densities for a=0.2, 0.5, 0.8 and 1
## and b=1 and 4.
plot_fTMRCA(x,n=c(2,5,10,15),a=1,b=c(1,4)) 
plot_fTMRCA(x,n=c(2,5,10,15),a=0.2,b=c(1,4))
plot_fTMRCA(x,n=c(2,5,10,15),a=0.5,b=c(1,4))
plot_fTMRCA(x,n=c(2,5,10,15),a=0.8,b=c(1,4))
