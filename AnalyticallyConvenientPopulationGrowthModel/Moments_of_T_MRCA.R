## Creating plots of the mean and variance of T_MRCA
## under an analytically convenient population growth model

##################################################
## The mean and variance of the time to the most
## recent common ancestor under a model with an
## analytically convenient population growth.
##################################################
## Name : mvTMRCA
## 
## Author: Jette Steinbach
## 
## Purpose : Compute the mean and variance 
##          of the time to the most recent 
##          common ancestor under an analytically
##          convenient population growth model.
##
## Input:
## - n : a natural number representing the sample
##       size. Must be greater than one. 
##       Defaults to 10.
## - b : A strictly positive real value less than 1/2 
##       that enters the exponential function.
##
## Output: 
## - The mean and variance of the time to the most
##   recent common ancestor for a population of size 
##   n
##
## Remark: Requires the package expm.

mvTMRCA <- function(n=10, b){
  ## Checking whether n is a positive natural number:
  if(n <=1 | abs(n - round(n)) > .Machine$double.eps^0.5) stop("The sample size must be a natural number greater than 1.")
  ## Checking whether b is a real number less than 1/2
  if(b >= 0.5) stop("b must be less than 1/2.")
  ## Checking whether b is strictly positive
  if(b <= 0) stop("b must be strictly positive.")
  
  if(n==2){
    ## For a sample size of n=2, there can only be one 
    ## single coalescent event. Hence, all vectors and matrices 
    ## are just numbers. The initial distribution is equal to 
    ## one and the sub-intensity rate matrix is equal to minus
    ## one.
    
    m <- 1/(1-b)
    v <- 1/(1-2*b)- m^2
    
  }else{
    ## Defining the initial distribution
    initDist <- c(1,replicate(n=n-2,0))
    ## Defining the sub-intensity rate matrix T
    Tmat <- diag(choose(n:3,2),nrow = n-2, ncol = n-2)
    Tmat <- cbind(0,Tmat)
    Tmat <- rbind(Tmat,0)
    diag(Tmat) <- -choose(n:2,2)
    
    m <- initDist%*%solve(diag(x=b, nrow = n-1)+Tmat)%*%rowSums(Tmat)
    v <- initDist%*%solve(diag(x=2*b, nrow = n-1)+Tmat)%*%rowSums(Tmat) - m^2
  }
  res <- c(m,v)
  names(res) <- c("mean", "variance")
  return(res)
}

##################################################
## Plotting the mean and variance of the time to 
## the most recent common ancestor under a model 
## with an analytically convenient population 
## growth.
##################################################
## Name : plot_mvTMRCA
## 
## Author: Jette Steinbach
## 
## Purpose : Plot the mean and variance 
##           of the time to the most recent
##           common ancestor under an analytically
##           convenient population growth model.
##
## Input:
## - n : a vector of natural number representing the 
##       sample sizes. All entries must be greater 
##       than one. 
## - b : A vector of strictly positive real values 
##       less than 1/2 and that enter the
##       exponential function.
##
## Output: 
## - A ggplot of the mean and variance of 
##   the time to the most recent common ancestor.
##
## Remark: Requires the package gglplot2, ggpubr,
## RColorBrewer and tidyr

plot_mvTMRCA <- function(n, b){
  
  for(j in b){
    
    ## Creating the data frame that will hold the 
    ## results
    Dx <- data.frame(n)
    
    ## Computing the mean and variance of the time 
    ## to the most recent common ancestor for all 
    ## sample sizes 
    Dx <- cbind(Dx, t(sapply(n, mvTMRCA, b=j)))
    ## Changing the format of the data frame to long
    Dx <- pivot_longer(Dx,cols = 2:3)
    
    ## Creating a simple ggplot
    p <- ggplot(Dx, aes(x=n, y = value, colour = name)) + 
      ggtitle(expression(paste("The mean and variance for ", tau["MRCA"])), 
              subtitle=paste("b =",j)) + 
      xlab("n") + 
      ylab(" ") +
      theme_minimal() + 
      geom_line() +
      scale_colour_brewer(name = " ", 
                          labels = c(expression(paste("E[", tau["MRCA"], "]")), 
                                     expression(paste("V[", tau["MRCA"], "]"))),
                          palette = "Dark2")
    
    print(p)
  }
}

plot_mvTMRCA(n=2:20,b=c(0.001,0.1,0.2,0.45))
