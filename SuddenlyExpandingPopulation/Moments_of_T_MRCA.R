## Creating plots of the mean and variance of T_MRCA
## under a suddenly expanding population growth model

##################################################
## The mean and variance of the time to the most
## recent common ancestor under a model with a 
## suddenly expanding population size
##################################################
## Name : mvTMRCA
## 
## Author: Jette Steinbach
## 
## Purpose : Compute the mean and variance 
##          of the time to the 
##          most recent common ancestor under a 
##          model with a suddenly expanding 
##          population size.
##
## Input:
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
## - The mean and variance of the time to the most
##   recent common ancestor for a population of size 
##   n
##
## Remark: Requires the package expm.

mvTMRCA <- function(n=10, a=1, b=0.05){
  ## Checking whether n is a positive natural number:
  if(n <=1 | abs(n - round(n)) > .Machine$double.eps^0.5) stop("The sample size must be a natural number greater than 1.")
  ## Checking whether a is a real number between 0 and 1
  if(a <= 0) stop("The factor a must be strictly positive.")
  if(a > 1) stop("The factor a must be less than or equal to one.")
  ## Checking whether b is a positive real number
  if(b <= 0) stop("b must be a positive real number.")
  
  
  if(n==2){
    ## For a sample size of n=2, there can only be one 
    ## single coalescent event. Hence, all vectors and matrices 
    ## are just numbers. The initial distribution is equal to 
    ## one and the sub-intensity rate matrix is equal to minus
    ## one.
    
    if(a==1){
      
      m <- v <- 1
    }else{
      
      m <- 1-(1-a)*exp(-b)
      v <- (2-2*b*(1-a)*exp(-b) - 2*(1-a^2)*exp(-b))-m^2
    }
  }else{
    ## Defining the initial distribution
    initDist <- c(1,replicate(n=n-2,0))
    ## Defining the sub-intensity rate matrix T
    Tmat <- diag(choose(n:3,2),nrow = n-2, ncol = n-2)
    Tmat <- cbind(0,Tmat)
    Tmat <- rbind(Tmat,0)
    diag(Tmat) <- -choose(n:2,2)
    
    if(a==1){
      
      m <- -initDist%*%rowSums(solve(Tmat))
      v <- 2*initDist%*%rowSums(solve(Tmat%^%2)) - m^2
      
    }else{
      
      m <- -initDist%*%rowSums(solve(Tmat)) + (1-a)*initDist%*%solve(Tmat%^%2)%*%
        expm(Tmat*b)%*%rowSums(Tmat)
      v <- 2*initDist%*%rowSums(solve(Tmat%^%2)) + 
        initDist%*%(2*b*(1-a)*solve(Tmat%^%2) - 2*(1-a^2)*solve(Tmat%^%3))%*%
        expm(Tmat*b)%*%rowSums(Tmat) - m^2
    }
    
  }
  res <- c(m,v)
  names(res) <- c("mean", "variance")
  return(res)
}

##################################################
## Plotting the mean and variance of the time to 
## the most recent common ancestor under a model 
## with a suddenly expanding population size
##################################################
## Name : plot_mvTMRCA
## 
## Author: Jette Steinbach
## 
## Purpose : Plot the mean and variance 
##           of the time to the 
##           most recent common ancestor under a 
##           model with a suddenly expanding 
##           population size.

## Input:
## - n : a vector of natural number representing the 
##       sample sizes. All entries must be greater 
##       than one. 
## - a : a number between 0 and 1 that represents
##       the factor by which the population size
##       decreases. Defaults to 1.
## - b : a vector of real numbers that represents the
##       generation scaled in units of N 
##       generations where the expansions occurred.
## Output: 
## - A ggplot of the mean and variance of 
##   the time to the most recent common ancestor.
##
## Remark: Requires the package gglplot2, ggpubr,
## RColorBrewer and tidyr

plot_mvTMRCA <- function(n, a=1, b){
  
  for(j in b){
    
    ## Creating the data frame that will hold the 
    ## results
    Dx <- data.frame(n)
    
    ## Computing the mean and variance of the time 
    ## to the most recent common ancestor for all 
    ## sample sizes 
    Dx <- cbind(Dx, t(sapply(n, mvTMRCA, a=a,b=j)))
    ## Changing the format of the data frame to long
    Dx <- pivot_longer(Dx,cols = 2:3)
    
    ## Creating a simple ggplot
    p <- ggplot(Dx, aes(x=n, y = value, colour = name)) + 
      ggtitle(expression(paste("The mean and variance for ", tau["MRCA"])), 
              subtitle=paste("b =",j, "and a =", a)) + 
      xlab("n") + 
      ylab(" ") +
      ylim(0,2) +
      theme_minimal() + 
      geom_line() +
      scale_colour_brewer(name = " ", 
                          labels = c(expression(paste("E[", tau["MRCA"], "]")), 
                                     expression(paste("V[", tau["MRCA"], "]"))),
                          palette = "Dark2")
    
    print(p)
  }
}

plot_mvTMRCA(n=2:20,a=1,b=1) 
plot_mvTMRCA(n=2:20,a=0.2,b=c(1,4))
plot_mvTMRCA(n=2:20,a=0.5,b=c(1,4))
plot_mvTMRCA(n=2:20,a=0.8,b=c(1,4))
