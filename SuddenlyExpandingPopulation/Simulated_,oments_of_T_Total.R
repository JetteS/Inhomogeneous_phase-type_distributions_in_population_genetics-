## Creating plots of the simulated mean and variance of T_Total
## under a suddenly expanding population growth model

##################################################
## Plotting the simulated mean and variance of the
## total branch length under a model with a
## suddenly expanding population size
##################################################
## Name : plot_mvTTotal
## 
## Author: Jette Steinbach
## 
## Purpose : Plotting the simulated mean and variance 
##           of the total branch length under  
##           a model with a suddenly expanding 
##           population size.
##
## Input:
## - n : a vector of natural numbers representing the 
##       sample sizes. All entries must be greater 
##       than one. 
## - a : a number between 0 and 1 that represents
##       the factor by which the population size
##       decreases. Defaults to 1.
## - b : The generation scaled in units of N 
##       generations where the expansion occured.
##
## Output: 
## - A ggplot of the mean and variance of 
##   the total branch length.
##
## Remark: Requires the package gglplot2,
## tidyr and Rfast

plot_mvTTotal <- function(n, a=1, b){
  
  ## Simulating the waiting times once for all sample sizes 
  Dx <- sapply(n, simT, out= 'wtimes', LambdaInv = LambdaInvSEP, a=a,b=b)[2,]
  
  ## Simulating the waiting times 10000 times more
  for(i in 1:10000){
    
    Dx <- rbind(Dx,sapply(n, simT, out= 'wtimes', LambdaInv = LambdaInvSEP, a=a,b=b)[2,])
  }
  
  ## Computing the means and variances of the columns
  m <- colMeans(Dx)
  v <- colVars(Dx)
  
  ## Creating a data frame 
  Dx <- data.frame(n=n,m=m,v=v)
  ## and changing the format of the data frame to long
  Dx <- pivot_longer(Dx,cols = 2:3)
  
  ## Creating a simple ggplot
  p <- ggplot(Dx, aes(x=n, y = value, colour = name)) + 
    ggtitle(expression(paste("The mean and variance for ", tau["Total"])), 
            subtitle=paste("b =",b, "and a =", a)) + 
    xlab("n") + 
    ylab(" ") +
    ylim(0,7) +
    theme_minimal() + 
    geom_line() +
    scale_colour_brewer(name = " ", 
                        labels = c(expression(paste("E[", tau["Total"], "]")), 
                                   expression(paste("V[", tau["Total"], "]"))),
                        palette = "Dark2")
  
  print(p)
}

plot_mvTTotal(n=2:20,a=1, b=1)
plot_mvTTotal(n=2:20,a=0.2,b=1)
plot_mvTTotal(n=2:20,a=0.2,b=3)
plot_mvTTotal(n=2:20,a=0.8,b=1)
plot_mvTTotal(n=2:20,a=0.8,b=3)
