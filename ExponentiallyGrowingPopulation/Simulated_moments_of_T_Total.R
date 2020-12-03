## Creating plots of the simulated mean and variance of T_Total
## under an exponentially growing population model

##################################################
## Plotting the simulated mean and variance of the
## total branch length under a model with an 
## exponentially growing population size
##################################################
## Name : plot_mvTTotal
## 
## Author: Jette Steinbach
## 
## Purpose : Plotting the simulated mean and variance 
##           of the total branch length under a model
##           with an exponentially growing population 
##           size.
##
## Input:
## - n : a vector of natural numbers representing the 
##       sample sizes. All entries must be greater 
##       than one. 
## - gamma : a number between 0 and N (the population 
##           size in the current generation). This 
##           number is the important part of the 
##           factor 1/(1-gamma/N) by which the 
##           population size increases. Defaults to 0.
##
## Output: 
## - A ggplot of the simulated mean and variance of 
##   the total branch length.
##
## Remark: Requires the package gglplot2,
## tidyr and Rfast

plot_mvTTotal <- function(n, gamma=0){
  
  ## Simulating the waiting times once for all sample sizes 
  Dx <- sapply(n, simT, out= 'wtimes', LambdaInv = LambdaInvEGP, gamma=gamma)[2,]
  
  ## Simulating the waiting times 10000 times more
  for(i in 1:10000){
    
    Dx <- rbind(Dx,sapply(n, simT, out= 'wtimes', LambdaInv = LambdaInvEGP, gamma=gamma)[2,])
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
            subtitle=TeX(paste("$\\gamma = ",gamma,"$"))) +
    xlab("n") + 
    ylab(" ") +
    ylim(0,8) +
    theme_minimal() + 
    geom_line() +
    scale_colour_brewer(name = " ", 
                        labels = c(expression(paste("E[", tau["Total"], "]")), 
                                   expression(paste("V[", tau["Total"], "]"))),
                        palette = "Dark2")
  
  print(p)
}

plot_mvTTotal(n=2:20, gamma=0)
plot_mvTTotal(n=2:20, gamma=0.1)
plot_mvTTotal(n=2:20, gamma=0.5)
plot_mvTTotal(n=2:20, gamma=1)
plot_mvTTotal(n=2:20, gamma=2)

