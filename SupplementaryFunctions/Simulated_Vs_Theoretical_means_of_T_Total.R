##################################################
## Plotting the simulated and theoretical 
## mean of the total branch length
## for different sample sizes under a model 
## with a suddenly expanding population size
##################################################
## Name : plot_meanTTotal
## 
## Author: Jette Steinbach
## 
## Purpose : Plot the simulated and theoretical 
##           mean of the total branch 
##           length for different sample sizes
##           under a model with a suddenly
##           expanding population size.
##
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
## - A ggplot of the theoretical and simulated 
##   mean of the total branch length.
##
## Remark: Requires the packages gglplot2, latex2exp 
## and RColorBrewer

plot_meanTTotal <- function(n, a=1, b){
  
  for(l in b){
    
    ## Creating the data frame that will hold the 
    ## results
    Dx <- data.frame(n)
    
    ## Computing the mean of the total branch length for all 
    ## sample sizes 
    Dx$means <- sapply(n,mTTotal,a=a,b=l)
    
    ## Simulating the means as well
    mTtotal <- replicate(length(n),NA)
    
    for (m in n) {
      
      ## Simulating the waiting times once
      Ttotal <- simT(n=m, out= 'wtimes', LambdaInv = LambdaInvSEP, a=a,b=b)[2]
      
      ## Simulating the holding times 10000 times more
      for(i in 1:10000){
        
        Ttotal <- c(Ttotal,simT(n=m, out= 'wtimes', LambdaInv = LambdaInvSEP, a=a,b=b)[2])
      }
      
      mTtotal[which(n==m)] <- mean(Ttotal)
    }
    
    Dx$SimMeans <- mTtotal
    
    ## Creating a simple ggplot
    p <- ggplot(Dx, aes(x=n)) + 
      ggtitle(expression(paste("The mean for ", tau["Total"])), 
              subtitle=paste("b =",l, "and a =", a)) + 
      xlab("n") + 
      ylab(TeX("$E\\[\\tau_{Total}\\]$")) +
      ylim(0,7) +
      theme_minimal() + 
      geom_line(aes(y=means), colour = brewer.pal(3,"Dark2")[1]) +
      geom_line((aes(y=SimMeans)), colour = "red")
    
    
    print(p)
  }
}

plot_meanTTotal(n=2:20,a=1,b=1)
plot_meanTTotal(n=2:20,a=0.2,b=1)
plot_meanTTotal(n=2:20,a=0.2,b=3)
plot_meanTTotal(n=2:20,a=0.5,b=1)
plot_meanTTotal(n=2:20,a=0.5,b=3)
plot_meanTTotal(n=2:20,a=0.8,b=1)
plot_meanTTotal(n=2:20,a=0.8,b=3)