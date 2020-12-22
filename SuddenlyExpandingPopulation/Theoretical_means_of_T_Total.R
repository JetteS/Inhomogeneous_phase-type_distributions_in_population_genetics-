## Computing the mean of the total branch length 
## T_Total for a sample of size n
## under a model with a suddenly expanding 
## population


##################################################
## Van Loan's method to evaluate specific mean 
## values
##################################################
## Name : VanLoan
## 
## Author: Jette Steinbach
## 
## Purpose : Evaluate mean values of the form
##           E[T_k 1(x_t=b)| x_0 = a]
##           using Van Loan's method.
##
## Input:
## - Q : A nxn-dimensional rate matrix, which satisfies
##       that all rows sum to zero and that
##       q_{ii} = - sum_{j \neq i} q_{ij}.
## - t : a positive real value representing the
##       coalescent time at which the process
##       is in state fstate
## - istate,fstate : natural numbers between one and n 
##          representing the initial and terminal
##          state, respectively. The initial state
##          defaults to one.
## - k : a natural number between one and n
##           that corresponds to the 
##           state for which the time T_k is 
##           measured.
##           
## Output: 
## - The mean value E[T_k 1(x_T=fstate)| x_0 = istate].
##
## Remark: Requires the package expm.

VanLoan <- function(Q,t,istate=1,fstate,k){
  
  ## Checking whether Q is a valid rate matrix
  if(sum(rowSums(Q)!= 0)>0) stop("The matrix Q is not a valid rate matrix. All rows must sum to zero")
  ## Checking whether t is a positive real number
  if(t <= 0) stop("The coalescent time t must be positive.")
  ## Checking whether istate is a natural number between one and n
  if(abs(istate - round(istate)) > .Machine$double.eps^0.5) stop("istate must be a natural number.")
  if(istate < 1 & istate > nrow(Q)) stop("istate must be a number between one and n.")
  ## Checking whether fstate is a natural number between one and n
  if(abs(fstate - round(fstate)) > .Machine$double.eps^0.5) stop("fstate must be a natural number.")
  if(fstate < 1 & fstate > nrow(Q)) stop("fstate must be a number between one and n.")
  ## Checking whether k is a natural number between one and n
  if(abs(k - round(k)) > .Machine$double.eps^0.5) stop("k must be a natural number.")
  if(k < 1 & k > nrow(Q)) stop("k must be a number between one and n.")
  
  ## Extracting the sample size
  n <- nrow(Q)
  ## Defining the matrix B, a matrix with the same dimensions as 
  ## Q, which entries are zero except for entry ((n+1)-k,(n+1)-k)
  ## (In this case entry (1,1) corresponds to state n).
  Bmat <- matrix(0,nrow = n, ncol = n)
  Bmat[(n+1)-k,(n+1)-k] <- 1
  
  ## Defining matrix A, that has a dimension which is 
  ## twice that of Q
  Amat <- cbind(Q,Bmat)
  Bmat[(n+1)-k,(n+1)-k] <- 0
  Amat <- rbind(Amat,cbind(Bmat,Q))
  
  ## Removing the matrix Bmat
  rm(Bmat)
  ## Taking the matrix exponential of A and returning entry
  ## (istate,fstate) of the matrix in the upper right corner
  Resmat <- expm(Amat*t)[1:n,(n+1):(2*n)]
  
  return(Resmat[(n+1)-istate, (n+1)-fstate])
}

##################################################
## The mean of the total branch length for a sample
## of size n under a model with a suddenly
## expanding population
##################################################
## Name : mTTotal
## 
## Author: Jette Steinbach
## 
## Purpose : Compute the mean of the total branch
##           length for a sample of size n under 
##           a model with a suddenly expanding 
##           population size.
##
## Input:
## - n : a natural number representing the sample
##       size. Must be greater than one. 
##       Defaults to 3.
## - a : a number between 0 and 1 that represents
##       the factor by which the population size
##       decreases. Defaults to 1.
## - b : The generation scaled in units of N 
##       generations where the expansion occurred.
##       Defaults to 1.
## Output: 
## - The mean of the total branch length
##
## Remark: Requires the package expm.

mTTotal <- function(n=3, a=1, b=1){
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
    ## one. Furthermore, T_Total = 2*T_2 = 2*T_MRCA, which means
    ## that we can compute the mean and variance of T_Total 
    ## using the mean and variance of T_MRCA.
    if(a==1){
      
      m <- 2
      
    }else{
      
      m <- 2*(1+(a-1)*exp(-b))
    }

  }else{ 
    
    if(a==1){
      
      ## For n>2 and a=1, we consider a constanr population
      ## size. In this situation, we can compute the mean
      ## of T_Total as the weighted sum of the means of 
      ## all holding times T_j ~ Exp(choose(j,2))
      m <- sum(2/(1:(n-1)))
      
    }else{
      
      ## For n>2 and a<1, we cannot use the mean and variance of 
      ## T_MRCA to compute the mean and variance of T_Total.
      ## Instead we have to use an approach that uses 
      ## Van Loan's method.
      ## In order to use this approach, we have to 
      ## define the rate matrix Q
      Q <- diag(choose(n:2,2),nrow = n-1, ncol = n-1)
      Q<- cbind(0,Q)
      Q <- rbind(Q,0)
      diag(Q) <- c(-choose(n:2,2),0)
    
      ## We initialize the means
      mVec <- replicate(n=n-1,NA)
    
      for(j in 2:n){
      
        FirstPart <- sum(sapply(1:n, VanLoan, Q=Q, t=b, istate=n, k=j))
        SecondPart <- sum(2*a/(j*(j-1))*expm(Q*b)[1,1:((n+1)-j)])
      
        mVec[j-1] <- (FirstPart+SecondPart)
      }
      ## The mean of T_Total can now be found as the weighted 
      ## sum of all entries in mVec
      m <- sum((2:n)*mVec)
    }
  }
  return(m)
}


##################################################
## Plotting the mean of the total branch length
## for different sample sizes under a model 
## with a suddenly expanding population size
##################################################
## Name : plot_mTTotal
## 
## Author: Jette Steinbach
## 
## Purpose : Plot the mean of the total branch 
##           length for different sample sizes
##           under a model with a suddenly
##           expanding population size.

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
## - A ggplot of the mean of the total branch length.
##
## Remark: Requires the packages gglplot2, latex2exp 
## and RColorBrewer

plot_mTTotal <- function(n, a=1, b){
  
  for(l in b){
    
    ## Creating the data frame that will hold the 
    ## results
    Dx <- data.frame(n)
    
    ## Computing the mean of the total branch length for all 
    ## sample sizes 
    Dx$means <- sapply(n,mTTotal,a=a,b=l)
    
    ## Creating a simple ggplot
    p <- ggplot(Dx, aes(x=n, y = means)) + 
      ggtitle(expression(paste("The mean for ", tau["Total"])), 
              subtitle=paste("b =",l, "and a =", a)) + 
      xlab("n") + 
      ylab(TeX("$E\\[\\tau_{Total}\\]$")) +
      ylim(0,7) +
      theme_minimal() + 
      geom_line(colour = brewer.pal(3,"Dark2")[1])

    
    print(p)
  }
}


plot_mTTotal(n=2:20,a=1,b=1)
plot_mTTotal(n=2:20,a=0.2,b=1)
plot_mTTotal(n=2:20,a=0.2,b=3)
plot_mTTotal(n=2:20,a=0.5,b=1)
plot_mTTotal(n=2:20,a=0.5,b=3)
plot_mTTotal(n=2:20,a=0.8,b=1)
plot_mTTotal(n=2:20,a=0.8,b=3)
