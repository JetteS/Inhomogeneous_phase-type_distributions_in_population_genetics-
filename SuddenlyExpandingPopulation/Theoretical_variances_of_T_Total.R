## Computing the variance of the total branch length 
## T_Total for a sample of size n
## under a model with a suddenly expanding 
## population


##################################################
## Van Loan's method to evaluate specific mean 
## values
##################################################
## Name : VanLoanExtended
## 
## Author: Jette Steinbach
## 
## Purpose : Evaluate mean values of the form
##           E[T_i T_j 1(x_t=b)| x_0 = n]
##           using Van Loan's method.
##
## Input:
## - Q : A nxn-dimensional rate matrix, which satisfies
##       that all rows sum to zero and that
##       q_{ii} = - sum_{j \neq i} q_{ij}.
## - t : a positive real value representing the
##       coalescent time at which the process
##       is in state fstate
## - fstate : a natural numbers between one and n 
##          representing the terminal state.
## - i,j : natural numbers between one and n
##         that correspond to the states for
##         which the times T_i and T_j are 
##         measured.
##           
## Output: 
## - The mean value E[T_i T_j 1(x_t=fstate)| x_0 = n].
##
## Remark: Requires the package expm.

VanLoanExtended <- function(Q,t,fstate,i,j){
  
  ## Checking whether Q is a valid rate matrix
  if(sum(rowSums(Q)!= 0)>0) stop("The matrix Q is not a valid rate matrix. All rows must sum to zero")
  ## Checking whether t is a positive real number
  if(t <= 0) stop("The coalescent time t must be positive.")
  ## Checking whether fstate is a natural number between one and n
  if(abs(fstate - round(fstate)) > .Machine$double.eps^0.5) stop("fstate must be a natural number.")
  if(fstate < 1 & fstate > nrow(Q)) stop("fstate must be a number between one and n.")
  ## Checking whether i is a natural number between one and n
  if(abs(i - round(i)) > .Machine$double.eps^0.5) stop("i must be a natural number.")
  if(i < 1 & i > nrow(Q)) stop("i must be a number between one and n.")
  ## Checking whether j is a natural number between one and n
  if(abs(j - round(j)) > .Machine$double.eps^0.5) stop("j must be a natural number.")
  if(j < 1 & j > nrow(Q)) stop("j must be a number between one and n.")
  
  ## Extracting the sample size
  n <- nrow(Q)
  
  ## Defining the matrices B1 and B2, two matrices of
  ## the same dimensions as Q which entries are zero except
  ## for entry ((n+1)-j,(n+1)-j) and ((n+1)-k,(n+1)-k), 
  ## respectively.
  ## (In this case entry (1,1) corresponds to state n).
  B1mat <- matrix(0,nrow = n, ncol = n)
  B1mat[(n+1)-fstate,(n+1)-j] <- 1
  
  B2mat <- matrix(0,nrow = n, ncol = n)
  B2mat[(n+1)-i,1] <- 1
  
  ## Defining matrix A, which is composed of Q,B1 and B2.
  ## Hence, A has a dimension which is three times that of Q.
  Amat <- cbind(Q,B1mat)
  B1mat[(n+1)-fstate,(n+1)-j] <- 0
  Amat <- cbind(Amat,B1mat)
  
  Asub <- cbind(cbind(B1mat,Q),B2mat)
  
  Amat <- rbind(Amat,Asub)
  
  Asub <- cbind(cbind(B1mat,B1mat),Q)
  
  Amat <- rbind(Amat,Asub)
  
  ## Removing the matrices Asub and Bmat
  rm(Asub,B1mat,B2mat)
  
  ## Taking the matrix exponential of A and returning entry
  ## (i,j) of the matrix in the upper right corner
  Resmat <- expm(Amat*t)[1:n,(2*n+1):(3*n)]
  res <- Resmat[(n+1)-i, (n+1)-j]
  ## Removing the matrix Amat
  rm(Amat,Resmat)
  
  ## And returning the mean value
  return(res)
}


##################################################
## The variance of the total branch length for a 
## sample of size n under a model with a suddenly
## expanding population
##################################################
## Name : vTTotal
## 
## Author: Jette Steinbach
## 
## Purpose : Compute the variance of the total branch
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
## - The variance of the total branch length
##
## Remark: Requires the package expm.

vTTotal <- function(n=3, a=1, b=1){
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
      
      v <- 4
    }else{
      
      m <- 2*(1+(a-1)*exp(-b))
      v <- 4*(2-(2*b*(1-a) + 2*(1-a^2))*exp(-b)-(m/2)^2)
    }
    
  }else{ 
    
    if(a==1){
      
      ## For n>2 and a=1, we consider a constanr population
      ## size. In this situation, we can compute the variance
      ## of T_Total as the weighted sum of the variances of 
      ## all holding times T_j ~ Exp(choose(j,2))
      
      v <- sum(4/((1:(n-1))^2))
      
    }else{
      
      ## For n>2 and a<1, we cannot use the mean and variance of 
      ## T_MRCA to compute the mean and variance of T_Total.
      ## Instead we have to use an approach that uses 
      ## Van Loan's method.
      ## In order to use this approach, we have to 
      ## define the full rate matrix Q
      Q <- diag(choose(n:2,2),nrow = n-1, ncol = n-1)
      Q<- cbind(0,Q)
      Q <- rbind(Q,0)
      diag(Q) <- c(-choose(n:2,2),0)
      
      ## In order to be able to compute the covariances, 
      ## we need to calculate the means of all holding 
      ## times T_i
      mVec <- replicate(n=n-1,NA)
      
      for(j in 2:n){
        
        FirstPart <- sum(sapply(1:n, VanLoan, Q=Q, t=b, istate=n, k=j))
        SecondPart <- sum(2*a/(j*(j-1))*expm(Q*b)[1,1:((n+1)-j)])
        
        mVec[j-1] <- (FirstPart+SecondPart)
      }
      
      rm(FirstPart,SecondPart,j)
      
      ## Now we initiatlize the matrix that will hold the 
      ## means of T_iT_j for all i,j=2,...,n
      mMat <- matrix(NA,nrow = n-1, ncol = n-1)
      
      for(i in 2:n){
        for(j in 2:n){
          
          FirstPart <- sum(sapply(1:n, VanLoanExtended, Q=Q, t=b, i=i,j=j))+
                       sum(sapply(1:n, VanLoanExtended, Q=Q, t=b, i=j,j=i))

          SecondPart <- sum(sapply(j:n, VanLoan, Q=Q, t=b, istate=n, k=i)*
                               2*a/(j*(j-1)))
          
          ThirdPart <- sum(sapply(i:n, VanLoan, Q=Q, t=b, istate=n, k=j)*
                               2*a/(i*(i-1)))
            
          FourthPart <- ifelse(i!=j, 4*a^2/(i*(i-1))*1/(j*(j-1))*sum(expm(Q*b)[1,1:(n-max(i,j))]),
                               8*a^2/(i^2*(i-1)^2)*sum(expm(Q*b)[1,1:(n-i)]))
            
          
          mMat[i-1,j-1] <- FirstPart+SecondPart+ThirdPart+FourthPart
        }
      }
      
      ## Now that we have computed all means, we are able to compute the 
      ## covariances
      covMat <- matrix(NA,nrow = n-1, ncol = n-1)
      
      for (i in 2:n) {
        for (j in 2:n) {
          
          covMat[i-1,j-1] <- mMat[i-1,j-1] - mVec[i-1]*mVec[j-1]
        }
      }
      
      ## Finally, we are able to compute the variance of
      ## the total branch length
      
      v <- 0
      
      for (i in 2:n){
        for (j in 2:n) {
          
          v <- v + i*j*covMat[i-1,j-1]
        }
      }
    }
  }
  return(v)
}


##################################################
## Plotting the variance of the total branch length
## for different sample sizes under a model 
## with a suddenly expanding population size
##################################################
## Name : plot_vTTotal
## 
## Author: Jette Steinbach
## 
## Purpose : Plot the variance of the total branch 
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
## - A ggplot of the variance of the total branch length.
##
## Remark: Requires the packages gglplot2, latex2exp 
## and RColorBrewer

plot_vTTotal <- function(n, a=1, b){
  
  for(l in b){
    
    ## Creating the data frame that will hold the 
    ## results
    Dx <- data.frame(n)
    
    ## Computing the variance of the total branch
    ## length for all sample sizes 
    Dx$variances <- sapply(n,vTTotal,a=a,b=l)
    
    ## Creating a simple ggplot
    p <- ggplot(Dx, aes(x=n, y = variances)) + 
      ggtitle(expression(paste("The variance for ", tau["Total"])), 
              subtitle=paste("b =",l, "and a =", a)) + 
      xlab("n") + 
      ylab(TeX("$V\\[\\tau_{Total}\\]$")) +
      theme_minimal() + 
      ylim(0,7)+
      geom_line(colour = brewer.pal(3,"Dark2")[2])
    
    
    print(p)
    ggsave(paste0("TheoreticalVarianceTTotala",a,"b",l,".pdf"), plot = p, device = "pdf",
           path = "D:/UNI/SpecialeE20/Graphics")
  }
}

plot_vTTotal(n=2:20,a=1,b=1)
plot_vTTotal(n=2:19,a=0.2,b=1)
plot_vTTotal(n=2:19,a=0.2,b=3)
plot_vTTotal(n=2:19,a=0.5,b=1)
plot_vTTotal(n=2:19,a=0.5,b=3)
plot_vTTotal(n=2:19,a=0.8,b=1)
plot_vTTotal(n=2:19,a=0.8,b=3)
