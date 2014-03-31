# implements an HMM model for CpG island finding
# arguments
#   ncg: vector of CG counts in segments (N_CG(s))
#   p: CG proportion in segments (p(s))
#   L: segment lengths
#   maxit: maximum number of iterations
#
# returns a list with following components
#   y: state assignements for each segment
#   pi: estimated state probability 
#   a1: O/E ratio for "high" CpG rate regions
#   a0: O/E ratio for "baseline" CpG rate regions
hmmCGI <- function(ncg, p, L, maxit = 100) {
  
  dampen <- function(a_in, p_in, tol=1e-6) pmax(tol, pmin(1-tol, p_in*a_in))
  loss_a <- function(a_in, p_in, ncg_in) {
    return(-sum(ncg_in*log(dampen(a_in, p_in)/(1-dampen(a_in,p_in))) + L*log(1-dampen(a_in,p_in))))
  }
  
  # initialize parameters
  len = length(ncg)
  
  a0 = 0.05
  a1 = 0.6
  pi = 0.2
  piMatrix = array(0.5, dim=c(2,2))
  y = as.numeric(pi*dbinom(ncg, L, p*a1)>(1-pi)*dbinom(ncg, L, p*a0))
  d = array(0, dim=c(2,len))
  d[1,1] = log(1-pi) + log(dbinom(ncg[1],L,p[1]*a0))
  d[2,1] = log(pi) + log(dbinom(ncg[1],L,p[1]*a1))
  path = array(0, dim=c(2,len-1))
  pre_y = y
  
  # repeat until convergence or maximum number of iterations is reached
  i = 1
  print(sprintf("=========== DEBUG ============"))  
  while (i<=maxit) {
    print(sprintf(">>> %d / %d", i, maxit))
    # update parameter estimates
    for (j in 1:(len-1)) {
      d_0from0 = log(piMatrix[1,1]) + log(dbinom(ncg[j+1], L, p[j+1]*a0)) + d[1,j]
      d_0from1 = log(piMatrix[1,2]) + log(dbinom(ncg[j+1], L, p[j+1]*a0)) + d[2,j]
      d[1,j+1] = max(d_0from0,d_0from1)
      path[1,j] = as.numeric(d_0from1>d_0from0)

      d_1from0 = log(piMatrix[2,1]) + log(dbinom(ncg[j+1], L, p[j+1]*a1)) + d[1,j]
      d_1from1 = log(piMatrix[2,2]) + log(dbinom(ncg[j+1], L, p[j+1]*a1)) + d[2,j]
      d[2,j+1] = max(d_1from0,d_1from1)
      path[2,j] = as.numeric(d_1from1>d_1from0)      
    }

    # assign variables y
    y[len] = as.numeric(d[2,len]>d[1,len])
    for (j in (len-1):1) {
      y[j] = path[y[j+1]+1,j]
    }
    print("Y:")
    print(summary(y))
    
    pi = mean(y)
    print("pi:")    
    print(pi)
    a0 = optimize(loss_a, c(0, 1), p=p[y==0], ncg=ncg[y==0])$minimum
    print("a0:") 
    print(a0)
    a1 = optimize(loss_a, c(0, 1), p=p[y==1], ncg=ncg[y==1])$minimum
    print("a1:") 
    print(a1)
    
    piMatrix = piMatrix*0
    for (j in 1:(len-1)) {
      if(y[j+1]==0 & y[j]==0) piMatrix[1,1] = piMatrix[1,1] + 1
      if(y[j+1]==0 & y[j]==1) piMatrix[1,2] = piMatrix[1,2] + 1
      if(y[j+1]==1 & y[j]==0) piMatrix[2,1] = piMatrix[2,1] + 1
      if(y[j+1]==1 & y[j]==1) piMatrix[2,2] = piMatrix[2,2] + 1
    }
    piMatrix = piMatrix/(len-1)
    print("pi matrix:")
    print(piMatrix)
    
    if (all(y==pre_y)) break
    
    i = i+1
    pre_y = y
  }
  # return final states y, and estimates pi, a1, a0 
  return(list(y=y, pi=pi, a1=a1, a0=a0))     
}