# implements a simple mixture model for CpG island finding
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
simpleMixtureCGI <- function(ncg, p, L, maxit = 100) {
  
  p=p^2
  
  dampen <- function(a_in, p_in, tol=1e-6) pmax(tol, pmin(1-tol, p_in*a_in))
  loss_a <- function(a_in, p_in, ncg_in, y_in) {
    return(-sum(ncg_in*log(dampen(a_in, p_in)/(1-dampen(a_in,p_in))) + L*log(1-dampen(a_in,p_in))))
  }
  
  # initialize parameters
  a0 = 0.05
  a1 = 0.6
  pi = 0.2
  y = vector(mode="numeric", length=length(ncg))
  pre_y = y
  
  # repeat until convergence or maximum number of iterations is reached
  i=1
  while (i<=maxit) {
    print(sprintf(">>> %d / %d", i, maxit))
    # assign variables y
    y = as.numeric(pi*dbinom(ncg, L, p*a1)>(1-pi)*dbinom(ncg, L, p*a0))
    print(summary(y))
    
    if (all(y==pre_y)) break
    
    # update parameter estimates
    pi = mean(y)
    print(pi)
    a0 = optimize(loss_a, c(0, 1), p=p[y==0], ncg=ncg[y==0])$minimum
    print(a0)
    a1 = optimize(loss_a, c(0, 1), p=p[y==1], ncg=ncg[y==1])$minimum
    print(a1)
    
    i = i+1
    pre_y = y
  }
  # return final states y, and estimates pi, a1, a0 
  return(list(y=y, pi=pi, a1=a1, a0=a0))
}