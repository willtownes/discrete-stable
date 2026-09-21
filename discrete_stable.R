# functions to compute discrete stable PMFs and other handy convenience functions

sigma2gamma<-function(sigma,alpha){
  #change from mixed Poisson-Stable(Nolan-1) to simpler scale params
  ifelse(alpha==1, -sigma*2/pi, sigma^alpha/cos(alpha*pi/2) )
}

gamma2sigma<-function(gamma,alpha){
  ifelse(alpha==1, -gamma*pi/2, (cos(alpha*pi/2)*gamma)^(1/alpha) )
}

validgamma<-function(gamma,alpha){
  if(alpha<=0){ #alpha must be positive
    return(FALSE)
  }
  if(gamma==0){ #if gamma=0, this is ordinary Poisson, so alpha=1 by convention
    return(alpha==1)
  }
  #hereafter assume gamma!=0
  gpos<- gamma>0
  if(alpha<1){
    return(gpos) #gamma should be positive
  } else if(alpha<=2){
    return(!gpos) #gamma should be negative
  } else { #alpha>2
    alpha_star<-ceiling(alpha)
    if(alpha==alpha_star){ #integer values of alpha>2 not allowed
      return(FALSE)
    } else if(alpha_star%%2==1){
      return(gpos) #gamma positive when ceil(alpha) is odd
    } else {
      return(!gpos)
    }
  }
}

check_delta<-function(delta,gamma,alpha){
  if(delta< -gamma*alpha){
    stop("Must have delta>= -gamma*alpha to ensure P(X=1) is nonnegative")
  } else if(delta< -gamma){
    stop("Must have delta>= -gamma to ensure P(X=0)<=1")
  }
  #if delta>= -alpha*gamma, all probabilities are valid for alpha<=2
  #if alpha>2, even if both of these conditions are met, it's still possible
  #to have invalid probabilities. These are caught dynamically by the 
  #PMF function since there doesn't appear to be a closed from constraint.
  #increasing delta or making |gamma| closer to zero will always
  #make the PMF more well behaved.
}

dstable_pmf_naive<-function(nmax,delta,gamma,alpha){
  #compute the discrete stable PMF for non-integer alpha
  #nmax is the largest nonnegative integer at which to evaluate the PMF
  #will return P(X=n) for n=0,1,...,nmax
  #delta is location parameter
  #gamma is scale parameter
  #alpha is shape (tail) parameter
  #all parameters can only be scalars
  #evaluates the PMF using a direct recursion method
  #may not have great numerical stability
  stopifnot(validgamma(gamma,alpha))
  #direct way (numerical underflow likely)
  #assume alpha is not an integer, for alpha=1,2 handle separately
  #nmax<-100
  #nvals<-0:nmax
  jvals<-2:nmax
  check_delta(delta,gamma,alpha)
  dag<-delta+alpha*gamma
  # alpha_non_int<- (alpha != round(alpha))
  if(alpha==1){
    weights<-c(dag, -gamma/(jvals-1)) #1/k for k=1...n-1
  } else if(alpha!=2){
    weights<-c(dag, gamma*jvals*choose(alpha,jvals)*(-1)^(jvals-1))
  }
  #if alpha=2 there are only two weights so compute directly in the loop
  probs<-vector("numeric",length=nmax+1)
  probs[1]<-ifelse(alpha==1,exp(-delta),exp(-delta-gamma)) #P(X=0)
  if(nmax==0){ #edge case where only want P(X=0)
    return(probs[1])
  } 
  for(n in 1:nmax){
    if(alpha!=2){
      probs[n+1]<- mean(probs[n:1]*weights[1:n]) #P(X=n), mean does 1/n part
    } else if(alpha==2){
      probs[n+1]<-(dag*probs[n]-2*gamma*probs[n-1])/n
    } else {
      stop("invalid alpha, must be 1,2, or a positive non-integer value")
    }
    if(probs[n+1]<0){
      stop("negative probability encountered, increase delta or decrease absolute value of gamma")
    }
  }
  return(probs)
}

log1mexp<-function(x){
  #slightly modified from original VGAM implementation
  stopifnot(all(x[!is.na(x)]>0))
  ifelse(x<=log(2), log(-expm1(-x)), log1p(-exp(-x)))
}

signed_logSumExp<-function(lx,nn){
  #lx is log|x| where x could be positive or negative
  #nn is a logical vector indicating which elements of x are nonnegative
  #ie, nn<- sign(x)>=0
  #computes log(sum(x)) with numerical stability
  #we require that the sum of x be positive otherwise it throws an error
  #the sum is never directly computed though
  #lx is allowed to be -Inf (ie, x can have exact zeros)
  A<-matrixStats::logSumExp(lx,nn)
  B<-matrixStats::logSumExp(lx,!nn)
  A+log1mexp(A-B)
}

dstable_pmf<-function(nmax,delta,gamma,alpha,log=FALSE){
  #compute the discrete stable PMF for non-integer alpha
  #nmax is the largest nonnegative integer at which to evaluate the PMF
  #will return P(X=n) for n=0,1,...,nmax
  #delta is location parameter
  #gamma is scale parameter
  #alpha is shape (tail) parameter
  #all parameters can only be scalars
  stopifnot(validgamma(gamma,alpha))
  #assume alpha is not an integer, for alpha=1,2 handle separately
  jvals<-2:nmax
  check_delta(delta,gamma,alpha)
  dag<-delta+alpha*gamma
  ldag<-log(dag)
  # if(alpha==1){
  #   weights<-c(dag, -gamma/(jvals-1)) #1/k for k=1...n-1
  # } else if(alpha!=2){
  #   weights<-c(dag, gamma*jvals*choose(alpha,jvals)*(-1)^(jvals-1))
  # }
  if(alpha==1){
    lwts<-c(ldag, log(-gamma)-log(jvals-1))
  } else if(alpha==2){
    lwts<-c(ldag, log(-2*gamma)) #note this is shorter than the others!
  } else {
    lwts<-c(ldag, log(abs(gamma))+log(jvals)+lchoose(alpha,jvals))
  } 
  lp<-vector("numeric",length=nmax+1)
  lp[1]<- ifelse(alpha==1, -delta, -delta-gamma) #logP(X=0)
  if(nmax==0){ #edge case where only want P(X=0)
    return(ifelse(log,lp[1],exp(lp[1])))
  }
  nn<-rep.int(TRUE,nmax) #1,2,3,...,nmax. Which weights are nonnegative?
  #negative weights can only occur with alpha>2
  if(alpha>2){
    if(gamma>0){ #this implies ceil(alpha) is odd
      #negative terms are j=2,4,...floor(alpha)
      negterms<-seq(from=2,to=floor(alpha),by=2)
    } else { #case where ceil(alpha) is even
      #negative terms are j=3,5,...,floor(alpha)
      negterms<-seq(from=3,to=floor(alpha),by=2)
    }
    nn[negterms]<-FALSE
  }
  for(n in 1:nmax){
    if(alpha!=2){
      terms<- lp[n:1] + lwts[1:n]
    } else if(n==1) { #alpha=2, n=1 has only one weight
      terms<- lp[1]+ldag
    } else { #alpha=2, n=2,3,... has two weights
      terms<- lp[c(n,n-1)] + lwts
    }
    if(alpha<=2){ #all terms are nonnegative
      lp[n+1]<-matrixStats::logSumExp(terms) - log(n)
    } else { #alpha>2, some terms negative
      lp[n+1]<- signed_logSumExp(terms,nn[1:n]) - log(n)
    }
  }
  if(log){
    return(lp)
  } else {
    return(exp(lp))
  }
}

find_gamma_limit_gt2<-function(alpha, delta, n=max(delta*10,100), tol=1e-12, gamma_start = 1e-3, gamma_max = 1e6){
  stopifnot(alpha > 2, delta >= 0, tol > 0)
  if(alpha==round(alpha)){ return(0)}
  # sign convention from ceiling(alpha) parity
  gamma_sign <- if (ceiling(alpha) %% 2 == 1) 1 else -1
  is_valid <- function(g){
    tryCatch({
      dstable_pmf(n=n,delta=delta,gamma=g,alpha=alpha)
      return(TRUE)
    }, error = function(e) FALSE)
  }
  # anchor: |gamma| = 0 treated as valid lower bound
  lo <- 0
  hi <- gamma_start
  # expand hi geometrically until we find an invalid magnitude
  repeat {
    if (!is_valid(gamma_sign * hi)) break
    lo <- hi
    hi <- hi * 2
    if (hi > gamma_max) {
      stop("No invalid gamma found up to gamma_max = ", gamma_max,
           "; increase gamma_max or check pmf_fun.")
    }
  }
  # bisection: lo is valid, hi is invalid, shrink until width < epsilon
  while ((hi - lo) > tol) {
    mid <- (lo + hi) / 2
    if (is_valid(gamma_sign * mid)) {
      lo <- mid
    } else {
      hi <- mid
    }
  }
  return(gamma_sign*lo)
}

find_gamma_limit<-function(alpha,delta,...){
  #for a given value of tail index (alpha) and location (delta), return the gamma value with largest absolute value
  #such that the parameter combination is still valid but it is not valid for any gamma with larger abs val. 
  #based on numerical pmf evaluation of P(X=0),P(X=1),...,P(X=n)
  #... additional args passed to find_gamma_boundary_gt2
  stopifnot(alpha>0)
  if(alpha<1){
    #mean is not finite, gamma>0 is only lower bounded by -delta/alpha in the weird case of delta<0
    return(Inf)
  } 
  stopifnot(delta>=0)
  if(alpha<=2){ #gamma<=0 for alpha in [1,2]
    return(-alpha*delta)
  }
  find_gamma_limit_gt2(alpha,delta,...)
}
