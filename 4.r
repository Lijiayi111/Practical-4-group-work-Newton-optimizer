


hb <- function(th,k=2) {
  h <- matrix(0,2,2)
  h[1,1] <- 2-k*2*(2*(th[2]-th[1]ˆ2) - 4*th[1]ˆ2)
  h[2,2] <- 2*k
  h[1,2] <- h[2,1] <- -4*k*th[1]
  h
}


grad <- function(theta,t,y) {
  ## grad of -ve log lik of Poisson AIDS early epidemic model
  alpha <- theta[1];beta <- theta[2] ## enhances readability
  ebt <- exp(beta*t) ## avoid computing twice
  -c(sum(y)/alpha - sum(ebt), ## -dl/dalpha
     sum(y*t) - alpha*sum(t*ebt)) ## -dl/dbeta
} ## gll



#finite difference
Hessian <- function(theta, f, y, t=t80,grad,eps){
  eps <- 1e-8 ## finite difference interval
  grad0 <- grad(th0,t=t80,y=y) ## grad of nll at th0 
  Hfd <- matrix(0L, nrow = length(theta), ncol = length(theta)) ## initilize finite diference Hessian
  for (i in 1:length(th0)) { ## loop over parameters
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps
    grad1 <- grad(th1,t=t80,y=y) ## compute resulting nll
    Hfd[i,] <- (grad1 - grad0)/eps ## approximate second derivs  #i, haishi ,i
    
  }
  Hfd  <- 0.5 * (t(Hfd ) + Hfd)
  return(Hfd)
}





newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  object <- func(theta,...)
  gradient <- grad(theta,...)
  
  # Hessian calculation by finite differencing if no hessian as input.
  if(any(is.null(hess))) hessian <- hess(theta, f,...)
  else hessian <- hess(theta,...) # Calling fd.Hessian defined above
  
  
  # Checking quantities needed for optimization are finite
  if (!is.finite(object) || any(!is.finite(gradient))|| any(!is.finite(hessian)))
    stop("The objective or derivatives are not finite at the initial theta")## waring 1 
  

  for(i in 1:maxit){
    if(max(abs(gradient))< tol*(abs(object)+fscale)){
      cat("Convergence is done")
      
      if (inherits(try(chol(h),silent=TRUE),"try-error")) {
        inver_hessian = NULL
        warning("The Hessian is not positive definite at convergence")
      }else{
        
        # inverse of hessian matrix. 
        inver_hessian <- chol2inv(chol(hessian)) 
      }
      return(list(f=object, theta=theta, iter=i, g=gradient, Hi=inver_hessian)) #when converse return
      
      break
    }else{
      number_perturb <- 0
      # Cholesky factor calculation (perturb H until this calculation is possible)
      while(inherits(try(chol_number_perturb <- chol(number_perturb), silent = TRUE), "try-error") == TRUE) {
        # Perturbing hessian to be positive definite (runs of if chol(H) fails)
        H <- H + diag(norm(hessian)*10^number_perturb*1e-6, nrow=nrow(hessian), ncol=ncol(hessian))
        number_perturb <- number_perturb + 1
      }
      # Calculate forward step towards optimum (minimum)
      Delta <- backsolve(H.chol, forwardsolve(t(H.chol), -gradient))
      
      
      # Ensuring steps are taken in right direction towards optimum (minimum)
      # a.k.a preventing optimizer from blowing up 
      number_halving = 0
      while ((func(theta + Delta, ...)> object)|!is.finite(func(theta + Delta, ...)){
           if (number_halving < max.half) {
             Delta <- Delta / 2
             number_halving <- number_halving + 1
           } else {
             stop(paste("The update step failed to reduce the objective
              after ", as.character(max.half), " halvings"))
           }
      }
     theta <- theta + Delta
     
     object<-func(theta,...)
     gradient<-grad(theta,...)
     # Hessian calculation by finite differencing if no hessian as input.
     if(any(is.null(hess))) hessian <- hess(theta, f,...)
     else hessian <- hess(theta,...) # Calling fd.Hessian defined above
             
    }
    
    
    # Checking whether convergence occurs when iter == maxit
    if (max(abs(gradient)) < (abs(f0)+fscale)*tol){
      cat("Converged")
      return(list(f0, theta, iter, gradient, Hi))
    } else {
      warning(paste("Newton optimizer failed to converage after
                  maxit = ", as.character(maxit), " iterations"))
    }
  }
  
}



