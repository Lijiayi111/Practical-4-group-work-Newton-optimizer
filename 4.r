## Jiayi Li, S2304671 :  commenting for the code, coding realization and correcting. 
## Rongkai Fang, s2310813 : commenting for the code, coding designing and guidance.
## Yifan Hu, s2287683 ： commenting for the code, coding realization and correcting. 
## every one contribute roughly equal on this project. 

## address of github: https://github.com/Lijiayi111/Practical-4-group-work-Newton-optimizer.git



## In this project, we write an R function, newt, implementing Newton’s method for
## minimization of functions. 


## finite difference of the gradient vector when hessian function is not given 
Hessian <- function(theta,grad,eps=1e-6,...){ 
## this function will return our approximate Hessian matrix and has the same 
## arguments as func. 
  grad0 <- grad(theta,...) # gradient vector at initial theta. 
  ## initialize the approximate Hessian
  Hfd <- matrix(0L, nrow = length(theta), ncol = length(theta)) 
  for (i in 1:length(theta)) { # loop over parameters
    th1 <- theta; th1[i] <- th1[i] + eps ## increase th0[i]/initial theta by eps
    grad1 <- grad(th1,...) ## compute resulting gradient vector
    Hfd[i,] <- (grad1 - grad0)/eps ## approximate second derivatives  
  }
  Hfd  <- 0.5 * (t(Hfd ) + Hfd) # generate a symmetric matrix
  return(Hfd) # return approximated hessian matrix. 
}



newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
## This newt function is to find the minimization of func using Newton's method
##Output: should return a list containing:
##      f: the minimum value of the objective function.
##      theta: the value of the parameters at the minimum.
##      iter: the number of iterations taken to reach the minimum.
##      g: the gradient vector at the minimum.
##      Hi: the inverse of the Hessian matrix at the minimum
##Input: theta is a vector of initial values for the optimization parameters.
##      func is the objective function to minimize.grad is the gradient function,
##      which has the same arguments as func but returns the gradient vector of 
##      the objective. hess is the Hessian matrix function, which has the same 
##      arguments as func but returns the Hessian matrix. 
##      tol is the convergence tolerance. fscale is  a rough estimate of the 
##      magnitude of func near the optimum - used in convergence testing.
##      maxit is the maximum number of Newton iterations to try before giving up.
##      max.half is the maximum number of times a step should be halved before 
##      concluding that the step has failed to improve the objective.
##      eps is the finite difference intervals to use when a Hessian function is not provided.
  
  
  object <- func(theta,...) # get objective value under initial theta.
  gradient <- grad(theta,...)# get the gradient vector under initial theta. 
  
  # if hess is not supplied as input
  if(any(is.null(hess))) hessian <- Hessian(theta,grad,...)  # calculation by finite differencing 
  # if hess function is defined, just call it 
  else hessian <- hess(theta,...) 
  
  
  # Checking objective or gradient or hessian is finite or not. 
  if (!is.finite(object) || any(!is.finite(gradient))|| any(!is.finite(hessian)))
    stop("The objective or gradient or hessian is not finite under the initial defined theta") ## stop and print warning 
  
  # loop over maxit of time to generate new theta in improving direction and find the minimum
  for(i in 1:maxit){
    # check convergence
    if(max(abs(gradient)) < tol*(abs(object)+fscale)){ 
      cat("Convergence is succeed\n") # convergence succeed 
      
      # Cholesky decomposition, a method of decomposing the positive-definite matrix
      # test whether hessian is positive-definite matrix
      if (inherits(try(chol(hessian),silent=TRUE),"try-error")) {
        inver_hessian = NULL # if not positive-definite matrix then print null 
        # print warning of not positive definite at convergence.
        warning("The Hessian is not a positive-definite matrix at convergence") 
      }else{
        # calculate the inverse of hessian matrix using Cholesky decomposition 
        inver_hessian <- chol2inv(chol(hessian)) 
      }
      # return converge objective function, theta during converge, number of 
      # iterations going through,corresponding gradient vector and inverse hessian
      return(list(f=object, theta=theta, iter=i, g=gradient, Hi=inver_hessian))  
      
    }else{
      number_perturb <- 0 # initialize the time of perturbing 
      test <- hessian # initialize the result of perturbing  
      
      # Cholesky factor calculation, this calculation is succeeded meaning it perturbs to positive definite. 
      while(inherits(try(chol(test), silent = TRUE), "try-error") == TRUE) {
        #the multiplier:
        #  a small multiple (e.g. 1e-6) of a matrix of the Hessian. 
        #  If that doesn't do it, repeatedly multiply the multiplier by 10 until positive definite matrix is got. 
        test <- hessian+ diag(norm(hessian)*10^number_perturb*1e-6, nrow=nrow(hessian), ncol=ncol(hessian))
        number_perturb <- number_perturb + 1 # count the number of perturbing time. 
      }
      # assign the perturbed and positive definite matrix into the hessian matrix
      hessian <- test  
      
      # Return the step towards optimum value 
      Delta <- backsolve(chol(hessian), forwardsolve(t(chol(hessian)), -gradient))
      
      # Test if the step really lead to improve objective value, if not, halving it for max of max.half times 
      number_halving = 0 # initialize the time of halving
      # when objective become worse or is not finite
      while ((func(theta + Delta, ...)> object)|!is.finite(func(theta + Delta, ...))){  
        if (number_halving < max.half) {  # if within the maximum halving time
          Delta <- Delta / 2 # delta halved 
          number_halving <- number_halving + 1 # halve time increase by 1 
        } else { # if after max.half time, still getting worse and infinite objective value. 
          stop(paste("The step fails to reduce the objective
          despite trying ", as.character(max.half), " halvings")) # stop and print warning. 
        }
      }
      
      theta <- theta + Delta # when delta was proven to be the finite and true improve direction. 
      object<-func(theta,...) # refresh the objective value for next iteration 
      gradient<-grad(theta,...) # refresh the gradient value for next iteration
  
      # get the refreshed hessian matrix under diverse condition
      if(any(is.null(hess))){ # if hessian function is not provided. 
        hessian <- Hessian(theta,grad,...) # use finite difference method we defined 
      } else {
        hessian <- hess(theta,...) # Calling function given 
      }
      
      # iterate again using new objective, hessian and gradient value. 
    }
  }
  
  # Checking the convergence when iteration time turns to maxit
  if (max(abs(gradient)) < (abs(object)+fscale)*tol){# converge?
    cat("Convergence is succeed\n") # if converge, then print information of succeed 
    
    # if converge, then return inverse of hessian 
    if (inherits(try(chol(hessian),silent=TRUE),"try-error")) { # test can be inverse or not 
      inver_hessian = NULL  # if cannot be inverse then print warning 
      warning("The Hessian is not positive definite at convergence")
    }else{ # if can be inverse then return inverse of hessian matrix (will be returned) 
      inver_hessian <- chol2inv(chol(hessian)) 
    }
    # return converge objective function, theta during converge, number of iteration going through,gradient and inverse hessian 
    return(list(f=object, theta=theta, iter=maxit, g=gradient, Hi=inver_hessian)) 
  } else {
    # if after maxit of time still not converge print warning message
    warning(paste("fail to converage after maxit = ", as.character(maxit), " iterations"))
  }
}




