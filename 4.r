## Jiayi Li, S2304671 :  commenting for the code, coding realization and correcting. 
## Rongkai Fang, s2310813 : commenting for the code, coding designing and guidance.
## Yifan Hu, s2287683 ： commenting for the code, coding realization and correcting. 
## every one contribute roughly equal on this project. 


#finite difference when hessian function is not given
#Input:
#theta, a vector of initial values for the optimization parameters.
#gradient, the gradient function.
#eps, finite difference interval
#Output:approximated hessian matrix.
Hessian <- function(theta,gradient,eps=1e-6,...){ # eps is finite difference interval
  grad0 <- grad(theta,...) # grad of function at initial theta. 
  Hfd <- matrix(0L, nrow = length(theta), ncol = length(theta)) ## initilize finite diference Hessian
  for (i in 1:length(theta)) { # loop over parameters
    th1 <- theta; th1[i] <- th1[i] + eps ## increase th0[i]/initial theta by eps
    grad1 <- grad(th1,...) ## compute resulting gradient
    Hfd[i,] <- (grad1 - grad0)/eps ## approximate second derivs  
  }
  Hfd  <- 0.5 * (t(Hfd ) + Hfd) # turn it to symmetric matrix
  return(Hfd) # return approximated hessian matrix. 
}


#This function implementing Newton’s method for minimization of functions
#Input:
#theta, a vector of initial values for the optimization parameters.
#func, the objective function to minimize. Its first argument is the vector of optimization parameters. Remaining
#arguments will be passed from newt using ‘...’.
#grad, the gradient function. It has the same arguments as func but returns the gradient vector of the objective
#hess, the Hessian matrix function. It has the same arguments as func but returns the Hessian matrix of the
#objective w.r.t. the elements of parameter vector. If not supplied then newt should obtain an approximation to the
#Hessian by finite differencing of the gradient vector.
#..., any arguments of func, grad and hess after the first (the parameter vector) are passed using this.
#tol, the convergence tolerance.
#fscale, a rough estimate of the magnitude of func near the optimum - used in convergence testing.
#maxit, the maximum number of Newton iterations to try before giving up.
#max.half, the maximum number of times a step should be halved before concluding that the step has failed to
#improve the objective.
#eps, the finite difference intervals to use when a Hessian function is not provided.
#Output:
#f, the value of the objective function at the minimum.
#theta, the value of the parameters at the minimum.
#iter, the number of iterations taken to reach the minimum.
#g, the gradient vector at the minimum (so the user can judge closeness to numerical zero).
#Hi, the inverse of the Hessian matrix at the minimum (useful if the objective is a negative log likelihood).


newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,max.half=20,eps=1e-6){
  object <- func(theta,...) # get objective value under initial theta.
  gradient <- grad(theta,...)# get the gradient under initial theta. 
  
  # if no hessian as input:
  if(any(is.null(hess))) hessian <- Hessian(theta,...)  # calculation by finite differencing  
  else hessian <- hess(theta,...) # if hess function is definded and input it, just call it 
  
  
  # Checking objective or gradient or hessian is finite or not. 
  if (!is.finite(object) || any(!is.finite(gradient))|| any(!is.finite(hessian)))
    stop("The objective or gradient or hessian is not finite under the initial defined theta") ## stop and print warning 
  
  
  for(i in 1:maxit){# loop over maxit of time to generate new theta in improving direction
    if(max(abs(gradient)) < tol*(abs(object)+fscale)){ # check overage or not 
      cat("Convergence is succeed") # covergence succeed 
      
      # Cholesky decomposition, a method of decomposing the positive-definite matrix
      if (inherits(try(chol(hessian),silent=TRUE),"try-error")) {# test whether hessian is positive-definite matrix
        inver_hessian = NULL # if not positive-definite matrix then print null 
        warning("The Hessian is not a positive-definite matrix at convergence")# print warning of not positive definite at convergence. 
        
      }else{
        
        # inverse of hessian matrix to be return  
        inver_hessian <- chol2inv(chol(hessian)) 
      }
      # return converge objective function, theta during converge, number of interation going through,gradient and inverse hessian
      return(list(f=object, theta=theta, iter=i, g=gradient, Hi=inver_hessian))  # return when converge 
      
      break # break when converge before the loop ends 
    }else{
      
      number_perturb <- 0 # initial the time of perturbing 
      test <- hessian # the result of perturbing after each interation of perturbing 
      
      # Cholesky factor calculation, this calculation is successed meaning it perturbs to positive definite. 
      while(inherits(try(chol(test), silent = TRUE), "try-error") == TRUE) {

        #the multiplier:
        #  a small multiple (e.g. 1e-6) of a matrix of the Hessian. 
        #  If that doesn't do it, repeatedly multiply the multiplier by 10 until positive definite matrix is got. 
        test <- hessian+ diag(norm(hessian)*10^number_perturb*1e-6, nrow=nrow(hessian), ncol=ncol(hessian))
        number_perturb <- number_perturb + 1 # count the number of perturbing time. 
        
      }
      hessian <- test # assign the perturbed and positive definite matrix into the hessian matrix 
      
      # Return the step towards optimum value 
      Delta <- backsolve(chol(hessian), forwardsolve(t(chol(hessian)), -gradient))
      
      
      # Test if the step really lead to improve objective value, if not, halving it for max of max.half times 
      number_halving = 0 # initialize the time of halving
      while ((func(theta + Delta, ...)> object)|!is.finite(func(theta + Delta, ...))){ #when objective become worse or is not finite 
        if (number_halving < max.half) {  # if within the maximum halving time
          Delta <- Delta / 2 # delta halved 
          number_halving <- number_halving + 1 # halve time increase by 1 
        } else { # if after max.half time, still getting worse and infinite objective value. 
          stop(paste("The step fails to reduce the objective
          despite trying ", as.character(max.half), " halvings")) # stop and print warning. 
        }
      }
      
      theta <- theta + Delta # when delta was proven to be the finite and true improve direction. 
      object<-func(theta,...) # refresh the objective value for next interation 
      gradient<-grad(theta,...) # refresh the gradient value for next interation
  
      # get the refreshed hession matrix under diverse condition
      if(any(is.null(hess))){ # if hession function is not provided. 
        hessian <- Hessian(theta, f,...) # use finite difference method we defined 
      } else {
        hessian <- hess(theta,...) # Calling function given 
      }
      
      # interate again using new objective, hessian and gradient value. 
    }
  }
  
  # Checking the convergence when iteration time turn to maxit
  if (max(abs(gradient)) < (abs(object)+fscale)*tol){# converge?
    cat("Convergence is succeed") # if converge, then print information of succeed 
    
    # if converge, then return inverse of hessian 
    if (inherits(try(chol(hessian),silent=TRUE),"try-error")) { # test can be inverse or not 
      inver_hessian = NULL  # if cannot be inverse then print warning 
      warning("The Hessian is not positive definite at convergence")
    }else{ # if can be inverse then return inverse of hessian matrix (will be returned) 
      inver_hessian <- chol2inv(chol(hessian)) 
    }
    # return converge objective function, theta during converge, number of interation going through,gradient and inverse hessian 
    return(list(f=object, theta=theta, iter=maxit, g=gradient, Hi=inver_hessian)) 
  } else {
    # if after maxit of time still not converge print warning message
    warning(paste("Newton optimizer failed to converage after
                  maxit = ", as.character(maxit), " times of iterations"))
  }
}



 



