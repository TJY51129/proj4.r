# ------- Group Member Names----------------------------------------------------
# Jiayi Tu(s2306508)
# Xiaoyu Zhu(s2296761)
# Jingyun Shan(s1943324)
#
# -------- GitHub Repo ---------------------------------------------------------
# https://github.com/TJY51129/proj4.r
#
#
# -------- Individual Contributions --------------------------------------------
# Our team discussed the overall idea of the project together and completed 
# the writing of the hess1 and newt functions together during the group 
# discussion. Besides, Xiaoyu Zhu and Jiayi Tu added all errors and warning
# statements in the code. And Jingyun Shan tested the code as a whole, 
# and gave some feedback to other team members.
# Finally, everyone modified the code together.
#
# Proportion: Jiayi Tu(33%),Xiaoyu Zhu(35%),Jingyun Shan(32%)
#
# -------- Finite difference Hessian function ----------------------------------
# In the condition of hess = null. The hess is not supplied, This function is 
# to calculate the new approximation to the Hessian by finite differencing of 
# the gradient vector.
#
############
hess1 <- function(theta, grad,...){
  eps=1e-6
  n = length(theta)
  grad0 <- grad(theta) # marked 
  fdH <- matrix(0,n,n)  #create a n*n matrix of 0 
  for (i in 1:n){     #loop over each value of the matrix
    
    #we need to get the new theta number by adding with 
    #the the finite difference intervals
    th <- theta; th[i] <- th[i] + eps 
    
    grad1 <- grad(th,...) #calculate the new theta value in the grad function
    
    fdH[i,] <- (grad1 - grad0)/eps #differencing the gradient vector
  }
  fdH=(fdH + t(fdH))/2
  return(fdH)
}

# ------- Outline of Project----------------------------------------------------
#The main function of the problem:
#This function is based on the Newton's method by minimizing successive
#quadratic approximations to the objective function.

#Firstly, in the case of null, we only need to use the hess1 function we just 
#calculated.

#Then, we check if the objective or derivatives are not finite at the initial
#theta.

#After that, we do the positive definate of the hessian which is used to find
#the minimum of the obective function.

#Keep iteration until the value of iteration reach the maxit.

#Finally, we check the convergence, if it converges to the minimum, return the 
#final result. Otherwise, return the error.

newt <- function(theta,func,grad,hess=NULL,...,tol=1e-8,fscale=1,maxit=100,
                 max.half=20,eps=1e-6){
  #the condition hess is null
  dummy <- 0
  if (any(is.null(hess))) {
    hess <- hess1 #set a new hess equal to hess1
    dummy <- 1
  }
  # If the objective or derivatives are not finite at the initial theta
  else if (is.infinite(func(theta))){ 
    stop("Infinite objective at the initial theta")  
  }
  else if (any(is.infinite(grad(theta)))){ 
    stop("Infinite objective at the initial theta")  
  }
  else if (any(is.infinite(hess(theta)))){ 
    stop("Infinite objective at the initial theta")  
  }
  
  
  iter = 0         #the initial number of iter is 0
  
  while (iter < maxit) {      #while loop for any iter smaller than 100
    iter = iter + 1
    
    grad1 = grad(theta,...)
    
    if (dummy == 1) {
      hess <- hess1
      hess2 <- hess(theta,grad)
    } else if (dummy == 0) {
      hess2 <- hess(theta)
    }
    

    
    #We need to do the positive definate of the hessian
    while (inherits(try(chol(hess2), silent=TRUE),"try-error")) { 
      hess2 <- hess2 + diag(length(theta))  #add the identity matrix till we can 
                                            #do the cholesky decomposition
    }
    H11 <- chol(hess2)                      #do the cholesky decomposition
    #set a identity matrix to compute the inverse of hessian
    identity <- diag(length(theta))            
    #the delta minimizing the quadratic (applies near D's minimum)
    step <- -1 * (backsolve(H11, forwardsolve(t(H11), identity)%*%grad1)) 
    
    
    oldf = func(theta,...)         #the original function of theta
    #the new function after using the Taylor's theorem
    newf <- oldf+ (t(step)%*%grad1)+(0.5*t(step)%*%hess2%*%step) 
    
    
    
    step_number = 0
    #while loop for the step until the new function smaller
    #than the original one
    while (newf > oldf) {          
      if (step_number > max.half){
        stop(paste("The update step failed to reduce the objective"))
      }
    print(step_number)
        
      step = step/2
      newf <- oldf+ (t(step)%*%grad1)+(0.5*t(step)%*%hess2%*%step)
      step_number = step_number + 1   #count the step number
    }
    
    theta = theta + step          #change the theta with the new step added to 
                                  #the old theta
    grad2 = grad(theta)           #take the new theta into the grad function
    
    #judge whether the convergence is reach
    if (max(abs(grad2)) < (abs(newf)+fscale)*tol){    
      if (dummy == 1){
        if(inherits(try(chol(hess(theta, grad)), silent=TRUE),"try-error")){
          warning("The Hessian is not positive definite at convergence")
        }else{
          cat("Converged")
          # Hii = chol(hess(theta,grad)) #do the decomposition for the hessian
          ifelse(dummy == 1,Hii <- chol(hess(theta,grad)),
                 Hii <- chol(hess(theta)))
          Hi = chol2inv(Hii)      #and then, we need to inverse the cholesky
          #return the final result
          return(list(f=newf, theta=theta, iter=iter, g=grad2,Hi=Hi)) 
        }
      } else if (dummy == 0) {
        if(inherits(try(chol(hess(theta)), silent=TRUE),"try-error")){
          warning("The Hessian is not positive definite at convergence")
        }else{
          cat("Converged")
          # Hii = chol(hess(theta,grad)) #do the decomposition for the hessian
          ifelse(dummy == 1,Hii <- chol(hess(theta,grad)),
                 Hii <- chol(hess(theta)))
          Hi = chol2inv(Hii)      #and then, we need to inverse the cholesky
          #return the final result
          return(list(f=newf, theta=theta, iter=iter, g=grad2,Hi=Hi))  
        }
      }
      
      } 
    if(iter > maxit){
      warning(paste("Newton optimizer failed to converge after
                  maxit = ", as.character(maxit), " iterations"))
    }
  }
}
