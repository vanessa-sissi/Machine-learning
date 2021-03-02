fit_logistic_lasso <- function(x, y, lambda, beta0 = NULL, 
                               eps = 0.0001, iter_max = 100){
  ## fit_logistic_lasso(x, y, lambda, beta0, eps, iter_max) fits a 
  ## lasso logistic regression. 
  ##
  ## Input:
  ## - x: matrix of predictors which does not include the intercept
  ## - y: vector of data -- response variable
  ## - lambda: penalty parameter
  ## - beta0: initial guess for optimization
  ## - eps: parameter for stopping criterion
  ## - iter_max: maximum number of iterations
  ##
  ## Output:
  ## - List of intercept, beta, lambda, factor levels, if converged
  ##   and the numebr of iteration when the training is finished ( 
  ##   converged or reach the max number of iterations)
  ##
  ## Example:
  ## library(tidyverse)
  ## library(tidymodels)
  ## n =  1000
  ## x = cbind(seq(-3,3, length.out = n), 3*cos(3*seq(-pi,pi, length.out = n)))
  ## y = rbinom(n,size = 1, prob = 1/(1 + exp(-x[, 2]+2*x[, 1])) )%>% as.numeric %>% factor
  ## fit = fit_logistic_lasso(x, y, 3)
  
  if(!require(tidymodels)) {
    stop("The tidymodels packages must be installed. Run install.packages(\"tidymodels\") and then try again.")
  }

  n <- dim(x)[1]
  P <- dim(x)[2]
  
  if(is.null((beta0))){
    beta0 <- rep(0, P)
  }
  names(beta0) <- colnames(x)
  fct_levels <- levels(y)
  y <- as.numeric(y) - 1
  
  beta <- beta0
  intercept0 <- 0
  x_beta0 <- (x %*% beta0) %>% as.numeric
  p <- 1/(1 + exp(-x_beta0-intercept0))
  intercept <- intercept0
  
  for(iter in 1:iter_max){
    w <- p * (1 - p)
    z <- x_beta0 + intercept0 + (y - p)/w
    # argmin_intercept_beta 1/2n * t(w)%*%(z - intercept - x%*%beta)^2 + 
    # lambda*sum(abs(beta))

    cd <- coordinate_descent(x, z, w, beta, intercept, lambda, P, eps)
    intercept <- cd$intercept
    beta <- cd$beta
 
    if(max(abs(c(beta-beta0, intercept-intercept0)))<eps){
      return(list(intercept = intercept, 
                  beta = beta, 
                  lambda = lambda,
                  fct_levels = fct_levels,
                  iter = iter,
                  converged = TRUE))
    }

    intercept0 <- intercept
    beta0 <- beta
    x_beta0 <- (x %*% beta0) %>% as.numeric
    p <- 1/(1 + exp(-x_beta0-intercept0))
  }
  
  warning(paste("Method did not converge in", iter_max, "iterations", sep = " "))
  return(list(intercept = intercept, 
              beta = beta, 
              lambda = lambda,
              fct_levels = fct_levels,
              iter = iter, converged = FALSE))
  

}

err <- function(x, z, w, beta, intercept, lambda){
  ## calculaste errors for iterative reweighted lesat square 
  ##
  ## Input:
  ## - x: matrix of predictors (not including the intercept)
  ## - z: z from IRLS
  ## - w: p*(1-p)
  ## - beta: cuurent step's coefficient for x 
  ## - intercept: current intercept
  ## - lambda: parameter for penalty term 
  ## Output:
  ## - error term (current IRLS step)
  ##
  ## Example:
  ## library(tidyverse)
  ## li  1000
  ## x = cbind(seq(-3,3, length.out = n), 3*cos(3*seq(-pi,pi, length.out = n)))
  ## y = rbinom(n,size = 1, prob = 1/(1 + exp(-x[, 2]+2*x[, 1])) )%>% as.numeric %>% factor
  ## beta0 <- rep(0, nrow(x))
  ## intercept0=0
  ## x_beta0 <- (x %*% beta0) %>% as.numeric
  ## p <- 1/(1 + exp(-x_beta0-intercept0))
  ## w <- p * (1 - p)
  ## z <- x_beta0 + intercept0 + (y - p)/w
  ## err(x, z, w, beta0, intercept0, 100)
  
  return( t(w)%*%(z - intercept - x%*%beta)^2 + lambda*sum(abs(beta)))
}

coordinate_descent<-function(x, z, w, beta, intercept, lambda, P, eps){
  ## Update  coordinates by coordinate descent  
  ##
  ## Input:
  ## - x: matrix of predictors (not including the intercept)
  ## - z: z from IRLS
  ## - w: p*(1-p)
  ## - beta: cuurent step's coefficient for x 
  ## - intercept: current intercept
  ## - lambda: parameter for penalty term 
  ## - P: numner rows of x
  ## - eps: tolerence 
  ## Output:
  ## - list of new beta and new intercept 
  ##
  ## Example:
  ## library(tidyverse)
  ## li  1000
  ## x = cbind(seq(-3,3, length.out = n), 3*cos(3*seq(-pi,pi, length.out = n)))
  ## y = rbinom(n,size = 1, prob = 1/(1 + exp(-x[, 2]+2*x[, 1])) )%>% as.numeric %>% factor
  ## beta0 <- rep(0, nrow(x))
  ## intercept0=0
  ## x_beta0 <- (x %*% beta0) %>% as.numeric
  ## p <- 1/(1 + exp(-x_beta0-intercept0))
  ## w <- p * (1 - p)
  ## z <- x_beta0 + intercept0 + (y - p)/w
  ## coordinate_descent(x, z, w, beta0, intercept0, 100, nrow(x), 1e-5)

  e0 <- 0
  e1 <- 100
  while(abs(e0 - e1)>eps){
    e0 <- e1
    for(j in 1:P){
      rj <- (z - x %*% matrix(beta)-intercept + x[,j] %*% matrix(beta[j]))
      beta[j] <- sign(t(w*2 * x[,j]) %*% rj)/sum(w*2 * x[,j]^2)*max(abs(t(w*2 * x[,j]) %*% rj)- lambda, 0)
    }
    rj <- z - x %*% beta
    intercept <- sum(w*rj) / sum(w)
    e1 <- err(x, z, w, beta, intercept, lambda)
  }
  return(list(beta=beta, intercept=intercept, error=e1))
}

predict_logistic_lasso <- function(object, new_x){
  ## predict_logistic_lasso(object, new_x) predicts values 
  ## for given input new_x based on the fitted lasso logistic
  ## model -- object.
  ##
  ## Input:
  ## - object: Output from fit_logistic_lasso -- List of intercept, 
  ##           beta, and lambda for the lasso logistic model
  ## - new_x: Data to predict at
  ## Output:
  ## - Predicted values
  ##
  ## Example:
  ## library(tidyverse)
  ## library(tidymodels)
  ## n =  1000
  ## x = cbind(seq(-3,3, length.out = n), 3*cos(3*seq(-pi,pi, length.out = n)))
  ## y = rbinom(n,size = 1, prob = 1/(1 + exp(-x[, 2]+2*x[, 1])) )%>% as.numeric %>% factor
  ## fit = fit_logistic_lasso(x, y, 3)
  ## pred = predict_logistic_lasso(fit, x)
  
  if(!require(tidymodels)) {
    stop("The tidymodels packages must be installed. Run install.packages(\"tidymodels\") and then try again.")
  }
  numeric_pred <- (new_x %*% object$beta + object$intercept >=0) %>% as.numeric
  # return(list(pred = object$fct_levels[numeric_pred + 1] %>% factor,
  #             intercept = object$intercept,
  #             beta = object$beta))
  return(object$fct_levels[numeric_pred + 1] %>% factor)
}
#set.seed(123)
#n =  1000
#x = cbind(seq(-3,3, length.out = n), 3*cos(3*seq(-pi,pi, length.out = n)))
#y = rbinom(n,size = 1, prob = 1/(1 + exp(-x[, 2]+2*x[, 1])) )%>% as.numeric %>% factor
#fit = fit_logistic_lasso(x, y, 0.3*2000)
#fit

#fit2 = glmnet::glmnet(x,y, family = "binomial", lambda = 0)
#fit$beta
#fit2$beta




