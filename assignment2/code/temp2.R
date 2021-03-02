set.seed(123)
n =  1000
x = cbind(seq(-3,3, length.out = n), 3*cos(3*seq(-pi,pi, length.out = n)))
y = rbinom(n,size = 1, prob = 1/(1 + exp(-x[, 2]+2*x[, 1])) )%>% as.numeric %>% factor

err <- function(x, z, w, beta, intercept, lambda){
  return(1/2 / length(w) * t(w)%*%(z - intercept - x%*%beta)^2 + lambda*sum(abs(beta)))
}

S <- function(z, lambda){
  return(sign(z)*max(abs(z)-lambda, 0))
}

update_b <- function(w, x, z, z_hat, lambda){
  temp = 0
  temp2 = 0
  n = length(w)
  for(i in 1:n){
    temp = temp + w[i]/n*x[i,j]*(z[i]-z_hat[i])
    temp2 = temp2 + w[i]/n*x[i,j]^2
  }
  # ans = S(t(w*x) %*% (z - z_hat), lambda)
  # ans = ans / (t(w) %*% x^2)
  ans = S(temp, lambda) / temp2
  
  return(ans)
}

get_z <- function(x, y, beta, intercept, p){
  z  = y
  for(i in 1:length(y)){
    z[i] = intercept + x[i,] %*% beta + (y[i]-p[i])/(p[i]*(1-p[i]))
  }
  return(z)
}

get_zh <-  function(x, y, beta, intercept, p, j){
  z  = y
  P = dim(x)[2]
  for(i in 1:length(y)){
    for(l in 1:P){
      if(l != j)
        z[i] = intercept + x[i,l] %*% beta[l]
    }
  }
  return(z)
}

cd <- function(x, y, lambda)


n <- dim(x)[1]
P <- dim(x)[2]

beta0 <- rep(0, P)

names(beta0) <- colnames(x)
fct_levels <- levels(y)
y <- as.numeric(y) - 1


fit = glmnet::glmnet(x,y, family = "binomial", lambda = lambda)

intercept0 <- 0

intercept0 <- fit$a0
beta0 <- fit$beta

beta <- beta0
x_beta0 <- (x %*% beta0) %>% as.numeric
p <- 1/(1 + exp(-x_beta0-intercept0))
intercept <- intercept0
lambda = 0.01



w <- p * (1 - p)
z <- x_beta0 + intercept0 + (y - p)/w
# z - get_z(x, y, beta, intercept, p)


j=1
z_hat <- intercept0 + x[,-j] %*% matrix(beta[-j])
# z_hat - get_zh(x, y, beta, intercept, p, j)
beta1 <- sign(sum(w/n * x[,j]*(z-z_hat))) / sum(w/n*x[,j]^2) * max( abs(sum(w/n * x[,j]*(z-z_hat)))-lambda, 0)
# beta1
beta2 <- update_b(w, x, z, z_hat, lambda)

rj <- z - x[,-j] %*% matrix(beta[-j]) - intercept
# max(abs(t(w * x[,j]) %*% rj)-n*lambda, 0)
beta[j] = update_b(w, x, z, z_hat, lambda)
# beta[j] <- sign(t(w/n * x[,j]) %*% rj)/sum(w/n * x[,j]^2) * max(abs(t(w/n * x[,j]) %*% rj)-lambda, 0)
# beta[j]

j=2
z_hat <- intercept0 + x[,-j] %*% matrix(beta[-j])
beta1 <- sign(sum(w/n * x[,j]*(z-z_hat))) / sum(w/n*x[,j]^2) * max( abs(sum(w/n * x[,j]*(z-z_hat)))-lambda, 0)
# beta1
rj <- (z - x[,-j] %*% matrix(beta[-j])-intercept)
# max(abs(t(w * x[,j]) %*% rj)-n*lambda, 0)
beta[j] = update_b(w, x, z, z_hat, lambda)
# beta[j] <- sign(t(w/n * x[,j]) %*% rj)/sum(w/n * x[,j]^2) * max(abs(t(w/n * x[,j]) %*% rj)-lambda, 0)
# beta[j]

rj <- z - x %*% beta
intercept <- sum(w*rj) / sum(w)
err(x, z, w, beta, intercept, lambda)
fit1 = lm(z ~ x, weights = w)
err(x, z, w, fit1$coef[2:3], fit1$coef[1], lambda)
fit2 = glmnet(x, z, lambda = lambda, weights = w/n, family = "gaussian", standardize = FALSE) 
err(x, z, w, fit2$beta, fit2$a0, lambda)

beta
fit2$beta
coordinate_descent(x, z, w, fit$beta%>%as.numeric, fit$a0, lambda, P, 1e-5)

err(x, z, w, fit$beta, fit$a0, lambda)

# grad <- t(x) %*% (y - p)
# sum(grad^2)
# max(abs(c(beta-beta0, intercept-intercept0)))


intercept0 <- intercept
beta0 <- beta
x_beta0 <- (x %*% beta0) %>% as.numeric
p <- 1/(1 + exp(-x_beta0-intercept0))


# fit_IRLS(x, y)
# fit$beta
# fit2$beta
