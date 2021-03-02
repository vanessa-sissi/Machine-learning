set.seed(123)
n =  1000
x = cbind(seq(-3,3, length.out = n), 3*cos(3*seq(-pi,pi, length.out = n)))
y = rbinom(n,size = 1, prob = 1/(1 + exp(-x[, 2]+2*x[, 1])) )%>% as.numeric %>% factor



n <- dim(x)[1]
P <- dim(x)[2]

beta0 <- rep(0, P)

names(beta0) <- colnames(x)
fct_levels <- levels(y)
y <- as.numeric(y) - 1

beta <- beta0
intercept0 <- 0
x_beta0 <- (x %*% beta0) %>% as.numeric
p <- 1/(1 + exp(-x_beta0-intercept0))
intercept <- intercept0
lambda = 0.3

fit2 = glmnet::glmnet(x,y, family = "binomial", lambda = lambda, alpha=1)


w <- p * (1 - p)
z <- x_beta0 + intercept0 + (y - p)/w

fit <- glmnet(x, z, lambda = lambda, alpha=1, weights = w) 
beta <- fit$beta
intercept <- fit$a0
intercept0 <- intercept

x_beta0 <- (x %*% beta) %>% as.numeric
p <- 1/(1 + exp(-x_beta0-intercept0))
beta
intercept

> beta
2 x 1 sparse Matrix of class "dgCMatrix"
s0
V1 -0.6282348
V2  0.1771945
> intercept
s0 
-0.01146842 
