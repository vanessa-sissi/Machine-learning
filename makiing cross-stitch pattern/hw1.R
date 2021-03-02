############################################################################
## problem 1
copula.hw1 <- read.csv("G:/Documents/Toronto/STA 410/copula-hw1.txt", sep="")

XX <- copula.hw1$X
YY <- copula.hw1$Y
eps <- 10^(-8)
# newton rapson 
LL <- function(x){
  mu1 <- x[1]
  mu2 <- x[2]
  theta <- x[3]
  L <- sum(log(dnorm(XX-mu1)*dnorm(YY - mu2)*theta*(1-exp(-theta))*exp(
    -theta*(pnorm(XX-mu1)+pnorm(YY-mu2)))/(1-exp(-theta)-(1-exp(
      -theta*pnorm(XX-mu1)))*(1-exp(-theta*pnorm(YY-mu2))))^2 ))
  return(L)
}
## in this case the second order gradient seems difficult to be calculated
## we use numerical approach
## dl/dmu1
l1 <- function(x){
  mu1 <- x[1]
  mu2 <- x[2]
  theta <- x[3]
  ll <- sum( (XX - mu1) + theta*dnorm(XX - mu1) - 
               2*( 1-exp( -theta*pnorm(YY - mu2) ) ) * 
               exp(-theta*pnorm(XX - mu1))*theta*dnorm(XX-mu1)/
               ( 1 - exp(-theta) - (1-exp(-theta*pnorm(XX-mu1)))*
                   (1 - exp(-theta*pnorm(YY-mu2))) ) )
  return(ll)
}
## debug
## (LL(c(1.00001,1,1)) - LL(c(1,1,1)))/0.00001
## l1(c(1,1,1))

## dl/dmu2
l2 <- function(x){
  mu1 <- x[1]
  mu2 <- x[2]
  theta <- x[3]
  ll <- sum( (YY - mu2) + theta*dnorm(YY - mu2) - 
               2*( 1-exp( -theta*pnorm(XX - mu1) ) ) * 
               exp(-theta*pnorm(YY - mu2))*theta*dnorm(YY-mu2)/
               ( 1 - exp(-theta) - (1-exp(-theta*pnorm(XX-mu1)))*
                   (1 - exp(-theta*pnorm(YY-mu2))) ) )
  return(ll)
}
## l2(c(1,1,1))
## (LL(c(1,1.00001,1)) - LL(c(1,1,1)))/0.00001

## dl/dtheta
l3 <-function(x){
  mu1 <- x[1]
  mu2 <- x[2]
  theta <- x[3]
  ll <- sum( 1/theta + exp(-theta)/(1 - exp(-theta)) - pnorm(XX-mu1)
             - pnorm(YY-mu2) - 2*(exp(-theta) - exp(-theta*pnorm(XX-mu1))*
                pnorm(XX-mu1) - exp(-theta*pnorm(YY-mu2))*
                  pnorm(YY-mu2) + exp(-theta*(pnorm(XX-mu1)+pnorm(YY-mu2)))*
                  (pnorm(XX-mu1) + pnorm(YY-mu2)) )/
               ( 1 - exp(-theta) - (1-exp(-theta*pnorm(XX-mu1)))*
                                        (1 - exp(-theta*pnorm(YY-mu2))) ))
  return(ll)
}
## l3(c(1,1,1))
## (LL(c(1,1,1.00001)) - LL(c(1,1,1)))/0.00001
l11 <- function(x0,x1){
  temp <- x1
  temp[1] <- x0[1]
  ll <- (l1(x1) - l1(temp))/(x1[1]-temp[1])
  return(ll)
}
l12 <- function(x0,x1){
  temp <- x1
  temp[2] <- x0[2]
  ll <- (l1(x1) - l1(temp))/(x1[2]-temp[2])
  return(ll)
}
l13 <- function(x0,x1){
  temp <- x1
  temp[3] <- x0[3]
  ll <- (l1(x1) - l1(temp))/(x1[3]-temp[3])
  return(ll)
}
l21 <- function(x0,x1){
  temp <- x1
  temp[1] <- x0[1]
  ll <- (l2(x1) - l2(temp))/(x1[1]-temp[1])
  return(ll)
}
l22 <- function(x0,x1){
  temp <- x1
  temp[2] <- x0[2]
  ll <- (l2(x1) - l2(temp))/(x1[2]-temp[2])
  return(ll)
}
l23 <- function(x0,x1){
  temp <- x1
  temp[3] <- x0[3]
  ll <- (l2(x1) - l2(temp))/(x1[3]-temp[3])
  return(ll)
}
l31 <- function(x0,x1){
  temp <- x1
  temp[1] <- x0[1]
  ll <- (l3(x1) - l3(temp))/(x1[1]-temp[1])
  return(ll)
}
l32 <- function(x0,x1){
  temp <- x1
  temp[2] <- x0[2]
  ll <- (l3(x1) - l3(temp))/(x1[2]-temp[2])
  return(ll)
}
l33 <- function(x0,x1){
  temp <- x1
  temp[3] <- x0[3]
  ll <- (l3(x1) - l3(temp))/(x1[3]-temp[3])
  return(ll)
}
s0 <- c(mean(XX)-0.0001,mean(YY)-0.0001,0.9999)
s1 <-c(mean(XX),mean(YY),1)
delta <- 1
while(delta>eps){
  dl <- c(l1(s1),l2(s1),l3(s1))
  ddl <- matrix(c(l11(s0,s1),l12(s0,s1),l13(s0,s1),
                  l21(s0,s1),l22(s0,s1),l23(s0,s1),
                  l31(s0,s1),l32(s0,s1),l33(s0,s1)),nrow = 3,byrow = T)
  s0 <- s1
  s1 <- s1 - solve(ddl)%*%dl
  delta <- sqrt(sum((s1-s0)^2))
  print(s1)
}
 nr<-function(x){
   ans <- x - solve(ddl)%*%dl
   return(ans)
 }





LL <- function(x){
  mu1 <- x[1]
  mu2 <- x[2]
  theta <- x[3]
  L <- -sum(log(dnorm(XX-mu1)*dnorm(YY - mu2)*theta*(1-exp(-theta))*exp(
    -theta*(pnorm(XX-mu1)+pnorm(YY-mu2)))/(1-exp(-theta)-(1-exp(
      -theta*pnorm(XX-mu1)))*(1-exp(-theta*pnorm(YY-mu2))))^2 ))
  return(L)
}
mle1 <- optim(c(mean(XX),mean(YY),0.5),LL,hessian = T)


############################################################################

############################################################################
#problem 2
# combinatorial optimization method based on local search
antithetic.boot <- read.table("G:/Documents/Toronto/STA 410/
                              antithetic-boot.txt", quote="\"", comment.char="")

YY2 <- antithetic.boot$V1

## maximize LL2
LL2 <- function(tao1,tao2){
  ll <- sum(YY2*YY2[tao1] + YY2*YY2[tao2] + YY2[tao1]*YY2[tao2])
  return(ll)
}


## local search
## initial solution
## below is an example of one neighbor.... 
tao1 <- sample(1:64,replace = F)
tao2 <- sample(1:64,replace = F)
tao_1 <- tao1
tao_2 <- tao2
## build neighborhood by swapping one element in tao1 and tao2
for(i in 1:63){
  for(j in (i+1):64){
    temp <- tao1
    temp[c(i,j)] <- tao1[c(j,i)]
    tao_1 <- cbind(tao_1,temp)
    temp <- tao2
    temp[c(i,j)]<- tao2[c(j,i)]
    tao_2 <- cbind(tao_2,temp)
  }
}
## calculate the objective function
## obtain the minimumal within the neighborhood
ll2<-rep(0,length(tao_1[1,]))
for(i in 1:length(tao_1[1,])){
  ll2[i]<-LL2(tao_1[,i],tao_2[,i])
}

min(ll2)
tao_1[,ll2==min(ll2)]
tao_2[,ll2==min(ll2)]


## improve by multi-start local search 

LocalSearch_0 <- function(tao1,tao2){
  tao_1 <- tao1
  tao_2 <- tao2
  ## build neighborhood by swapping one element in tao1 and tao2
  for(i in 1:63){
    for(j in (i+1):64){
      temp <- tao1
      temp[c(i,j)] <- tao1[c(j,i)]
      tao_1 <- cbind(tao_1,temp)
      temp <- tao2
      temp[c(i,j)]<- tao2[c(j,i)]
      tao_2 <- cbind(tao_2,temp)
    }
  }
  ## calculate the objective function
  ## obtain the minimumal within the neighborhood
  ll2<-rep(0,length(tao_1[1,]))
  for(i in 1:length(tao_1[1,])){
    ll2[i]<-LL2(tao_1[,i],tao_2[,i])
  }
  
  min(ll2)
  x1<-tao_1[,ll2==min(ll2)]
  x2<-tao_2[,ll2==min(ll2)]
  return(cbind(x1,x2))
}

LocalSearch <- function(tao1,tao2){
  tao12 <- tao1
  tao22 <- tao2
  ll1 <- 1
  ll2 <- 0
  while(ll1 > ll2){
    tao11 <- tao12
    tao21 <- tao22
    temp <- LocalSearch_0(tao11,tao21)
    tao12 <- temp[,1]
    tao22 <- temp[,2]
    ll1 <- LL2(tao11,tao21)
    ll2 <- LL2(tao12,tao22)
  }
  return(cbind(tao11,tao21))
}

iter <- 20
ll2<-c(0)
t1<-rep(0,64)
t2<-rep(0,64)
for(i in 1 : iter){
  tao1 <- sample(1:64,replace = F)
  tao2 <- sample(1:64,replace = F)
  tt <- LocalSearch(tao1,tao2)
  ll2[i] <- LL2(tt[,1],tt[,2])
  t1<- cbind(t1,tt[,1])
  t2<- cbind(t2,tt[,2])
  print(i)
}
t1 <- t1[,-1]
t2 <- t2[,-1]

min(ll2)
t1[,ll2==min(ll2)]
t2[,ll2==min(ll2)]

> min(ll2)
[1] 45.98373
> t1[,ll2==min(ll2)]
[1] 13 45 22 24 32 18 59 44 31 19 50 23 56 42 52 14 25 33 39 26  6 40  9 16  2 58 43
[28]  5 29 53 30 54 36 55 10 64 51 46  3 38  8 28 48 20 21 41 35  7 17 60 62 27 47 37
[55] 63 15 57 34 12  1 61 11 49  4
> t2[,ll2==min(ll2)]
[1] 50 51 47 30 20 33 12 11 46 10 23 17 25 57 21 63 27 28 41 55 36  6 62 43  9 61  1
[28] 53 56  5  3 45 59 38  2 34 58  4 32 52 15 29 40 13 39 14 42 26 19 24 35 60  8 44
[55] 37 54 22 48  7 31 18 16 49 64

############################################################################
# problem 3 
wine <- read.csv("G:/Documents/Toronto/STA 410/wine.dat", sep="")
XX3 <- wine[,-1]

ff <- function(ind1,ind2,ind3){
  ans <-  sum( sum((XX3[ind1,]-mean(XX3[ind1,]))^2) + 
                 sum((XX3[ind2,]-mean(XX3[ind2,]))^2) + 
                 sum((XX3[ind3,]-mean(XX3[ind3,]))^2) )
  return(ans)
}

#ff1 <- function(ind1,ind2,ind3){
#  temp1 <- 0
#  temp2 <- 0
#  temp3 <- 0
#  for(i in ind1){
#    temp1 <- temp1 + sum((XX3[i,]-colMeans(XX3[ind1,]))^2)
#  }
#  for(i in ind2){
#    temp2 <- temp2 + sum((XX3[i,]-colMeans(XX3[ind2,]))^2)
#  }
#  for(i in ind3){
#    temp3 <- temp3 + sum((XX3[i,]-colMeans(XX3[ind3,]))^2)
#  }
#  ans <- c(temp1,temp2,temp3)
#  return(ans)
#}

#ff1<-function(ind1,ind2,ind3){
#  col1 <- colMeans(XX3[ind1,])
#  col2 <- colMeans(XX3[ind2,])
#  col3 <- colMeans(XX3[ind3,])
#  fun1 <- function(i){
#    return(sum((XX3[i,]-col1)^2))
#  }
#  fun2 <- function(i){
#    return(sum((XX3[i,]-col2)^2))
#  }
#  fun3 <- function(i){
#    return(sum((XX3[i,]-col3)^2))
#  }
#  return(c(sum(sapply(ind1,fun1)),sum(sapply(ind2,fun2)),sum(sapply(ind3,fun3))))
#}

ff2<-function(ind){
  set <- XX3[ind,]
  ave <- colMeans(set)
  result <- colSums((t(set)-ave)^2)
  return(sum(result))
}
ff1<-function(ind1,ind2,ind3){
  return(sapply(list(ind1,ind2,ind3),ff2))
}

LocalSearch_03 <- function(ind1,ind2,ind3){
  n1 <- length(ind1)
  n2 <- length(ind2)
  n3 <- length(ind3)
  indx <- list(ind1,ind2,ind3)
  names(indx)<- c('ind1','ind2','ind3')
  for (i  in 1:(n1+n2+n3)){
    for (j in 1:3){
      if(i<=n1&j>1){
        temp1 <- ind1[-i]
        if(j==2){
          temp2 <- c(ind2,ind1[i])
          temp3 <- ind3
        }
        if(j>2){
          temp3 <- c(ind3,ind1[i])
          temp2 <- ind2
        }
        temp <- list(temp1,temp2,temp3)
        names(temp) <- c('ind1','ind2','ind3')
        indx<-cbind(indx,temp)
      }
      if(i>n1&i<=(n1+n2)&j!=2){
        temp2 <- ind2[-(i-n1)]
        if(j<=1){
          temp1 <- c(ind1,ind2[(i-n1)])
          temp3 <- ind3
        }
        if(j>2){
          temp3 <- c(ind3,ind2[(i-n1)])
          temp1 <-ind1
        }
        temp <- list(temp1,temp2,temp3)
        names(temp) <- c('ind1','ind2','ind3')
        indx<-cbind(indx,temp)
      }
      if(i>(n1+n2)&j<=2){
        temp3 <- ind3[-(i-n1-n2)]
        if(j<=1){
          temp1<-cbind(ind1,ind3[(i-n1-n2)])
          temp2 <-ind2
        }
        if(j==2){
          temp2<-cbind(ind2,ind3[i-n1-n2])
          temp1 <-ind1
        }
        temp <- list(temp1,temp2,temp3)
        names(temp) <- c('ind1','ind2','ind3')
        indx<-cbind(indx,temp)
      }
    }
  }
  ll3<-rep(0,dim(indx)[2])
  #for (i in 1:dim(indx)[2]){
    #ll3[i] <- ff(indx[,i]$ind1,indx[,i]$ind2,indx[,i]$ind3)
  #  ll3[i] <- sum(ff1(indx[,i]$ind1,indx[,i]$ind2,indx[,i]$ind3))
  #}
  fun4 <- function(i){
    sum(ff1(indx[,i]$ind1,indx[,i]$ind2,indx[,i]$ind3))
  }
  ll3<-sapply(1:dim(indx)[2],fun4)
  min(ll3)
  return(indx[,which.min(ll3)])
}

LocalSearch3 <- function(ind1,ind2,ind3){
  ind12 <- ind1
  ind22 <- ind2
  ind32 <- ind3
  ll1 <- 1
  ll2 <- 0
  while(ll1 > ll2){
    ind11 <- ind12
    ind21 <- ind22
    ind31 <- ind32
    temp <- LocalSearch_03(ind11,ind21,ind31)
    ind12 <- temp$ind1
    ind22 <- temp$ind2
    ind32 <- temp$ind3
    ll1 <- sum(ff1(ind11,ind21,ind31))
    ll2 <- sum(ff1(ind12,ind22,ind32))
    #ll1 <- ff(ind11,ind21,ind31)
    #ll2 <- ff(ind12,ind22,ind32)
    print(ll1)
  }
  indx <- list(ind11,ind21,ind31)
  names(indx) <- c('ind1','ind2','ind3')
  return(indx)
}


ans <- LocalSearch3(1:59,60:130,131:178)
ind1 <- ans$ind1
ind2 <- ans$ind2
ind3 <- ans$ind3
ff1(ind1,ind2,ind3)
## k-mean cluster
fit <- kmeans(XX3,3)
fit$withinss
sum(fit$withinss)
aaa<-1:178
ind1 <- aaa[fit$cluster==3]
ind2 <- aaa[fit$cluster==2]
ind3 <- aaa[fit$cluster==1]
#ff(ind1,ind2,ind3)
ff1(ind1,ind2,ind3)



> ans <- LocalSearch3(1:59,60:130,131:178)
[1] 1292.847
[1] 1286.508
[1] 1279.635
[1] 1273.675
[1] 1272.269
[1] 1271.405
[1] 1270.939
> fit$withinss
[1] 326.4147 558.8050 385.7191
> sum(fit$withinss)
[1] 1270.939
############################################################################

############
AA=0.5
G=function(xx=1,AA=1)
{
  rez=1628*(AA+(1-AA)*exp(-xx)) / (1013*AA+4075*(1-AA)*exp(-xx))
  return(rez)
}
G.f=expression(1628*(AA+(1-AA)*exp(-xx)) / (1013*AA+4075*(1-AA)*exp(-xx)))
dd=deriv(G.f,"xx",func=T)
lmt=10
a=seq(-lmt,lmt,0.1)
deriv=a
for(i in 1:length(a))
{
  m=dd(a[i])
  deriv[i]=unlist(attributes(m), use.names=FALSE)
}
plot(a,deriv,type="l")




