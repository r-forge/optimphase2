### Name: Test2stage
### Title: Statistical test for two-stage designs from function OptimDes
### Aliases: Test2stage
### Keywords: htest optimize

### ** Examples

B.init <- c(1, 2, 3, 4, 5)
m.init <- c(15, 20, 25, 20, 15)
alpha <- 0.05
beta <- 0.1
param <- c(1, 1.09, 2, 1.40)
x <- 1

# H0: S0=0.40 H1: S1=0.60

shape0 <- param[1]
scale0 <- param[2]
shape1 <- param[3]
scale1 <- param[4]

object1 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="EDA")
n <- object1$n[2]
t1 <- object1$stageTime[1]
C1 <- object1$boundary[1]
C2 <- object1$boundary[2]
b <- length(B.init)
l <- rank(c(cumsum(m.init),n),ties.method="min")[b+1]
mda <- ifelse(l>1,B.init[l-1]+(B.init[l]-B.init[l-1])*(n-sum(m.init[1:(l-1)]))/m.init[l],B.init[l]*(n/m.init[l]))

### set up values to create a stepwise uniform distribution for accrual
B <- B.init[1:l]
B[l] <- mda
xv <- c(0,B)
M <- m.init[1:l]
M[l] <- ifelse(l>1,n-sum(m.init[1:(l-1)]),n)
yv <- c(0,M/(diff(xv)*n),0)

# density function of accrual 
dens.Y <- stepfun(xv,yv,f=1,right=TRUE)
# pool of time points to be simulated from
t.Y <- seq(0,mda,by=0.01)

# simulate study times of length n
sample.Y <- sample(t.Y,n,replace=TRUE,prob=dens.Y(t.Y))

# simulate failure times of length n under the alternative hypothesis
sample.T <- rweibull(n,shape=shape1,scale=scale1)

Y1 <- sample.Y[sample.Y<=t1]
T1 <- sample.T[sample.Y<=t1]
Y2 <- sample.Y[sample.Y>t1]
T2 <- sample.T[sample.Y>t1]

# event rate under null hypothesis
p0<-pweibull(x,shape=shape0,scale=scale0)

# interim analysis
Test2stage(Y1, T1, Y2 = NULL, T2 = NULL, p0, x, C1, C2, t1, MTSL = NULL)

# final analysis if the study continues
Test2stage(Y1, T1, Y2, T2, p0, x, C1, C2, t1, mda+x)

# simulate failure times of length n under the null hypothesis
sample.T <- rweibull(n,shape=shape0,scale=scale0)

Y1 <- sample.Y[sample.Y<=t1]
T1 <- sample.T[sample.Y<=t1]
Y2 <- sample.Y[sample.Y>t1]
T2 <- sample.T[sample.Y>t1]

# interim analysis
Test2stage(Y1, T1, Y2 = NULL, T2 = NULL, p0, x, C1, C2, t1, MTSL = NULL)

# final analysis if the study continues
Test2stage(Y1, T1, Y2, T2, p0, x, C1, C2, t1, mda+x)




