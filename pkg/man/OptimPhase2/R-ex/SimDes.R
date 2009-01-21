### Name: SimDes
### Title: Simulation Studies for Two-Stage Designs with Time to Event
###   Endpoints from function OptimDes
### Aliases: SimDes
### Keywords: iteration optimize

### ** Examples

B.init <- c(1, 2, 3, 4, 5)
m.init <- c(15, 20, 25, 20, 15)
alpha <- 0.05
beta <- 0.1
param <- c(1, 1.09, 2, 1.40)
x <- 1

# H0: S0=0.40 H1: S1=0.60

object1 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="EDA")

SimDes(object1)

SimDes(object1,interimRule='t1')

### accrual rates differ from planned

SimDes(object1,m.init=c(5,5,25,25,25))



