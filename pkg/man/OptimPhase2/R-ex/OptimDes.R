### Name: OptimDes
### Title: Construct an Optimal Two-stage Design for a Time to Event
###   Endpoint Versus a Known Control Standard
### Aliases: OptimDes
### Keywords: design optimize

### ** Examples

B.init <- c(1, 2, 3, 4, 5)
m.init <- c(15, 20, 25, 20, 15)
alpha <- 0.05
beta <- 0.1
param <- c(1, 1.09, 2, 1.40)
x <- 1

# H0: S0=0.40 H1: S1=0.60

object1 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="EDA")
print(object1)

object2 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="ETSL")
plot(object2)

m.init <- c(25,20,15,20,25)
object3 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="ETSL")




