### Name: np.OptimDes
### Title: Optimal Two-Stage Design with User-specified Combined Sample
###   Size
### Aliases: np.OptimDes
### Keywords: design optimize

### ** Examples

B.init <- c(1, 2, 3, 4, 5)
m.init <- c(15, 20, 25, 20, 15)
alpha <- 0.05
beta <- 0.1
param <- c(1, 1.09, 2, 1.40)
x <- 1

# H0: S0=0.40 H1: S1=0.60

np.OptimDes(B.init,m.init,alpha,beta,param,x,pt=1.1,target="ETSL")



