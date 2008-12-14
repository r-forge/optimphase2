\name{Test2stage}
\alias{Test2stage}

\title{Statistical test for two-stage designs from function OptimDes}

\description{
  This function performs the hypothesis tests for the two-stage designs
with event-free endpoint from \code{\link{OptimDes}}.
}
\usage{
Test2stage(Y1, T1, Y2 = NULL, T2 = NULL, p0, x, C1, C2, t1, MTSL =
NULL, printTest = TRUE, cen1=rep(1,length(T1)),cen2=rep(1,length(T2)))
}

\arguments{
  \item{Y1}{A vector containing the study start times (measured from
the beginning of the study) of first-stage patients.}
  \item{T1}{A vector containing the event times of first-stage patients.}
  \item{Y2}{A vector containing the study start times (measured from
the beginning of the study) of second-stage patients.}
  \item{T2}{A vector containing the event times of second-stage
    patients.}
  \item{p0}{The event rate under the null hypothesis.}
  \item{x}{Pre-specified time for the event-free endpoint (e.g., 1 year).}
  \item{C1}{The two-stage study is terminated for futility after the
first stage if the Z-statistic is < C1}
  \item{C2}{Reject the null hypothesis for the one-sided alternative
hypothesis at the second stage if the Z-statistic>C2.}
  \item{t1}{Time for the interim analysis (actual time, which may
differ from the planned time).}
  \item{MTSL}{Maximum total study length, (actual, not planned)}
  \item{printTest}{If TRUE (default), the result of the test and the
interim decision is printed.}
  \item{cen1}{The times in \code{T1} are regarded as events unless
they are set to censored by setting the corresponding value in
\code{cen1} to zero.}
  \item{cen2}{The times in \code{T2} are regarded as events unless
they are set to censored by setting the corresponding value in
\code{cen1} to zero.}
}

\details{
  The hypothesis tests are performed in two stages as described in Case
  and Morgan (2003).

  Stage 1. Accrue patients between time 0 and time \code{t1}. Each
  patient will be followed until failure, or for \code{x} years or until
  time \code{t1}, whichever is less. Calculate the normalized interim test
  statistic \code{Z1}. If \code{Z1<C1}, stop the study for futility;
  otherwise, continue to the next stage.

  Stage 2. Accrue patients between \code{t1} and \code{MDA}. Follow all
  patients until failure or for \code{x} years, then calculate the
  normalized final test statistic \code{Z2}, and reject H0 if
  \code{Z2>C2}.

  The test statistic is based on the Nelson-Aalen estimator of the
  cumulative hazard function.
}

\value{
  A vector containing results for the interim analysis if \code{T2} and
\code{Y2} are NULL, and the 
final analysis if \code{T2} or \code{Y2} are non-NULL:

  \item{z}{The test statistic}
  \item{se}{Standard error of the cummulative hazard (not log
cummulative hazard) at time \code{x}.}
  \item{cumL}{Cummulative hazard estimator at time \code{x}.}
}

\author{Bo Huang \email{<bo.huang@pfizer.com>} and Neal Thomas \email{<neal.thomas@pfizer.com>}}

\references{
  Case M. D. and Morgan T. M. (2003) Design of Phase II cancer trials
  evaluating survival probabilities. \emph{BMC Medical Research
  Methodology}, \bold{3}, 7.

  Lin D. Y., Shen L., Ying Z. and Breslow N. E. (1996) Group seqential
  designs for monitoring survival probabilities. \emph{Biometrics},
  \bold{52}, 1033--1042.

  Simon R. (1989) Optimal two-stage designs for phase II clinical
  trials. \emph{Controlled Clinical Trials}, \bold{10}, 1--10.
  
}
  
\seealso{\code{\link{OptimDes}}, \code{\link{SimDes}}}

\examples{
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


}
\keyword{htest}
\keyword{optimize}
