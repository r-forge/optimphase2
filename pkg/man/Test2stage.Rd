\name{Test2stage}
\alias{Test2stage}

\title{Statistical test for two-stage designs from function OptimDes}

\description{
  This function performs the hypothesis tests for the two-stage designs
with event-free endpoint from \code{\link{OptimDes}}.
}
\usage{
Test2stage(x,C1,C1U,C2,tan,num.arm,Y11,T11,p0=NULL,Y21=NULL,T21=NULL,Y10=NULL,
           T10=NULL,Y20=NULL,T20=NULL,printTest=TRUE,cen11=rep(1,length(T11)),
           cen21=rep(1,length(T21)),cen10=rep(1,length(T10)),cen20=rep(1,length(T20)))
}

\arguments{
  \item{x}{Pre-specified time for the event-free endpoint (e.g., 1 year).}
  \item{C1}{The two-stage study is terminated for futility after the
first stage if the Z-statistic is < C1.}
  \item{C1U}{The two-stage study is terminated for efficacy after the
first stage if the Z-statistic is > C1U.}
  \item{C2}{Reject the null hypothesis for the one-sided alternative
hypothesis at the second stage if the Z-statistic is > C2.}
  \item{tan}{Interim study time.}
  \item{num.arm}{Number of treatment arms. \code{num.arm}=1 for single-arm trial assuming a known
  standard control. \code{num.arm}=2 for two-arm randomized trial with a comparative control arm.}
  \item{Y11}{A vector containing the study start times (measured from
the beginning of the study) of first-stage patients in the treatment arm.}
  \item{T11}{A vector containing the event times of first-stage patients in the treatment arm.}
  \item{p0}{The event rate under the null hypothesis.}
  \item{Y21}{A vector containing the study start times (measured from
the beginning of the study) of second-stage patients in the treatment arm.}
  \item{T21}{A vector containing the event times of second-stage
    patients in the treatment arm.}

  \item{Y10}{A vector containing the study start times (measured from
the beginning of the study) of first-stage patients in the control arm when \code{num.arm}=2.}
  \item{T10}{A vector containing the event times of first-stage patients 
  in the control arm when \code{num.arm}=2.}
  \item{Y20}{A vector containing the study start times (measured from
the beginning of the study) of second-stage patients in the control arm when \code{num.arm}=2.}
  \item{T20}{A vector containing the event times of second-stage
    patients in the treatment arm when \code{num.arm}=2.}
  \item{printTest}{If TRUE (default), the result of the test and the
interim decision is printed.}
  \item{cen11}{The times in \code{T11} are regarded as events unless
they are set to censored by setting the corresponding value in
\code{cen11} to zero.}
  \item{cen21}{The times in \code{T21} are regarded as events unless
they are set to censored by setting the corresponding value in
\code{cen21} to zero.}
  \item{cen10}{The times in \code{T10} are regarded as events unless
they are set to censored by setting the corresponding value in
\code{cen10} to zero.}
  \item{cen20}{The times in \code{T20} are regarded as events unless
they are set to censored by setting the corresponding value in
\code{cen20} to zero.}
}

\details{
  The hypothesis tests are performed in two stages as described in Huang, Talukder
  and Thomas (2010) and Case and Morgan (2003) for single-arm designs, and 
  extended to the randomized two-arm designs.

  Stage 1. Accrue patients between time 0 and time \code{t1}. Each
  patient will be followed until failure, or for \code{x} years or until
  time \code{t1}, whichever is less. Calculate the normalized interim test
  statistic \code{Z1}. If \code{Z1<C1}, stop the study for futility; For
  randomized two-arm trials, if \code{Z1>C1U}, stop the study for efficacy;
  otherwise, continue to the next stage.

  Stage 2. Accrue patients between \code{t1} and \code{MDA}. Follow all
  patients until failure or for \code{x} years, then calculate the
  normalized final test statistic \code{Z2}, and reject H0 if
  \code{Z2>C2}.

  The test statistic is based on the Nelson-Aalen estimator of the
  cumulative hazard function.
}

\value{
  A vector containing results for the interim analysis if \code{T21}, \code{T20},
\code{Y21} and \code{Y20} are NULL, and the 
final analysis if \code{T21}, \code{T20}, \code{Y21} and \code{Y20} are non-NULL:

  \item{z}{The test statistic}
  \item{se}{Standard error of sum of the cummulative hazards (not log
cummulative hazards) at time \code{x}.}
  \item{cumL}{A two-element vector of cummulative hazard estimators at time \code{x}.}
}

\author{Bo Huang \email{<bo.huang@pfizer.com>} and Neal Thomas \email{<neal.thomas@pfizer.com>}}

\references{

  Huang B., Talukder E. and Thomas N. Optimal two-stage Phase II designs with long-term
  endpoints. \emph{Statistics in Biopharmaceutical Research}, \bold{2(1)}, 51--61. 
  
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

### single arm trial

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

object1 <- OptimDes(B.init,m.init,alpha,beta,param,x,target="EDA",sf="futility",num.arm=1)
n <- object1$n[2]
t1 <- object1$stageTime[1]
C1 <- object1$boundary[1]
C1U <- object1$boundary[2]
C2 <- object1$boundary[3]
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
Test2stage(x, C1, C1U, C2, t1,1,Y1, T1, p0=p0)

# final analysis if the study continues
Test2stage(x, C1, C1U, C2, t1,1, Y1, T1, p0=p0)

# simulate failure times of length n under the null hypothesis
sample.T <- rweibull(n,shape=shape0,scale=scale0)

Y1 <- sample.Y[sample.Y<=t1]
T1 <- sample.T[sample.Y<=t1]
Y2 <- sample.Y[sample.Y>t1]
T2 <- sample.T[sample.Y>t1]

# interim analysis
Test2stage(x, C1, C1U, C2, t1,1,Y1, T1, p0=p0)

# final analysis if the study continues
Test2stage(x, C1, C1U, C2, mda+x,1,Y1, T1, p0=p0,Y21=Y2,T21=T2)


}
\keyword{htest}
\keyword{optimize}
