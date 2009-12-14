"FixDes" <-
function(B.init,m.init,alpha,beta,param,x)
{  
  shape0 <- param[1]
  scale0 <- param[2]
  shape1 <- param[3]
  scale1 <- param[4]

  if(length(B.init)!=length(m.init))
  {
    stop("projected patient times and numbers should be of equal length")
  }

  if(x>=max(B.init))
  {
    stop("the survival time of interest is beyond the projected accrual time")
  }
  
  if(qnorm(1-alpha)+qnorm(1-beta)<=0)stop("Stopped because alpha, beta are invalid")

  p0<-s(shape0,scale0,x)
  p1<-s(shape1,scale1,x)

  if(p0>=p1)stop("The null event-free rate exceeds the alternative rate")


  ### sample size and study times based on normal approximation to log
  ### hazard function
  sig21 <- sqrt((1-p1)/p1)  ### delta method applied to binomial variance (based on Taylor expansion) to estimate sig21 in equation 9
  n0 <-(sig21*(qnorm(1-alpha)+qnorm(1-beta))/((log(lambda(shape0,scale0,x))-
        log(lambda(shape1,scale1,x)))*lambda(shape1,scale1,x)))^2
  n0 <- ceiling(n0)
  if(n0>floor(sum(m.init)))warning("Sample size exceeds specified accrual rates/time")

  da<-compMDA(B.init,m.init,n0)[1]

  ### sample size and study times based on exact binomial calculations

  n0E<-single.exact(n0,alpha,beta,p0,p1)
  daE<-compMDA(B.init,m.init,n0E)[1]

  return(list(n0=n0,DA=da,SL=da+x,n0E=n0E,DAE=daE,SLE=daE+x,C=qnorm(1-alpha)))
}  
