np.OptimDes <-
function(B.init,m.init,alpha,beta,param,x,n=NULL,pn=NULL,pt=NULL,
         target=c("EDA","ETSL","ES"),sf=c("futility","OF","Pocock"),
         num.arm,r=0.5,recover=TRUE,
         control=OptimDesControl(), CMadj=F, ...)
{
  
  target <- match.arg(target)
  sf <- match.arg(sf)

  shape0 <- param[1]
  scale0 <- param[2]
  shape1 <- param[3]
  scale1 <- param[4]

  trace <- control$trace
  tol <- control$tol
  conv <- control$conv
  rho.int <- control$rho.int
  
  if(num.arm==1&sf!="futility")
  {
    stop("only futility interim analysis for single-arm trial")
  }
  
  if(num.arm!=1&num.arm!=2)
  {
    stop("The number of treatment arms must be either one or two")
  }   
 
  if(num.arm!=1 & CMadj==T)
  {
    stop("The CM adjustment is for single-group trials only")
  }   


  if(num.arm==1&r!=0.5)
  {
    warning("randomization ratio not equal to 0.5 for a single-arm study")
  }
  
  if(num.arm==2&r*(1-r)<=0)
  {
    stop("invalid randomization ratio r")
  }
  
  if(length(B.init)!=length(m.init))
  {
    stop("projected patient times and numbers should be of equal length")
  }

  if(x>=max(B.init))
  {
    stop("the survival time of interest is beyond the projected accrual time")
  }

  if(target!="EDA"&target!="ETSL" & target!="ES")
    stop("unknown target value")
  
  if(trace!=FALSE&trace!=TRUE)
    stop("trace should be a logical value.")

  if(rho.int[1]<0|rho.int[2]>1)
    stop("the endpoints should be between 0 and 1")
    
  if(s(shape1,scale1,x)<=s(shape0,scale0,x))
  {
    stop("the alternative survial rate should be greater than the null survival rate")
  }

  if(length(c(n,pn,pt))!=1)
    stop("There should be one input from n, pn or pt")

  b <- length(B.init)
  # single stage sample size, duration of accrual and study length
  fix.d <- FixDes(B.init,m.init,alpha,beta,param,x,num.arm,r)
  n0 <- fix.d$n0
  da <- fix.d$DA
  sl <- fix.d$SL

  n0E<- fix.d$n0E
  DAE<-fix.d$DAE
  SLE<-fix.d$SLE

  if(!CMadj){
     nfix <- n0
     sl <- fix.d$SL
  }else{
     nfix <- n0E
     sl <- fix.d$SLE
  }

  if(!is.null(pn))
    n <- ceiling(pn*nfix)
  if(!is.null(pt))
  {
    t <- pt*sl-x
    n<-ceiling(compTime(B.init,m.init,t))
  }  
  
  optout <- optimize(f.Des1,interval=rho.int,tol=tol,B.init=B.init,
                     m.init=m.init,alpha=alpha,beta=beta,param=param,
                     x=x,n=n,num.arm=num.arm,r=r,target=target,sf=sf,conv=conv,recover=recover)
  outcome <- optout$objective
  rho1 <- optout$minimum
  f.out <- f.Des(B.init,m.init,alpha,beta,param,x,n,rho1,num.arm,r,target,sf,conv,recover)

  EDA <- f.out$EDA
  ETSL <- f.out$ETSL
  ES <- f.out$ES
  mda <- f.out$mda
  t1.last <- f.out$t1
  C1.last <- f.out$C1
  C1U.last <- f.out$C1U
  C2.last <- f.out$C2
  n.last <- n
  n1 <- ceiling(f.out$n1)
  se<-f.out$se    
  u <-f.out$u
  EW<-truncC(B.init,m.init,n.last,x,t1.last) #potential exposure at t1
  EW<-n1*EW
  tadj<-(n0E/n0)*t1.last
  if(ceiling((n0E/n0)*n.last)>floor(sum(m.init))) 
  {
    stop("projected patient sample size is below the minimum requirement for normal approximation adjustment")
  }  
  EWadj<-truncC(B.init,m.init,ceiling((n0E/n0)*n.last),x,tadj)
  EWadj<-n1*(n0E/n0)*EWadj
  EW<-c(EW,EWadj)

  res<-structure(list(target=target,sf=sf,test=c(alpha=alpha,beta=beta,param=param,
      x=x,recover=recover), design=c(num.arm=num.arm,r=r),
      accrual=list(B.init=B.init,m.init=m.init),
      result=c(EDA=EDA,ETSL=ETSL,ES=ES),n=c(n1=n1,n.last=n.last),
      stageTime=c(t1=t1.last,MTSL=mda+x),boundary=c(C1L=C1.last,
      C1U=C1U.last,C2=C2.last),se=se,u=u,exposure=EW,all.info=NULL,
      single.stageTime=c(n0=n0,DA=da,SL=da+x,n0E=n0E,DAE=DAE,SLE=SLE)),
      class="OptimDes")
 
  return(res)
  
    
}
