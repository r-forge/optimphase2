"SimDes" <- 
function(object,B.init,m.init,weib0,weib1,interimRule='e1',
         sim.n=1000,e1conv=1/365,CMadj=F,attainI=1,attainT=1,FixDes="F",
         Rseed)
{  

  interimRule<-match.arg(interimRule,c('e1','n1','t1'))
  if(interimRule!='n1' & interimRule!='t1' & interimRule!='e1')stop('Invalid interimRule')


  if(missing(B.init))B.init<-object$accrual$B.init
  if(missing(m.init))m.init<-object$accrual$m.init
  if(length(B.init)!=length(m.init))stop(
     "The length of the accrual rates and intervals must be equal")
  if(missing(weib0))weib0<-c(object$test["param1"],object$test["param2"])
  if(missing(weib1))weib1<-c(object$test["param3"],object$test["param4"])

  if(!missing(Rseed))set.seed(Rseed)

  if(FixDes=="E" | FixDes=="N")interimRule='t1'

  cat("\n The interimRule is ",interimRule,"\n\n")

  x<-object$test["x"]
  n1<-object$n[1]
  n<- object$n[2]
  t1 <- object$stageTime[1]
  C1 <- object$boundary[1]
  C2 <- object$boundary[2]

  if(length(B.init)!=length(m.init))
  {
    stop("projected patient times and numbers should be of equal length")
  }

  if(x>=max(B.init))
  {
    stop("the survival time of interest is beyond the projected accrual time")
  } 

  ### apply CM correction factor

  Stime<-object$single.stageTime

  ### apply adjustment to exact binomial fixed design size
  ### per CM 2003
  if (CMadj){
        CMadjfac=Stime["n0E"]/Stime["n0"]
        exposure<-object$exposure[2]
  }else{
        CMadjfac<-1
        exposure<-object$exposure[1]
  }

  if(FixDes=="E"){n<-Stime["n0E"]
  }else if(FixDes=="N"){n<-Stime["n0"]
  }else{
      n<-ceiling(CMadjfac*attainT*n)
      ### the attained levels are applied at the beginning of each
      ### iteration because they may change based on interimRule
      n1<-ceiling(CMadjfac*n1)
      t1<-CMadjfac*t1
      if(ceiling(CMadjfac*n1*attainI)>=n)stop("Single-stage design is optimal")
  }


  shape0 <- weib0[1]
  scale0 <- weib0[2]
  shape1 <- weib1[1]
  scale1 <- weib1[2]  

  ### event rate under null hypothesis (not event-free as in OptimDes)
  p0<-pweibull(x,shape=shape0,scale=scale0)
  
  ### cutoff for number of events if exact test is used at second stage
  cexact<-qbinom(object$test["alpha"],n,p0)
  if(pbinom(cexact,n,p0)>object$test["alpha"])cexact<-cexact-1


  if(n>sum(m.init))stop("sample size n greater than total potential accrual")

  ### compute mda and distribution restricted to (0,mda)
  ### the CMadjfac is applied to mda through the adjustment to n
  mdaout<-compMDA(B.init,m.init,n)
  mda<-mdaout[1]
  l<-mdaout[2]

  B <- B.init[1:l]  # "l" not "one"
  B[l]<-mda
  M <- m.init[1:l]
  M[l]<- ifelse(l>1, n-sum(m.init[1:(l-1)]) , n )
  Pi<-M/sum(M)  # proportion in each accrual interval
  B<-c(0,B)     # add B0 for diff calculation


  ### compute truncated exposure conditional on Y<=t1

  EW<-exposure
  tlow<-0
  thigh<-max(B.init)+x

  
  alphaNorm <- 0
  alphaExact<- 0
  powNorm   <- 0
  powExact  <- 0
  eda<-0
  etsl<-0
  es<-0
  aveE<-0
  intvarNull<-0
  finvarNull<-0
  intvarAlt<-0
  finvarAlt<-0
  n1int  <- 0
  t1int <- 0
  phatKl<-1
  phatKh<-0
  phatRl<-1
  phatRh<-0
  lambda0<-lambda(shape0,scale0,x)
  pstopNull<-0
  pstopAlt<-0

  for(i in 1:sim.n)
  {
   
    randI <- sample(1:l,n,replace=TRUE,prob=Pi) #sample accrual intervals
    u<-runif(n)
    sample.Y <- sort( B[randI] + u*(B[randI+1]-B[randI]) )
    
    ### compute time for interim to match planned exposure
    if(interimRule=='e1'){
        
       t1.root <- uniroot(ft,c(tlow,thigh),tol=e1conv,y=sample.Y,x=x,EW=EW)
       if(suppressWarnings(warning()!=""))
          stop("The algorithm of uniroot() did not converge for exposure time calculation")
       else{
          t1<-t1.root$root
          }
       t1<-t1*attainI
       if(t1<x+sample.Y[1]){
          warning("Matched exposure yielded interim time less than x from enrollment")
          t1<-x+sample.Y[1]
       }
    }else if(interimRule=='t1'){
       t1<-object$stageTime[1]*CMadjfac*attainI
       if(t1<x+sample.Y[1]){
          warning("Adjusted t1 occurred before time x from enrollment")
          t1<=x+sample.Y[1]
       }
    }else n1<-ceiling(object$n[1]*CMadjfac*attainI)
       
    ### alpha level
    sample.T <- rweibull(n,shape=shape0,scale=scale0) #null survival times
    if(interimRule=='t1' | interimRule=='e1'){
       n1<-sum(sample.Y<=t1)
       Y1 <- sample.Y[sample.Y<=t1]
       T1 <- sample.T[sample.Y<=t1]
       Y2 <- sample.Y[sample.Y>t1]
       T2 <- sample.T[sample.Y>t1]
    }else{
       if(x>sample.Y[n1]-sample.Y[1]){
          t1<-x+sample.Y[1]
          n1<-min(n,sum(sample.Y<t1)+1)
          warning("n1 patients were accrued before time x from enrollment")
       }else  t1<-sample.Y[n1]
       Y1 <- sample.Y[1:n1]
       T1 <- sample.T[1:n1]
       if(n1<n){
          Y2 <- sample.Y[(n1+1):n]
          T2 <- sample.T[(n1+1):n]
       }else{
          Y2<-numeric(0)
          T2<-numeric(0)
       }
    }

    
    ### compute interim and z-statistics
    ### the final z-statistic is computed even if the interim
    ### is less than c1.  the stopping rule is imposed below
    intout <- Test2stage(Y1, T1, Y2 = NULL, T2 = NULL, p0, x, C1,
                        C2, t1, MTSL = NULL, printTest=FALSE)

    finout <- Test2stage(Y1, T1, Y2, T2, p0, x, C1, C2, t1, mda+x,
                        printTest=FALSE)

    r.int<-intout["z"]
    r.fin<-finout["z"]
    yobs<-sum(c(T1,T2)<=x)

    if(FixDes=="E" | FixDes=="N"){
       alphaExact<- alphaExact + 1*(yobs<=cexact)
       alphaNorm<-  alphaNorm  + 1*(r.fin>C2)
    }else{
        intY<-(r.int<=C1)
        eda<-eda + (!intY)*max(sample.Y) + intY*max(Y1)
        etsl<-etsl + (!intY)*(max(sample.Y)+x) + intY*t1
        es<-es + intY*n1 + (!intY)*n
        aveE<-aveE + sum(pmin(x,pmax(0,t1-Y1)))
        pstopNull<-pstopNull + 1*intY*(n1<n)
         
        phat1<-exp(-intout["cumL"])
        phat2<-exp(-finout["cumL"])
  
        if(r.int>C1)alphaExact<-alphaExact + 1*(yobs<=cexact)
  
        if(r.int>C1&(r.fin>C2)){
          alphaNorm <- alphaNorm + 1
          phatRl<-min(phatRl,phat2)
          phatKl<-min(phatKl,phat1)
        }else if(r.int>C1){
          phatRh<-max(phatRh,phat2)
          phatKl<-min(phatKl,phat1)
        }else{
           phatKh<-max(phatKh,phat1)
        }
        ### proportion of info at interim
        intvarNull<-intvarNull+(intout["se"])^2
        finvarNull<-finvarNull+(finout["se"])^2
  
        n1int<-n1int+n1
        t1int<-t1int+t1
     } 

    ### power
    ### re-use patient start times and associated calculations
    sample.T <- rweibull(n,shape=shape1,scale=scale1)
    
    if(interimRule=='t1' | interimRule=='e1'){
       T1 <- sample.T[sample.Y<=t1]
       T2 <- sample.T[sample.Y>t1]
    }else{
       T1 <- sample.T[1:n1]
       if(n1<n){
          T2 <- sample.T[(n1+1):n]
       }else{
          T2<-numeric(0)
       }
    }

    intout <- Test2stage(Y1, T1, Y2 = NULL, T2 = NULL, p0, x, C1,
                        C2, t1, MTSL = NULL, printTest=FALSE)

    finout <- Test2stage(Y1, T1, Y2, T2, p0, x, C1, C2, t1, mda+x,
                        printTest=FALSE)

    r.int<-intout["z"]
    r.fin<-finout["z"]
    yobs<-sum(c(T1,T2)<=x)

    if(FixDes=="E" | FixDes=="N"){
       powExact<- powExact + 1*(yobs<=cexact)
       powNorm<-  powNorm  + 1*(r.fin>C2)
    }else{
       phat1<-exp(-intout["cumL"])
       phat2<-exp(-finout["cumL"])
       pstopAlt<-pstopAlt+1*(r.int<=C1)*(n1<n)
       if(r.int>C1)powExact<-powExact + 1*(yobs<=cexact)
   
       if(r.int>C1&(r.fin>C2)){
         powNorm <- powNorm + 1
         phatRl<-min(phatRl,phat2)
         phatKl<-min(phatKl,phat1)
       }else if(r.int>C1){
         phatRh<-max(phatRh,phat2)
         phatKl<-min(phatKl,phat1)
       }else{
          phatKh<-max(phatKh,phat1)
       }
       ### proportion of info at interim
       intvarAlt<-intvarAlt + (intout["se"])^2
       finvarAlt<-finvarAlt + (finout["se"])^2
    }
  }

  alphaExact <- alphaExact/sim.n
  alphaNorm  <- alphaNorm/sim.n
  powExact   <- powExact/sim.n
  powNorm    <- powNorm/sim.n
  eda   <- eda/sim.n
  etsl  <- etsl/sim.n
  es    <- es/sim.n
  aveE  <- aveE/sim.n
  n1int <- n1int/sim.n
  t1int <- t1int/sim.n
  intvarNull<-intvarNull/sim.n
  finvarNull<-finvarNull/sim.n
  pinfoNull <- finvarNull/intvarNull
  pstopNull <- pstopNull/sim.n
  intvarAlt<-intvarAlt/sim.n
  finvarAlt<-finvarAlt/sim.n
  pinfoAlt  <- finvarAlt/intvarAlt
  pstopAlt  <- pstopAlt/sim.n
  names(t1int)<-""
  names(n1int)<-""
  names(eda)<-""
  names(etsl)<-""
  names(pinfoNull)<-""
  names(pinfoAlt)<-""
  names(pstopNull)<-""
  names(pstopAlt)<-""
  names(es)<-""
  names(alphaExact)<-""
  names(alphaNorm)<-""
  names(powExact)<-""
  names(powNorm)<-""


  return(c(alphaExact=alphaExact,alphaNorm=alphaNorm,
           powerExact=powExact,powerNorm=powNorm,
           eda=eda,etsl=etsl,
           es=es,pstopNull=pstopNull,pstopAlt=pstopAlt,aveE=aveE,
           pinfoNull=pinfoNull,pinfoAlt=pinfoAlt,
           n1=n1int,t1=t1int,
           phatKl=phatKl,phatKh=phatKh,
           phatRl=phatRl,phatRh=phatRh))

}


ft<-function(ti,y,x,EW){
  tymin<-pmax(0,ti-y)
  return( sum(pmin(x,tymin))-EW )
}

