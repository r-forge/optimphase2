OptimDes<-
function(B.init,m.init,alpha,beta,param,x,target=c("EDA","ETSL","ES"),
         sf=c("futility","OF","Pocock"),num.arm,r=0.5,recover=TRUE,
         control=OptimDesControl(),...)
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
  aboveMin<- control$aboveMin

  if(num.arm==1&sf!="futility")
  {
    stop("only futility interim analysis for single-arm trial")
  }
  
  if(num.arm!=1&num.arm!=2)
  {
    stop("the number of treatment arms must be either one or two")
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
    stop("Projected patient times and numbers should be of equal length")
  }

  if(x>=max(B.init))
  {
    stop("The survival time of interest is beyond the projected accrual time")
  }

  if(target!="EDA"&target!="ETSL" & target!="ES")
    stop("unknown target value")
  
  if(trace!=FALSE&trace!=TRUE)
    stop("trace should be a logical value.")

  if(rho.int[1]<0|rho.int[2]>1)
    stop("the endpoints should be between 0 and 1")
    
  ## Step 1, input n, rho1
  b <- length(B.init)

  p0<-s(shape0,scale0,x)
  p1<-s(shape1,scale1,x)
  if(p1<=p0)
  {
    stop("the alternative survival rate should be greater than the null survival rate")
  }

  # single stage sample size, duration of accrual and study length
  fix.d <- FixDes(B.init,m.init,alpha,beta,param,x,num.arm,r)
  n0 <- fix.d$n0
  da <- fix.d$DA
  sl <- fix.d$SL

  ## compute fixed sample size based on exact distributions for
  ## use in optional proportional adjustment of sample sizes and times
  n0E<- fix.d$n0E
  DAE<-fix.d$DAE
  SLE<-fix.d$SLE

  if(n0>floor(sum(m.init)))
  {
    stop(paste("projected patient sample size is below the minimum requirement of ", 
                n0))
  }

  # store  information
  nrho.res <- matrix(0,nrow=sum(m.init)-n0+1,ncol=2)
  nrho.res[,1] <- c(n0:sum(m.init))
  n.res <- nrho.res[,1]
  nvec<-length(n.res)
  # set up vector to hold the minimum EDA, ES or ETSL for each n
  min.res <- numeric(nvec)
  t1.res <- numeric(nvec)
  C1.res <- numeric(nvec)
  C1U.res <- numeric(nvec)
  C2.res <- numeric(nvec)
  MTSL.res <- numeric(nvec)
  MDA.res <- numeric(nvec)
  ETSL.res <- numeric(nvec)
  EDA.res <- numeric(nvec)
  ES.res <- numeric(nvec)
  n1.res <- numeric(nvec)
#  if(num.arm==1)
#  { u.res <-  numeric(nvec) }
#  if(num.arm==2)
  u.res <-  matrix(numeric(2*nvec),ncol=2)       
  se.res <- matrix(numeric(4*nvec),ncol=4)
  
  currentMin<-Inf
  for(n in n0:sum(m.init))
  {
    i<-n-n0+1
    optout <- optimize(f.Des1,interval=rho.int,tol=tol,B.init=B.init,
            m.init=m.init,alpha=alpha,beta=beta,param=param,x=x,n=n,
            num.arm=num.arm,r=r,target=target,sf=sf,conv=conv,recover=recover)
   
    min.res[i] <- optout$objective
    nrho.res[i,2] <- optout$minimum
    
    if(optout$objective==1e10) next
    
    f.out <- f.Des(B.init,m.init,alpha,beta,param,x,n,
                   optout$minimum,num.arm,r,target,sf,conv,recover)
    t1.res[i] <- f.out$t1
    C1.res[i] <- f.out$C1
    C1U.res[i] <- f.out$C1U
    C2.res[i] <- f.out$C2
    MDA.res[i] <- f.out$mda
    MTSL.res[i] <- f.out$mda+x
    se.res[i,] <- f.out$se
    u.res[i,] <- f.out$u
    EDA.res[i]<-f.out$EDA
    ETSL.res[i] <- f.out$ETSL
    ES.res[i] <- f.out$ES
    n1.res[i] <- ceiling(f.out$n1)
  

    if(trace==TRUE)
      cat("n=",nrho.res[i,1],"optimal rho=",nrho.res[i,2],
           target,"=",min.res[i],"\n")
    if(min.res[i]<currentMin)currentMin<-min.res[i]
    else if(min.res[i]>aboveMin*currentMin)break
  }

  ## truncate result vectors at last sample size evaluated
  nrho.res <- nrho.res[1:i,]
  min.res <-  min.res[1:i]
  n.res <-    n.res[1:i]
  t1.res <-   t1.res[1:i]
  C1.res <-   C1.res[1:i]
  C1U.res <- C1U.res[1:i]
  C2.res <-   C2.res[1:i]
  MTSL.res <- MTSL.res[1:i]
  MDA.res <-  MDA.res[1:i]
  ETSL.res <- ETSL.res[1:i]
  EDA.res <-  EDA.res[1:i]
  ES.res <-   ES.res[1:i]
  n1.res <-   n1.res[1:i]
  se.res <-   se.res[1:i,]
  u.res <-  u.res[1:i,]
  

  ## select the optimal n
  order.min <- order(min.res)[1]
  n.last <- nrho.res[order.min,1]
  outcome <- min.res[order.min]
 
  EDA<-EDA.res[order.min]
  ETSL<-ETSL.res[order.min]
  ES <- ES.res[order.min]
  mda <- MDA.res[order.min]
  t1.last <- t1.res[order.min]
  C1.last <- C1.res[order.min]
  C1U.last <- C1U.res[order.min]
  C2.last <- C2.res[order.min]
  n1 <- n1.res[order.min]
  se <- se.res[order.min,]
  u  <- u.res[order.min,]
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


  res <-structure(list(target=target,sf=sf,test=c(alpha=alpha,beta=beta,
             param=param,x=x,recover=recover), design=c(num.arm=num.arm,r=r),
             accrual=list(B.init=B.init,m.init=m.init),
             result=c(EDA=EDA,ETSL=ETSL,ES=ES),n=c(n1=n1,n.last=n.last),
             stageTime=c(t1=t1.last,MTSL=mda+x),boundary=c(C1L=C1.last,
             C1U=C1U.last,C2=C2.last),se=se,u=u,exposure=EW,
             all.info=data.frame(n=n.res,t1=t1.res,C1L=C1.res,C1U=C1U.res,C2=C2.res,
             MDA=MDA.res,MTSL=MTSL.res,EDA=EDA.res,ETSL=ETSL.res,ES=ES.res),
             single.stageTime=c(n0=n0,DA=da,SL=da+x,n0E=n0E,DAE=DAE,SLE=SLE)),
             class="OptimDes")
    return(res)        
}

##############################################################################################  

  ### functions for the weibull distribution
  ### h and s use R definitions of weibull parameters

  # hazard function of Weibull distribution 
  h <- function(shape,scale,x)
  {
    (shape/scale)*(x/scale)^(shape-1)
  }

  # survival function of Weibull distribution
  s <- function(shape,scale,x)
  {
    exp(-(x/scale)^shape)
  }

  # cummulative hazard function of weibull distribution
  lambda <- function(shape,scale,x)
  {
    (x/scale)^shape
  }
  
##############################################################################################

  ### accrual distribution function

  # function for Fn(t)
  fnt <- function(B.init,m.init,n,t,mda,l)
  {
    B <- B.init[1:l]   # "L" not one
    B[l] <- mda
    M <- m.init[1:l]
    M[l] <- ifelse(l>1,n-sum(m.init[1:(l-1)]),n)
    ll <- length(t)
    m1 <- matrix(rep(B,ll),nrow=l,byrow=FALSE)
    m2 <- matrix(0,nrow=l+1,ncol=ll)
    m2[1:l,] <- m1
    m2[l+1,] <- t
    # the interval that t falls in    
    count <- apply(m2,2,indX)
    M.sum <- cumsum(M)
    count.1 <- ifelse(count==1,1,count-1) # avoid error message in the following ifelse function
    pr <-ifelse(count==1,(M[1]/n)*(t/B[1]),ifelse(count<=l,M.sum[count.1]/n+
         (M[count]/n)*((t-B[count.1])/(B[count]-B[count.1])),1))
    pr<-ifelse(t>=mda,1,pr)
    return(pr)
      
  }    


  ### function to find smallest index for entry exceeding an input

  indX<-function(x) {
    rank(x,ties.method="min")[length(x)]
  }

  ### 
    compMDA<-function(B.init,m.init,n){
       l<-indX(c(cumsum(m.init),n))
       mda <- ifelse(l>1,B.init[l-1]+(B.init[l]-B.init[l-1])*(n-
            sum(m.init[1:(l-1)]))/m.init[l],B.init[l]*(n/m.init[l]))
       return(c(mda,l))
    }
    compTime<-function(B.init,m.init,ti){
       l1<-indX(c(B.init,ti))
       n1 <- ifelse(l1>1,sum(m.init[1:(l1-1)])+ m.init[l1]*(ti-
              B.init[l1-1])/(B.init[l1]-B.init[l1-1]),
              m.init[1]*ti/B.init[1])    
       return(n1)
    }

  ### compute the single stage design sample size using exact binomial for one-arm study
  ## one arm
  single.exact<-function(n0,alpha,beta,p0,p1){
     if(p0>=p1)stop("The null event-free rate exceeds the alternative rate")
     fixN<-F
     n0E<-floor(n0/2)
     while(!fixN){

        bx<-qbinom(alpha,n0E,p0,lower.tail=F)
        if( pbinom(bx,n0E,p1)<=beta){
           fixN<-T
        }else n0E<-n0E+1

     }
  return(n0E)
  }


    
  
  
  

truncC<-function(B.init,m.init,n,x,t1){
  ### Form distribution of Y given Y<=t1
  ### then return the expectation of  the conditional
  ### exposure time

  ### compute mda
  mdaout<-compMDA(B.init,m.init,n)
  mda<-mdaout[1]

  if(t1>mda)t1<-mda
  
  ### compute conditional density
  k<-1
  while(B.init[k]<t1)k<-k+1
  m<-m.init[1:k]
  B<-c(0,B.init[1:k])  
  m[k]<-m.init[k]*(t1-B[k])/(B[k+1]-B[k])
  B[k+1]<-t1

  Pi<-m/sum(m)
  Di<-Pi/diff(B)

  ### compute conditional expectation
  EW<-truncM(B,Di,x,t1)
  return(EW)
}



truncM<-function(B,Di,x,t1){
  ### Compute the expectation of exposure over the conditional distribution
  trM<-0
  b<-length(B)
  for (i in 2:b){

     if(t1-x>=B[i]){
         trM<-trM + x*Di[i-1]*(B[i]-B[i-1])
     }else if(t1-x>=B[i-1]){
         tmp<-x*(t1-x-B[i-1]) + t1*(B[i]-t1+x) - 0.5*(B[i]^2-(t1-x)^2)
         trM<-trM+Di[i-1]*tmp
     }else{
         tmp<- t1*(B[i]-B[i-1]) - 0.5*(B[i]^2-B[i-1]^2)
         trM<- trM+Di[i-1]*tmp
     }
  }
  return(trM)
}


  # function for equation 1
  sig2 <- function(B.init,m.init,shape,scale,x,n,t,mda,l)
  {
    ### function to be integrated in function sig2
    f.int <- function(a){
      h(shape,scale,a)/(s(shape,scale,a)*fnt(B.init,m.init,n,t-a,mda,l))
}

    sig2.int <- integrate(f.int,lower=0,upper=x)
    if(sig2.int$message!="OK")
      stop("integration for sig2 cannot be completed")
    else
      sig2.int$value
  }

##############################################################################################

  # for each n and rho1, obtain C1L, C1U and C2 by the type I&II error constraints
  f.Des <- function(B.init,m.init,alpha,beta,param,x,n,rho1,num.arm,r,
                    target=c("EDA","ETSL","ES"),sf=c("futility","OF","Pocock"),
                    conv,recover)
  {

    target <- match.arg(target)
    sf <- match.arg(sf)

    shape0 <- param[1]
    scale0 <- param[2]
    shape1 <- param[3]
    scale1 <- param[4]    

    lam0<-lambda(shape0,scale0,x)
    lam1<-lambda(shape1,scale1,x)

    smallC<-1e-5       ### small offset to avoid division by 0

    if(rho1>=1-smallC)return(list(min=1e10)) ### avoid numerical boundary problems


    ## Step 2, Calculate MDA by n and the accrual function
    mdaout<-compMDA(B.init,m.init,n)
    mda<-mdaout[1]
    l<-mdaout[2]

    ## Step 3
    sig21 <- sqrt(sig2(B.init,m.init,shape1,scale1,x,n,mda+x,mda,l))
    sig20 <- sqrt(sig2(B.init,m.init,shape0,scale0,x,n,mda+x,mda,l))
 
    ## Steps 4 and 5
    if(num.arm==1){
        sig11 <- sig21/rho1

        ### find design associated with rho  
        ### function supplied to root finder for step 5 
        f5 <- function(t){
               sig2(B.init,m.init,shape1,scale1,x,n,t,mda,l)-sig11^2
        }
    }else{
        ## find design associated with rho1
        f5 <-  function(t) {
              sqrt(sig2(B.init,m.init,shape0,scale0,x,n,t,mda,l)/((1-r)*(lam0^2))+
              sig2(B.init,m.init,shape1,scale1,x,n,t,mda,l)/(r*(lam1^2)))-
              (1/rho1)*sqrt((sig20^2)/((1-r)*(lam0^2))+(sig21^2)/(r*(lam1^2)))
        }
    }

    t1.root <- uniroot(f5,c(x+smallC,mda+x))
    if(suppressWarnings(warning()!=""))
      stop("the algorithm of uniroot() does not converge in 'maxiter' steps") else{ 
         t1 <- t1.root$root } 

    sig10 <- sqrt(sig2(B.init,m.init,shape0,scale0,x,n,t1,mda,l))
    if(num.arm==2)sig11 <- sqrt(sig2(B.init,m.init,shape1,scale1,x,n,t1,mda,l))

    ## Step 6
    t2 <- mda-t1

    ## Step 7
    rho0 <- sig20/sig10

    ## Step 8
    if(num.arm==1){
    u2 <- sqrt(n)*(log(lam0)-log(lam1))*lam1/sig21
    u1 <- rho1*u2
    u <- c(u1,u2)
    }else{
        v10 <- (sig10^2)/((1-r)*(lam0^2))
        v11 <- (sig11^2)/(r*(lam1^2))
        v20 <- (sig20^2)/((1-r)*(lam0^2))
        v21 <- (sig21^2)/(r*(lam1^2))
        u1 <- sqrt(n)*(log(lam0)-log(lam1))/sqrt(v10+v11)
        u2 <- sqrt(n)*(log(lam0)-log(lam1))/sqrt(v20+v21) 
        
        u <- c(u1,u2)
    }
    # the number of patients accrued at the time of interim analysis
    n1<-compTime(B.init,m.init,min(t1,mda))

    ## Step 9
    sigma0 <- diag(2)
    sigma1 <- diag(2)
    sigma0[1,2] <- rho0
    sigma0[2,1] <- rho0
    sigma1[1,2] <- rho1
    sigma1[2,1] <- rho1

    C1U <- Inf   
    if(num.arm==2)
        C1U <- qnorm(1-spfun(rho0,alpha,sf))
  

    if(recover){  ### recover alpha from interim stops
	    ### bivariate normal integrands for C1s and C2
        if(num.arm==1){
            bvne <- function(z)
            {
              (pmvnorm(lower=c(z[1],z[2]),upper=c(Inf,Inf),sigma=sigma0)-alpha)^2+
              (pmvnorm(lower=c((z[1]-u1),(z[2]-u2)),upper=c(Inf,Inf),sigma=sigma1)-
                                                        (1-beta))^2
            }
        }else{
          
            bvne <- function(z)
            {
              (pmvnorm(lower=c(min(z[1],C1U),z[2]),upper=c(C1U,Inf),sigma=sigma0)-(alpha-spfun(rho0,alpha,sf)))^2+
              (1-pnorm(C1U-u1)+pmvnorm(lower=c(min((z[1]-u1),(C1U-u1)),(z[2]-u2)),upper=c((C1U-u1),Inf),sigma=sigma1)-
                                                       (1-beta))^2
            }
        }

	    opt <- optim(c(1,1),bvne,method="Nelder-Mead")
	    if(opt$value<conv&opt$convergence==0)
	    {
	      C1 <- opt$par[1]
	      C2 <- opt$par[2]
	    }else{
	      return(list(min=1e10))} #when there is no solution

	}else{ ### do not recover alpha from interim stops

        
        if(num.arm==1){
            C2<- qnorm(1-alpha)
            C1low<- -3   ### essentially no stopping at interim

            C1up<-qnorm(1-beta)+u2
            powC1<-function(C1){
                pmvnorm(lower=c((C1-u1),(C2-u2)),upper=c(Inf,Inf),sigma=sigma1)-
                                                                    (1-beta)
            }
        }else{
        powC2<-function(C2){
            1-pnorm(C2)- pmvnorm(lower=c(C1U,C2),upper=c(Inf,Inf),sigma=sigma0)-
                                                 (alpha-spfun(rho0,alpha,sf))
        }
        C2up <- qnorm(1-(alpha-spfun(rho0,alpha,sf)))
        C2root <- uniroot(powC2,c(-3,C2up))
        if(suppressWarnings(warning()!="")){
           return(list(min=1e10))  #when there is no solution
        }else{
           C2 <- C2root$root
        }        
               
        
        C1up<-C1U
            powC1<-function(C1){  # C1 is C1L in the two-arm
                1-pnorm(C1U-u1)+pmvnorm(lower=c((C1-u1),(C2-u2)),upper=c((C1U-u1),Inf),sigma=sigma1)
                                                                  -(1-beta)
            }
        }

        C1root <- uniroot(powC1,c(C1low,C1up))  ###C1 is C1L in the two-arm
        if(suppressWarnings(warning()!="")){
           return(list(min=1e10))  #when there is no solution
        }else{
           C1 <- C1root$root
        }
    }

    EDA<- ifelse(t1<mda,t1+(pnorm(C1U)-pnorm(C1))*t2,mda)
    ETSL<- t1+(pnorm(C1U)-pnorm(C1))*(t2+x)
    ES<- n1+(pnorm(C1U)-pnorm(C1))*(n-n1)
    if(target=="EDA")
      outcome1 <- EDA
    if(target=="ETSL")
      outcome1 <- ETSL
    if(target=="ES")
      outcome1<- ES
    return(list(min=outcome1,EDA=EDA,ETSL=ETSL,ES=ES,n1=n1,t1=t1,mda=mda,
            C1=C1,C1U=C1U,C2=C2,se=c(sig10,sig20,sig11,sig21),u=u))

  }

##############################################################################################



##############################################################################################

  ### wrapper function for f.Des so it returns the limited output 
  ### expected by function Optimize

  f.Des1 <- function(B.init,m.init,alpha,beta,param,x,n,rho1,num.arm,r,
                     target=c("EDA","ETSL","ES"),sf=c("futility","OF","Pocock"),
                     conv,recover)
  {
    target <- match.arg(target)
    sf <- match.arg(sf)
    out <- f.Des(B.init,m.init,alpha,beta,param,x,n,rho1,num.arm,r,target,sf,conv,recover)
    return(out$min)
  }
  

  ### numerical control values that can be changed by a user
  OptimDesControl<-function(trace=FALSE,tol=0.01,conv=1e-6,rho.int=c(0,1),aboveMin=1.05)
  {
  list(trace=trace,tol=tol,conv=conv,rho.int=rho.int,aboveMin=aboveMin)
}

  ### alpha spending function
  spfun <- function(rho,alpha,sf=c("futility","OF","Pocock"))
  {
    sf <- match.arg(sf)
    if(sf=="futility")
      alpha.rho <- ifelse(rho<1&rho>=0,0,1)
    
    if(sf=="OF")
      alpha.rho <- min(2-2*pnorm(qnorm(1-alpha/2)/sqrt(rho)), alpha)
    
    if(sf=="Pocock")
      alpha.rho <- min(alpha*log(1+(exp(1)-1)*rho),alpha)
  
    return(alpha.rho)
  }