################################################################
######## Multiple Scaled Contaminated Asymmetric Laplace #######
###############################################################
## Tortora C., Franczak B.C., Bagnato L., Punzo A.
## A Laplace-based model with flexible tail behavior
##



library(gtools)
library(mvtnorm)
library(MASS)
library(Matrix)


################GAMMA#############
###tuning gamma with Monte Carlo search

###input data nxp matrix
### n.trials number of gamma martices used in the search, 2000 by default

MSCAL <- function(data,n.trials=2000){
  G=1
  n=nrow(data)
  p=ncol(data)
  gamMain=array(1,dim=c(p,p,1))
  
  d        <- ncol(data)          # Data dimension
  n.par    <- d*(d-1)/2  # Numer of Gamma parameters
  
  S          <- cov(data)
  temp       <- eigen(S)
  Gamma.init <- temp$vectors
  temp2      <- QPl(Gamma.init)
  P          <- temp2$P
  PAR.init   <- temp2$l
  
  PAR        <- array(NA,c(n.trials,n.par),dimnames = list(paste("par",1:n.trials,sep=" "),paste("dim",1:n.par,sep=" ")))
  Gamma      <- list()     # Gamma matrices to be tried
  ll.grid    <- numeric(n.trials) ##stores the lik
  PAR[1,]    <- PAR.init
  Gamma[[1]] <- Gamma.init
  cont <- 0
  for(m in 1:n.trials){
    #  cat("*")
    cont <- m/n.trials*100
    # cat(paste(cont,"%",sep=""))
    #cat("\n")
    PAR[m,]     <- PAR.init +runif(n.par,-0.1,0.1)  # rnorm(n.par,mean=0,sd=.01)      
    Gamma[[m]]  <- PlQ(P = P, l = PAR[m,])
    gamMain[,,1]=as.matrix(Gamma[[m]])
    out         <- mainMSCMSAL(data,gamMain,G=G)
    ll.grid[m]  <- max(out$ll)
  }
  
  Best.pos   <- which.max(ll.grid)
  Best.Gamma <- Gamma[[Best.pos]]
  Best.PAR   <- PAR[Best.pos,]
  Best.ll    <- max(ll.grid)
  gamMain[,,1]<-  as.matrix(Best.Gamma)
  #Best.out   <- mainMSCMSAL(data,gamMain,G=G)
  temp2_2      <- QPl(Best.Gamma)
  
  plot(ll.grid, type="b",cex=0.6,xlab="trial",ylab="Obs. Log-Lik")
  
  n.trials2 <- 5
  PAR2      <- array(NA,c(n.trials2,n.par),dimnames = list(paste("par",1:n.trials2,sep=" "),paste("dim",1:n.par,sep=" ")))
  Gamma2    <- list()     # Gamma matrices to be tried
  ll.grid2  <- numeric(n.trials2) ##stores the lik
  PAR2[1,]    <- Best.PAR
  Gamma2[[1]] <- Best.Gamma
  ll.grid2[1] <- Best.ll 
  cont <- 0
  for(m in 2:n.trials2){
    #  cat("*")
    cont <- m/n.trials2*100
    # cat(paste(cont,"%",sep=""))
    #cat("\n")
    PAR2[m,]    <- Best.PAR  + jitter(rep(0,n.par))  
    Gamma2[[m]]  <- PlQ(P = P, l = PAR2[m,])
    gamMain[,,1]=as.matrix(Gamma2[[m]])
    out2         <- mainMSCMSAL(data,gamMain,G=G)
    ll.grid2[m]  <- max(out2$ll)
  }
  
  points(ll.grid2, type="b",cex=0.9,col='red')
  Best.pos2   <- which.max(ll.grid2)
  Best.Gamma2 <- Gamma2[[Best.pos2]]
  gamMain[,,1] = as.matrix(Best.Gamma2)
  Best.out2   <- mainMSCMSAL(data,gamMain,G=G)
  
  plot(Best.out2$ll,type="l",main='Loglikelihood per iteration')
  plot(data,col=Best.out2$cluster,main='Clustering results with outlier detection')
  points(data[which(Best.out2$detect[,1]==1),],col="blue",pch=3)
  points(data[which(Best.out2$detect[,2]==1),],col="green",pch=3)
  
  return(list(out=Best.out2, ll.grid=ll.grid,ll.grid2=ll.grid2 ))
  
}



###############################Random generation###############################

rSAL <-function(n,p,skew,sigma,mu){
  sig=sigma
  alpha=skew
  Y <- matrix(mvrnorm(n,rep(0,p),sig),nrow=n,ncol=p,byrow=TRUE)
  w <- rexp(n)
  alpha <- matrix(alpha,nrow=n,ncol=p,byrow=TRUE)
  X.1 <- w*alpha + sqrt(w)*Y
  mu <- matrix(mu,ncol=p,nrow=n,byrow=TRUE)
  X <- X.1 + mu
  return(X)
}
rSL <-function(n,p,sigma,mu){
  sig=sigma
  alpha=rep(0,p)
  Y <- mvrnorm(n,rep(0,p),Sigma=sig)
  w <- rexp(n)
  alpha <- matrix(alpha,nrow=n,ncol=p,byrow=TRUE)
  X.1 <- w*alpha + sqrt(w)*Y
  mu <- matrix(mu,ncol=p,nrow=n,byrow=TRUE)
  X <- X.1 + mu
  return(X)
}
rMSSAL <- function(n,p,mu=rep(0,p),skew=rep(0,p),sigma=diag(p),gam=NULL,phi=NULL) {
  
  x = matrix(0, nrow=n, ncol=p)
  y = matrix(0, nrow=n, ncol=p)
  if(is.null(gam)){ phi=diag(p)
  diag(phi)=(eigen( sigma)$values)
  gam=eigen( sigma)$vectors}
  skew_y = as.vector(t(gam)%*%skew)
  mu_y =as.vector(t(gam)%*%mu)
  
  for (i in 1:n) {
    ww=matrix(0,p,1)
    for(j in 1:p){
      w = rexp(n = 1, 1)
      ww[j,1]=w}
    wd=diag(p)
    diag(wd)=ww
    #y[i,] = t((wd%*%sqrt(phi))%*%t(skew_y))+(rnorm(p)%*%(sqrt(phi)%*%sqrt(wd)))#wd%*%
    y[i,] = t(wd%*%skew_y)+(rnorm(p)%*%(sqrt(phi)%*%sqrt(wd)))#wd%*%
    
  }
  
  
  yn=sweep(y,2,mu_y,FUN="+")
  x =  yn%*% t(gam)
  val = x
  return(val)
}
rCSAL <- function(n,p,mu=rep(0,p),skew=rep(0,p),rho=0.5,sigma=diag(p),eta=1.01) {
  good=rbinom(n,1,rho)
  X<-matrix(0,n,p)
  for(i in 1:n){
    # Delta[i,,]     <- diag((good[i,]+(1-good[i,])*rho*phi)^(-1))
    sigstar  <- (good[i]+(1-good[i])/eta)^(-1)*sigma
    skewstar<- (good[i]+(1-good[i])/sqrt(eta))^(-1)*(skew)
    X[i,]<- rSAL(1,p,mu=mu, sigma=sigstar,skew=skewstar)
  }
  
  return(X)
}
rMSCAL <- function(n,p,mu=rep(0,p),skew=rep(0,p),rho=rep(0.5,p),sigma=diag(p),eta=rep(1.01,p),gam=NULL, phi=NULL) {
  if(is.null(gam)){
    gam=eigen(sigma)$vectors
    phi=eigen(sigma)$values
  }else if(is.null(phi)){cat("phi cannot be NULL")}
  
  good=matrix(0,n,p)
  for(i in 1:p){
    good[,i]=rbinom(n,1,rho[i])}
  SigmaCont <- matrix(0,n,p)
  SkewCont<-matrix(0,n,p)
  X<-matrix(0,n,p)
  for(i in 1:n){
    SigmaCont[i,]   <- (good[i,]*phi+(1-good[i,])*eta*phi)
    SkewCont[i,]<- good[i,]*(skew)%*%(gam)+(1-good[i,])*sqrt(eta)*(skew)%*%(gam)
    X[i,]<- rMSSAL(1,p, gam=diag(p),phi= diag(p)*SigmaCont[i,],skew=SkewCont[i,])
  }
  X=sweep(X%*%t(gam),2,mu,FUN="+")
  
  return(X)
}




###############################Densities###############################

dMSCAL <- function(data,p,mu=rep(0,p),rho=rep(0.5,p),sigma=diag(p),eta=rep(1.01,p),skew=rep(0,p),gam=NULL,phi=NULL,rotated=F) {
  if(!is.matrix(data))data=matrix(data,length(data)/p,p)
  if(is.null(gam)){gam=eigen(sigma)$vectors
  phi=eigen(sigma)$values
  }else if(is.null(phi)){cat("phi cannot be NULL")}
  
  dg= dMSSALdim(data,p,mu=mu,sigma=sigma,skew=skew,gam=gam,phi=phi,rotated = rotated)
  #dd=sweep(data%*%gam,2,mu%*%gam+sqrt(rho)*skew%*%gam,FUN="-")
  if(rotated==F){db= dMSSALdim(data%*%gam,p,gam=diag(p),phi=eta*phi,skew=sqrt(eta)*(skew%*%gam),mu=mu%*%gam) }
  else{db= dMSSALdim(data,p,gam=diag(p),phi=eta*phi,skew=sqrt(eta)*(skew),mu=t(mu))}
  den11=sweep(dg,2,rho,FUN="*")
  den12=sweep(db,2,(1-rho),FUN="*")
  den1=den11+den12
  den = apply(den1, 1, prod)
  return(den)
}


dSAL <- function(data,p,skew,sigma,mu){
  sig=sigma
  alpha=skew
  x <- as.matrix(data)
  talpha <- as.matrix(alpha,nrow=1,ncol=p);
  alpha <- t(talpha);
  n <- nrow(x);
  nu <- (2-p)/2;
  if(p>1){ log.det.sig <- log(det(sig))
  }else{log.det.sig=log(sig)}
  if(log.det.sig == "NaN") stop('Degenerate Solution Reached - The Determinate of this matrix is Zero');
  ty.t <- sweep(x,2,mu,"-");
  inv.sig <- solve(sig)
  t1.num <- sweep(ty.t%*%inv.sig%*%talpha,2,log(2),"+");
  t1.den <- (p/2)*log(2*pi)+0.5*log.det.sig;
  t1 <- sweep(t1.num,2,t1.den,"-");
  maha <- mahalanobis(x=x,center=mu,cov=inv.sig,inverted=TRUE); 
  ata <- 2 + alpha%*%inv.sig%*%talpha;
  l.t2.num <- log(maha)
  l.t2.den <- log(ata)
  l.t2.den <- rep(l.t2.den,n)
  t2 <- 0.5*nu*(l.t2.num-l.t2.den)
  u <- exp( 0.5*(l.t2.den + l.t2.num) )
  t3 <- log(besselK(u,nu,expon.scaled=TRUE)) - u
  val1 <- t1+t2+t3
  val <- exp(val1)
  return(c(val))
}

dSL <- function(data,p,sigma,mu){
  sig=sigma
  alpha=rep(0,p)
  x <- as.matrix(data)
  talpha <- as.matrix(alpha,nrow=1,ncol=2);
  alpha <- t(talpha);
  n <- nrow(x);
  nu <- (2-p)/2;
  if(p>1){ log.det.sig <- log(det(sig))
  }else{log.det.sig=log(sig)}
  if(log.det.sig == "NaN") stop('Degenerate Solution Reached - The Determinate of this matrix is Zero');
  ty.t <- sweep(x,2,mu,"-");
  inv.sig <- solve(sig)
  t1.num <- sweep(ty.t%*%inv.sig%*%talpha,2,log(2),"+");
  t1.den <- (p/2)*log(2*pi)+0.5*log.det.sig;
  t1 <- sweep(t1.num,2,t1.den,"-");
  maha <- mahalanobis(x=x,center=mu,cov=inv.sig,inverted=TRUE); 
  ata <- 2 + alpha%*%inv.sig%*%talpha;
  l.t2.num <- log(maha)
  l.t2.den <- log(ata)
  l.t2.den <- rep(l.t2.den,n)
  t2 <- 0.5*nu*(l.t2.num-l.t2.den)
  u <- exp( 0.5*(l.t2.den + l.t2.num) )
  t3 <- log(besselK(u,nu,expon.scaled=TRUE)) - u
  val1 <- t1+t2+t3
  val <- exp(val1)
  return(c(val))
}

dMSSAL <- function(data,p,mu=rep(0,p),sigma=diag(p),skew=rep(0,p),gam=NULL,phi=NULL,rotated=F,PARrotated=F) {
  if(!is.matrix(data))data=matrix(data,length(data)/p,p)
  x=data
  ###if rotated=F the data, mu, and skew are in the original space
  ### x is the data set
  ##par is a list with all the parameter
  ## this is the density for 1 component so g is fixed
  # x is a n x p matrix
  if(is.null(gam)){gam=eigen(sigma)$vectors
  phi=eigen(sigma)$values
  }else if(is.null(phi)){cat("phi cannot be NULL")}
  
  if(rotated==FALSE){
    x=x%*%gam}
  if(PARrotated==FALSE) {skew=t(gam)%*%skew
  mu=t(gam)%*%mu
  }####t(gam )%*% par per ruotare parametri
  den=rep(1,nrow(x))
  for(i in 1:p){
    den=den*dSAL(x[,i],1,skew[i],phi[i],mu[i])
  }
  
  return(den)
}

dMSSALdim <- function(data,p,mu=rep(0,p),sigma=diag(p),skew=rep(0,p),gam=NULL,phi=NULL,rotated=F) {
  ###density of an MSSAL per dimension, it does not multiply for all dimensions
  
  if(!is.matrix(data))data=matrix(data,length(data)/p,p)
  y=data
  ### y is the data set
  ##par is a list with all the parameter
  ## this is the density for 1 component so g is fixed
  # x is a n x p matrix
  if(is.null(gam)){gam=eigen(sigma)$vectors
  phi=eigen(sigma)$values
  }else if(is.null(phi)){cat("phi cannot be NULL")}
  if(rotated==F){
    x     = y %*% (gam)  
    # mu    = par$mu; phi = par$phi;
    # alpha = par$alpha;
    skew =skew%*%gam#t(skew)# 
    mu =mu%*%gam}
  else{x=y}
  den=matrix(NA,nrow(x),ncol(x))
  for(j in 1:p){
    den[,j]=dSAL(x[,j],1,skew[j],phi[j],mu[j])
  }  
  return(den)
}

dCSAL <- function(data,p,mu=rep(0,p),rho=0.5,sigma=diag(p),eta=1.01,skew=rep(0,p)) {
  if(!is.matrix(data))data=matrix(data,length(data)/p,p)
  
  dg=dSAL(data,p,skew=skew,sigma=sigma,mu=mu)  #dd=sweep(data%*%gam,2,mu%*%gam+sqrt(rho)*skew%*%gam,FUN="-")
  db=dSAL(data,p,skew=sqrt(eta)*skew,sigma=eta*sigma,mu=mu) 
  den = rho*dg+(1-rho)*db
  return(den)
}

dMSCMSSAL <- function(data,p,mu=rep(0,p),rho=rep(0.5,p),sigma=diag(p),eta=rep(1.01,p),skew=rep(0,p),gam=NULL,phi=NULL,rotated=F) {
  if(!is.matrix(data))data=matrix(data,length(data)/p,p)
  if(is.null(gam)){gam=eigen(sigma)$vectors
  phi=eigen(sigma)$values
  }else if(is.null(phi)){cat("phi cannot be NULL")}
  
  dg= dMSSALdim(data,p,mu=mu,sigma=sigma,skew=skew,gam=gam,phi=phi,rotated = rotated)
  #dd=sweep(data%*%gam,2,mu%*%gam+sqrt(rho)*skew%*%gam,FUN="-")
  if(rotated==F){db= dMSSALdim(data%*%gam,p,gam=diag(p),phi=eta*phi,skew=sqrt(eta)*(skew%*%gam),mu=mu%*%gam) }
  else{db= dMSSALdim(data,p,gam=diag(p),phi=eta*phi,skew=sqrt(eta)*(skew),mu=t(mu))}
  den11=sweep(dg,2,rho,FUN="*")
  den12=sweep(db,2,(1-rho),FUN="*")
  den1=den11+den12
  den = apply(den1, 1, prod)
  return(den)
}


dMSCMSSAL <- function(data,p,mu=rep(0,p),rho=rep(0.5,p),sigma=diag(p),eta=rep(1.01,p),skew=rep(0,p),gam=NULL,phi=NULL,rotated=F) {
  if(!is.matrix(data))data=matrix(data,length(data)/p,p)
  if(is.null(gam)){gam=eigen(sigma)$vectors
  phi=eigen(sigma)$values
  }else if(is.null(phi)){cat("phi cannot be NULL")}
  
  dg= dMSSALdim(data,p,mu=mu,sigma=sigma,skew=skew,gam=gam,phi=phi,rotated = rotated)
  #dd=sweep(data%*%gam,2,mu%*%gam+sqrt(rho)*skew%*%gam,FUN="-")
  if(rotated==F){db= dMSSALdim(data%*%gam,p,gam=diag(p),phi=eta*phi,skew=sqrt(eta)*(skew%*%gam),mu=mu%*%gam) }
  else{db= dMSSALdim(data,p,gam=diag(p),phi=eta*phi,skew=sqrt(eta)*(skew),mu=t(mu))}
  den11=sweep(dg,2,rho,FUN="*")
  den12=sweep(db,2,(1-rho),FUN="*")
  den1=den11+den12
  den = apply(den1, 1, prod)
  return(den)
}





dMSCSAL <- function(data,p,mu=rep(0,p),sigma=diag(p),skew=rep(0,p),rho=(0.5),eta=rep(1.01,p),gam=NULL,phi=NULL,rotated=F) {
  if(!is.matrix(data))data=matrix(data,length(data)/p,p)
  y=data
  
  n=nrow(y)
  if(is.null(n))n=1
  
  
  ### y is the data set
  ## this is the density for 1 component so g is fixed
  # x is a n x p matrix
  
  if(is.null(gam)){gam=eigen(sigma)$vectors
  phi=eigen(sigma)$values
  }else if(is.null(phi)){cat("phi cannot be NULL")}
  if(rotated==F){
    y = data%*%gam
    skew =t(skew%*%gam)#t(skew)# 
    mu =mu%*%gam
  }
  vs=c(0,1)
  mvs= permutations(n=2,r=p,v=vs,repeats.allowed=T)
  dskew=diag(p)
  dvar=diag(p)
  den=numeric(n)
  for(i in 1:n){
    deni=0
    for(j in 1:nrow(mvs)){
      vvs=mvs[j,]
      diag(dskew)=(vvs+(1-vvs)/sqrt(eta))^(-1)
      diag(dvar)=(vvs+(1-vvs)/eta)^(-1)
      molt=prod(rho^vvs*(1-rho)^(1-vvs))
      sigstar=gam%*%dvar%*%(diag(p)*phi)%*%t(gam)
      #skewstar=as.vector(gam%*%dskew%*%skew)
      skewstar=as.vector(dskew%*%skew)
      densal=dSAL(t(y[i,]),p,mu=mu,skew=skewstar,sigma=sigstar)
      deni=deni+molt*densal
    }
    den[i]=deni}
  
  
  return(den)
}

Lbar <- function(n){
  
  Lb      <- matrix(0,(n*(n-1)/2),n^2)
  indcol  <- sapply(0:(n-2), function(i) 
    sapply((i*n+(i+2)):(n+i*n), function(j) j )
  )
  
  c1      <- 1:(n*(n-1)/2)
  c2      <- unlist(indcol)
  ind     <- cbind(c1,c2)
  Lb[ind] <- 1
  
  return(Lb)
  
}

#-------------------------------------#
# Stricly half-vectorization operator #
#-------------------------------------#

bvech <- function(A){   # A: square marrix
  
  n   <- ncol(A)
  res <- Lbar(n)%*%c(A)
  
  return(res)
  
}

#--------------#
# From Q to Pl #
#--------------#

QPl <- function(Q   # Orthogonal matrix
){
  
  luQ    <- Matrix::lu(Q)
  eluQ   <- Matrix::expand(luQ)
  l      <- bvech(as.matrix(eluQ$L))
  U      <- eluQ$U
  P      <- eluQ$P
  return(list(      # It returns:
    P=P,  # Permutation matrix
    l=l   # Vector of n(n-1)/2 elements
  )
  )
  
}


#--------------#
# From Q to PL #
#--------------#

QPL <- function(Q   # Orthogonal matrix
){
  
  luQ    <- lu(Q)
  eluQ   <- expand(luQ)
  L      <- eluQ$L
  U      <- eluQ$U
  P      <- eluQ$P
  return(list(      # It returns:
    P=P,  # Permutation matrix
    L=L   # Unit lower triangular matrix
  )
  )
  
}

#--------------#
# From Pl to Q #
#--------------#

PlQ <- function(P,  # Permutation matrix
                l   # Vector of n(n-1)/2 elements
){
  
  n  <- (1+sqrt(1+8*length(l)))/2
  L  <- diag(n) + matrix(t(Lbar(n))%*%l,n,n) 
  
  Q <- P %*% L %*% solve(qr.R(qr(P%*%L)), tol = 1e-50) 
  
  return(Q)         # It returns an orthogonal matrix
  
}

#--------------#
# From PL to Q #
#--------------#

PLQ <- function(P,  # Permutation matrix
                L   # Unit lower triangular matrix
){
  
  Q <- P %*% L %*% solve(qr.R(qr(P%*%L)), tol = 1e-50) 
  
  return(Q)         # It returns an orthogonal matrix
  
}

#----------#
# Examples #
#----------#

#  d=3

# d <- 3
# 
# X <- matrix(rnorm(d*d),d,d)
# Q <- qr.Q(qr(X))
# 
# P <- QPl(Q)$P
# L <- QPl(Q)$l
# 
# P <- QPL(Q)$P
# L <- QPL(Q)$L
# 
# Q1 <- PLQ(P,L)
# 
# sum(Q-Q1)



library(MixGHD)

##################Densities
#source("~/Dropbox/Multiple Scaled SAL/R/FINAL CODE DEC 22/densities.r")
#source("~/Dropbox/Multiple Scaled SAL/R/FINAL CODE DEC 22/PLRfactorization.r")
#source("~/Dropbox/Multiple Scaled SAL/R/FINAL CODE DEC 22/ALSM_EM.R")
##################################### Initialization ##################################
#######Parameters
## n              number of units
## p              number of variables
##mu              mean p
##skew            skewness p
##Sigma           scale matrix pxp
## gamma          eigenvectors pxp
## phi            eignevalues p
## rho            degree of contamination  p
## eta            contamination  p

iniMSCMSSAL<-function(data,G,gam){
  p=ncol(data)
  l=kmeans(data,G)
  pi=rep(1/G,G)
  mu=l$centers
  ini=l$cluster
  z=matrix(0,length(ini),max(ini))
  for(i in 1:length(ini)){
    z[i,ini[i]]=1}
  rho=matrix(0.9,G,p)
  eta=matrix(2,G,p)
  skew=matrix(1,G,p)
  
  sigma=array(diag(p),dim=c(p,p,G))
  # gam=array(1,dim=c(p,p,G))##array(eigen(sigma[,,1])$vectors,dim=c(p,p,G))
  phi=matrix(1,G,p,1)#eigen(sigmaO)$values
  for(g in 1:G){
    dataS=data[z[,g]==1,]
    #gam[,,g]=eigen(cov(dataS))$vectors
    # phi[g,]=eigen(cov(data))$values
    mu[g,]=mu[g,]%*%gam[,,g]
    skew[g,]=skew[g,]%*%gam[,,g]
  }
  
  par=list(mu=mu,skew=skew,gam=gam, phi=phi,sigma=sigma,rho=rho,eta=eta,pi=pi)
  return(par)
}


##################################### Convergence ##################################

getall <- function(loglik) {
  if (length(loglik) <3) stop("must have at least 3 likelihood values")
  n = length(loglik)
  lm1 = loglik[n]
  lm  = loglik[(n-1)]
  lm_1  = loglik[(n-2)]
  am = (lm1 - lm)/(lm - lm_1)
  lm1.Inf = lm + (lm1 - lm)/(1-am)
  val = lm1.Inf - lm
  if (is.nan(val)) val=0
  if (val < 0) val= 1
  return( val )
}

##################################### Log likelihood ##################################

llikMSCMSS <- function(data=NULL, gpar=NULL, label=NULL) {
  fx = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
  p=ncol(data)
  for (k in 1:length(gpar$pi) ) fx[,k] = dMSCMSSAL(data=data%*%gpar$gam[,,k], p=p,mu=gpar$mu[k,],eta=gpar$eta[k,],rho=gpar$rho[k,],skew=gpar$skew[k,],gam=gpar$gam[,,k],phi=gpar$phi[k,],rotated = T)
  
  lablik = matrix(1, nrow=nrow(data), ncol=length(gpar$pi))
  if (!is.null(label)) lablik = combinewk(lablik, label= label)
  
  val = apply(fx*lablik, 1, weighted.sum, wt=gpar$pi)
  val = sum(log(val))
  
  return(val)
}

combinewk <- function(weights=NULL, label=NULL)	{
  # known is a numeric with 
  # 0 if unknown group membership 
  # 1,2,3,.. for label of known group
  if (is.null(label)) stop('label is null')
  kw     = label !=0
  for (j in 1:ncol(weights)) weights[kw,j] = (label == j)[kw]
  return(weights)	
}

weighted.sum <- function(z,wt,...) return( sum(z*wt,...) )

######################################EM########################

EstepZMSCMSSAL <- function(data=NULL, gpar=NULL) {
  G = length(gpar$pi)
  p=ncol(data)
  if (G > 1) {
    fx = matrix(0, nrow=nrow(data), ncol=length(gpar$pi))
    for (k in 1:G ) fx[,k] = dMSCMSSAL(data=data%*%gpar$gam[,,k], p=p,mu=gpar$mu[k,],eta=gpar$eta[k,],rho=gpar$rho[k,],skew=gpar$skew[k,],gam=gpar$gam[,,k],phi=gpar$phi[k,],rotated=T)
  }
  x=sweep(fx,2,gpar$pi,FUN="*")
  xs=apply(x,1,sum)
  x=x/xs
  return(x)
}

EstepVMSCMSSAL<-function(data=NULL,gpar=NULL){
  G = length(gpar$pi)
  n=nrow(data)
  p=ncol(data)
  rho=gpar$rho
  v=array(0,dim=c(n,p,G))
  gam=gpar$gam
  for(k in 1:G){
    dd=data%*%gam[,,k]
    for(j in 1:p){
      num <- dMSSAL(data=dd[,j], p=1,mu=gpar$mu[k,j],skew=gpar$skew[k,j],gam=gam[j,j,k],phi=gpar$phi[k,j],rotated=T,PARrotated=T)
      den <- dMSCMSSAL(data=dd[,j], p=1,mu=gpar$mu[k,j],eta=gpar$eta[k,j],rho=gpar$rho[k,j],skew=gpar$skew[k,j],gam=gam[j,j,k],phi=gpar$phi[k,j],rotated=T)
      v[,j,k] <- (rho[k,j]*num)/den
    }
  }
  
  v[is.nan(v)] <- 0 # to avoid 0/0 in v
  return(v)
}



gigS <- function(a=NULL,b=NULL) {
  p=ncol(b)
  n=nrow(b)
  av=a
  bv=b
  val1=matrix(0,n,p)
  val2=matrix(0,n,p)
  for(j in 1:p){
    a=av[j]
    b=bv[,j]
    # returns a matrix with dim length(a) x 3
    t1 = exp( (log(a) + log(b))/2 );
    kv1 = log(besselK( t1, nu=1.5, expon.scaled =TRUE))-t1
    kv  =log( besselK( t1, nu=0.5, expon.scaled =TRUE))-t1
    kv12 = kv1-kv
    
    sb.a = (log(b)-log(a))/2;
    v1=sb.a+kv12
    val1[,j]=exp(v1)
    val2[,j]=exp( -1*sb.a+kv12 ) - sign(.5)*exp( log(2)+log(abs(.5)) - log(b));
    # logw = log(sb.a) + grad( logbesselKvFA, x=rep(v,length(sab)), y=sab, method="Richardson",  method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
  }
  
  val =list(w= val1,invw= val2)#,logw)
  return(val)
}
######GIG MSSAL
gigSAL <- function(x=NULL, mu=NULL,skew=NULL,phi=NULL) {
  # returns the same as gig
  
  a1 = 2 + (skew^2)*(1/phi) 
  xmu= 	sweep(x,2,mu, "-")
  b1=sweep(xmu^2,2,1/phi, FUN="*")
  
  b1[b1<0.00000000000000000001]=0.00000000000000000001
  
  val = gigS(b=b1,a=a1)
  return(val)
}
gigSALC <- function(x=NULL, mu=NULL,skew=NULL,phi=NULL,eta=NULL) {
  # returns the same as gig
  
  
  
  a1 = 2 + (skew^2/phi )
  xmu= 	sweep(x,2,mu, "-")
  b1=sweep(xmu^2,2,1/(phi*eta), FUN="*")
  b1[b1<0.00000000000000000001]=0.00000000000000000001
  
  val = gigS(b=b1,a=a1)
  return(val)
}

###### core EM

EMMSCMSAL<-function(data,gpar){
  
  n=nrow(data)
  p=ncol(data)
  pi=gpar$pi
  G=length(pi)
  mu=gpar$mu
  skew=gpar$skew
  gam=gpar$gam
  phi=gpar$phi
  rho=gpar$rho
  eta=gpar$eta
  
  ######E step#######
  
  if(G>1){z=EstepZMSCMSSAL(data,gpar)
  }else{z=matrix(1,n,1)}
  v=EstepVMSCMSSAL(data,gpar)
  for(g in 1:G){
    
    x=data%*%gam[,,g]
    expected=gigSAL(x,mu=mu[g,],skew=skew[g,],phi=phi[g,])
    expectedC=gigSALC(x,mu=mu[g,],skew=skew[g,],phi=phi[g,],eta=eta[g,])
    wijg=expected$w
    invwijg=expected$invw
    wijgB=expectedC$w
    invwijgB=expectedC$invw
    
    aig=sweep(v[,,g],1,z[,g],FUN="*")
    big=aig*invwijg
    cig=sweep((1-v[,,g]),1,z[,g],FUN="*")
    dig=cig*invwijgB
    fig=aig*wijg
    hig=cig*wijgB
    
    A=apply(fig+hig,2,sum)
    dovereta=sweep(dig,2,eta[g,],FUN="/")
    coversqrteta=sweep(cig,2,sqrt(eta[g,]),FUN="/")
    B=apply(big+ dovereta,2,sum)
    C=apply(aig+coversqrteta,2,sum)
    
    ##### M-step
    
    pi.new = mean(z[,g])
    rho.new = pmax(apply(aig,2,sum)/sum(z[,g]),0.5)
    
    p1=apply(x*(big+dovereta),2,sum)
    p2=apply(x*(aig+coversqrteta),2,sum)
    
    ### mu
    
    mu.new=(A*p1-C*p2)/(A*B-C*C)
    
    ### skew
    
    skew.new=(B*p2-C*p1)/(A*B-C*C)
    
    ### phi
    
    xmu= 	sweep(x,2,mu.new, "-")
    xmus= 	sweep(xmu,2,skew.new, "*")
    p1=apply(xmu^2*(big+dovereta),2,sum)
    p2=2*apply(xmus*(aig+coversqrteta),2,sum)
    p3=skew.new^2*A
    phi.new=(p1-p2+p3)/sum(z[,g])
    
    ##### eta 
    
    a=apply(sweep(cig,2,phi.new,FUN="*"),2,sum)
    b=apply(xmus*cig,2,sum)
    c=-apply(dig*xmu^2,2,sum)
    
    y1=((-b+sqrt(b^2-4*a*c))/(2*a))^2
    y2=((-b-sqrt(b^2-4*a*c))/(2*a))^2
    
    gpar$pi[g]=pi.new
    gpar$mu[g,]=mu.new
    gpar$skew[g,]=skew.new
    gpar$phi[g,]=phi.new
    
    gpar$rho[g,]=rho.new
    eta.new=numeric(p)
    for(j in 1:p){
      gpar$eta[g,j]=y1[j]
      l1=llikMSCMSS(data,gpar)
      gpar$eta[g,j]=y2[j]
      l2=llikMSCMSS(data,gpar)
      if(l1>l2){eta.new[j]=y1[j]}
      else if(l1==l2){eta.new[j]=min(y1[j],y2[j])}
      else{eta.new[j]=y2[j]}
      eta.new=pmax(eta.new,1.01)
    }
    gpar$eta[g,]=eta.new
    
    #######old gamma##########3
    #   #  gam.new=updategamF1(gam0=gam[,,g], data=data,  phig=phi[g,], skewg=skew[g,], mug=mu[g,], zg=z[,g],  big,dig,aig,cig,rho[g,]) 
    #    #  gam.new=updategamF2(gam0=gam.new, data=data,  phig=phi[g,], skewg=skew[g,], mug=mu[g,], zg=z[,g],  big,dig,aig,cig,rho[g,]) 
    #   # gam[,,g]=gam.new
    
  }
  
  return(gpar)
}

###### Main function given gamma

mainMSCMSAL<-function(data,gam,G=2,max.iter=1000,eps=0.1){
  
  data=as.matrix(data)
  gpar=iniMSCMSSAL(data,G,gam)
  
  n=nrow(data)
  p=ncol(data)
  
  loglik=NULL
  for(i in 1:3){
    gpar=EMMSCMSAL(data, gpar)
    
    loglik=c(loglik,llikMSCMSS(data,gpar))
    if(i>2){if(loglik[i]<loglik[i-1])cat('Decreasing llik, difference between llik',loglik[i-1]- loglik[i])}
  }
  while(( getall(loglik[1:i]) > eps) & i<max.iter){
    i=i+1
    gpar=EMMSCMSAL(data, gpar)
    
    loglik=c(loglik,llikMSCMSS(data,gpar))
    if(i>2){if(loglik[i]<loglik[i-1])cat('Decreasing llik, difference between llik',loglik[i-1]- loglik[i])}
  }
  if(G>1){z=EstepZMSCMSSAL(data,gpar)
  }else{z=matrix(1,n,1)}
  v=EstepVMSCMSSAL(data,gpar)
  
  realp=gpar
  
  for(g in 1:G){
    # ord=order(gpar$phi[g,])
    # realp$phi[g,]=realp$phi[g,ord]
    # realp$gam[,,g]=realp$gam[,ord,g]
    # realp$rho[g,]=realp$rho[g,ord]
    # realp$eta[g,]=realp$eta[g,ord]
    # realp$mu[g,]=realp$mu[g,ord]
    # realp$skew[g,]=realp$skew[g,ord]
    realp$mu[g,]=realp$mu[g,]%*%t(realp$gam[,,g])
    realp$skew[g,]=realp$skew[g,]%*%t(realp$gam[,,g])
    realp$sigma[,,g]=realp$gam[,,g]%*%(diag(p)*realp$phi[g,])%*%t(realp$gam[,,g])
    realp$rho[g,which(realp$eta[g,]<1.02)]=1
  }
  cluster <- apply(z,1,which.max)
  detect <- array(1/n,c(n,p),dimnames=list(1:n,paste("dim.",1:p,sep="")))
  for(h in 1:p){
    for(i in 1:n){
      detect[i,h] <- ifelse(v[i,h,cluster[i]] > 1/2,0,1)
    }
  }
  BIC=2 * max(loglik) - ((G-1)+G*(4*p+p*(p+1)/2)) * log(n)
  return(list(detect=detect,gpar=gpar,realp=realp,ll=loglik,BIC=BIC))
  
}



