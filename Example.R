################################################################
######## Multiple Scaled Contaminated Asymmetric Laplace #######
###############################################################
## Tortora C., Franczak B.C., Bagnato L., Punzo A.
## A Laplace-based model with flexible tail behavior
##

## Example

source("MSCAL.R")

### Parameter generation

eta=c(6,10)
rho=c(.7,.9)
skew=c(-2,6)
n=500
sigma=matrix(-0.8,2,2)
diag(sigma)=c(4,1)
ed=eigen(sigma)
gam=ed$vectors
phi=ed$values

## Data generation

data=rMSCAL(n,p=2,skew=skew,rho=rho,eta=eta,gam=gam, phi=phi)
plot(data)

##Running the algorithm using the MCEM algorithm
##n.trailas set to 100 for speed, please consider higher n.trials, e.g. 2000
res_Mon=MSCAL(data=data,n.trials=100)

### BIC
res_Mon$out$BIC

###parameters
res_Mon$out$realp

