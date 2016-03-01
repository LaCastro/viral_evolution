##################################################################################
# An R script to perform a stochastic epidemic simulation of
# an SIR model using Stochastic Differential Equations
#
# Author:  Sherry Towers
#          admin@sherrytowers.com
# Created: Dec 4, 2012
#
# Copyright Sherry Towers, 2012
#
# This script is not guaranteed to be free of bugs and/or errors
#
# This script can be freely used and shared as long as the author and
# copyright information in this header remain intact.
##################################################################################
rm(list = ls(all = TRUE))  # resets R to fresh

set.seed(166829)
require("sfsmisc")
require("deSolve")

niter = 50
ldynamic_time_step = 0
##################################################################################
##################################################################################
##################################################################################
SIRfunc=function(t, x, vparameters){
  S = x[1]  # the value of S at time t
  I = x[2]  # the value of I at time t
  R = x[3]  # the value of R at time t
  
  if (I<0) I=0 # this is a cross check to ensure that we always have sensical values of I
  
  with(as.list(vparameters),{
    npop = S+I+R   # the population size is always S+I+R because there are no births or deaths in the model
    dS = -beta*S*I/npop            # the derivative of S wrt time
    dI = +beta*S*I/npop - gamma*I  # the derivative of I wrt time
    dR = +gamma*I                  # the derivative of R wrt time
    out = c(dS,dI,dR)
    list(out)
  })
}


##################################################################################
# set up the model parameters, and some initial conditions at time t=0
##################################################################################
gamma = 1/3         # recovery period in days^{-1}
R0    = 1.5         # R0 of R of the disease
beta = gamma*R0

N_0 = 1000  # population size
I_0 = 25  # number intially infected people in the population
S_0 = N_0-I_0
R_0 = 0

vt = seq(0,1000,1)
vparameters=c(gamma=gamma,beta=beta)
inits=c(S=S_0,I=I_0,R=R_0)
sirmodel = as.data.frame(lsoda(inits, vt, SIRfunc, vparameters))

##################################################################################
##################################################################################
# Perform niter realisations of the stochastic model
# The list object, zstate, will have all the S,I, and R values over
# time for each realisation, along with the time, and realisation number
# (such that each element of the list is a vector of length 5)
#
# Using a list to store this information, instead of appending to a data frame at
# each step, is computationally more expedient
##################################################################################

delta_t = 0.1
zstate = list()
i = 1
for (iter in 1:niter){
  time = 0
  vstate = c(S_0,I_0,R_0)
  
  K = length(vstate)  # number of compartments
  J = 2               # number of possible state changes
  lambda = matrix(0,nrow=J,ncol=length(vstate))
  lambda[1,] = c(-1,1,0)
  lambda[2,] = c(0,-1,1)
  
  while(vstate[2]>0&vstate[1]>0){
    zstate[[i]] = c(vstate,time,iter)
    S = vstate[1]
    I = vstate[2]
    R = vstate[3]
    N = S+I+R
    vec_p = c(beta*S*I/N
              ,gamma*I)
    
    scale = 10
    if (ldynamic_time_step) delta_t = scale/sum(vec_p)
    
    G = lambda
    for (irow in 1:nrow(G)){
      G[irow,] = G[irow,]*sqrt(vec_p[irow])
    }
    G = t(G)
    
    W = rnorm(ncol(G))
    delta_S = delta_t*(-beta*S*I/N)         + sqrt(delta_t)*sum(G[1,]*W)
    delta_I = delta_t*(+beta*S*I/N-gamma*I) + sqrt(delta_t)*sum(G[2,]*W)
    delta_R = delta_t*(+gamma*I)            + sqrt(delta_t)*sum(G[3,]*W)
    
    vstate[1] = vstate[1] + delta_S 
    vstate[2] = vstate[2] + delta_I 
    vstate[3] = vstate[3] + delta_R 
    
    
    vstate[vstate<0] = 0
    i = i+1
    time = time + delta_t
  }
  cat("Doing realisation:",iter,niter," ",time,vstate,"\n")
}

##################################################################################
##################################################################################
# the sapply functions here extract the elements of the zstate list into vectors
##################################################################################
par(mfrow=c(2,2))
vS = sapply(zstate, "[[", 1)
vI = sapply(zstate, "[[", 2)
vR = sapply(zstate, "[[", 3)
vtime = sapply(zstate, "[[", 4)
viter = sapply(zstate, "[[", 5)

##################################################################################
##################################################################################
# now lets plot the results for each iteration, all overlaid on the same
# plot for S, I and R
##################################################################################
mult.fig(4)
for (iter in 1:niter){
  l = which(viter==iter)
  if (iter==1){
    plot(vtime[l],vS[l]/N_0,type="l",xlab="Time, in Days",ylab="S",main=paste("N=",N_0," and I_0=",vI[1],sep=""),ylim=c(0,1),xlim=c(0,max(vtime)))
  }else{
    lines(vtime[l],vS[l]/N_0,type="l")
  }
}
lines(sirmodel$time,sirmodel$S/N_0,col=2,lwd=2)
legend("bottomleft",legend=c("SDE","Deterministic"),col=c(1,2),lwd=4,bty="n",cex=0.8)

for (iter in 1:niter){
  l = which(viter==iter)
  if (iter==1){
    plot(vtime[l],vI[l]/N_0,xlab="Time, in Days",ylab="I",main=paste("N=",N_0," and I_0=",vI[1],sep=""),ylim=c(0,1.1*max(vI/N_0)),xlim=c(0,max(vtime)),type="l")
  }else{
    lines(vtime[l],vI[l]/N_0,type="l")
  }
}
lines(sirmodel$time,sirmodel$I/N_0,col=2,lwd=2)

##################################################################################
##################################################################################
# Let's also histogram the final size from the realisations.
# In this particular model, we can assess that through the
# maximum value of R for each realisation
##################################################################################
vfinal = numeric(0)
for (iter in 1:niter){
  l = which(viter==iter)
  vfinal = append(vfinal,max(vR[l]/N_0))
  if (iter==1){
    plot(vtime[l],vR[l]/N_0,xlab="Time, in Days",ylab="R",main=paste("N=",N_0," and I_0=",vI[1],sep=""),ylim=c(0,1.0),type="l",xlim=c(0,max(vtime)))
  }else{
    lines(vtime[l],vR[l]/N_0,type="l")
  }
}
lines(sirmodel$time,sirmodel$R/N_0,col=2,lwd=2)

hist(vfinal,breaks=seq(0,1,0.05),xlab="Final size",main="Final size")
Rmax = max(sirmodel$R)/N_0
lines(c(Rmax,Rmax),c(0,1e6),col=2,lty=3,lwd=3)
legend("topleft",legend=c("SDE realisations","Deterministic model"),col=c(1,2),lwd=4,bty="n",cex=0.65)

##################################################################################
##################################################################################
# See section 3.6.1 of Mathematical Epidemiology, edited by
# Brauer, van den Driessche, and Wu
##################################################################################
cat("The expected probability of outbreak is = ",1-(1/R0)^I_0,"\n")







