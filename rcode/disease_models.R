# Disease Model Scripts
require(deSolve)


#Basic SIR model 
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta*S*I
    dI <- beta*S*I-gamma*I
    #dcumI <- beta*S*I
    dR <- gamma*I
    
    return(list(c(dS, dI, dR)))
  })
}


##################################################################################
##################################################################################
# this is a function which, given a value of S,I and R at time t
# calculates the time derivatives of S I and R
# vparameters contains the parameters of the model, like the
# recovery period, gamma, and the transmission rate, beta
# this function gets passed to the deSolve package
##################################################################################

SIRfunc=function(t, x, vparameters){
  S = x[1]  
  I = x[2]  
  R = x[3]  
  if (I<0) I=0 
  
  with(as.list(vparameters),{
    npop = S+I+R   
    dS = -beta*S*I/npop            
    dI = +beta*S*I/npop - gamma*I  
    dR = +gamma*I                  
    out = c(dS,dI,dR)
    list(out)
  })
}




# basic SIR model with vaccination 

sir_vacc <- function(time, state, parameters) {
  with(as.list(c(state,parameters)), {
    dS <- -beta*S*I-nu*S
    dI <- beta*S*I-gamma*I
    #dcumI <- beta*S*I
    dR <- gamma*I+nu*S
    
    return(list(c(dS, dI, dR)))
  })
}