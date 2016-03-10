require(deSolve)
require(tictoc)

#basic SIR 
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    # dcumI <- beta*S*I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}


#basic SIR with vaccination
sir_vacc <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I - nu*S
    dI <- beta * S * I - gamma * I
    # dcumI <- beta*S*I
    dR <- gamma * I+nu*S
    
    return(list(c(dS, dI, dR)))
  })
}