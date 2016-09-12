##################################################################################
##################################################################################
# An R script to perform a stochastic epidemic simulation using an Agent Based model
# (with homogenous mixing equivalent to an SIR model)
#
# Author:  Sherry Towers
#          smtowers@asu.edu
# Created: Feb 13, 2016
#
# Copyright Sherry Towers, 2016
#
# This script is not guaranteed to be free of bugs and/or errors
#
# This script can be freely used and shared as long as the author and 
# copyright information in this header remain intact.
##################################################################################

##################################################################################
##################################################################################
##################################################################################
SIR_agent = function(N         # population size
                     ,I_0       # initial number infected
                     ,S_0       # initial number susceptible
                     ,gamma     # recovery rate in days^{-1}
                     ,R0        # reproduction number
                     ,tbeg      # begin time of simulation
                     ,tend      # end time of simulation
                     ,delta_t=0 # time step (if 0, then dynamic time step is implemented)
){
  
  ########################################################################
  # begin by reverse engineering the transmission rate, beta, from
  # R0 and gamma
  #
  # beta is the number of contacts sufficient to transmit infection per unit time
  #
  # Thus, in time delta_t, a susceptible person will contact a Poisson random
  # number of infected people, with mean beta*I*delta_t/N
  #
  # The probability an infected person will recover in time
  # delta_t is p=1-exp(-gamma*delta_t)
  # (note that p does not depend on t!  This is due to the 
  #  memoryless nature of the Exponential distribution)
  ########################################################################
  beta = R0*gamma 
  
  ########################################################################
  # now set up the state vector of the population
  # vstate = 0   means susceptible
  # vstate = 1   means infected    
  # vstate = 2   means recovered   
  # lstate = types of strains 
  # randomly pick I_0 of the people to be infected
  ########################################################################
  vstate = rep(0,N)
  index_inf = sample(N,I_0) # randomly pick I_0 people from N
  vstate[index_inf] = 1     # these are the infected people

  ########################################################################
  # now begin the time steps
  ########################################################################
  t     = tbeg
  S     = S_0
  I     = I_0
  vS    = numeric(0)  # this will be filled with the number of susceptibles in the population over time
  vI    = numeric(0)  # this will be filled with the number of infectives in the population over time
  vtime = numeric(0)  # this will be filled with the time steps
  #recover_prob = (1-exp(-gamma*delta_t)) # the probability an infective recovers in this time step, # this probabaly doesn't have to be set each time? 
  
  while (t<tend&I>0){ # continue the simulation until we have no more infectious people, or t>=tend
    S = length(vstate[vstate==0])  # count the number of susceptibles, based on the state vector 
    I = length(vstate[vstate==1])  # count the number of infectives, based on the state vector
    
    vS = append(vS,S)  # append this info to vectors that we will return from the function
    vI = append(vI,I)
    
    vtime = append(vtime,t)
    
    #cat(t,S,I,"\n")
    
    deltat=delta_t
    if (delta_t==0){ # this is the calculation of a dynamic time step
      deltat = 1/(beta*I*S/N + gamma*I)
    }
  
    recover_prob = (1-exp(-gamma*deltat))
    
    # sample Poisson random numbers of infected people contacted by each person
    avg_num_infected_people_contacted = beta*I*deltat/N
    vnum_infected_people_contacted = rpois(N,avg_num_infected_people_contacted) 
    vprob = runif(N)   # sample uniform random numbers
    vnewstate = vstate # copy the state vector to a temporary vector used for calculations
    
    # Infected people recover if the sampled uniform random
    # number is less than the recovery probability
    vnewstate[vstate==1&vprob<recover_prob] = 2              
    
    # If a susceptible contacted at least one infective, they are infected
    vnewstate[vstate==0&vnum_infected_people_contacted>0] = 1 
    
    vstate = vnewstate # update the state vector
    
    t = t + deltat     # update the time
  } 
  final = 0
  if (length(vS)>0) final = 1-min(vS)/N
  
  return(list(time=vtime,I=vI,S=vS,final_size=final))
  
}


##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
SIR_agent_erlang = function(N         # population size
                            ,I_0       # initial number infected 
                            ,S_0       # initial number susceptible
                            ,gamma     # recovery rate in days^{-1}
                            ,k         # shape parameter of the erlang distribution
                            ,R0        # reproduction number
                            ,tbeg      # begin time of simulation
                            ,tend      # end time of simulation
                            ,delta_t=0 # time step (if 0, then dynamic time step is implemented)
){
  
  cat(N,I_0,S_0,gamma,k,R0,tbeg,tend,delta_t,"\n")
  ########################################################################
  # begin by reverse engineering the transmission rate, beta, from
  # R0 and gamma
  #
  # beta is the number of contacts sufficient to transmit infection per unit time
  #
  # Thus, in time delta_t, a susceptible person will contact a Poisson random
  # number of infected people, with mean beta*I*delta_t/N
  #
  # The probability an infected person will recover in time
  # delta_t is p=1-exp(-gamma*delta_t)
  # (note that p does not depend on t!  This is due to the 
  #  memoryless nature of the Exponential distribution)
  ########################################################################
  beta = R0*gamma 
  
  ########################################################################
  # now set up the state vector of the population
  # vstate = 0   means susceptible
  # vstate = 1   means infected    
  # vstate = 2   means recovered   
  # randomly pick I_0 of the people to be infected
  ########################################################################
  vstate = rep(0,N)
  index_inf = sample(N,sum(I_0)) # randomly pick I_0 people from N
  vstate[index_inf] = 1     # these are the infected people
  
  ########################################################################
  # randomly set up the sojourn time spent so far in each state.
  # when an individual leaves a state, the sojourn time gets reset to 0
  ########################################################################
  vsojourn = rep(0,N)    
  
  ########################################################################
  # now begin the time steps
  ########################################################################
  t     = tbeg
  S     = S_0
  I     = sum(I_0)
  vS    = numeric(0)  # this will be filled with the number of susceptibles in the population over time
  vI    = numeric(0)  # this will be filled with the number of infectives in the population over time
  vtime = numeric(0)  # this will be filled with the time steps
  
  while (t<tend&I>0){ # continue the simulation until we have no more infectious people, or t>=tend
    S = length(vstate[vstate==0])  # count the number of susceptibles, based on the state vector 
    I = length(vstate[vstate==1])  # count the number of infectives, based on the state vector
    
    vS = append(vS,S)  # append this info to vectors that we will return from the function
    vI = append(vI,I)
    
    vtime = append(vtime,t)
    
    #cat(t,S,I,"\n")
    
    deltat=delta_t
    if (delta_t==0){ # this is the calculation of a dynamic time step
      deltat = 1/(beta*I*S/N + gamma*I)
    }
    
    # the probability an infective recovers in this time step
    theta = 1/(gamma*k) # scale factor for the gamma distribution, this never changes 
    a = 1-pgamma(vsojourn,shape=k,scale=theta) # the probability that they are still in the compartment
    
    recover_prob = rep(1,N) # just setting up the vector to remember 
    l = which(a>1e-4) # not sure if this ever changes?
    
    recover_prob[l] = (pgamma(vsojourn[l]+deltat,shape=k,scale=theta)-pgamma(vsojourn[l],shape=k,scale=theta))/a[l]
      #getting the probability of delta t

    
    # sample Poisson random numbers of infected people contacted by each person
    avg_num_infected_people_contacted = beta*I*deltat/N
    vnum_infected_people_contacted = rpois(N,avg_num_infected_people_contacted) 
    vprob = runif(N)   # sample uniform random numbers
    vnewstate = vstate # copy the state vector to a temporary vector used for calculations
    
    vsojourn = vsojourn+deltat # adding the time step to track how long each person has been in a sate 
    # Infected people recover if the sampled uniform random number is less than the recovery probability
    vnewstate[vstate==1&vprob<recover_prob] = 2              
    # reset the time they have spent in their new state to zero!
    vsojourn[vstate==1&vprob<recover_prob] = 0              
    
    # If a susceptible contacted at least one infective, they are infected
    vnewstate[vstate==0&vnum_infected_people_contacted>0] = 1 
    # reset the time they have spent in their new state to zero!
    vsojourn[vstate==0&vnum_infected_people_contacted>0] = 0 
    
    vstate = vnewstate # update the state vector
    
    t = t + deltat     # update the time
  } 
  final = 0
  if (length(vS)>0) final = 1-min(vS)/N
  
  return(list(time=vtime,I=vI,S=vS,final_size=final))
  
}



