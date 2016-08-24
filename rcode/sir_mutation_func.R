
##################################################################################
##################################################################################
sir_mutation_agent = function(N         # population size
                              ,I_0       # initial number infected
                              ,S_0       # initial number susceptible
                              ,gamma     # recovery rate in days^{-1}
                              ,R0        # reproduction number
                              ,tbeg      # begin time of simulation
                              ,tend      # end time of simulation
                              ,delta_t=0 # time step (if 0, then dynamic time step is implemented)
                              ,seq_len   # the number of bases in the sequence
                              ,alphabet  # DNA code 
                              ,mut_rate # mutation rate for time step 
){
  
  ########################################################################
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
  # set up the state vector of the population
  # vstate = 0   means susceptible
  # vstate = 1   means infected    
  # vstate = 2   means recovered   
  # randomly pick I_0 of the people to be infected
  ########################################################################
  
  vstate = rep(0,N) # keeps track of the state of the individuals 
  index.inf = sample(N,I_0) # randomly pick I_0 people from N
  vstate[index.inf] = 1     # these are the infected people
  
  ########################################################################
  # set up parameters for the mutation model and initalize with the same strain
  #strain.state = as.data.frame(matrix(nrow = N, ncol = seq_len, dimnames = list(c(seq(1:N)), paste("base", seq(1:seq_len),sep = "_"))))
  strain.state = numeric(length = N)
  strain.state[index.inf] <- 1
  
  base.haplotype = as.vector(rep(1, seq_len))
  #strain.state[index.inf] <- base.haplotype
  
  
  ########################################################################
  ## Creating dataframes to hold individual trial analysis 
  diversity <- numeric(0)
  divergence <- numeric(0)
  cum.strains <- numeric(0)# 
  circulating.strains <- numeric(0)
  number.mutations <- numeric(0)
  
  #population <- init_population(times = tend/delta_t, base.haplotype = base.haplotype)
  #population <- init_population(base.haplotype)
  
  population.strains <- init_population_strain(base.haplotype)
  population.frequency <- init_population_freqdf()
  
  ########################################################################
  # now begin the time steps
  ########################################################################
  t     = tbeg
  S     = S_0
  I     = I_0
  vS    = numeric(0)  # this will be filled with the number of susceptibles in the population over time
  vI    = numeric(0)  # this will be filled with the number of infectives in the population over time
  vtime = numeric(0)  # this will be filled with the time steps
  #recover.prob = (1-exp(-gamma*delta_t)) # the probability an infective recovers in this time step, # this probabaly doesn't have to be set each time? 
  
  
  # Main Simulation 
  
  while (t<tend&I>0) { # continue the simulation until we have no more infectious people, or t>=tend
    
    ############################################################
    # Book-keeping 
    ############################################################
    S = length(vstate[vstate==0])  # count the number of susceptibles, based on the state vector 
    I = length(vstate[vstate==1])  # count the number of infectives, based on the state vector
    
    vS = append(vS,S)  # append this info to vectors that we will return from the function
    vI = append(vI,I)
    
    vtime = append(vtime, t)
    
    # Mutation book-keeping 
    #current.haplotypes <- get_current(population, t+1) # Strains that have above 0 frequency 
    current.haplotypes <- get_current(population.strains, population.frequency,timestep =  length(vtime))
    diversity.time  =  get_diversity(current.haplotypes)
    divergence.time  = get_divergence(current.haplotypes, base.haplotype) 
    
    circulating.strains = append(circulating.strains, length(current.haplotypes$hindex))
    cum.strains = append(cum.strains, length(population.strains))
    diversity = append(diversity, diversity.time)
    divergence = append(divergence, divergence.time)
    
    deltat=delta_t
    if (delta_t==0){ # this is the calculation of a dynamic time step
      deltat = 1/(beta*I*S/N + gamma*I)
    }
    
    recover.prob = (1-exp(-gamma*deltat))
    
    
    ############################################################
    # Mutation 
    ############################################################
    
    # Current Infection Index
    current.inf.index = which(vstate == 1)
    
    #get expected number of mutations assuming one sequence per individual
    
    mean.mutation <- mut_rate*seq_len*I
    num.mutations <- rpois(1, mean.mutation)
    number.mutations <- append(number.mutations, num.mutations)
    
    if (num.mutations > 0) {    
      # choose indices from Infecteds
      
      mut.current.index = sample(current.inf.index, num.mutations, replace = FALSE) #Same strain can't mutate twice  
      
      for (j in 1:num.mutations) { # for each of the mutations
        #infected_strain = strain.state[mut.current.index[j],]
        #new.strain = mutate.strain(strain = infected_strain)
        #strain.state[mut.current.index[j],] = new.strain
        
        original.strain.index = strain.state[mut.current.index[j]] # Get hindex of the original strain
        original.strain = population.strains[[original.strain.index]]  # Get actual strain
        new.strain = mutate.strain(strain = original.strain)
        
        new.strain.index <- length(population.strains)+1
        #print(paste("new.strain.index", new.strain.index, sep = "."))
        
        population.frequency[new.strain.index] <- numeric()
        # print(head(population.frequency))
        
        #Adding the strain to the population 
        population.strains[[new.strain.index]] <- new.strain
        #population.strains
        names(population.strains) <- paste("s", seq_along(population.strains), sep = ".") #setting up all the names everytime, not efficient, can do this at the end 
        #update strain index of population
        strain.state[mut.current.index[j]] <- new.strain.index
      }
    } 
    
    
    ############################################################
    # Replication  
    ############################################################
    # sample Poisson random numbers of infected people contacted by each person
    avg.num.infected.people.contacted = beta*I*deltat/N
    vnum.infected.people.contacted = rpois(N,avg.num.infected.people.contacted) 
    vprob = runif(N)   # sample uniform random numbers
    vnewstate = vstate # copy the state vector to a temporary vector used for calculations
    
    # Infected people recover if the sampled uniform random
    # number is less than the recovery probability
    vnewstate[vstate==1&vprob<recover.prob] = 2   
    recovered_ind = which(vstate==1&vprob<recover.prob)
    
    # If a susceptible contacted at least one infective, they are infected
    vnewstate[vstate==0&vnum.infected.people.contacted>0] = 1 
    new.infects.index = which(vstate==0&vnum.infected.people.contacted>0)
    num.new.infects = length(new.infects.index)
    
    vstate = vnewstate # update the state vector
    
    # Assign New Strains    
    for(k in 1:num.new.infects) {
      responsible.individual = sample(current.inf.index, 1)
      responsible.strain <- strain.state[responsible.individual] 
      strain.state[new.infects.index[k]] <- responsible.strain
    }
    
    ## Recovered Individual Strains are a dead-end
    strain.state[recovered_ind] <- 0
    
    I = length(vstate[vstate==1])
    if (I == 0) {
      return
    } else {
      ## Add any new strains to population and calculate overall frequencies
      unique.strains <- unique(strain.state[strain.state > 0])
      
      current.freq =  data.frame(aaply(.data = unique.strains, .margins = 1, function(x) 
        freq = length(which(strain.state == x))/length(which(strain.state > 0)))); colnames(current.freq) = "frequency"
      
      current.freq = cbind(unique.strains, current.freq)
      
      for (n in 1:nrow(current.freq)) {
        population.frequency[length(vtime)+1, current.freq[n,1]] = current.freq$frequency[n]
      }
      
      t = t + deltat # update the time
    }
  } 
  
  
  final = 0
  if (length(vS)>0) final = 1-min(vS)/N
  
  return(list(time=vtime,I=vI,S=vS,final_size=final, strain.freq = population.frequency, divergence = divergence, diversity = diversity, 
              ciruclating.strains = circulating.strains, cumulative.strains = cum.strains, number.mutations= number.mutations, population.strains = population.strains))
}

