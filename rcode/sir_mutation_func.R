#N = 1000
#I_0 = .01*N
#S_0 = N-I_0
#delta_t = 1 
#tbeg = 1
#tend = 250
#gamma = 1/3
#R0 = 1.5
#beta = R0*gamma
#contact.per.day = 4
#seq_len = 100
#alphabet = c(1, 2, 3, 4)
#year_mut_rate = 1.8*10^-2 #1.8*10^-3
#mut_rate = year_mut_rate / (365 / delta_t) 

##################################################################################
##################################################################################
sir_mutation_agent = function(params) {
  with(params, {
    
    ########################################################################
    # beta is the number of contacts sufficient to transmit infection per unit time
    beta = R0*gamma*delta_t
    infectivity = beta/(contact.per.day*delta_t) # contact rate already in terms of time step
    
    recover.prob = (1-exp(-gamma*delta_t)) # Probability of Recovery given still infectious
    ########################################################################
    # set up the state vector of the population
    # vstate = 0   means susceptible
    # vstate = 1   means infected    
    # vstate = 2   means recovered   
    vstate = rep(0,N) 
    
    # randomly pick I_0 people from N
    index.inf = sample(N,I_0)
    vstate[index.inf] = 1     
    
    # Initialize infected individuals with the same strain
    strain.state = vector("numeric", length = N)
    strain.state[index.inf] <- 1
    base.haplotype = as.vector(rep(1, seq_len))
    
    population.strains <- init_population_strain(base.haplotype) # Creates first entry for keeping track of strains
    population.frequency <- init_population_freqdf() # 100% frequency to start 
    
    ########################################################################
    # now begin the time steps and Main Simulation
    ########################################################################
    t     = tbeg
    S     = S_0
    I     = I_0
    
    vtime = vector("numeric")  # this will be filled with the time steps
    vtime = append(vtime, 0)
    
    # Time Recrod
    
    time_record = data.frame(matrix(data=0, ncol = 8))
    colnames(time_record) = c('vS', 'vI', 'cir.strains', 'cum.strains',
                              'diverge', 'diversity', 'entropy', 'num.mutations')
    time_record$vS = S
    time_record$vI = I
    time_record$cir.strains=1
    time_record$cum.strains=1
    
    while (t<tend&I>0) { 
      # continue the simulation until we have no more infectious people or t>=tend
      
      # Not Sure yet if I still need this
      deltat=delta_t
      
      #if (delta_t==0){ # this is the calculation of a dynamic time step
      #  deltat = 1/(beta*I*S/N + gamma*I)
      #}
      # Count number of individuals in each state 
      current.inf.index = which(vstate == 1) # Current Infection Index
      I = length(current.inf.index) 
      
      
      ############################################################
      # Mutation 
      ############################################################
 
      # Expected number of mutations Assuming one sequence per individual
      mean.mutation <- mut_rate*seq_len*I
      num.mutations <- rpois(1, mean.mutation)
      
      if (num.mutations > 0) {    
        ## possible add a section if there is only one infected 
        
        if (length(current.inf.index) == 1) {
          for(i in 1:num.mutations) {
            original.strain.index = strain.state[current.inf.index] 
            
            original.strain = population.strains[[original.strain.index]]  
            
            new.strain = mutate.strain(original.strain, alphabet, seq_len)
            new.strain.index <- length(population.strains)+1
            strain.state[current.inf.index] <- new.strain.index  # update strain index of population
            
            population.frequency[new.strain.index] <- numeric()  # Add frequency marker to tell if eventually replicates
            population.strains[[new.strain.index]] <- new.strain
          }
        } else {
          # choose indices from Current Infecteds
          mut.current.index = sample(x = current.inf.index,size =  num.mutations, 
                                     replace = FALSE) # Same strain can't mutate twice  
          
          
          for (j in 1:num.mutations) { 
            
            # For each of the mutation get original strain
            # Mutate strain and update strain state
            # Add new strain to the population 
            
            original.strain.index = strain.state[mut.current.index[j]] 
            
            original.strain = population.strains[[original.strain.index]]  
            
            new.strain = mutate.strain(original.strain, alphabet, seq_len)
            new.strain.index <- length(population.strains)+1
            strain.state[mut.current.index[j]] <- new.strain.index  # update strain index of population
            
            population.frequency[new.strain.index] <- numeric()  # Add frequency marker to tell if eventually replicates
            population.strains[[new.strain.index]] <- new.strain
          }
        }
      } 
      
      ############################################################
      # Replication/Infection 
      ############################################################
      
      # sample Poisson random numbers of infected people contacted by each person
      contacts.by.infected = rpois(n = current.inf.index, lambda = contact.per.day)
      
      vnewstate = vstate
      
      for (i in 1:length(current.inf.index)) {
        responsible.ind = current.inf.index[i]
        infected.contacts = sample(c(1:(responsible.ind-1), (responsible.ind+1):N), contacts.by.infected[i])
        
        contact.disease.state = vstate[infected.contacts]
        contact.strian.state = strain.state[infected.contacts]
        contact.prob = runif(contacts.by.infected[i])
        
        new.infected =  which(contact.disease.state == 0 & contact.prob < infectivity)
        
        vnewstate[infected.contacts[new.infected]] = 1
       
        strain.state[infected.contacts[new.infected]] = strain.state[responsible.ind]
      }
      
      
      #contact.strain.state = new.s.state
      #vnewstate[vstate==0&vprob>infectivity] = 1 
      #  responsible.individual = sample(current.inf.index, 1)
      #  responsible.strain <- new.s.state[responsible.individual]
      #  new.s.state[new.infects.index[k]] <- responsible.strain
      #}
      
      #avg.num.infected.people.contacted = 1-exp(-beta*I/N*deltat)
      #vnum.infected.people.contacted = rpois(N,avg.num.infected.people.contacted)
      
      vprob = runif(N)   # sample uniform random numbers
      #vnewstate = vstate # copy the state vector to a temporary vector used for calculations
      #new.s.state = strain.state # copy the strain vector to a temporary vector used for calculations
      
      # Infected people recover if the sampled uniform random
      # number is less than the recovery probability
      vnewstate[vstate==1&vprob<recover.prob] = 2   
      recovered_ind = which(vstate==1&vprob<recover.prob)
      num.recovered = length(recovered_ind)
      strain.state[recovered_ind] = 0
      
      # If a susceptible contacted at least one infective, they are infected
      #vnewstate[vstate==0&vnum.infected.people.contacted>0] = 1 
      #new.infects.index = which(vstate==0&vnum.infected.people.contacted>0)
      #num.new.infects = length(new.infects.index)
      #infects.per.person = num.new.infects/I
      
      
      vstate = vnewstate # update the state vector
      
      #if (num.new.infects > 0) {
      #  # Assign strains to newly infected:
      # if (I == 1) { 
      
      # if only one strain available all new strains get that one
      #    new.s.state[new.infects.index] <- new.s.state[current.inf.index]
      #  } else { 
      # if more than one present infected, then sample 
      #   for(k in 1:num.new.infects) {
      
      #      responsible.individual = sample(current.inf.index, 1)
      #      responsible.strain <- new.s.state[responsible.individual]
      #      new.s.state[new.infects.index[k]] <- responsible.strain
      #    }
      #  }
      #} 
      
      ## Recovered Individual Strains are a dead-end
      #new.s.state[recovered_ind] <- 0
      #strain.state <- new.s.state
      
      
      # update I for while condition and frequency calculations 
      I = length(vstate[vstate==1])
      sstate  = length(which(strain.state > 0))
      
      if (I == 0) {
        return
      } else {
        ## Add any new strains to population and calculate overall frequencies
        unique.strains <- unique(strain.state[strain.state > 0])
        
        current.freq =  data.frame(aaply(.data = unique.strains, .margins = 1, function(x) 
          freq = length(which(strain.state == x))/length(which(strain.state > 0))))
        colnames(current.freq) = "frequency"
        current.freq = cbind(unique.strains, current.freq)
        
        for (n in 1:nrow(current.freq)) {
          # Add frequencies to population frequency storage 
          population.frequency[length(vtime)+1, current.freq[n,1]] = current.freq$frequency[n]
        }
        
        #################### Book-Keeping After Replication
        S = length(vstate[vstate==0]) 
        
        # Calculate mutation metrics 
        current.haplotypes <- get_current(population.strains, population.frequency, timestep =  length(vtime)+1)
        circulating.strains <- length(current.haplotypes$hindex)
        diversity.time  =  get_diversity(current.haplotypes)
        divergence.time  = get_divergence(current.haplotypes, base.haplotype) 
        entropy.time = get_s.entropy(current.haplotypes)
        
        cum.strains = length(population.strains)
        
        time_record <- rbind(time_record, c(S,I, circulating.strains, cum.strains,  
                                            divergence.time, diversity.time, entropy.time,
                                            num.mutations))
        vtime = append(vtime, t)
        t = t + deltat 
        #print(paste0("Finished Replication", t-1))
      }
    } 
    
    # Calculate final size of each epidemic 
    final = 0
    if (length(time_record$vS)>0) final = 1-min(time_record$vS)/N
    
    # Naming all the population strains for easier analsis 
    names(population.strains) <- paste("s", seq_along(population.strains), sep = ".") 
    colnames(population.frequency) <- paste("s", seq_along(population.frequency), sep = ".")
    time_record <- cbind(as.data.frame(vtime), time_record)
    
    return(list(time_record = time_record, final_size = final, strain.freq = population.frequency,
                population.strains = population.strains))
  })
}

