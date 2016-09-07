
##################################################################################
##################################################################################
sir_mutation_agent = function(params) {
  with(params, {
    
    ########################################################################
    # beta is the number of contacts sufficient to transmit infection per unit time
    beta = R0*gamma 
    recover.prob = (1-exp(-gamma*delta_t))
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
    
    population.strains <- init_population_strain(base.haplotype)
    population.frequency <- init_population_freqdf()
    
    ########################################################################
    # now begin the time steps and Main Simulation
    ########################################################################
    t     = tbeg
    S     = S_0
    I     = I_0
    
    vtime = vector("numeric")  # this will be filled with the time steps
    vtime = append(vtime, 0)
    
    # Time Recrod
    time_record = data.frame(matrix(data=0, ncol = 7))
    colnames(time_record) = c('vS', 'vI', 'cir.strains', 'cum.strains',
                              'diverge', 'diversity', 'num.mutations')
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
      
      ############################################################
      # Mutation 
      ############################################################
      # Count number of individuals in each state 
      current.inf.index = which(vstate == 1) # Current Infection Index
      I = length(current.inf.index) 
      
      # Expected number of mutations Assuming one sequence per individual
      mean.mutation <- mut_rate*seq_len*I
      num.mutations <- rpois(1, mean.mutation)
      
      if (num.mutations > 0) {    
        # choose indices from Current Infecteds
        mut.current.index = sample(current.inf.index, num.mutations, 
                                   replace = FALSE) # Same strain can't mutate twice  
        
        ### Trying SAPPLY 
        undebug(mutation.addition)
        mutation.addition <- function(index, strain.state, population.strains) {
          original.strain.index = strain.state[index]
          original.strain = population.strains[[original.strain.index]]  
          new.strain = mutate.strain(strain = original.strain)
          new.strain.index <- length(population.strains)+1
          
          population.frequency[new.strain.index] <- numeric()  # Add frequency marker to tell if eventually replicates
          population.strains[[new.strain.index]] <- new.strain
          return(list(c(unlist(strain.state),population.frequency, population.strains)))
        }
        
        
        mutation.addition(index = mut.current.index, strain.state = strain.state, population.strains = population.strains)
        
        
        sapply(mut.current.index, mutation.addition, strain.state = strain.state, population.strains = population.strains) 
        
        for (j in seq_along(mut.current.index)) {
          original.strain.index = strain.state[mut.current.index[j]] 
          original.strain = population.strains[[original.strain.index]]  
          new.strain = mutate.strain(strain = original.strain)
          new.strain.index <- length(population.strains)+1
          strain.state[mut.current.index[j]] <- new.strain.index  # update strain index of population
          
          population.frequency[new.strain.index] <- numeric()  # Add frequency marker to tell if eventually replicates
          population.strains[[new.strain.index]] <- new.strain
        })
        
        # APLY HERE 
        for (j in seq_along(mut.current.index)) { 
          # For each of the mutation get original strain
          # Mutate strain and update strain state
          # Add new strain to the population 
          
          original.strain.index = strain.state[mut.current.index[j]] 
          original.strain = population.strains[[original.strain.index]]  
          new.strain = mutate.strain(strain = original.strain)
          new.strain.index <- length(population.strains)+1
          strain.state[mut.current.index[j]] <- new.strain.index  # update strain index of population
          
          population.frequency[new.strain.index] <- numeric()  # Add frequency marker to tell if eventually replicates
          population.strains[[new.strain.index]] <- new.strain
        }
      } 
      
      ############################################################
      # Replication/Infection 
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
      
      # Assign strains to newly infected  
     ## APLY HERE 
       for(k in seq_along(num.new.infects)) {
        responsible.individual = sample(current.inf.index, 1)
        responsible.strain <- strain.state[responsible.individual] 
        strain.state[new.infects.index[k]] <- responsible.strain
      }
      
      ## Recovered Individual Strains are a dead-end
      strain.state[recovered_ind] <- 0
      
      # update I for while condition and frequency calculations 
      I = length(vstate[vstate==1])
      if (I == 0) {
        return
      } else {
        ## Add any new strains to population and calculate overall frequencies
        unique.strains <- unique(strain.state[strain.state > 0])
        
        current.freq =  data.frame(aaply(.data = unique.strains, .margins = 1, function(x) 
          freq = length(which(strain.state == x))/length(which(strain.state > 0))))
        colnames(current.freq) = "frequency"
        
        current.freq = cbind(unique.strains, current.freq)
        
        # APPLY HERE 
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
        cum.strains = length(population.strains)
        
        time_record <- rbind(time_record, c(S,I, circulating.strains, cum.strains,  
                                            divergence.time, diversity.time, 
                                            num.mutations))
        vtime = append(vtime, t)
        t = t + deltat 
      }
    } 
    
    # Calculate final size of each epidemic 
    final = 0
    if (length(time_record$vS)>0) final = 1-min(time_record$vS)/N
    
    # Naming all the population strains for easier analsis 
    names(population.strains) <- paste("s", seq_along(population.strains), sep = ".") 
    colnames(population.frequency) <- paste("s", seq_along(population.frequency), sep = ".")
    colnames(strain.freq) <- paste("s", seq_along(strain.freq), sep = ".") 
    time_record <- cbind(as.data.frame(vtime), time_record)
    
    return(list(time_record = time_record, final_size = final, strain.freq = population.frequency,
                population.strains = population.strains))
  })
}

