library(deSolve)
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dcumI <- beta*S*I
    dR <- gamma * I
    
    return(list(c(dS, dI, dcumI,dR)))
  })
}



# Need to convert cumulative to incidence 

init <- c(S = 1-1e-6, I = 1e-6, 0.0) #1e-6
parameters <- c(beta = 1.4247, gamma = 0.14286)
times <- seq(1, 70, by = 1)
out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
population.out <- 10^6*out #population x proportion to get numbers 
population.out$time <- NULL


seq_len <- 100
alphabet = c(1, 2, 3, 4)
mut_rate = 0.0001 # per gen per individual per site

#For testing
infected <- population.out$I[1:6]
times <- times[1:6]

mutation.drift <- function(times, infected, seq_len, mut_rate) {
  #initalize diveristy and divergence vectors
  diversity <- rep(0, length(times))
  divergence <- rep(0, length(times))
  
  
  #initialize sequence 
  base_haplotype = rep(1, seq_len)

  #Setting up the population list to store haplotypes through time 
    frequency.df <- data.frame(rep(0, length(times))) #storing the frequencies of each strain
    colnames(frequency.df) <- paste("strain",1, sep = ".")
  
    #storing number of unique strains in the population 
    hindex <- c()
    hindex[1] <- 1 #initializing the first 
   
    #List of strains 
    strain <- list()
    strain[[1]] <- base_haplotype
  
    population <- list(hindex = hindex, frequency = frequency.df, strain = strain )
    population$frequency$strain.1[1] <- 1 #When there is only one strain introduced 


    
    #Beginning Mutation and Drift 
    for (i in 1:(length(times))) {
      
    #marking the time step 
    timestep <- times[i]
    
    
    #calculate diversity and divergence for the timestep before the mutation and replication event
    current.haplotypes <- get_current(population, timestep)  
    num.current.strains <- length(current.haplotypes$hindex)
    diversity[timestep] <- get_diversity(current.haplotypes)
    divergence[timestep] <- get_divergence(current.haplotypes, base_haplotype)
    
    if (timestep == length(timestep)) {
      return
    } else {
    #get expected number of mutations assuming one sequence per individual
    mean.mutation <- mut_rate*seq_len*infected[timestep]
    num.mutations <- rpois(1, mean.mutation)
  
    #If there are mutations....   
    if(num.mutations > 0){
      
      #choose equal number of haplotypes as number of mutations from infected by their frequency in the population
      chosen.strains <- as.vector(rmultinom(n = 1, size = num.mutations, prob = current.haplotypes$frequencies))
      original.current.strains <- length(chosen.strains) #how many strains are there in the population to mutate
      
      for (j in 1:original.current.strains) { #for each of the original strains in the population
        num.current.strains <- length(current.haplotypes$hindex) # how many are there-this will be added to 
        num.times.specific.strain <- chosen.strains[j] # how many times did we sample this particular strains
       
        #if we sampled this particular strain at least once 
        if (num.times.specific.strain > 0) {
        for (m in 1:num.times.specific.strain) { #for however many times we did 
          
          # Mutate the strain 
          new.strain <- mutate.strain(strain = current.haplotypes$strain[[j]]) 
          
          #Index the next 
          num.current.strains <- num.current.strains+1
          
          #Check if new strain is the same as any previous strains #probably don't need this here 
          #for (i in 1:length(current.haplotypes)) {
            # Is it the same as any of the current haplotypes? 
           # if(length(which(new.strain == current.haplotypes[[i]]$strain)) == 100) {
             
              
              #update the frequency of the strain for the time step 
              #Need to change this from current haplotypes
              #current.haplotypes$frequencies[i] <- (current.haplotypes$frequencies[i]*infected[timestep]+1)/infected[timestep]
            #} else {
          current.haplotypes$strain[[num.current.strains]] <- new.strain
          current.haplotypes$hindex[num.current.strains] <- "temporary" #Good through here 
          
          # Calculate and change the frequencies - Messed up here 
          #calulate.mutated.frequencies(current.haplotypes, infected, og.strain.index, strain.index, timestep)
          new.strain.frequency <- 1/infected[timestep]
          old.strain.frequency <- (current.haplotypes$frequencies[j]*infected[timestep]-1)/infected[timestep] # new frequency of old strain
          current.haplotypes$frequencies[j] <- old.strain.frequency
          current.haplotypes$frequencies[num.current.strains] <- new.strain.frequency
          print(sum(current.haplotypes$frequencies))
              }
            }
          } 
        }
  
      
      
    #replicate based on number of infections in next timestep
    if(timestep == length(infected)) {
      return
    } else {
      num.new.infecteds <- ceiling(infected[timestep+1])
      new.strains <- rmultinom(n = 1, size = num.new.infecteds, prob = current.haplotypes$frequencies)
     
      #list of new haplotypes
      nex.gen.hap <- list(hindex = current.haplotypes$hindex, frequencies = as.vector(new.strains)/num.new.infecteds, strain = current.haplotypes$strain)
      print(sum(nex.gen.hap$frequencies))
      
      #check if any have temporary have 0 frequency - mutated didn't replicate 
      for (i in 1:length(nex.gen.hap$hindex)) {
        if(length(which(nex.gen.hap$hindex == "temporary")) > 0) {
          temporary.strains <- which(nex.gen.hap$hindex == "temporary")
          for(j in 1:length(temporary.strains)) {
            if (nex.gen.hap$frequencies[temporary.strains[j]] == 0) {
          nex.gen.hap$hindex <- nex.gen.hap$hindex[-temporary.strains]
          nex.gen.hap$frequencies <-nex.gen.hap$frequencies[-temporary.strains]
          nex.gen.hap$strain <- nex.gen.hap$strain[-temporary.strains]
            }
          }
        }
      }      
      
      #update historical record 
      #update_population <- function(nex.gen.hap, population, timestep) {
        new.hindex <- as.character(nex.gen.hap$hindex)
        pop.hindex <- as.character(population$hindex)
        same.hindex <- which(new.hindex%in%pop.hindex)
        different.hindex <- which(!new.hindex%in%pop.hindex)
        older.strains <- which(!pop.hindex%in%new.hindex) #which hindex were once present but not now
        
        #For those that are the same, update their frequencies in next time step
        for(i in 1:length(same.hindex)) {
        population$frequency[timestep+1,same.hindex[i]] <- nex.gen.hap$frequencies[same.hindex[i]]
        }
        
        #If temporary, check if strain matches previous that have gone extinct
        #Loops through all those strains that weren't in the previous generation
        if (length(different.hindex) > 0 ) {
        for (m in 1:length(different.hindex)) {
          temp_strain <- nex.gen.hap$strain[[different.hindex[m]]]
          
          #Loops through all the strains that have been in the population but weren't in the previous, if not any, add new 
          if (length(older.strains) > 0){
          for(i in 1:length(population$strain[[older.strains[i]]])) {
            ancestor_strain <- population$strain[[older.strains[i]]]
            match = 0
            #checks the distance between the two, if its 0 adds the frequency 
            if(get_distance(temp_strain, ancestor_strain) == 0) {
              population$frequency[timestep+1,older.strains[i]] <- population$frequency[timestep+1,older.strains[i]] + nex.gen.hap$frequency[m]
              match = 1
              return
              }
              if(match == 1) return
              } 
            } else {
          
          #Adding to population
          new.strain.hindex <- length(population$hindex)+1
          population$hindex[new.strain.hindex] <- new.strain.hindex
          population$strain[[new.strain.hindex]] <- temp_strain
          new.strain.name <- paste("strain", new.strain.hindex, sep = ".")
          population$frequency$new.strain.hindex <- rep(0, length(times))
          names(population$frequency)[names(population$frequency) == 'new.strain.hindex'] <- new.strain.name
          population$frequency[timestep+1,new.strain.hindex] <- nex.gen.hap$frequencies[different.hindex[m]]    
              }
            }
          }
        }
      }
    }
 
                 
      
matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SIR Model", lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(40, 0.7, c("Susceptibles", "Infecteds", "Recovereds"), pch = 1, col = 2:4)



