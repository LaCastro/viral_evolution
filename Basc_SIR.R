library(deSolve)
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I
    dI <- beta * S * I - gamma * I
    dR <- gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

init <- c(S = 1-1e-6, I = 1e-6, 0.0) #1e-6
parameters <- c(beta = 1.4247, gamma = 0.14286)
times <- seq(1, 70, by = 1)
out <- as.data.frame(ode(y = init, times = times, func = sir, parms = parameters))
head(out)
population.out <- 10^6*out #population x proportion to get numbers 
population.out$time <- NULL
head(population.out)


seq_len <- 100
alphabet = c(1, 2, 3, 4)
mut_rate = 0.0001 # per gen per individual per site
infected <- population.out$I


mutation.drift <- function(times, infected, seq_len, mut_rate) {
  #initialize sequence 
  base_haplotype = rep(1, seq_len)

  #Setting up list to store the frequencies for each strain, keeping track of when its in the population 
  frequency.df <- data.frame(cbind(times, 0))
  colnames(frequency.df) <- c("time", "frequency")
  
  hap1 <- list(historical.index = 1, frequency = frequency.df, strain = base_haplotype)
  hap1$frequency[1,2] <- 1 #When there is only one strain in the population
  
  #just for testing
  base_haplotype2 = rep(2, seq_len)
  hap1$frequency[1,2] <- 0.7 #When there is only one strain in the population
  hap2 <- list(historical.index = 2, frequency = frequency.df, strain = base_haplotype2)
  hap2$frequency[1,2] <- 0.3 #When there is only one strain in the population
  
  hap.pop <- list(hap1, hap2)
  
  #Beginning Mutation and Drift 
  for (i in 1:length(infected)) {
    
    #marking the time step 
    timestep <- times[i]
    
    #get expected number of mutations assuming one sequence per individual
    mean.mutation <- mut_rate*seq_len*infected[i]
    num.mutations <- rpois(1, mean.mutation)
    
    current.haplotypes <- get_haplotypes(hap.pop, timestep)  
  
    #If there are mutations....   
    if(num.mutations > 0){
      
      #choose equal number of haplotypes as number of mutations from infected by their frequency in the population
      chosen.strains <- rmultinom(n = 1, size = num.mutations, prob = current.haplotypes$frequencies)
      
      for (j in 1:length(chosen.strains)) {
        
        num.times.specific.strain <- chosen.strains[j,1] # how many times did we sample the same strain
        og.strain.index <- j
        
        for (m in 1:num.times.specific.strain) {
          
          # Mutate the strain 
          new.strain <- mutate.strain(strain = current.haplotypes$strain[[og.strain.index]])
          
          #Index the next 
          strain.index <- length(current.haplotypes$strain)+1
          
          #Check if new strain is the same as any previous strains 
          for (i in 1:length(current.haplotypes)) {
            # Is it the same as any of the current haplotypes? 
            if(length(which(new.strain == current.haplotypes[[i]]$strain)) == 100) {
             
              
              #update the frequency of the strain for the timte step 
              #Need to change this from current haplotypes
              current.haplotypes$frequencies[i] <- current.haplotypes$frequencies[i]*infected[timestep]+1)/infected[timestep]
            }
          else {
          current.haplotypes$strain[[strain.index]] <- new.strain
          current.haplotypes$hindex[strain.index] <- "temporary" #Good through here 
          
          # Calculate and change the frequencies - Messed up here 
          #calulate.mutated.frequencies(current.haplotypes, infected, og.strain.index, strain.index, timestep)
          new.strain.frequency <- 1/infected[timestep]
          old.strain.frequency <- (current.haplotypes$frequencies[og.strain.index]*infected[timestep]-1)/infected[timestep] # new frequency of old strain
          current.haplotypes$frequencies[og.strain.index] <- old.strain.frequency
          current.haplotypes$frequencies[strain.index] <- new.strain.frequency
        }
      }
    } 
  }
}
      
    #replicate based on number of infections in next timestep
   
    if(infected[timestep] == length(infected)) {
      return
    }
    else {
      num.new.infecteds <- infected[timestep+1]
      new.strains <- rmultinom(n = 1, size = num.new.infecteds, prob = current.haplotypes$frequencies)
     
      #list of new haplotypes - #need to change this so that frequencies aren't in a column
      nex.gen.hap <- list(hindex = current.haplotypes$hindex, frequencies = as.vector(new.strains)/num.new.infecteds, strain = current.haplotypes$strain)
      
      #check if any have temporary have above 0 frequency 
      for (i in 1:length(nex.gen.hap$hindex) {
        
      }      
      #update historical record 
      
      
    }
   
    else{}
            
        }
    
      }
      
      
      #Mutation
      
      #get random integer in sequence
      
      #new base pair 
      
      
      
       
       
      
        
       #Function to randomly sample a haplotype from the present population and mutate it   
       #First step, get the total population present 
       #Choose one based on its frequency 
  
          #choose the haplotype
          num.sequence <- seq(1:length(present.sequences))
          chosen <- rmultinom(, size = length(present.sequences), prob = present.frequences)
          
          
        
        



matplot(times, out, type = "l", xlab = "Time", ylab = "Susceptibles and Recovereds", main = "SIR Model", lwd = 1, lty = 1, bty = "l", col = 2:4)
legend(40, 0.7, c("Susceptibles", "Infecteds", "Recovereds"), pch = 1, col = 2:4)



