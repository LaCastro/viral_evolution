get_haplotypes <- function(hap.pop,timestep) {
  
  #get the population of present strains and their frequencies 
  present.sequences <- list()
  present.frequences <- c()
  present.hindex <- c()
  
  total.present <- 0
  for (m in 1:length(hap.pop)) {
    strain <- hap.pop[[m]]
    
    #Is this sequence present in the present timestep? If it is, record it 
    if (strain$frequency[timestep,2] > 0) {
      total.present <- total.present+ 1
      present.sequences[[total.present]] <- strain$strain
      present.frequences[total.present] <- strain$frequency[timestep,2]
      present.hindex[total.present] <- strain$historical.index
    }
    
  } 
  current.haplotypes <- list(hindex = present.hindex,frequencies =  present.frequences, strain = present.sequences)
  return(current.haplotypes)
}



#Mutate Strain and Recalculate frequencies 
mutate.strain <- function(strain)  {
  strain <- current.haplotypes$strain[[j]] #Get selected strain from list of current haplotypes
  baseindex <- round(runif(1, min = 1, max = seq_len)) # select an index to mutate
  nucleotide <- strain[baseindex]
  possible.mutations <- alphabet[which(alphabet != nucleotide)]
  strain[baseindex] <- sample(possible.mutations, 1) #Mutate the strain to any of the other options
  return(strain)
}
  

#calculate frequencies -can't get this to work yet 
#calulate.mutated.frequencies <- function(current.haplotypes, og.strain.index, infected, strain.index, timestep) {
  
#  new.strain.frequency <- 1/infected[timestep]
#  old.strain.frequency <- (current.haplotypes$frequencies[og.strain.index]*infected[timestep]-1)/infected[timestep] # new frequency of old strain
#  current.haplotypes$frequencies[og.strain.index] <- old.strain.frequency
#  current.haplotypes$frequencies[strain.index] <- new.strain.frequency
}





