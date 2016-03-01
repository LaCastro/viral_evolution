get_haplotypes <- function(hap.pop) {
  
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



 <- get_haplotypes(hap.pop)
