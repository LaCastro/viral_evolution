# Analysis Functions

epi_size_all <- function(trial) {
  ## function for collecting the final sizes of all epidemics
  ## combines in a dataframe
  
    ldply(.data = trial, .fun = function(x) {
    size <- x$final_size
    return(size)
  })
}
  

population_strains_all <- function(trial) {
  ## function for collecting all the strains from multiple trials
  ## each entry of a list is a different trial
  llply(.data = trial, .fun = function(x) {
    strains <- x$population.strains
    return(strains)
  })
}


time_records_all <- function(trial) {
  ## function to combine records 
  ## returns a list 
  llply(.data = trial, .fun = function(x) {
    records <- x$time_record
    return(records)
  })
}
  

combine_time_records <- function(time.records.all) {
  # Takes a list of time records
  # adds an iteration column to each entry
  # combines into a dataframe 
  for (iter in 1:length(time.records.all)) {
    time.records.all[[iter]] <- cbind(iter, time.records.all[[iter]])
  }
  time.records.all  <- do.call("rbind", time.records.all) 
  return(time.records.all)
}


strain_freq_all <- function(trial) {
  ## function to combine records 
  strain.records <-llply(.data = trial, .fun = function(x) {
    s.records  <- x$strain.freq
    return(s.records)
  })
  for (iter in 1:length(strain.records)) {
    strain.records[[iter]] <- cbind(iter, strain.records[[iter]])
  }
  return(strain.records)
}

###################  Finding Maximum Functions

## Find the max diversity
max_diversity <- function(x) {
  return(max(x[,"diversity"]))
}
all_max_diversity <- function(x) {
  return(unlist(laply(x,max_diversity)))
}

# Finding divergence
max_divergence <- function(x) {
  return(max(x[,"diverge"]))
}
all_max_divergence <- function(x) {
  return(unlist(laply(x,max_divergence)))
}


# Finding Max Infected 
max_infected <- function(x) {
  return(max(x[,"vI"]))
}
all_max_infected <- function(x) {
  return(unlist(laply(x,max_infected)))
}


#Get Time of Max Infection
time_max_infected <- function(x) {
  return(x[which.max(x[,"vI"]), "vtime"])
}
all_time_max_infected <- function(x) {
  return(unlist(laply(x,time_max_infected)))
}

#Get Time of Max Diversity 
time_max_diversity <- function(x) {
  return(x[which.max(x[,"diversity"]), "vtime"])
}
all_time_max_diversity <- function(x) {
  return(unlist(laply(x,time_max_diversity)))
}

#Get Time of Max Divergence 
time_max_divergence <- function(x) {
  return(x[which.max(x[,"diverge"]), "vtime"])
}
all_time_max_diverge <- function(x) {
  return(unlist(laply(x, time_max_divergence)))
}


# Get Divergence/Diversity Estimates when 
# 10% of the population is still infected 

threshold_metrics <- function(threshold, N, records.list) {
  metrics.all <- ldply(records.list, .fun = function(x) {
    number.individuals <- threshold * N
    rows <- which(x[, "vI"] < number.individuals)
    diff.rows <- diff(rows)
    index = which.max(diff.rows)+1 ## Check to see if this is working 
    
    metrics <- cbind(x[rows[index], "diverge"], x[rows[index], "diversity"], x[rows[index], "vtime"])
    return(metrics)
  })
  colnames(metrics.all) <- c("diverge", "diversity", "time")
  return(metrics.all)
}





