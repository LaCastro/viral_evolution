population_strains_all <- function(trial) {
  ## function for collecting all the strains from multiple trials
  ## of the same parameter set 
  ## each entry of a list is a different trial
  llply(.data = trial, .fun = function(x) {
    strains <- x$population.strains
    return(strains)
  })
}


time_records_all <- function(trial) {
  ## function to combine records 
  ## of the same parameter set 
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
  # useful for plotting 
  for (iter in 1:length(time.records.all)) {
    time.records.all[[iter]] <- cbind(iter, time.records.all[[iter]])
  }
  time.records.all  <- do.call("rbind", time.records.all) 
  return(time.records.all)
}


strain_freq_all <- function(trial) {
  # function to combine records 
  # of the same parameter set 
  # can't rbind at the end because there are differ columns
  strain.records <-llply(.data = trial, .fun = function(x) {
    s.records  <- x$strain.freq
    return(s.records)
  })
  for (iter in 1:length(strain.records)) {
    strain.records[[iter]] <- cbind(iter, strain.records[[iter]])
  }
  return(strain.records)
}

## Combine all num of circulating strains for
## a set of parameters 
combine_mutations <- function(time.records) {
  # Takes a list of time records
  # subsets to just time and cum strains
  # combines into a dataframe 
  num.mutations.master <- ldply(.data = time.records, function(x) {
    time.cum.mutations <- data.frame(cbind(x[,"vtime"], x[,"cum.strains"]))
    return(time.cum.mutations)
  })
  #time.records.all  <- do.call("rbind", time.records.all) 
  colnames(num.mutations.master) <- c("vtime", "cum.strains")
  return(num.mutations.master)
}



################### Functions for getting info from files 
get_vec_of_files <- function(dir_path, type, r0_seq, N){
  # creates a dataframe of desired paths based on 
  # r0 and N
  # to data files so can read in 
  # for analysis 
  data_files <- c()
  for(r_not in r0_seq){
    #pattern <- paste0("*_", "t" paste(r_not,  disc_prob, intro_rate, sep="_"), ".Rdata")  
    pattern  <- paste0(type, "_rnott", r_not, "_N")
    data_files <- c(data_files, list.files(path=paste0(dir_path,type), pattern=pattern, full.names=T, recursive=FALSE))
  }
  data_files
}


get_params <- function(x) {
  # extracts the R0 and N infor
  # and keeps track so easier for plotting
  rnott = as.numeric(sub("_.*", "", x$combos))
  pop.size = as.numeric(sub(".*_", "", x$combos))
  trial.params = data.frame(cbind(rnott, pop.size))
  return(trial.params)
}

######################  Finding Maximum Functions

## Find the Max Diversity
max_diversity <- function(x) {
  return(max(x[,"diversity"]))
}
all_max_diversity <- function(x) {
  return(unlist(laply(x,max_diversity)))
}

# Finding Max Point Divergence
max_divergence <- function(x) {
  return(max(x[,"diverge"]))
}
all_max_divergence <- function(x) {
  return(unlist(laply(x,max_divergence)))
}

# Finding Max S.Entropy
max_sentropy <- function(x) {
  return(max(x[,"entropy"]))
}
all_max_sentropy <- function(x) {
  return(unlist(laply(x,max_sentropy)))
}

# Combines max divergence/diversity for a set of 
# simulations with the same parameters 
get_max_genetic_metrics <- function(time.records, rnott) {
  # Finds Max Diversity/Divergence
  # Combines into a data-frame
  max.divergence <- all_max_divergence(time.records)
  max.diversity <- all_max_diversity(time.records)
  max.entropy <- all_max_sentropy(time.records)
  data <- cbind(max.divergence, max.diversity, max.entropy)
  data <- data.frame(cbind(rep(rnott, nrow(data)), data))
  colnames(data)[1] <- "rnott"
  return(data)
}

# Finding Max Infected 
max_infected <- function(x) {
  return(max(x[,"vI"]))
}
all_max_infected <- function(x) {
  return(unlist(laply(x,max_infected)))
}


#Finding the Time of Max Infection
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

#Get Time of Max Entropy
time_max_entropy <- function(x) {
  return(x[which.max(x[,"entropy"]), "vtime"])
}
all_time_max_entropy <- function(x) {
  return(unlist(laply(x, time_max_entropy)))
}

# Get Max Time of Strain 1/wildtype
time_life_wildtype <- function(x) {
  # how long is the wildtype strain around
  dead.end <- which(is.na(x[,"s.1"]))
  if (length(dead.end) != 0) {
    life.span = dead.end[1]  
  } else {
    life.span  = nrow(x) # strain lasted whole simulation 
  }
  return(life.span)
}
all_time_life_wildtype <- function(x) {
  return(unlist(laply(x, time_life_wildtype)))
}  


## Get Total Time of epidemic 
epi_time <- function(x) {
  return(x[which.max(x[,"vtime"]), "vtime"])
}
all_epi_time <- function(x) {
  return(unlist(laply(x, epi_time)))
}


############################## 
# Functions Based on a Threshold
###############################
# Get Divergence/Diversity Estimates when 
# 10% of the population is still infected 
end_threshold_metrics <- function(threshold, N, records.list, rnott) {
  # Given a threshold of individuals still infected
  # what are the metrics at this point
  metrics.all <- ldply(records.list, .fun = function(x) {
    number.individuals <- threshold * N
    # first check to see if it ever reaches that level of individuals 
    max.infected = max_infected(x = x)
    if (max.infected > number.individuals) {
     
      diff.rows <- diff(rows)
      index = which.max(diff.rows)+1 ## Check to see if this is working 
      
      metrics <- cbind(rnott, x[rows[index], "diverge"], x[rows[index], "diversity"], x[rows[index], "vtime"], x[rows[index], "entropy"])
      colnames(metrics) <- c("rnott", "diverge", "diversity", "time", "entropy")
      return(metrics)
    } else {
      return()
    }
    
  })
  return(metrics.all)
}
 
# Get Divergence/Diversity Estimates when 
# 10% of the population is infected 


beginning_threshold_metrics <- function(threshold, N, records.list, rnott) {
  # Given a threshold of individuals infected (not cumulative)
  # what are the metrics at this point
  metrics.all <- ldply(records.list, .fun = function(x) {
    number.individuals <- threshold * N
    time.of.max <- time_max_infected(x = x)
    below.threshold <- filter(x,vI < number.individuals, vtime < time.of.max)
    #browser()
    if (nrow(below.threshold) == 0) {
      index = 0 
      return()
    } 
    else if (nrow(below.threshold) == 1) {
     # browser()
      index = below.threshold$vtime + 1
    } 
    else {
      # should just be the last row in this 
      #browser()
      index = which.max(below.threshold$vtime)
    }
    #browser()
    metrics <- cbind(rnott, below.threshold[index, "diverge"], 
                     below.threshold[index, "diversity"], below.threshold[index, "entropy"], below.threshold[index, "vtime"])
    colnames(metrics) <- c("rnott", "diverge", "diversity", "entropy", "time")
    return(metrics)
  })
  return(metrics.all)
}



###################################
# MISC Functions 
###################################

# Life span data of each strain
mutant_lifespan <- function(s.record) {
  # Collects lifespan of each mutant strain
  # In each trial
  # Combines all trials 
  if (ncol(s.record) == 2) {
    # The wildtype is the only strain present 
    return()
  }
  s.record <- data.frame(s.record[,3:ncol(s.record)]) # getting the mutants 
  
  if (ncol(s.record) > 1) {
    # More than one mutant in trial
    m.life.span <- aaply(.data = s.record, .margins = 2, .expand = F, function(x) {
      alive = length(which(x > 0)) ## How many days is it alive 
      return(alive)
    })
  } else {
    # Only one mutant in trial
    m.life.span = length(which(s.record > 0)) 
  }
  return(unname(m.life.span))
}

all_mutant_lifespan <- function(strain.records) {
  # Combines mutant lifespans from all runs
  # of same parameter set 
  all.life.span <- llply(strain.records, mutant_lifespan)
  all.life.span = unlist(all.life.span)
  return(all.life.span)
}


# Final Size Incted 
epi_size_all <- function(trial) {
  # function for collecting the final sizes of all epidemics
  # combines in a dataframe
  # only works if you have the full trial data 
  ldply(.data = trial, .fun = function(x) {
    size <- x$final_size
    return(size)
  })
}




