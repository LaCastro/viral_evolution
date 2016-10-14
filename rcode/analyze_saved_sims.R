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
  browser()
  for (iter in 1:length(time.records.all)) {
 
    time.records.all[[iter]] <- cbind(iter, time.records.all[[iter]])
  }
  browser()
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
    time.cum.mutations <- data.frame(cbind(x[,"vtime"], x[,"shifted.time"],x[,"cum.strains"]))
    return(time.cum.mutations)
  })
  #time.records.all  <- do.call("rbind", time.records.all) 
  colnames(num.mutations.master) <- c("vtime", "shifted.time","cum.strains")
  return(num.mutations.master)
}

combine_circulating <- function(time.records) {
  num.circulating <- ldply(.data = time.records, function(x) {
    time.cum.mutations <- data.frame(cbind(x[,"vtime"], x[,"cir.strains"]))
    return(time.cum.mutations)
  })
  colnames(num.circulating) <- c("vtime", "cir.strains")
  return(num.circulating)
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

# Get Time of Max Strains Circulating
time_max_strains <- function(x) {
  return(x[which.max(x[,"cir.strains"]), "vtime"])
}
all_time_max_strains <- function(x) {
  return(unlist(laply(x,time_max_strains)))
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



## Get Total Size of Infected 
infected_trial <- function(time.record, trial.N) {
  remaining.susceptible <- min(time.record$vS)
  total.infected <- trial.N - remaining.susceptible
  total.infected.prop <- total.infected/trial.N
  return(total.infected.prop)
}

all_total_infected <- function(time.records, trial.N) {
  all.infected <- laply(.data = time.records, .fun = function(x) {
    final.infected = infected_trial(x, trial.N)
  })
 return(unlist(all.infected))
}

#Get Metrics When Max is infected
metrics_at_max_infect <- function(time.records.a) {
  # Has to be a shifted time series
  metrics <- ldply(.data = time.records.a, function(x) {
    index = which(x$shifted.time == 0)
    diversity = x[index, "diversity"]
    entropy = x[index, "entropy"]
    num.cir.strains = x[index, "cir.strains"]
    prop.strains.est = x[index, "cum.strains"]/max(x$cum.strains)
    return(cbind(diversity, entropy, num.cir.strains, prop.strains.est))
  })
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


## For getting the results at multiple thresholds 
beg_threshold_metrics <- function(data.files, threshold) { 
  beg.threshold.master <-  adply(.data = data.files, .margins = 1, .id=NULL, .expand = F, function(x) {
    load(as.character(x$file.list))
    trial.params <- get_params(x)
    threshold.metrics = beginning_threshold_metrics(threshold = threshold,trial.params$pop.size, records.list = time.records, rnott = trial.params$rnott)
    threshold.metrics = cbind(rep(trial.params$pop.size, nrow(threshold.metrics)), threshold.metrics)
    colnames(threshold.metrics)[1] <- "pop.size"
    return(threshold.metrics)
  })
  return(beg.threshold.master)
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


align_time_series <- function(trial) {
  # function for assigning peak time as 0
  # time before becomes negative time,
  # adds separate column
  trial$shifted.time <- NA
  max.time <- time_max_infected(trial)
  trial.before <- filter(trial, vtime < max.time)
  if (nrow(trial.before) > 0) {
    trial.before$shifted.time <- -rev(seq(1:nrow(trial.before)))
    trial.after <- filter(trial, vtime >= max.time)
    trial.after$shifted.time <- seq(from = 0, to = (nrow(trial.after) - 1), by = 1)
    trial.shifted <- rbind(trial.before, trial.after)
  } else {
    trial.shifted <- trial
    trial.shifted$shifted.time <- trial.shifted$vtime
  }
  return(trial.shifted)
}

align_time_series_all <- function(time.records) {
  
  time.records.a <- llply(.data = time.records, function(x) {
    trial <- align_time_series(x)
    return(trial)
  })
  return(time.records.a)
}

#Function to identify entries that are epidemics and ones that are not 
get_epidemic_index <- function(time.records, threshold.prev, threshold.prop, trial.N) {
  final.infected.trial <- data.frame(all_total_infected(time.records = time.records, trial.N = trial.N))
  epidemics.prop <- which(final.infected.trial > threshold.prop)
  max.infected <- all_max_infected(time.records)
  epidemics.prev <- which(max.infected > (threshold.prev*trial.N))
  #browser()
  epidemic.trials <- intersect(epidemics.prop, epidemics.prev)
  return(epidemic.trials)
}

set_epidemic_criteria <- function(time.records, threshold.prev, threshold.prop) {
  epidemic.trials <- alply(.data = data.files, .margins = 1, function(x) {
    load(as.character(x$file.list))
    trial.params <- get_params(x)
    # Set thresholds here 
    epidemic.trials <- get_epidemic_index(time.records, trial.N = trial.params$pop.size, 
                                          threshold.prev = threshold.prev, threshold.prop = threshold.prop)
    if (length(epidemic.trials) == 0) epidemic.trials <- 0
    epidemic.trials <- cbind(trial.params, epidemic.trials)
    return(epidemic.trials)
  })
  return(epidemic.trials)
}
