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
  ## can make the second step faster, do this! 
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

end_threshold_metrics <- function(threshold, N, records.list, rnott) {
  metrics.all <- ldply(records.list, .fun = function(x) {
    number.individuals <- threshold * N
    # first check to see if it ever reaches that level of individuals 
    max.infected = max_infected(x = x)
    if (max.infected > number.individuals) {
     
      diff.rows <- diff(rows)
      index = which.max(diff.rows)+1 ## Check to see if this is working 
      
      metrics <- cbind(rnott, x[rows[index], "diverge"], x[rows[index], "diversity"], x[rows[index], "vtime"])
      colnames(metrics) <- c("rnott", "diverge", "diversity", "time")
      return(metrics)
    } else {
      return()
    }
    
  })
  return(metrics.all)
}
 

beginning_threshold_metrics <- function(threshold, N, records.list, rnott) {
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
                     below.threshold[index, "diversity"], below.threshold[index, "vtime"])
    colnames(metrics) <- c("rnott", "diverge", "diversity", "time")
    return(metrics)
  })
  return(metrics.all)
}


# Need to function to combine max divergence/diversity of each and returning 
get_max_genetic_metrics <- function(time.records, rnott) {
  max.divergence <- all_max_divergence(time.records)
  max.diversity <- all_max_diversity(time.records)
  data <- cbind(max.divergence, max.diversity)
  data <- data.frame(cbind(rep(rnott, nrow(data)), data))
  colnames(data)[1] <- "rnott"
  return(data)
}



get_vec_of_files <- function(dir_path, r0_seq, N){
  data_files <- c()
  for(r_not in r0_seq){
        #pattern <- paste0("*_", "t" paste(r_not,  disc_prob, intro_rate, sep="_"), ".Rdata")  
      pattern  <- paste0("trial_rnott", r_not, "_N")
      data_files <- c(data_files, list.files(path=dir_path, pattern=pattern, full.names=T, recursive=FALSE))
  }
  data_files
}



## combine all num of circulating strains for time-series 
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


