# Sampling Functions for Post-analysis of strain records

##  Uniform Cases 
sample_uniform<- function(time.record, n.samples, low.detection.threshold) {
  ##  Calculates a Time-Series of detected cases
  ##  Assumes that the number of samples in each time step is the same
  ##  There is a low detection threshold where a proption of samples are collected 
  ##  When infections are not high 
  ## Works on one time record
  cases.sample <- adply(.data = time.record, .margins = 1,.id=NULL, .expand = F, function(x) {
    detected.time <- x$vtime
    infected <- x$vI
    if (infected < n.samples) {
      detected.cases = round(low.detection.threshold*infected,0)
    } else {
      detected.cases = n.samples
    } 
    detected = cbind(detected.time, infected, detected.cases)
    return(detected)
  })
  return(cases.sample)
}

#######  
sample_uniform_all <- function(time.records, ...) {
  # Returns list of all time series
  # Of detected cases for uniform
  sample.cases.all <- llply(.data = time.records, .fun = function(x) {
    cases = sample_uniform(x, ...)
  }) 
  return(sample.cases.all)
}



#### 
sample_proportion <- function(time.record, prop) {
  # Cases are detected by a constant proportion
  # as a function of Infected Individuals
  # Returns results for one time record 
  cases.sample <- adply(.data = time.record, .margins = 1,.id=NULL, .expand = F, function(x) {
    detected.time <- x$vtime
    infected <- x$vI
    detected.cases <- round(infected*prop, 0)
    detected = cbind(detected.time, infected, detected.cases)
    return(detected)
  })
  return(cases.sample)
}

sample_prop_all <- function(time.records, ...) {
  # Returns list of all time series
  # Of detected cases for proportional sampler 
  sample.cases.all <- llply(.data = time.records, .fun = function(x) {
    cases = sample_proportion(x, ...)
  }) 
  return(sample.cases.all)
}


sampler_strains <- function(strain.record, cases.sample) {
  ## Strain Sampler, works with any sampler strategy
  ## Takes in a list of detected cases 
  ## Samples strains according to the number of detected cases
  ## And frequency in each time step 
  ## Function for one strain recorod and one set of detected cases
  strain.combined <- cbind(cases.sample, strain.record)
  strains.sample <- adply(.data = strain.combined, .margins = 1,.id=NULL, .expand = F, function(x) {
    time <- x$detected.time
    probs <- x[,5:ncol(x)] # Probability For Each Strain
    probs[is.na(probs)] <- 0 # Replace NA for 0
    
    samples <- data.frame(rmultinom(n = 1, size = x$detected.cases, prob = probs)) # how many of each strain, by the number of detected
    strain.names <- rownames(samples)
    samples <- cbind(time, strain.names, samples, x$detected.cases) #combine the time and the number of time each strain was sampled 
    return(samples)
  })
  colnames(strains.sample)[3:4] <- c("num.strain.samples", "detected.cases")
  return(strains.sample)
}

sampler_strains_all <- function(strain.records, cases.sample.all) {
  # Returns N iterations of sampled strains 
  strains.sampled.all <- Map(sampler_strains, strain.record = strain.records, cases.sample = cases.sample.all)
  return(strains.sampled.all)
}


###################################################################################
# Functions for Analyzing Sampled Strains
###################################################################################

unique_samples_record <- function(sample.strains.record) {
  unique.strains <- daply(sample.strains.record, .(time), function(x) {
    sampled <-  filter(x, num.strain.samples >=1)
    unique <- nrow(sampled)
    return(unique)
  })
  return(unique.strains)
}

get_all_unique_samples <- function(samples.strains.all) {
  ## Count up number of unique strains
  ## sampled in a time step
  ## Return in a list 
  unique.strains.list <- lapply(samples.strains.all, function(x) {
    unique.strains <- data.frame(unique_samples_record(x))
    time.steps <- seq(1:nrow(unique.strains))
    unique.strains <- cbind(time.steps, unique.strains)
    return(unique.strains)
  })
  unique.strains.all <- rbindlist(unique.strains.list, idcol = "iter")
  colnames(unique.strains.all)[3] <- "unique.strains"
  return(unique.strains.all)
}


