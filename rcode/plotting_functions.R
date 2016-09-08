plot_final_sizes <- function(epi.size.all) {
  qplot(epi.size.all, geom = "histogram", bins = 20) +
    scale_y_continuous(expand=c(0,0)) +
    scale_x_continuous(expand=c(0,0)) +
    labs(x= "Final Epidemic Size", y ="Frequency")
}



plot_max_diversity <- function(time.records.all) {
 max.diversity <- all_max_diversity(time.records.all)
 qplot(max.diversity, geom = "histogram", bins = 25) +
   scale_y_continuous(expand = c(0,0)) +
   scale_x_continuous(expand = c(0,0)) + 
   labs(x = "Max Diversity")
}


plot_max_divergence <- function(time.records.all) {
  max.divergence <- all_max_divergence(time.records.all)
  qplot(max.divergence, geom = "histogram", bins = 25) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) + 
    labs(x = "Max Divergence")
}



plot_max_times <- function(time.records.all, type) {  
  time.max.infected <- all_time_max_infected(time.records.all)
  if (type == "diversity") {
    time.max.variable <- all_time_max_diversity(time.records.all)
  } else { 
    time.max.variable <- all_time_max_diverge(time.records.all)
  }
  data <- data.frame(cbind(time.max.variable, time.max.infected))
  plot <- ggplot(data, aes(time.max.variable, time.max.infected)) +
    geom_point(color = "black", size = 3) + 
    geom_abline(intercept = 0, slope = 1) +
    labs(x = ifelse(type == "diversity", "Day of Maximum Strain Diversity", "Day of Maximum Strain Divergence"), y = "Day of Maximum Infected Cases")
  return(plot)
}


