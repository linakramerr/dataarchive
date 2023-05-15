
# content for table 2

simulationconditions <- data.frame("True proportion of measurement error" = "0, 0.2, 0.4, 0.6", "Sample size" = "300, 10000", "Number of timepoints" = "3, 8", "Cross-lagged effects" = ".2, .2; .2, 0; -.2, -.2", "Proportions of between:within variances" = "70:30, 50:50, 30:70")

simulationconditions <- plyr::ldply(simulationconditions, cbind) 
colnames(simulationconditions) <- c("Factor", "Values")

simulationconditions %>% knitr::kable(caption = "Simulation study factors") %>%
  kable_styling(latex_options = c("striped", "hold_position"))
