library(tidyverse)
library(viridis)

# Parameters
juvenile_survival_rate <- 1 - ((0.429 * 2) / 2)
subadult_survival_rate <- 1 - ((0.429 * 1.5) / 2)
adult_survival_rate <- 1 - (0.429 / 2)

juvenile_to_subadult_rate <- 1
subadult_to_adult_rate <- 0.5

reproduction_rate <- 2
reproduction_rate_2 <- 1.5

carrying_capacity <- 68
time_steps <- 500

# Parameter ranges
fishing_effort_values <- seq(0, 0.5, by = 0.05)
restocking_values <- c(0, 1, 3, 5) # Restocking scenarios


run_simulation <- function(F_adults, F_juveniles, restocked_juveniles, burn_in_time = 100) {
  juvenile_population <- 10
  subadult_population <- 10
  adult_population <- 10
  
  # Store population over time
  population_over_time <- numeric(time_steps)
  
  for (t in 1:time_steps) {
    if (t <= burn_in_time) {
      # Burn-in period with zero fishing and zero restocking
      F_adults_current <- 0
      F_juveniles_current <- 0
      restocking <- 0
    } else {
      # Post burn-in period with specified fishing and restocking
      F_adults_current <- F_adults
      F_juveniles_current <- F_juveniles
      restocking <- restocked_juveniles
    }
    
    if (t %% 2 == 0) {
      current_reproduction_rate <- reproduction_rate
    } else {
      current_reproduction_rate <- reproduction_rate_2
    }
    
    total_population <- juvenile_population + subadult_population + adult_population
    new_juveniles <- adult_population * current_reproduction_rate * (1 - total_population / carrying_capacity)
    surviving_juveniles <- juvenile_population * juvenile_survival_rate * (1 - F_juveniles_current)
    new_subadults <- juvenile_population * juvenile_to_subadult_rate
    surviving_subadults <- subadult_population * subadult_survival_rate
    new_adults <- subadult_population * subadult_to_adult_rate
    surviving_adults <- adult_population * adult_survival_rate * (1 - F_adults_current)
    
    juvenile_population <- max(0, new_juveniles + surviving_juveniles - new_subadults + restocking)
    subadult_population <- max(0, new_subadults + surviving_subadults - new_adults)
    adult_population <- max(0, new_adults + surviving_adults)
    
    # Store total population at each timestep
    population_over_time[t] <- juvenile_population + subadult_population + adult_population
  }
  
  # Return the average of the last 20 timesteps
  avg_population_last_20 <- mean(population_over_time[(time_steps - 19):time_steps])
  
  return(avg_population_last_20)
}



# Sensitivity analysis loop with burn-in period
all_results <- data.frame()


for (restocked_juveniles in restocking_values) {
  # Create a dataframe to store results for this restocking value
  sensitivity_results <- expand.grid(
    F_adults = fishing_effort_values,
    F_juveniles = fishing_effort_values
  ) %>%
    mutate(
      Avg_Total_Population = 0
    )
  
  # Reference population with zero fishing effort
  reference_population <- run_simulation(
    F_adults = 0,
    F_juveniles = 0,
    restocked_juveniles = 0,
    burn_in_time = 100
  )
  
  # Sensitivity analysis loop
  for (i in seq_len(nrow(sensitivity_results))) {
    sensitivity_results$Avg_Total_Population[i] <- run_simulation(
      F_adults = sensitivity_results$F_adults[i],
      F_juveniles = sensitivity_results$F_juveniles[i],
      restocked_juveniles = restocked_juveniles,
      burn_in_time = 100
    )
  }
  
  # Normalize by the reference population
  sensitivity_results <- sensitivity_results %>%
    mutate(
      Relative_Population = Avg_Total_Population / reference_population,
      Restocking = paste0(restocked_juveniles, " g/m²")
    )
  
  # Combine results for all restocking values
  all_results <- bind_rows(all_results, sensitivity_results)
}



# Plot results as heatmaps
ggplot(all_results, aes(x = F_juveniles, y = F_adults, fill = Relative_Population)) +
  geom_tile() +
  facet_wrap(~Restocking, labeller = label_both) +
  scale_fill_gradientn(
    colors = c("#d73027", "#fee08b", "#1a9850"), # Color-blind friendly palette
    name = bquote("Relative biomass")
  ) +
  labs(
    x = "Fishing effort on mañahak",
    y = "Fishing effort on hiteng kahlao",
    title = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


