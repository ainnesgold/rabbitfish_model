library(tidyverse)
library(viridis)

# Parameters
burn_in_time <- 100  # Define burn-in period of 100 timesteps
time_steps <- 500 

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

# Define the range of restocking mortality rates
restocking_mortality_rates <- seq(0, 0.8, by = 0.1)

#restocking scenario
restocked_juveniles = 1

#Current fishing for burn in
F_current_instantaneous <- 0.09
F_current_discrete <- (1 - exp(-F_current_instantaneous)) / 2
F_current_discrete_juv <- 0.01




# Update run_simulation to include restocking mortality
run_simulation <- function(F_adults, F_juveniles, restocked_juveniles, restocking_mortality_rate, burn_in_time = 100) {
  juvenile_population <- 10
  subadult_population <- 10
  adult_population <- 10
  
  # Store population over time
  population_over_time <- numeric(time_steps)
  
  for (t in 1:time_steps) {
    if (t <= burn_in_time) {
      # Burn-in period with burn in fishing and no restocking
      F_adults_current <- F_current_discrete
      F_juveniles_current <- F_current_discrete_juv
      restocking <- 0
    } else {
      # Post burn-in period with specified fishing and restocking
      F_adults_current <- F_adults
      F_juveniles_current <- F_juveniles
      restocking <- restocked_juveniles * (1 - restocking_mortality_rate)
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
    
    if (t %% 2 == 0) {
      subadult_population <- max(0, new_subadults + surviving_subadults - new_adults + restocking)
    } else {
      subadult_population <- max(0, new_subadults + surviving_subadults - new_adults)
    }
    
    juvenile_population <- max(0, new_juveniles + surviving_juveniles - new_subadults)
    adult_population <- max(0, new_adults + surviving_adults)
    # Store total population at each timestep
    population_over_time[t] <- juvenile_population + subadult_population + adult_population
  }
  
  # Return the average of the last 20 timesteps
  avg_population_last_20 <- mean(population_over_time[(time_steps - 19):time_steps])
  
  return(avg_population_last_20)
}

# Sensitivity analysis loop
all_results <- data.frame()

for (restocking_mortality_rate in restocking_mortality_rates) {
  # Create a dataframe to store results for this restocking mortality rate
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
    restocked_juveniles = 1,
    restocking_mortality_rate = restocking_mortality_rate,
    burn_in_time = 100
  )
  
  # Sensitivity analysis loop
  for (i in seq_len(nrow(sensitivity_results))) {
    sensitivity_results$Avg_Total_Population[i] <- run_simulation(
      F_adults = sensitivity_results$F_adults[i],
      F_juveniles = sensitivity_results$F_juveniles[i],
      restocked_juveniles = 1,
      restocking_mortality_rate = restocking_mortality_rate,
      burn_in_time = 100
    )
  }
  
  # Normalize by the reference population
  sensitivity_results <- sensitivity_results %>%
    mutate(
      Relative_Population = Avg_Total_Population / reference_population,
      Restocking_Mortality = paste0(restocking_mortality_rate * 100, "%")
    )
  
  # Combine results for all restocking mortality rates
  all_results <- bind_rows(all_results, sensitivity_results)
}

# Custom labeller function for restocking mortality
restocking_labeller <- function(values) {
  paste("Restocked M =", values)
}

restocked_juveniles / reference_population


# Update the code to use this custom labeller
ggplot(all_results, aes(x = F_juveniles, y = F_adults, fill = Relative_Population)) +
  geom_tile() +
  geom_contour(
    aes(z = Relative_Population), 
    breaks = 0.5, 
    color = "black", 
    linetype = "dashed", 
    linewidth = 1
  ) +
  facet_wrap(
    ~Restocking_Mortality, 
    labeller = labeller(Restocking_Mortality = function(values) restocking_labeller(values))
  ) +
  scale_fill_gradientn(
    colors = c("#d73027", "#fee08b", "#1a9850"), # Color-blind friendly palette
    name = bquote("Biomass relative to"~B[0]),
  ) +
  labs(
    x = "Fishing effort on mañahak",
    y = "Fishing effort on hiteng kahlao",
    caption = bquote("Restocking = 2% of"~B[0]),
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )








