library(tidyverse)

#Parameters
# Define parameters
juvenile_survival_rate <- 1 - ((0.429*2) / 2)   # Probability of a juvenile surviving to the next time step
subadult_survival_rate <- 1 - ((0.429*1.5) / 2)   # Probability of a subadult surviving to the next time step
adult_survival_rate <-  1 - (0.429 / 2)  # 0.429 is the annual mortality, divided in 2 for 6 months time steps - Probability of an adult surviving to the next time step

juvenile_to_subadult_rate <- 1 #0.2 # Probability of a juvenile transitioning to subadult
subadult_to_adult_rate <- 0.5 #0.2   # Probability of a subadult transitioning to adult

reproduction_rate <- 2 #0.19
reproduction_rate_2 <- 1.5 #0.03

carrying_capacity <- 68 #g/m2, in range of Friedlander and Sandin studies

#restocking
restocked_juveniles <- 0

# Number of time steps to simulate
time_steps <- 500


# Parameter ranges
fishing_effort_values <- seq(0, 0.5, by = 0.05)

# Create a dataframe to store results
sensitivity_results <- expand.grid(
  F_adults = fishing_effort_values,
  F_juveniles = fishing_effort_values
) %>%
  mutate(
    Final_Total_Population = 0
  )

# Function to run the simulation and return the final total population
run_simulation <- function(F_adults, F_juveniles) {
  juvenile_population <- 10
  subadult_population <- 10
  adult_population <- 10
  
  for (t in 2:time_steps) {
    if (t %% 2 == 0) {
      current_reproduction_rate <- reproduction_rate
    } else {
      current_reproduction_rate <- reproduction_rate_2
    }
    
    total_population <- juvenile_population + subadult_population + adult_population
    new_juveniles <- adult_population * current_reproduction_rate * (1 - total_population / carrying_capacity)
    surviving_juveniles <- juvenile_population * juvenile_survival_rate * (1 - F_juveniles)
    new_subadults <- juvenile_population * juvenile_to_subadult_rate
    surviving_subadults <- subadult_population * subadult_survival_rate
    new_adults <- subadult_population * subadult_to_adult_rate
    surviving_adults <- adult_population * adult_survival_rate * (1 - F_adults)
    
    if (t %% 2 == 0) {
      juvenile_population <- max(0, new_juveniles + surviving_juveniles - new_subadults + restocked_juveniles)
    } else {
      juvenile_population <- max(0, new_juveniles + surviving_juveniles - new_subadults)
    }
    subadult_population <- max(0, new_subadults + surviving_subadults - new_adults)
    adult_population <- max(0, new_adults + surviving_adults)
  }
  
  return(juvenile_population + subadult_population + adult_population)
}

# Reference population with zero fishing effort
reference_population <- run_simulation(F_adults = 0, F_juveniles = 0)

# Sensitivity analysis loop
for (i in seq_len(nrow(sensitivity_results))) {
  sensitivity_results$Final_Total_Population[i] <- run_simulation(
    F_adults = sensitivity_results$F_adults[i],
    F_juveniles = sensitivity_results$F_juveniles[i]
  )
}

# Normalize by the reference population
sensitivity_results <- sensitivity_results %>%
  mutate(Relative_Population = Final_Total_Population / reference_population)

# Plot results as a heatmap
ggplot(sensitivity_results, aes(x = F_juveniles, y = F_adults, fill = Relative_Population)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("red", "yellow", "green"),
    name = bquote("Relative Population")
  ) +
  labs(
    x = "Juvenile Fishing Effort",
    y = "Adult Fishing Effort",
    title = "Relative Population Sensitivity Analysis"
  ) +
  theme_minimal()



