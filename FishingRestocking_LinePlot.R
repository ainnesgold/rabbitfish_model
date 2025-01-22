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
time_steps <- 200

# Parameter ranges
fishing_effort_values <- seq(0, 0.5, by = 0.05)
restocking_values <- seq(0, 5, by = 0.5) # Restocking scenarios

#Current fishing for burn in
F_current_instantaneous <- 0.09 
F_current_discrete <- (1 - exp(-F_current_instantaneous)) / 2
F_current_discrete

F_current_discrete_juv <- 0.01


run_simulation <- function(F_adults, F_juveniles, restocked_juveniles, burn_in_time = 100) {
  juvenile_population <- 10
  subadult_population <- 10
  adult_population <- 10
  
  # Store population over time
  population_over_time <- numeric(time_steps)
  
  for (t in 1:time_steps) {
    if (t <= burn_in_time) {
      # Burn-in period with zero fishing and zero restocking
      # Try burn in with "current" fishing
      F_adults_current <- F_current_discrete #0
      F_juveniles_current <- F_current_discrete_juv #0
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
  
  # Normalize by the reference population and calculate restocking as a percentage
  sensitivity_results <- sensitivity_results %>%
    mutate(
      Relative_Population = Avg_Total_Population / reference_population,
      Restocking_Percent = round((restocked_juveniles / reference_population) * 100, 0) # Convert to percentage
    )
  
  # Combine results for all restocking values
  all_results <- bind_rows(all_results, sensitivity_results)
}


line_data <- all_results %>%
  group_by(Restocking_Percent, F_adults, F_juveniles) %>%
  summarize(
    Avg_Relative_Population = mean(Relative_Population),
    .groups = "drop"
  ) %>%
  mutate(Fishing_Scenario = paste0(F_adults, ", ", F_juveniles))


subset_fishing_scenarios <- c(
  "0, 0",
  "0.1, 0.05",
  "0.2, 0.1",
  "0.3, 0.15",
  "0.4, 0.2",
  "0.5, 0.25",
  "0.5, 0.5"
)

# Prepare data for line plot, filtered by selected fishing scenarios
line_data_subset <- line_data %>%
  filter(Fishing_Scenario %in% subset_fishing_scenarios)

# Create the plot for the subset
ggplot(line_data_subset, aes(x = Restocking_Percent, y = Avg_Relative_Population, color = Fishing_Scenario, group = Fishing_Scenario)) +
  geom_line(size = 1) +
  geom_point(data = filter(line_data_subset, Avg_Relative_Population == 0.5), aes(x = Restocking_Percent, y = Avg_Relative_Population), 
             color = "black", shape = 21, fill = "yellow", size = 3) +
  scale_color_viridis_d(name = "Fishing scenario \n(F hiteng kahlao, F mañahak)") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    x = bquote("Restocking (% of "~B[0]~")"),
    y = bquote("Biomass relative to"~B[0]),
    title = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )



#adjusting so that labels are annual F instead of /6 months
# Modify the fishing scenario labels by multiplying each number by 2
subset_fishing_scenarios <- c(
  "0, 0",
  "0.1, 0.05",
  "0.2, 0.1",
  "0.3, 0.15",
  "0.4, 0.2",
  "0.5, 0.25",
  "0.5, 0.5"
)

# Convert the labels into numeric pairs, multiply by 2, and then convert back to string
subset_fishing_scenarios_modified <- sapply(subset_fishing_scenarios, function(x) {
  # Split the string into numbers, multiply by 2, and then rejoin as a string
  nums <- as.numeric(strsplit(x, ",")[[1]]) * 2
  paste(nums, collapse = ", ")
})

# Prepare data for line plot, filtered by selected fishing scenarios
line_data_subset <- line_data %>%
  filter(Fishing_Scenario %in% subset_fishing_scenarios)

# Create the plot for the subset
ggplot(line_data_subset, aes(x = Restocking_Percent, y = Avg_Relative_Population, color = Fishing_Scenario, group = Fishing_Scenario)) +
  geom_line(size = 1) +
  geom_point(data = filter(line_data_subset, Avg_Relative_Population == 0.5), aes(x = Restocking_Percent, y = Avg_Relative_Population), 
             color = "black", shape = 21, fill = "yellow", size = 3) +
  scale_color_manual(values = viridis::viridis(length(subset_fishing_scenarios)), 
                     labels = subset_fishing_scenarios_modified,
                     name = "Annual fishing effort \n(F hiteng kahlao, F mañahak)") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    x = bquote("Restocking (% of "~B[0]~")"),
    y = bquote("Biomass relative to"~B[0]),
    title = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )



##alternative plot
# Modify the fishing scenario labels by multiplying each number by 2
subset_fishing_scenarios <- c(
  "0, 0",
  "0.1, 0.05",
  "0.1, 0.1",
  "0.2, 0.1",
  "0.2, 0.2",
  "0.3, 0.15",
  "0.3, 0.3",
  "0.4, 0.2",
  "0.4, 0.4",
  "0.5, 0.25",
  "0.5, 0.5"
)

# Convert the labels into numeric pairs, multiply by 2, and then convert back to string
subset_fishing_scenarios_modified <- sapply(subset_fishing_scenarios, function(x) {
  # Split the string into numbers, multiply by 2, and then rejoin as a string
  nums <- as.numeric(strsplit(x, ",")[[1]]) * 2
  paste(nums, collapse = ", ")
})

# Prepare data for line plot, filtered by selected fishing scenarios
line_data_subset <- line_data %>%
  filter(Fishing_Scenario %in% subset_fishing_scenarios)

# Create the plot for the subset
ggplot(line_data_subset, aes(x = Restocking_Percent, y = Avg_Relative_Population, color = Fishing_Scenario, group = Fishing_Scenario)) +
  geom_line(size = 1) +
  geom_point(data = filter(line_data_subset, Avg_Relative_Population == 0.5), aes(x = Restocking_Percent, y = Avg_Relative_Population), 
             color = "black", shape = 21, fill = "yellow", size = 3) +
  scale_color_manual(values = viridis::viridis(length(subset_fishing_scenarios)), 
                     labels = subset_fishing_scenarios_modified,
                     name = "Annual fishing effort on \nhiteng kahlao and mañahak") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.8) +
  labs(
    x = bquote("Restocking (% of "~B[0]~")"),
    y = bquote("Biomass relative to"~B[0]),
    title = ""
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )











#Quantifying differences in fishing effort
# Calculate sustainable fishing effort for adults
sustainable_fishing_adults <- all_results %>%
  filter(Relative_Population >= 0.5) %>%
  group_by(Restocking_Percent) %>%
  summarize(
    Max_F_Adults = max(F_adults),
    Max_F_Juveniles = max(F_juveniles), # Maximum sustainable fishing effort for adults
    .groups = "drop"
  )

# Extract the baseline sustainable adult fishing effort when restocking is zero
baseline_fishing_adults <- sustainable_fishing_adults %>%
  filter(Restocking_Percent == 0) %>%
  pull(Max_F_Adults)

# Calculate the percentage increase in adult fishing effort for each restocking level
sustainable_fishing_adults <- sustainable_fishing_adults %>%
  mutate(
    Percent_Increase_From_Zero = ifelse(
      Restocking_Percent == 0, 
      0, # No increase for restocking = 0
      ((Max_F_Adults - baseline_fishing_adults) / baseline_fishing_adults) * 100
    )
  )

# Print the table
print(sustainable_fishing_adults)







