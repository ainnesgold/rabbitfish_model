library(tidyverse)
library(viridis)
library(ggpubr)

# Parameters

juvenile_mortality <- 0.9 / 2
subadult_mortality <- 0.54 / 2
adult_mortality <- 0.33 / 2
juvenile_survival_rate <- 1 - juvenile_mortality
subadult_survival_rate <- 1 - subadult_mortality
adult_survival_rate <-  1 - adult_mortality

juvenile_to_subadult_rate <- 1
subadult_to_adult_rate <- 0.5

reproduction_rate <- 2
reproduction_rate_2 <- 1.5

carrying_capacity <- 68
time_steps <- 200

# Parameter ranges
fishing_effort_values <- seq(0, 0.5, by = 0.05)
restocking_values <- seq(0, 5.5, by = 0.5) # Restocking scenarios

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
  restocked_dagge_values <- numeric(time_steps)  # To store restocked dagge values
  
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
      juvenile_population <- max(0, new_juveniles + surviving_juveniles - new_subadults - restocking)
      restocked_dagge <- restocking + (2.35 * (1 - restocking / carrying_capacity))
      
      restocked_dagge_values[t] <- restocked_dagge  # Store restocked dagge value
      subadult_population <- max(0, new_subadults + surviving_subadults - new_adults + restocked_dagge)
    } else {
      juvenile_population <- max(0, new_juveniles + surviving_juveniles - new_subadults)
      subadult_population <- max(0, new_subadults + surviving_subadults - new_adults)
    }
    
    adult_population <- max(0, new_adults + surviving_adults)
    
    # Store total population at each timestep
    total_population <- juvenile_population + subadult_population + adult_population
    population_over_time[t] <- total_population 
    
  }
  
  # Calculate the average restocked dagge for the last 20 timesteps
  # Return the average of the last 20 timesteps and the average restocked dagge
  return(list(
    avg_population_last_20 = mean(population_over_time[(time_steps - 19):time_steps]),
    avg_restocked_dagge_last_20 = mean(restocked_dagge_values[(time_steps - 19):time_steps])
  ))
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
      Avg_Total_Population = 0,
      Avg_Restocked_Dagge = 0  # New column for average restocked dagge
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
    sim_result <- run_simulation(
      F_adults = sensitivity_results$F_adults[i],
      F_juveniles = sensitivity_results$F_juveniles[i],
      restocked_juveniles = restocked_juveniles,
      burn_in_time = 100
    )
    
    sensitivity_results$Avg_Total_Population[i] <- sim_result$avg_population_last_20
    sensitivity_results$Avg_Restocked_Dagge[i] <- sim_result$avg_restocked_dagge_last_20*2  #times 2 because it averages it with 0 in the 6 months when there is no restocking
  }
  
  # Normalize by the reference population and calculate restocking as a percentage
  sensitivity_results <- sensitivity_results %>%
    mutate(
      Relative_Population = Avg_Total_Population / reference_population$avg_population_last_20,
      Restocking_Percent = round((restocked_juveniles / reference_population$avg_population_last_20) * 100, 0) # Convert to percentage
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




#alternative plot
#subset_fishing_scenarios <- c(
#  "0, 0",
#  "0, 0.5",
#  "0.1, 0.25",
#  "0.15, 0.3",
#  "0.2, 0.2",
#  "0.4, 0.4",
#  "0.5, 0"
#)


subset_fishing_scenarios <- c(
  "0, 0",
  "0, 0.5",
  "0.15, 0.3",
  "0.2, 0.2",
  "0.4, 0.4",
  "0.5, 0",
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
figure4a <- ggplot(line_data_subset, aes(x = Restocking_Percent, y = Avg_Relative_Population, 
                                         color = Fishing_Scenario, linetype = Fishing_Scenario, 
                                         group = Fishing_Scenario)) +
  geom_line(size = 1.5) +
  scale_color_manual(values = viridis::viridis(length(subset_fishing_scenarios)), 
                     labels = subset_fishing_scenarios_modified,
                     name = "Annual fishing effort on \nhiteng kahlao and mañahak") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "twodash", "longdash", "F1"), 
                        labels = subset_fishing_scenarios_modified,  # Use same labels as color
                        name = "Annual fishing effort on \nhiteng kahlao and mañahak") +  
  geom_hline(yintercept = 0.5, linetype = "solid", color = "red", linewidth = 1, alpha = 0.5) +  
  labs(
    x = bquote("Restocking (% of "~B[0]~")"),
    y = bquote("Biomass relative to"~B[0]),
    title = ""
  ) +
  ggtitle("A.") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 24),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key.width = unit(1.5, "cm")  # Adjust legend width if needed
  )
figure4a

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



source('Feb24_Eigenvalue_Analysis.R')

figure4<-ggarrange(figure4a, figure4b, common.legend=TRUE)

ggsave("~/Desktop/rabbitfish_figure4.png", figure4, width=8, height=6, bg="transparent")




