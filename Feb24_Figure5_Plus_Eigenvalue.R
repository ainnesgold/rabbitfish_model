library(tidyverse)
library(viridis)

# Parameters
burn_in_time <- 100  # Define burn-in period of 100 timesteps
time_steps <- 200 
juvenile_mortality <- 0.9 / 2
subadult_mortality <- 0.54 / 2
adult_mortality <- 0.33 / 2
juvenile_survival_rate <- 1 - juvenile_mortality
subadult_survival_rate <- 1 - subadult_mortality
adult_survival_rate <-  1 - adult_mortality

juvenile_to_subadult_rate <- 1
subadult_to_adult_rate <- 0.5

#reproduction_rate <- 2
#reproduction_rate_2 <- 1.5

carrying_capacity <- 68
time_steps <- 500

# Parameter ranges
fishing_effort_values <- seq(0, 0.5, by = 0.05)

# Define the range of restocking mortality rates
restocking_mortality_rate <- 0.2
reproduction_rates <- seq(1, 3.5, by = 0.5)
#reproduction_rates_2 <- reproduction_rates * 0.75

#restocking scenario
restocked_juveniles = 1

#Current fishing for burn in
F_current_instantaneous <- 0.09
F_current_discrete <- (1 - exp(-F_current_instantaneous))/2
F_current_discrete_juv <- 0.01





run_simulation <- function(F_adults, F_juveniles, restocked_juveniles, reproduction_rate, burn_in_time = 100) {
  juvenile_population <- 10
  subadult_population <- 10
  adult_population <- 10
  
  # Store population over time
  population_over_time <- numeric(time_steps)
  eigenvalues_over_time <- numeric(time_steps) 
  
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
      current_reproduction_rate <- reproduction_rate * 0.75
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
      
      #restocked_dagge_values[t] <- restocked_dagge  # Store restocked dagge value
      subadult_population <- max(0, new_subadults + surviving_subadults - new_adults + restocked_dagge)
    } else {
      juvenile_population <- max(0, new_juveniles + surviving_juveniles - new_subadults)
      subadult_population <- max(0, new_subadults + surviving_subadults - new_adults)
    }
    
    adult_population <- max(0, new_adults + surviving_adults)
    
    # Store total population at each timestep
    total_population <- juvenile_population + subadult_population + adult_population
    population_over_time[t] <- total_population 
    # Construct the population transition matrix
    transition_matrix <- matrix(c(
      0, 0, current_reproduction_rate - (restocking / carrying_capacity),  # Juvenile production
      juvenile_survival_rate * (1 - F_juveniles_current) * juvenile_to_subadult_rate + (restocking + (2.35 * (1 - restocking / carrying_capacity))) / carrying_capacity, 0, 0,  # Juvenile survival to subadult
      0, subadult_survival_rate * subadult_to_adult_rate, adult_survival_rate * (1 - F_adults_current)  # Subadult survival & transition
    ), nrow = 3, byrow = TRUE)
    
    # Compute the dominant eigenvalue (growth rate)
    eigenvalues <- eigen(transition_matrix)$values
    dominant_eigenvalue <- max(Re(eigenvalues))  # Take the real part
    eigenvalues_over_time[t] <- dominant_eigenvalue
    
  }
  
  return(list(
    avg_population_last_20 = mean(population_over_time[(time_steps - 19):time_steps]),
    avg_eigenvalue_last_20 = mean(eigenvalues_over_time[(time_steps - 19):time_steps])
  ))
}






# Sensitivity analysis loop
all_results <- data.frame()

for (reproduction_rate in reproduction_rates) {
  # Create a dataframe to store results for this restocking mortality rate
  sensitivity_results <- expand.grid(
    F_adults = fishing_effort_values,
    F_juveniles = fishing_effort_values
  ) %>%
    mutate(
      Avg_Total_Population = 0,
      Avg_Eigenvalue = 0  # New column for eigenvalues
    )
  
  # Reference population with zero fishing effort
  reference_population <- run_simulation(
    F_adults = 0,
    F_juveniles = 0,
    restocked_juveniles = 1,
    reproduction_rate = reproduction_rate,
    burn_in_time = 100
  )
  
  for (i in seq_len(nrow(sensitivity_results))) {
    sim_result <- run_simulation(
      F_adults = sensitivity_results$F_adults[i],
      F_juveniles = sensitivity_results$F_juveniles[i],
      restocked_juveniles = restocked_juveniles,
      reproduction_rate = reproduction_rate,
      burn_in_time = 100
    )
    
    sensitivity_results$Avg_Total_Population[i] <- sim_result$avg_population_last_20
    sensitivity_results$Avg_Eigenvalue[i] <- sim_result$avg_eigenvalue_last_20
  }
  
  sensitivity_results <- sensitivity_results %>%
    mutate(
      Relative_Population = Avg_Total_Population / reference_population$avg_population_last_20,
      Reproduction_Rate = reproduction_rate
    )
  
  all_results <- bind_rows(all_results, sensitivity_results)
}

# Custom labeller function for restocking mortality
reproduction_labeller <- function(values) {
  paste("Reproduction =", values)
}

#restocked_juveniles / reference_population

# Update the code to use this custom labeller
figure5<-ggplot(all_results, aes(x = F_juveniles*2, y = F_adults*2, fill = Relative_Population)) +
  geom_tile() +
  geom_contour(
    aes(z = Relative_Population), 
    breaks = 0.5, 
    color = "black", 
    linetype = "dashed", 
    linewidth = 1
  ) +
  facet_wrap(
    ~Reproduction_Rate, 
    labeller = labeller(Reproduction_Rate = function(values) reproduction_labeller(values))
  ) +
  scale_fill_gradientn(
    colors = c("#d73027", "#fee08b", "#1a9850"), # Color-blind friendly palette
    name = bquote("Biomass relative to"~B[0]),
  ) +
  labs(
    x = "Annual fishing effort on mañahak",
    y = "Annual fishing effort on hiteng kahlao",
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
figure5

#ggsave("~/Desktop/rabbitfish_figure5.png", figure5, width=10, height=8, bg="transparent")



figure5_eigen<-ggplot(all_results, aes(x = F_juveniles*2, y = F_adults*2, fill = Avg_Eigenvalue)) +
  geom_tile() +
  geom_contour(
    aes(z = Avg_Eigenvalue), 
    breaks = 1, 
    color = "black", 
    linetype = "dashed", 
    linewidth = 1
  ) +
  facet_wrap(
    ~Reproduction_Rate, 
    labeller = labeller(Reproduction_Rate = function(values) reproduction_labeller(values))
  ) +
  scale_fill_viridis_c(name = "Dominant Eigenvalue") +
  labs(
    x = "Annual fishing effort on mañahak",
    y = "Annual fishing effort on hiteng kahlao",
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

ggsave("~/Desktop/rabbitfish_figure5_eigenversion.png", figure5_eigen, width=10, height=8, bg="transparent")


