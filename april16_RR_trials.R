library(ggplot2)

# Parameters
juvenile_mortality <- 0.9 / 2
subadult_mortality <- 0.54 / 2
adult_mortality <- 0.33 / 2

juvenile_survival_rate <- 1 - juvenile_mortality
subadult_survival_rate <- 1 - subadult_mortality
adult_survival_rate <- 1 - adult_mortality

juvenile_to_subadult_rate <- 1
subadult_to_adult_rate <- 0.5

carrying_capacity <- 68
time_steps <- 200
burn_in_time <- 100

# Burn-in fishing values
F_current_instantaneous <- 0.09 
F_current_discrete <- (1 - exp(-F_current_instantaneous)) / 2
F_current_discrete_juv <- 0.01

# Parameter ranges
reproduction_rate_values <- seq(0, 3, by = 0.5)
restocking_values <- seq(0, 5.5, by = 0.5)
fishing_scenarios <- expand.grid(F_simulation_adult = c(0, 0.5, 1),
                                 F_simulation_juv = c(0, 0.5, 1))

# Store results
all_results <- data.frame()

# Loop over fishing scenarios
for (fs in 1:nrow(fishing_scenarios)) {
  F_simulation_adult <- fishing_scenarios$F_simulation_adult[fs]
  F_simulation_juv <- fishing_scenarios$F_simulation_juv[fs]
  
  results <- expand.grid(reproduction_rate = reproduction_rate_values,
                         restocking = restocking_values)
  results$biomass <- NA
  results$eigenvalue <- NA
  results$F_adult <- F_simulation_adult
  results$F_juv <- F_simulation_juv
  
  for (i in 1:nrow(results)) {
    reproduction_rate <- results$reproduction_rate[i]
    reproduction_rate_2 <- 0.75 * reproduction_rate
    restocking <- results$restocking[i]
    
    juvenile <- numeric(time_steps)
    subadult <- numeric(time_steps)
    adult <- numeric(time_steps)
    
    juvenile[1] <- 10
    subadult[1] <- 10
    adult[1] <- 10
    
    for (t in 2:time_steps) {
      total_population <- juvenile[t-1] + subadult[t-1] + adult[t-1]
      
      if (t <= burn_in_time) {
        F_adults_current <- F_current_discrete
        F_juveniles_current <- F_current_discrete_juv
        restocked_fish <- 0
      } else {
        F_adults_current <- F_simulation_adult
        F_juveniles_current <- F_simulation_juv
        restocked_fish <- if (t %% 2 == 0) restocking else 0
      }
      
      if (total_population > 0) {
        reproduced_juveniles <- if (t %% 2 == 0) {
          adult[t-1] * reproduction_rate * (1 - total_population / carrying_capacity)
        } else {
          adult[t-1] * reproduction_rate_2 * (1 - total_population / carrying_capacity)
        }
      } else {
        reproduced_juveniles <- 0
      }
      
      juvenile[t] <- reproduced_juveniles +
        juvenile[t-1] * juvenile_survival_rate * (1 - F_juveniles_current) -
        juvenile_to_subadult_rate * juvenile[t-1]
      juvenile[t] <- max(juvenile[t], 0)
      
      restocked_subadult <- restocked_fish + (2.35 * (1 - restocked_fish / carrying_capacity))
      subadult[t] <- subadult[t-1] * subadult_survival_rate * (1 - F_adults_current) +
        juvenile_to_subadult_rate * juvenile[t-1] -
        subadult_to_adult_rate * subadult[t-1] +
        restocked_subadult
      subadult[t] <- max(subadult[t], 0)
      
      adult[t] <- adult[t-1] * adult_survival_rate * (1 - F_adults_current) +
        subadult_to_adult_rate * subadult[t-1]
      adult[t] <- max(adult[t], 0)
    }
    
    total_biomass <- juvenile + subadult + adult
    
    transition_matrix <- matrix(c(
      - (restocking / carrying_capacity), 0, reproduction_rate,
      juvenile_survival_rate * (1 - F_current_discrete_juv) * juvenile_to_subadult_rate + 
        (restocking + (2.35 * (1 - restocking / carrying_capacity))) / carrying_capacity, 0, 0,
      0, subadult_survival_rate * subadult_to_adult_rate, adult_survival_rate * (1 - F_current_discrete)
    ), nrow = 3, byrow = TRUE)
    
    eigenvalues <- eigen(transition_matrix)$values
    dominant_eigenvalue <- max(Re(eigenvalues))
    
    results$biomass[i] <- mean(total_biomass[(time_steps - 19):time_steps])
    results$eigenvalue[i] <- dominant_eigenvalue
  }
  
  all_results <- rbind(all_results, results)
}

# Plot: Biomass faceted by fishing effort
ggplot(all_results, aes(x = reproduction_rate, y = restocking, fill = biomass)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  facet_grid(F_juv ~ F_adult, labeller = label_both) +
  labs(title = "Biomass under Different Fishing Scenarios",
       x = "Reproduction Rate",
       y = "Restocking",
       fill = "Biomass") +
  theme_minimal()

# Plot: Population growth rate (lambda) faceted by fishing effort
ggplot(all_results, aes(x = reproduction_rate, y = restocking, fill = eigenvalue)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c() +
  facet_grid(F_juv ~ F_adult, labeller = label_both) +
  labs(title = "Population Growth Rate (λ) under Different Fishing Scenarios",
       x = "Reproduction Rate",
       y = "Restocking",
       fill = "λ") +
  theme_minimal()



