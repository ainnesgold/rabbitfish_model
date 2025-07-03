# Define parameters
juvenile_survival_rate <- 0.5  # Probability of a juvenile surviving to the next time step
adult_survival_rate <- 0.8     # Probability of an adult surviving to the next time step
juvenile_to_adult_rate <- 0.3  # Probability of a juvenile transitioning to adult
reproduction_rate <- 2.0       # Number of juveniles produced by each adult
carrying_capacity <- 500       # Carrying capacity of the environment
fishing_mortality_rate <- 0.2  # Proportion of adult population caught by fishing each time step
restocked_juveniles <- 20

# Initial population sizes
juvenile_population <- 100
adult_population <- 50

# Number of time steps to simulate
time_steps <- 20

# Vectors to store population sizes over time
juvenile_population_over_time <- numeric(time_steps)
adult_population_over_time <- numeric(time_steps)
harvest_over_time <- numeric(time_steps)


# Initialize population sizes
juvenile_population_over_time[1] <- juvenile_population
adult_population_over_time[1] <- adult_population

# Simulation loop
for (t in 2:time_steps) {
  # Calculate the total population size
  total_population <- juvenile_population + adult_population
  
  # Calculate the density-dependent factor
  density_factor <- 1 - total_population / carrying_capacity
  if (density_factor < 0) {
    density_factor <- 0
  }
  
  # Calculate the number of juveniles and adults in the next time step
  new_juveniles <- adult_population * reproduction_rate * density_factor
  surviving_juveniles <- juvenile_population * juvenile_survival_rate #* density_factor
  new_adults <- juvenile_population * juvenile_to_adult_rate #* density_factor
  surviving_adults <- adult_population * adult_survival_rate #* density_factor
  
  # Apply fishing mortality
  caught_adults <- adult_population * fishing_mortality_rate
  surviving_adults <- surviving_adults - caught_adults
  
  # Update population sizes
  juvenile_population <- new_juveniles + surviving_juveniles - new_adults + restocked_juveniles
  adult_population <- new_adults + surviving_adults
  
  # Store population sizes
  juvenile_population_over_time[t] <- juvenile_population
  adult_population_over_time[t] <- adult_population
  harvest_over_time[t] <- caught_adults
}

# Plot the results
plot(1:time_steps, juvenile_population_over_time, type = "o", col = "blue", ylim = c(0, max(juvenile_population_over_time, adult_population_over_time)), xlab = "Time", ylab = "Population Size", main = "Stage-Structured Fish Population Model with Carrying Capacity and Fishing")
lines(1:time_steps, adult_population_over_time, type = "o", col = "red")
lines(1:time_steps, harvest_over_time, type = "o", col = "black")
legend("topright", legend = c("Juveniles", "Adults", "Harvest"), col = c("blue", "red", "black"), lty = 1)



