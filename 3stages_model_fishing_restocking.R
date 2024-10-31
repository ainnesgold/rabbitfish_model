# Define parameters
juvenile_survival_rate <- 1 - (0.429 / 2)   # Probability of a juvenile surviving to the next time step
subadult_survival_rate <- 1 - (0.429 / 2)   # Probability of a subadult surviving to the next time step
adult_survival_rate <-  1 - (0.429 / 2)  # 0.429 is the annual mortality, divided in 2 for 6 months time steps - Probability of an adult surviving to the next time step

juvenile_to_subadult_rate <- 1 #0.2 # Probability of a juvenile transitioning to subadult
subadult_to_adult_rate <- 1 #0.2   # Probability of a subadult transitioning to adult

reproduction_rate <- 0.6 #0.19
reproduction_rate_2 <- 0.3 #0.03

carrying_capacity <- 50 #g/m2

#restocking
restocked_juveniles <- 0

# Fishing rates
juvenile_fishing_rate <- 0    # Fraction of juvenile population caught by fishing per time step
adult_fishing_rate <- 0      # Fraction of adult population caught by fishing per time step

# Initial population sizes in g/m2 - just made up values
juvenile_population <- 20
subadult_population <- 10
adult_population <- 5

# Number of time steps to simulate
time_steps <- 200

# Vectors to store population sizes over time
juvenile_population_over_time <- numeric(time_steps)
subadult_population_over_time <- numeric(time_steps)
adult_population_over_time <- numeric(time_steps)

# Initialize population sizes
juvenile_population_over_time[1] <- juvenile_population
subadult_population_over_time[1] <- subadult_population
adult_population_over_time[1] <- adult_population

# Simulation loop
for (t in 2:time_steps) {
  # Alternate reproduction rate based on the timestep (even or odd)
  if (t %% 2 == 0) {
    current_reproduction_rate <- reproduction_rate
  } else {
    current_reproduction_rate <- reproduction_rate_2
  }
  
  # Calculate the number of juveniles, subadults, and adults in the next time step
  total_population <- juvenile_population + subadult_population + adult_population
  new_juveniles <- adult_population * current_reproduction_rate * (1 - total_population / carrying_capacity)
  surviving_juveniles <- juvenile_population * juvenile_survival_rate * (1 - juvenile_fishing_rate)
  new_subadults <- juvenile_population * juvenile_to_subadult_rate
  surviving_subadults <- subadult_population * subadult_survival_rate
  new_adults <- subadult_population * subadult_to_adult_rate
  surviving_adults <- adult_population * adult_survival_rate * (1 - adult_fishing_rate)
  
  # Update population sizes
  juvenile_population <- new_juveniles + surviving_juveniles - new_subadults + restocked_juveniles
  subadult_population <- new_subadults + surviving_subadults - new_adults
  adult_population <- new_adults + surviving_adults
  
  # Store population sizes
  juvenile_population_over_time[t] <- juvenile_population
  subadult_population_over_time[t] <- subadult_population
  adult_population_over_time[t] <- adult_population
}

# Plot the results
plot(1:time_steps, juvenile_population_over_time, type = "o", col = "blue", 
     ylim = c(0, max(juvenile_population_over_time, subadult_population_over_time, adult_population_over_time)), 
     xlab = "Time", ylab = "Population Size") #main = "Stage-Structured Fish Population Model with Fishing & Restocking"
lines(1:time_steps, subadult_population_over_time, type = "o", col = "green")
lines(1:time_steps, adult_population_over_time, type = "o", col = "red")
legend("topright", legend = c("Juveniles", "Subadults", "Adults"), col = c("blue", "green", "red"), lty = 1)

























#Version with restocking once a year (instead of every 6 month time step)

# Vectors to store population sizes over time
juvenile_population_over_time <- numeric(time_steps)
subadult_population_over_time <- numeric(time_steps)
adult_population_over_time <- numeric(time_steps)

# Initialize population sizes
juvenile_population_over_time[1] <- juvenile_population
subadult_population_over_time[1] <- subadult_population
adult_population_over_time[1] <- adult_population

# Simulation loop
for (t in 2:time_steps) {
  # Alternate reproduction rate based on the timestep (even or odd)
  if (t %% 2 == 0) {
    current_reproduction_rate <- reproduction_rate
  } else {
    current_reproduction_rate <- reproduction_rate_2
  }
  
  # Calculate the number of juveniles, subadults, and adults in the next time step
  total_population <- juvenile_population + subadult_population + adult_population
  new_juveniles <- adult_population * current_reproduction_rate * (1 - total_population / carrying_capacity)
  surviving_juveniles <- juvenile_population * juvenile_survival_rate * (1 - juvenile_fishing_rate)
  new_subadults <- juvenile_population * juvenile_to_subadult_rate
  surviving_subadults <- subadult_population * subadult_survival_rate
  new_adults <- subadult_population * subadult_to_adult_rate
  surviving_adults <- adult_population * adult_survival_rate * (1 - adult_fishing_rate)
  
  # Update population sizes
  if (t %% 2 == 0) {
    # Add restocked juveniles every other timestep (on even timesteps)
    juvenile_population <- new_juveniles + surviving_juveniles - new_subadults + restocked_juveniles
  } else {
    # No restocked juveniles on odd timesteps
    juvenile_population <- new_juveniles + surviving_juveniles - new_subadults
  }
  
  subadult_population <- new_subadults + surviving_subadults - new_adults
  adult_population <- new_adults + surviving_adults
  
  # Store population sizes
  juvenile_population_over_time[t] <- juvenile_population
  subadult_population_over_time[t] <- subadult_population
  adult_population_over_time[t] <- adult_population
}


plot(1:time_steps, juvenile_population_over_time, type = "o", col = "blue", 
     ylim = c(0, max(juvenile_population_over_time, subadult_population_over_time, adult_population_over_time)), 
     xlab = "Time", ylab = "Population Size", main = "Stage-Structured Fish Population Model with Fishing & Restocking")
lines(1:time_steps, subadult_population_over_time, type = "o", col = "green")
lines(1:time_steps, adult_population_over_time, type = "o", col = "red")
legend("topright", legend = c("Juveniles", "Subadults", "Adults"), col = c("blue", "green", "red"), lty = 1)
