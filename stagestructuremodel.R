# Define parameters
juvenile_survival_rate <- 0.5  # Probability of a juvenile surviving to the next time step
adult_survival_rate <- 0.8     # Probability of an adult surviving to the next time step
juvenile_to_adult_rate <- 0.1  # Probability of a juvenile transitioning to adult
reproduction_rate <- 2.0       # Number of juveniles produced by each adult
carrying_capacity <- 500

# Initial population sizes
juvenile_population <- 100
adult_population <- 50

# Number of time steps to simulate
time_steps <- 100

# Vectors to store population sizes over time
juvenile_population_over_time <- numeric(time_steps)
adult_population_over_time <- numeric(time_steps)

# Initialize population sizes
juvenile_population_over_time[1] <- juvenile_population
adult_population_over_time[1] <- adult_population

# Simulation loop
for (t in 2:time_steps) {
  # Calculate the number of juveniles and adults in the next time step
  total_population <- juvenile_population + adult_population
  new_juveniles <- adult_population * reproduction_rate * (1 - total_population / carrying_capacity)
  surviving_juveniles <- juvenile_population * juvenile_survival_rate 
  new_adults <- juvenile_population * juvenile_to_adult_rate
  surviving_adults <- adult_population * adult_survival_rate
  
  # Update population sizes
  juvenile_population <- new_juveniles + surviving_juveniles - new_adults
  adult_population <- new_adults + surviving_adults
  
  # Store population sizes
  juvenile_population_over_time[t] <- juvenile_population
  adult_population_over_time[t] <- adult_population
}

# Plot the results
plot(1:time_steps, juvenile_population_over_time, type = "o", col = "blue", ylim = c(0, max(juvenile_population_over_time, adult_population_over_time)), xlab = "Time", ylab = "Population Size", main = "Stage-Structured Fish Population Model")
lines(1:time_steps, adult_population_over_time, type = "o", col = "red")
legend("topright", legend = c("Juveniles", "Adults"), col = c("blue", "red"), lty = 1)

