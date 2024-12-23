library(tidyverse)

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

# Fishing rates
F_instantaneous <- 0.09
F_discrete <- 1 - exp(-F_instantaneous)
F_discrete

juvenile_fishing_rate <- 0.01 #(0.09/4)/2   # Fraction of juvenile population caught by fishing per time step
adult_fishing_rate <- F_discrete/2     # Fraction of adult population caught by fishing per time step

# Initial population sizes in g/m2 - made up values
juvenile_population <- 10
subadult_population <- 10
adult_population <- 10
total_population <- juvenile_population + subadult_population + adult_population

# Number of time steps to simulate
time_steps <- 100

# Vectors to store population sizes over time
juvenile_population_over_time <- numeric(time_steps)
subadult_population_over_time <- numeric(time_steps)
adult_population_over_time <- numeric(time_steps)
total_population_over_time <- numeric(time_steps)

# Initialize population sizes
juvenile_population_over_time[1] <- juvenile_population
subadult_population_over_time[1] <- subadult_population
adult_population_over_time[1] <- adult_population
total_population_over_time[1] <- total_population


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
    subadult_population <- max(0, new_subadults + surviving_subadults - new_adults + restocked_juveniles)
  } else {
    subadult_population <- max(0, new_subadults + surviving_subadults - new_adults)
  }
  
  juvenile_population <- max(0, new_juveniles + surviving_juveniles - new_subadults)
  adult_population <- max(0, new_adults + surviving_adults)
  #Store
  juvenile_population_over_time[t] <- juvenile_population
  subadult_population_over_time[t] <- subadult_population
  adult_population_over_time[t] <- adult_population
  total_population_over_time[t] <- total_population
}

# Plot the results
#plot(1:time_steps, juvenile_population_over_time, type = "o", col = "blue", 
#     ylim = c(0, max(juvenile_population_over_time, subadult_population_over_time, adult_population_over_time, total_population_over_time)), 
#     xlab = "Time", ylab = "Population Size") #main = "Stage-Structured Fish Population Model with Fishing & Restocking"
#lines(1:time_steps, subadult_population_over_time, type = "o", col = "green")
#lines(1:time_steps, adult_population_over_time, type = "o", col = "red")
#lines(1:time_steps, total_population_over_time, type = "o", col = "black")
#legend("topright", legend = c("Juveniles", "Subadults", "Adults", "Total"), col = c("blue", "green", "red", "black"), lty = 1)


population_data <- data.frame(
  Time = rep(1:time_steps, 4),
  Population = c(juvenile_population_over_time, subadult_population_over_time, adult_population_over_time, total_population_over_time),
  Stage = rep(c("Juveniles", "Subadults", "Adults", "Total"), each = time_steps)
)


ggplot(population_data, aes(x = Time/2, y = Population, color = Stage)) +
  geom_line(alpha = 0.3) +
  geom_smooth(se = FALSE, size = 1) +
  labs(
    x = "Years",
    y = bquote("Fish density"~(g/m^2)),
    color = "Stage",
    title = "Burn in period: No restocking, \nF hiteng kahlao = 0.086/year, F mañahak = 0.02/year"
  ) +
  scale_color_manual(
    values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), # Okabe-Ito palette
    labels = c("Hiteng kahlao (adults)", "Mañahak (juveniles)", "Dagge (subadults)", "Total")
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )




#Relative to total
# Preprocess the data: add a fraction column
population_data_relative <- population_data %>%
  group_by(Time) %>%
  mutate(Fraction = Population / Population[Stage == "Total"])

# Create the plot
ggplot(population_data_relative, aes(x = Time/2, y = Fraction, color = Stage)) +
  # Transparent lines for raw data
  geom_line(alpha = 0.3) +
  # Solid best fit lines for each stage
  geom_smooth(se = FALSE, size = 1) +
  labs(
    x = "Years",
    y = "Fish density (fraction of total population)",
    color = "Stage",
    title = "Burn in period: No restocking, \nF hiteng kahlao = 0.086/year, F mañahak = 0.02/year"
  ) +
  scale_color_manual(
    values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442"), # Okabe-Ito palette
    labels = c("Hiteng kahlao (adults)", "Mañahak (juveniles)", "Dagge (subadults)", "Total")
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 18),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


