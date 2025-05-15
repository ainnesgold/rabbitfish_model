library(tidyverse)
library(viridis)

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
fishing_effort_values <- seq(0, 0.5, by = 0.01)
restocking_values <- c(0, 1, 3, 5.5) # Restocking scenarios

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

ggplot(all_results, aes(x = F_juveniles, y = F_adults, fill = Relative_Population)) +
  geom_tile() +
  facet_wrap(~Restocking_Percent)


# Set Restocking_Percent as a factor with the desired order
#all_results$Restocking_Percent <- factor(all_results$Restocking_Percent, 
 #                                        levels = c(0, 2, 5, 10))

all_results$Restocking_Label <- dplyr::case_when(
  all_results$Restocking_Percent == 0 ~ "A: 0*'%'~of~B[0]~(0~g/m^2)",
  all_results$Restocking_Percent == 2 ~ "B: 2*'%'~of~B[0]~(1~g/m^2)",
  all_results$Restocking_Percent == 5 ~ "C: 5*'%'~of~B[0]~(3~g/m^2)",
  all_results$Restocking_Percent == 10 ~ "D: 10*'%'~of~B[0]~(5.5~g/m^2)"
)

# Convert Restocking_Label to a factor and ensure its order aligns with Restocking_Percent
all_results$Restocking_Label <- factor(all_results$Restocking_Label,
                                       levels = c("A: 0*'%'~of~B[0]~(0~g/m^2)", 
                                                  "B: 2*'%'~of~B[0]~(1~g/m^2)", 
                                                  "C: 5*'%'~of~B[0]~(3~g/m^2)",
                                                  "D: 10*'%'~of~B[0]~(5.5~g/m^2)"))



ggplot(all_results, aes(x = F_juveniles*2, y = F_adults*2, fill = Relative_Population)) +
  geom_tile() +
  geom_contour(
    aes(z = Relative_Population), 
    breaks = 0.5, 
    color = "black", 
    linetype = "dashed", 
    linewidth = 1
  ) +
  facet_wrap(~Restocking_Label, labeller = label_parsed) +
  scale_fill_gradientn(
    colors = c("#d73027", "#fee08b", "#1a9850"),
    name = bquote("Biomass relative to"~B[0])
  ) +
  labs(
    x = "Annual fishing effort on mañahak",
    y = "Annual fishing effort on hiteng kahlao",
    title = "Restocking amount"
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


#Zooming in

figure4<-ggplot(all_results %>% filter(F_adults > 0.3), aes(x = F_juveniles*2, y = F_adults*2, fill = Relative_Population)) +
  geom_tile() +
  geom_contour(
    aes(z = Relative_Population), 
    breaks = 0.5, 
    color = "black", 
    linetype = "dashed", 
    linewidth = 1
  ) +
  facet_wrap(~Restocking_Label, labeller = label_parsed) +
  scale_fill_gradientn(
    colors = c("#d73027", "#fee08b", "#1a9850"),
    name = bquote("Biomass relative to"~B[0])
  ) +
  labs(
    x = " Mañahak annual fishing mortality",
    y = "Hiteng kahlao annual fishing mortality",
    title = "Restocking amount"
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


ggsave("~/Desktop/Fig4_Heatmap.png", figure4, width=10, height=8, bg="transparent")




#trying a different version of the figure


# Filter to avoid zeroes that might confuse contouring
#contour_data <- all_results %>% 
#  filter(Relative_Population >= 0.01)

#figure3_alt <- ggplot(contour_data, aes(x = F_juveniles * 2, y = F_adults * 2)) +
#  geom_contour(
#    aes(z = Relative_Population, color = factor(Restocking_Percent), group = Restocking_Percent),
#    breaks = c(0.5),  # or use more like: breaks = c(0.25, 0.5, 0.75)
#    linewidth = 1
#  ) +
#  scale_color_viridis_d(name = bquote("Restocking ("~"% B"["0"]~")")) +
#  labs(
#    x = "Mañahak annual fishing mortality",
#    y = "Hiteng kahlao annual fishing mortality",
#    title = bquote("Maximum fishing mortality that sustains 50%"~B[0])
#  ) +
#  theme_minimal() +
#  theme(
#    text = element_text(size = 16),
#    plot.title = element_text(size = 16, hjust = 0.5),
#    axis.title = element_text(size = 16),
#    axis.text = element_text(size = 14),
#    legend.title = element_text(size = 14),
#    legend.text = element_text(size = 12)
#  )

#ggsave("~/Desktop/rabbitfish_figure3_alt.png", figure3_alt, width=7, height=6, bg="transparent")

