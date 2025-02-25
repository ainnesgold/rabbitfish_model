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
fishing_effort_values <- seq(0, 0.5, by = 0.05)
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
  
  # Store population over time and eigenvalues
  population_over_time <- numeric(time_steps)
  eigenvalues_over_time <- numeric(time_steps) 
  
  for (t in 1:time_steps) {
    if (t <= burn_in_time) {
      F_adults_current <- F_current_discrete
      F_juveniles_current <- F_current_discrete_juv
      restocking <- 0
    } else {
      F_adults_current <- F_adults
      F_juveniles_current <- F_juveniles
      restocking <- restocked_juveniles
    }
    
    current_reproduction_rate <- ifelse(t %% 2 == 0, reproduction_rate, reproduction_rate_2)
    
    total_population <- juvenile_population + subadult_population + adult_population
    new_juveniles <- adult_population * current_reproduction_rate * (1 - total_population / carrying_capacity)
    surviving_juveniles <- juvenile_population * juvenile_survival_rate * (1 - F_juveniles_current)
    new_subadults <- juvenile_population * juvenile_to_subadult_rate
    surviving_subadults <- subadult_population * subadult_survival_rate
    new_adults <- subadult_population * subadult_to_adult_rate
    surviving_adults <- adult_population * adult_survival_rate * (1 - F_adults_current)
    
    juvenile_population <- max(0, new_juveniles + surviving_juveniles - new_subadults)
    subadult_population <- max(0, new_subadults + surviving_subadults - new_adults)
    adult_population <- max(0, new_adults + surviving_adults)
    
    # Store total population at each timestep
    total_population <- juvenile_population + subadult_population + adult_population
    population_over_time[t] <- total_population 
    
    # Construct the population transition matrix
    transition_matrix <- matrix(c(
      0, 0, reproduction_rate + (restocking / carrying_capacity),  # Juvenile production
      juvenile_survival_rate * (1 - F_juveniles_current) * juvenile_to_subadult_rate, 0, 0,  # Juvenile survival to subadult
      0, subadult_survival_rate * subadult_to_adult_rate, adult_survival_rate * (1 - F_adults_current)  # Subadult survival & transition
    ), nrow = 3, byrow = TRUE)
    
    # Compute the dominant eigenvalue (growth rate)
    eigenvalues <- eigen(transition_matrix)$values
    dominant_eigenvalue <- max(Re(eigenvalues))  # Take the real part
    eigenvalues_over_time[t] <- dominant_eigenvalue
  }
  
  return(list(
    avg_population_last_20 = mean(population_over_time[(time_steps - 19):time_steps]),
    avg_eigenvalue_last_20 = mean(eigenvalues_over_time[(time_steps - 19):time_steps])  # Average last 20 timesteps
  ))
}



all_results <- data.frame()

for (restocked_juveniles in restocking_values) {
  sensitivity_results <- expand.grid(
    F_adults = fishing_effort_values,
    F_juveniles = fishing_effort_values
  ) %>%
    mutate(
      Avg_Total_Population = 0,
      Avg_Eigenvalue = 0  # New column for eigenvalues
    )
  
  reference_population <- run_simulation(
    F_adults = 0,
    F_juveniles = 0,
    restocked_juveniles = 0,
    burn_in_time = 100
  )
  
  for (i in seq_len(nrow(sensitivity_results))) {
    sim_result <- run_simulation(
      F_adults = sensitivity_results$F_adults[i],
      F_juveniles = sensitivity_results$F_juveniles[i],
      restocked_juveniles = restocked_juveniles,
      burn_in_time = 100
    )
    
    sensitivity_results$Avg_Total_Population[i] <- sim_result$avg_population_last_20
    sensitivity_results$Avg_Eigenvalue[i] <- sim_result$avg_eigenvalue_last_20
  }
  
  sensitivity_results <- sensitivity_results %>%
    mutate(
      Relative_Population = Avg_Total_Population / reference_population$avg_population_last_20,
      Restocking_Percent = round((restocked_juveniles / reference_population$avg_population_last_20) * 100, 0)
    )
  
  all_results <- bind_rows(all_results, sensitivity_results)
}

# Create labels with the correct expressions for each Restocking_Percent
all_results$Restocking_Label <- dplyr::case_when(
  all_results$Restocking_Percent == 0 ~ "0*'%'~of~B[0]~(0~g/m^2)",
  all_results$Restocking_Percent == 2 ~ "2*'%'~of~B[0]~(1~g/m^2)",
  all_results$Restocking_Percent == 6 ~ "5*'%'~of~B[0]~(3~g/m^2)",
  all_results$Restocking_Percent == 10 ~ "10*'%'~of~B[0]~(5.5~g/m^2)"
)

# Convert Restocking_Label to a factor and ensure its order aligns with Restocking_Percent
all_results$Restocking_Label <- factor(all_results$Restocking_Label,
                                       levels = c("0*'%'~of~B[0]~(0~g/m^2)", 
                                                  "2*'%'~of~B[0]~(1~g/m^2)", 
                                                  "5*'%'~of~B[0]~(3~g/m^2)",
                                                  "10*'%'~of~B[0]~(5.5~g/m^2)"))




eigen_heatmap<-ggplot(all_results %>% filter(F_adults <=0.3), aes(x = F_juveniles, y = F_adults, fill = Avg_Eigenvalue)) +
  geom_tile() +
  geom_contour(
    aes(z = Avg_Eigenvalue), 
    breaks = 1, 
    color = "black", 
    linetype = "dashed", 
    linewidth = 1
  ) +
  facet_wrap(~Restocking_Label, labeller = label_parsed) +
  scale_fill_viridis_c(name = "Dominant Eigenvalue") +
  labs(
    x = "Fishing effort on mañahak",
    y = "Fishing effort on hiteng kahlao",
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

ggsave("~/Desktop/rabbitfish_eigenheatmap.png", eigen_heatmap, width=6, height=6, bg="transparent")


#line plot version?
# Modify the fishing scenario labels by multiplying each number by 2
line_data <- all_results %>%
  group_by(Restocking_Percent, F_adults, F_juveniles) %>%
  summarize(
    Avg_Eigenvalue = mean(Avg_Eigenvalue),
    .groups = "drop"
  ) %>%
  mutate(Fishing_Scenario = paste0(F_adults, ", ", F_juveniles))


subset_fishing_scenarios <- c(
  "0, 0",
  "0, 0.5",
  "0.1, 0.25",
  "0.15, 0.3",
  "0.2, 0.2",
  "0.4, 0.4",
  "0.5, 0"
)

subset_fishing_scenarios_modified <- sapply(subset_fishing_scenarios, function(x) {
  nums <- as.numeric(strsplit(x, ",")[[1]]) * 2
  paste(nums, collapse = ", ")
})

# Prepare data for line plot, filtered by selected fishing scenarios
line_data_subset <- line_data %>%
  filter(Fishing_Scenario %in% subset_fishing_scenarios)

figure4b<-ggplot(line_data_subset, aes(x = Restocking_Percent, y = Avg_Eigenvalue, 
                             color = Fishing_Scenario, linetype = Fishing_Scenario, 
                             group = Fishing_Scenario)) +
  geom_line(size = 1.5) +  # Increase line thickness
  scale_color_manual(values = viridis::viridis(length(subset_fishing_scenarios)), 
                     labels = subset_fishing_scenarios_modified,
                     name = "Annual fishing effort on \nhiteng kahlao and mañahak") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "twodash", "longdash", "F1"), 
                        labels = subset_fishing_scenarios_modified,  # Use same labels as color
                        name = "Annual fishing effort on \nhiteng kahlao and mañahak") +  
  geom_hline(yintercept = 1, linetype = "solid", color = "red", linewidth = 1, alpha = 0.5) +  
  labs(
    x = bquote("Restocking (% of "~B[0]~")"),
    y = "Eigenvalue",
    title = ""
  ) +
  ggtitle("B.") +
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

figure4b




