library(tidyverse)
library(viridis)
library(purrr)
library(dplyr)

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
fishing_effort_values <- c(0, 0.25, 0.5)
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
  eigenvectors_over_time <- vector("list", time_steps)
  restocked_dagge_values <- numeric(time_steps) 
  
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
    
    
    # Construct the population transition matrix
    transition_matrix <- matrix(c(
      #0, 0, current_reproduction_rate - (restocking / carrying_capacity),  # Juvenile production
      - (restocking / carrying_capacity), 0, current_reproduction_rate,  # Juvenile production
      juvenile_survival_rate * (1 - F_juveniles_current) * juvenile_to_subadult_rate + (restocking + (2.35 * (1 - restocking / carrying_capacity))) / carrying_capacity, 0, 0,  # Juvenile survival to subadult
      0, subadult_survival_rate * subadult_to_adult_rate, adult_survival_rate * (1 - F_adults_current)  # Subadult survival & transition
    ), nrow = 3, byrow = TRUE)
    
    # Compute the dominant eigenvalue (growth rate)
    eigenvalues <- eigen(transition_matrix)$values
    dominant_eigenvalue <- max(Re(eigenvalues))  # Take the real part
    eigenvalues_over_time[t] <- dominant_eigenvalue
    
    
    #Second matrix for reproductive contribution
    #transition matrix for reproductive contribution
    transition_matrix_r <- matrix(c(
      - (restocking / carrying_capacity), juvenile_survival_rate * (1 - F_juveniles_current) * juvenile_to_subadult_rate + 
        (restocking + (2.35 * (1 - restocking / carrying_capacity))) / carrying_capacity, 0,
      0, 0, subadult_survival_rate * subadult_to_adult_rate,
      current_reproduction_rate, 0, adult_survival_rate * (1 - F_adults_current)
    ), nrow=3, byrow = TRUE)
    
    # Compute the dominant eigenvalue (growth rate)
    eigenvectors <- eigen(transition_matrix_r)$vectors
    principle_eigenvector <- eigenvectors[,1]
    eigenvectors_over_time[[t]] <- principle_eigenvector
    
    
  }
  
  return(list(
    avg_population_last_20 = mean(population_over_time[(time_steps - 19):time_steps]),
    avg_eigenvalue_last_20 = mean(eigenvalues_over_time[(time_steps - 19):time_steps]),  # Average last 20 timesteps
    end_eigenvector = eigenvectors_over_time[[time_steps]]
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
      Avg_Eigenvalue = 0,
      Eigenvector = vector("list", n())
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
    sensitivity_results$Eigenvector[[i]] <- sim_result$end_eigenvector
  }
  
  sensitivity_results <- sensitivity_results %>%
    mutate(
      Relative_Population = Avg_Total_Population / reference_population$avg_population_last_20,
      Restocking_Percent = round((restocked_juveniles / reference_population$avg_population_last_20) * 100, 0)
    )
  
  all_results <- bind_rows(all_results, sensitivity_results)
}









###Visualizing eigenvectors

all_results_long <- all_results %>%
  mutate(
    juvenile_contrib = map_dbl(Eigenvector, ~ Re(.x[1])),
    subadult_contrib = map_dbl(Eigenvector, ~ Re(.x[2])),
    adult_contrib = map_dbl(Eigenvector, ~ Re(.x[3]))
  )




all_results_long <- all_results_long %>%
  rowwise() %>%
  mutate(
    total = juvenile_contrib + subadult_contrib + adult_contrib,
    juvenile_prop = juvenile_contrib / total,
    subadult_prop = subadult_contrib / total,
    adult_prop = adult_contrib / total
  ) %>%
  ungroup()


# Create custom facet labels using expressions
all_results_long$F_adults <- factor(all_results_long$F_adults,
                                         levels = c(0, 0.25, 0.5),
                                         labels = c("italic(f)[H]==0", "italic(f)[H]==0.5", "italic(f)[H]==1")
)

all_results_long$F_juveniles <- factor(all_results_long$F_juveniles,
                                       levels = c(0, 0.25, 0.5),
                                       labels = c("italic(f)[M]==0", "italic(f)[M]==0.5", "italic(f)[M]==1")
)


# Plot
stabledist <- all_results_long %>%
  pivot_longer(cols = c(juvenile_prop, subadult_prop, adult_prop),
               names_to = "stage", values_to = "proportion") %>%
  mutate(
    stage = factor(stage, levels = c("juvenile_prop", "subadult_prop", "adult_prop"),
                   labels = c("Mañahak", "Dagge", "Hiteng kahlao"))
  ) %>%
  ggplot(aes(x = factor(Restocking_Percent), y = proportion, fill = stage)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = round(proportion, 2)),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  facet_grid(F_juveniles ~ F_adults, labeller = label_parsed) +
  labs(x = bquote("Restocking (% of "~B[0]~")"), y = "Stable population convergence",
       fill = "Life stage") +
  theme_minimal() +
  theme(strip.text = element_text(face = "italic")) +
  scale_fill_brewer(palette = "Set2")



all_results_long %>%
  group_by(F_adults, F_juveniles, Restocking_Percent) %>%
  summarise(across(ends_with("_prop"), mean))



#ggsave("~/Desktop/stabledist.png", stabledist, width=8, height=6, bg="transparent")












