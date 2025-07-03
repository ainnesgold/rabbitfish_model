# Parameter ranges
restocking_values <- seq(0, 5.5, by = 0.5)
juv_to_sub_values <- seq(0.1, 1, by = 0.1)
sub_to_adult_values <- seq(0.1, 1, by = 0.1)

# Fixed reproduction rate (if needed)
reproduction_rate <- 1
reproduction_rate_2 <- 0.75 * reproduction_rate

# Fishing levels to iterate over
fishing_levels <- expand.grid(F_simulation_adult = c(0, 0.25, 0.5),
                              F_simulation_juv = c(0, 0.25, 0.5))

# Store results
all_results <- data.frame()

for (f in 1:nrow(fishing_levels)) {
  F_simulation_adult <- fishing_levels$F_simulation_adult[f]
  F_simulation_juv <- fishing_levels$F_simulation_juv[f]
  
  results <- expand.grid(restocking = restocking_values,
                         juvenile_to_subadult_rate = juv_to_sub_values,
                         subadult_to_adult_rate = sub_to_adult_values)
  results$biomass <- NA
  results$eigenvalue <- NA
  
  for (i in 1:nrow(results)) {
    restocking <- results$restocking[i]
    juvenile_to_subadult_rate <- results$juvenile_to_subadult_rate[i]
    subadult_to_adult_rate <- results$subadult_to_adult_rate[i]
    
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
    
    # transition matrix for population stability
    transition_matrix <- matrix(c(
      0, 0, reproduction_rate,
      juvenile_survival_rate * (1 - F_simulation_juv) * juvenile_to_subadult_rate, 0, 0,
      0, subadult_survival_rate * subadult_to_adult_rate, adult_survival_rate * (1 - F_simulation_adult)
    ), nrow = 3, byrow = TRUE)
    
    eigenvalues <- eigen(transition_matrix)$values
    dominant_eigenvalue <- max(Re(eigenvalues))
    
    results$biomass[i] <- mean(total_biomass[(time_steps - 19):time_steps])
    results$eigenvalue[i] <- dominant_eigenvalue
  }
  
  results$F_simulation_adult <- F_simulation_adult
  results$F_simulation_juv <- F_simulation_juv
  all_results <- rbind(all_results, results)
}





# Plots: raw biomass
ggplot(all_results, aes(x = juvenile_to_subadult_rate, y = restocking, fill = biomass)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "Juvenile → Subadult Rate",
       y = expression("Restocking (g/m"^2*")"),
       fill = "Biomass") +
  facet_grid(F_simulation_juv ~ F_simulation_adult, labeller = label_parsed) +
  theme_minimal()


ggplot(all_results, aes(x = subadult_to_adult_rate, y = restocking, fill = biomass)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "Subadult → Adult Rate",
       y = expression("Restocking (g/m"^2*")"),
       fill = "Biomass") +
  facet_grid(F_simulation_juv ~ F_simulation_adult, labeller = label_parsed) +
  theme_minimal()


#Plots: eigenvalue
ggplot(all_results, aes(x = juvenile_to_subadult_rate, y = restocking, fill = eigenvalue)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "Juvenile → Subadult Rate",
       y = expression("Restocking (g/m"^2*")"),
       fill = "Eigenvalue") +
  facet_grid(F_simulation_juv ~ F_simulation_adult, labeller = label_parsed) +
  theme_minimal()


ggplot(all_results, aes(x = subadult_to_adult_rate, y = restocking, fill = eigenvalue)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "Subadult → Adult Rate",
       y = expression("Restocking (g/m"^2*")"),
       fill = "Eigenvalue") +
  facet_grid(F_simulation_juv ~ F_simulation_adult, labeller = label_parsed) +
  theme_minimal()


### Relative biomass versions

##version 1 - relative to bottom of each column

relative_results_1 <- all_results %>%
  group_by(F_simulation_adult, F_simulation_juv,
           juvenile_to_subadult_rate, subadult_to_adult_rate) %>%
  mutate(baseline_biomass = biomass[restocking == 0],
         rel_biomass = biomass / baseline_biomass) %>%
  ungroup()

ggplot(relative_results_1, aes(x = juvenile_to_subadult_rate, y = restocking, fill = rel_biomass)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Relative Biomass", option = "viridis", limits = c(0, NA)) +
  labs(x = "Juvenile → Subadult Rate",
       y = expression("Restocking (g/m"^2*")")) +
  facet_grid(F_simulation_juv ~ F_simulation_adult, labeller = label_parsed) +
  theme_minimal()


ggplot(relative_results_1, aes(x = subadult_to_adult_rate, y = restocking, fill = rel_biomass)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Relative Biomass", option = "viridis", limits = c(0, NA)) +
  labs(x = "Subadult → Adult Rate",
       y = expression("Restocking (g/m"^2*")")) +
  facet_grid(F_simulation_juv ~ F_simulation_adult, labeller = label_parsed) +
  theme_minimal()





###Relative biomass version 2 - relative to bottom left square

baseline_lookup <- all_results %>%
  filter(restocking == 0,
         juvenile_to_subadult_rate == min(juvenile_to_subadult_rate),
         subadult_to_adult_rate == min(subadult_to_adult_rate)) %>%
  rename(baseline_biomass = biomass) %>%
  select(F_simulation_adult, F_simulation_juv, baseline_biomass)

# Join with full results
relative_results_2 <- all_results %>%
  left_join(baseline_lookup,
            by = c("F_simulation_adult", "F_simulation_juv")) %>%
  mutate(rel_biomass = biomass / baseline_biomass)

ggplot(relative_results_2, aes(x = juvenile_to_subadult_rate, y = restocking, fill = rel_biomass)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Relative Biomass", option = "viridis", limits = c(0, NA)) +
  labs(x = "Juvenile → Subadult Rate",
       y = expression("Restocking (g/m"^2*")")) +
  facet_grid(F_simulation_juv ~ F_simulation_adult, labeller = label_parsed) +
  theme_minimal()


ggplot(relative_results_2, aes(x = subadult_to_adult_rate, y = restocking, fill = rel_biomass)) +
  geom_tile() +
  scale_fill_viridis_c(name = "Relative Biomass", option = "viridis", limits = c(0, NA)) +
  labs(x = "Subadult → Adult Rate",
       y = expression("Restocking (g/m"^2*")")) +
  facet_grid(F_simulation_juv ~ F_simulation_adult, labeller = label_parsed) +
  theme_minimal()



