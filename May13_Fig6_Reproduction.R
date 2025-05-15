library(ggplot2)
library(dplyr)

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
reproduction_rate_values <- seq(0, 3, by = 0.05)
restocking_values <- seq(0, 5.5, by = 0.5)

# Fishing levels to iterate over
fishing_levels <- expand.grid(F_simulation_adult = c(0, 0.25, 0.5),
                              F_simulation_juv = c(0, 0.25, 0.5))

# Store results for all simulations
all_results <- data.frame()

# Loop over fishing combinations
for (f in 1:nrow(fishing_levels)) {
  F_simulation_adult <- fishing_levels$F_simulation_adult[f]
  F_simulation_juv <- fishing_levels$F_simulation_juv[f]
  
  results <- expand.grid(reproduction_rate = reproduction_rate_values,
                         restocking = restocking_values)
  results$biomass <- NA
  results$eigenvalue <- NA

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
    
    #transition matrix for population stability
    transition_matrix <- matrix(c(
      - (restocking / carrying_capacity), 0, reproduction_rate,
      juvenile_survival_rate * (1 - F_simulation_juv) * juvenile_to_subadult_rate + 
        (restocking + (2.35 * (1 - restocking / carrying_capacity))) / carrying_capacity, 0, 0,
      0, subadult_survival_rate * subadult_to_adult_rate, adult_survival_rate * (1 - F_simulation_adult)
    ), nrow = 3, byrow = TRUE)
    
    eigenvalues <- eigen(transition_matrix)$values
    dominant_eigenvalue <- max(Re(eigenvalues))
    
    
    #Save results
    results$biomass[i] <- mean(total_biomass[(time_steps - 19):time_steps])
    results$eigenvalue[i] <- dominant_eigenvalue

    
  }
  
  results$F_simulation_adult <- F_simulation_adult
  results$F_simulation_juv <- F_simulation_juv
  all_results <- rbind(all_results, results)
}

# Convert variables to factors with parsed expression labels
all_results$F_simulation_adult <- factor(all_results$F_simulation_adult,
                                         levels = c(0, 0.25, 0.5),
                                         labels = c("italic(f)[H]==0", "italic(f)[H]==0.5", "italic(f)[H]==1")
)

all_results$F_simulation_juv <- factor(all_results$F_simulation_juv,
                                       levels = c(0, 0.25, 0.5),
                                       labels = c("italic(f)[M]==0", "italic(f)[M]==0.5", "italic(f)[M]==1")
)

rr_biomass <- ggplot(all_results, aes(x = reproduction_rate, y = restocking, fill = biomass)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "Reproduction rate",
    y = expression("Restocking (g/m"^2*")"),
    fill = expression("Biomass (g/m"^2*")")) +
  facet_grid(F_simulation_juv ~ F_simulation_adult, labeller = label_parsed) +
  theme_minimal()


# #relative version - relative to same repduction rate, no restocking
# 
# # Step 1: Join each row with the biomass at restocking = 0 for the same grouping
# all_rel <- all_results %>%
#   group_by(reproduction_rate, F_simulation_juv, F_simulation_adult) %>%
#   mutate(baseline_biomass = biomass[restocking == 0],
#          relative_biomass = biomass / baseline_biomass) %>%
#   ungroup()
# 
# # Step 2: Plot the relative biomass
# rr_biomass_rel1 <- ggplot(all_rel, aes(x = reproduction_rate, y = restocking, fill = relative_biomass)) +
#   geom_tile() +
#   scale_fill_viridis_c() +
#   labs(x = "Reproduction rate",
#        y = expression("Restocking (g/m"^2*")"),
#        fill = "Relative biomass",
#        caption = "Biomass relative to restocking = 0 for the corresponding reproduction rate"
#   ) +
#   facet_grid(F_simulation_juv ~ F_simulation_adult, labeller = label_parsed) +
#   theme_minimal()
# 
# #ggsave("~/Desktop/rr_biomass_rel1.png", rr_biomass_rel1, width=6, height=4, bg="transparent")
# 
# 
# #other relative version where they are rall relative to restocking = 0 and reproduction = 0
# 
# # Step 1: Get baseline where restocking == 0 and reproduction_rate == 0
# baseline <- all_results %>%
#   filter(restocking == 0, reproduction_rate == 0) %>%
#   select(F_simulation_juv, F_simulation_adult, baseline_biomass = biomass)
# 
# # Step 2: Join baseline to all data based on facet groupings
# all_rel_base00 <- all_results %>%
#   left_join(baseline, by = c("F_simulation_juv", "F_simulation_adult")) %>%
#   mutate(relative_biomass = biomass / baseline_biomass)
# 
# # Step 3: Plot
# rr_biomass_rel2<- ggplot(all_rel_base00, aes(x = reproduction_rate, y = restocking, fill = relative_biomass)) +
#   geom_tile() +
#   scale_fill_viridis_c() +
#   labs(x = "Reproduction rate",
#        y = expression("Restocking (g/m"^2*")"),
#        fill = "Relative biomass",
#        caption = "Biomass relative to restocking = 0, reproduction = 0") +
#   facet_grid(F_simulation_juv ~ F_simulation_adult, labeller = label_parsed) +
#   theme_minimal()

#ggsave("~/Desktop/rr_biomass_rel2.png", rr_biomass_rel2, width=6, height=4, bg="transparent")






# Eigen plot - stability
rr_eigen <- ggplot(all_results, aes(x = reproduction_rate, y = restocking, fill = eigenvalue)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(x = "Reproduction rate",
       y = expression("Restocking (g/m"^2*")"),
       fill = "λ") +
  facet_grid(F_simulation_juv ~ F_simulation_adult, labeller = label_parsed) +
  theme_minimal()



#Saving
ggsave("~/Desktop/Fig6_rr_biomass.png", rr_biomass, width=6, height=4, bg="transparent")
ggsave("~/Desktop/FigS2_rr_eigen.png", rr_eigen, width=6, height=4, bg="transparent")



