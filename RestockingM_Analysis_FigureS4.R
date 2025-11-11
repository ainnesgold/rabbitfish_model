library(tidyverse)
library(viridis)
library(ggpubr)

# Parameters
burn_in_time <- 100 
time_steps <- 200 

#mortality
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
time_steps <- 500

# Parameter ranges
fishing_effort_values <- seq(0, 0.5, by = 0.01)
restocking_mortality_rates <- c(0, 0.4, 0.8, 1)

#restocking scenario
restocked_juveniles = 2.5

#Current fishing for burn in
F_current_instantaneous <- 0.09
F_current_discrete <- (1 - exp(-F_current_instantaneous)) / 2
F_current_discrete_juv <- 0.01



run_simulation <- function(F_adults, F_juveniles, restocked_juveniles, restocking_mortality_rate, burn_in_time = 100) {
  juvenile_population <- 10
  subadult_population <- 10
  adult_population <- 10
  

  population_over_time <- numeric(time_steps)
  eigenvalues_over_time <- numeric(time_steps) 
  
  for (t in 1:time_steps) {
    if (t <= burn_in_time) {
      # Burn-in period with zero fishing and zero restocking
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
      subadult_population <- max(0, new_subadults + surviving_subadults - new_adults + (restocked_dagge* (1 - restocking_mortality_rate)))
    } else {
      juvenile_population <- max(0, new_juveniles + surviving_juveniles - new_subadults)
      subadult_population <- max(0, new_subadults + surviving_subadults - new_adults)
    }
    
    adult_population <- max(0, new_adults + surviving_adults)
    

    total_population <- juvenile_population + subadult_population + adult_population
    population_over_time[t] <- total_population 

    transition_matrix <- matrix(c(
      - (restocking / carrying_capacity), 0, current_reproduction_rate,  # Juvenile production
      juvenile_survival_rate * (1 - F_juveniles_current) * juvenile_to_subadult_rate + (restocking + (2.35 * (1 - restocking / carrying_capacity))) / carrying_capacity, 0, 0,  # Juvenile survival to subadult
      0, subadult_survival_rate * subadult_to_adult_rate, adult_survival_rate * (1 - F_adults_current)  # Subadult survival & transition
    ), nrow = 3, byrow = TRUE)
    

    eigenvalues <- eigen(transition_matrix)$values
    dominant_eigenvalue <- max(Re(eigenvalues)) 
    eigenvalues_over_time[t] <- dominant_eigenvalue
    
  }
  
  return(list(
    avg_population_last_20 = mean(population_over_time[(time_steps - 19):time_steps]),
    avg_eigenvalue_last_20 = mean(eigenvalues_over_time[(time_steps - 19):time_steps])
  ))
}
















# Sensitivity analysis loop
all_results <- data.frame()

for (restocking_mortality_rate in restocking_mortality_rates) {

  sensitivity_results <- expand.grid(
    F_adults = fishing_effort_values,
    F_juveniles = fishing_effort_values
  ) %>%
    mutate(
      Avg_Total_Population = 0,
      Avg_Eigenvalue = 0
    )
  

  reference_population <- run_simulation(
    F_adults = 0,
    F_juveniles = 0,
    restocked_juveniles = 1,
    restocking_mortality_rate = restocking_mortality_rate,
    burn_in_time = 100
  )
  
  for (i in seq_len(nrow(sensitivity_results))) {
    sim_result <- run_simulation(
      F_adults = sensitivity_results$F_adults[i],
      F_juveniles = sensitivity_results$F_juveniles[i],
      restocked_juveniles = restocked_juveniles,
      restocking_mortality_rate = restocking_mortality_rate,
      burn_in_time = 100
    )
    
    sensitivity_results$Avg_Total_Population[i] <- sim_result$avg_population_last_20
    sensitivity_results$Avg_Eigenvalue[i] <- sim_result$avg_eigenvalue_last_20
  }
  
  sensitivity_results <- sensitivity_results %>%
    mutate(
      Relative_Population = Avg_Total_Population / reference_population$avg_population_last_20,
      Restocking_Mortality = restocking_mortality_rate
    )
  
  all_results <- bind_rows(all_results, sensitivity_results)
}

restocking_labeller <- function(values) {
  paste("Restocked M =", values)
}


figure6<-ggplot(all_results, aes(x = F_juveniles*2, y = F_adults*2, fill = Relative_Population)) +
  geom_tile() +
  geom_contour(
    aes(z = Relative_Population), 
    breaks = 0.5, 
    color = "black", 
    linetype = "dashed", 
    linewidth = 1
  ) +
  facet_wrap(
    ~Restocking_Mortality, 
    labeller = labeller(Restocking_Mortality = function(values) restocking_labeller(values))
  ) +
  scale_fill_gradientn(
    colors = c("#d73027", "#fee08b", "#1a9850"), 
    name = bquote("Biomass relative to"~B[0]),
  ) +
  labs(
    x = "Annual fishing effort on mañahak",
    y = "Annual fishing effort on hiteng kahlao",
    caption = bquote("Restocking = 2% of"~B[0]),
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


figure6

#ggsave("~/Desktop/rabbitfish_figure6.png", figure6, width=8, height=6, bg="transparent")


fig6eigen<-ggplot(all_results, aes(x = F_juveniles*2, y = F_adults*2, fill = Avg_Eigenvalue)) +
  geom_tile() +
  geom_contour(
    aes(z = Avg_Eigenvalue), 
    breaks = 1, 
    color = "black", 
    linetype = "dashed", 
    linewidth = 1
  ) +
  facet_wrap(
    ~Restocking_Mortality, 
    labeller = labeller(Restocking_Mortality = function(values) restocking_labeller(values))
  ) +
  scale_fill_viridis_c(name = "Dominant Eigenvalue") +
  labs(
    x = "Annual fishing effort on mañahak",
    y = "Annual fishing effort on hiteng kahlao",
    caption = bquote("Restocking = 2% of"~B[0]),
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )


#ggsave("~/Desktop/rabbitfish_figure6_eigenversion.png", fig6eigen, width=8, height=6, bg="transparent")



#COMBO version

M_biomass<-ggplot(all_results, aes(x = F_juveniles*2, y = F_adults*2, fill = Relative_Population)) +
  geom_tile() +
  geom_contour(
    aes(z = Relative_Population), 
    breaks = 0.5, 
    color = "black", 
    linetype = "dashed", 
    linewidth = 1
  ) +
  facet_wrap(
    ~Restocking_Mortality, 
    labeller = labeller(Restocking_Mortality = function(values) restocking_labeller(values)),
    nrow=1
  ) +
  scale_fill_gradientn(
    colors = c("#d73027", "#fee08b", "#1a9850"), # Color-blind friendly palette
    name = bquote("Biomass relative to"~B[0]),
  ) +
  labs(
    x = "Mañahak annual fishing mortality",
    y = "Hiteng kahlao annual fishing mortality",
    title = "A."
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

M_biomass

M_eigen<-ggplot(all_results, aes(x = F_juveniles*2, y = F_adults*2, fill = Avg_Eigenvalue)) +
  geom_tile() +
  geom_contour(
    aes(z = Avg_Eigenvalue), 
    breaks = 1, 
    color = "black", 
    linetype = "dashed", 
    linewidth = 1
  ) +
  facet_wrap(
    ~Restocking_Mortality, 
    labeller = labeller(Restocking_Mortality = function(values) restocking_labeller(values)),
    nrow=1
  ) +
  scale_fill_viridis_c(name = "Dominant eigenvalue") +
  labs(
    x = "Mañahak annual fishing mortality",
    y = "Hiteng kahlao annual fishing mortality",
    caption = bquote("Restocking = 5% of"~B[0]),
    title = "B."
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    plot.title = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_blank() 
  )

M_eigen

figures4combo<-ggarrange(M_biomass + rremove("xlab") + rremove("ylab"), M_eigen + rremove("ylab"), nrow=2)

figures4combo <- annotate_figure(
  figures4combo,
  left = text_grob("Hiteng kahlao annual fishing mortality", 
                   rot = 90, 
                   size = 16)
)


figures4combo

ggsave("~/Desktop/figures4.png", figures4combo, width=12, height=8, bg="transparent")





