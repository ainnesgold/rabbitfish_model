library(tidyverse)
library(viridis)
library(purrr)
library(dplyr)
library(patchwork)

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
fishing_effort_values <- c(0, 0.5)
restocking_values <- c(0, 5.5) # Restocking scenarios

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
  eigenvectors_over_time_r <- vector("list", time_steps)
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
    eigenvectors <- eigen(transition_matrix)$vectors
    principle_eigenvector <- eigenvectors[,1]
    eigenvectors_over_time[[t]] <- principle_eigenvector
    
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
    eigenvectors_r <- eigen(transition_matrix_r)$vectors
    principle_eigenvector_r <- eigenvectors_r[,1]
    eigenvectors_over_time_r[[t]] <- principle_eigenvector_r
    
    #transition sensitivity
    
  }
  
  return(list(
    avg_population_last_20 = mean(population_over_time[(time_steps - 19):time_steps]),
    avg_eigenvalue_last_20 = mean(eigenvalues_over_time[(time_steps - 19):time_steps]),  # Average last 20 timesteps
    end_eigenvector = eigenvectors_over_time[[time_steps]],
    end_eigenvector_r = eigenvectors_over_time_r[[time_steps]]
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
      Eigenvector = vector("list", n()),
      Eigenvector_r = vector("list", n()),
      #Restocking_Percent = round((restocked_juveniles / reference_population$avg_population_last_20) * 100, 0),
      #scenario = paste0(
      #  "Restocking~'(% B0)':~", Restocking_Percent,
      #  "*','~F[H]==", F_adults*2,
      #  "*','~F[M]==", F_juveniles*2
      #)
      scenario = paste0(
        "f[H]==", F_adults * 2,
        "*','~f[M]==", F_juveniles * 2
      )
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
    sensitivity_results$Eigenvector_r[[i]] <- sim_result$end_eigenvector_r
  }
  
  sensitivity_results <- sensitivity_results %>%
    mutate(
      Relative_Population = Avg_Total_Population / reference_population$avg_population_last_20,
      Restocking_Percent = round((restocked_juveniles / reference_population$avg_population_last_20) * 100, 0)
    )
  
  all_results <- bind_rows(all_results, sensitivity_results)
}









###Visualizing population structure

all_results_long <- all_results %>%
  mutate(
    juvenile_contrib = map_dbl(Eigenvector_r, ~ Re(.x[1])),
    subadult_contrib = map_dbl(Eigenvector_r, ~ Re(.x[2])),
    adult_contrib = map_dbl(Eigenvector_r, ~ Re(.x[3]))
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

stabledist



#ggsave("~/Desktop/stabledist.png", stabledist, width=8, height=6, bg="transparent")










#################### New stuff - transition sensitivity #############################

all_results_long2 <- all_results_long

##Example looking at one matrix
# Get the first row's eigenvectors
v1 <- unlist(all_results_long2$Eigenvector[[1]])
v2 <- unlist(all_results_long2$Eigenvector_r[[1]])

#real part
v1_real <- Re(v1)
v2_real <- Re(v2)

# Outer product (matrix multiplication)
mat <- v1_real %*% t(v2_real)

# Normalize the matrix so it sums to 1
mat_norm <- mat / sum(mat)

# View the result
mat_norm

##### Do for all entries

all_results_long2$outer_products <- map2(all_results_long2$Eigenvector, all_results_long2$Eigenvector_r, ~ {
  v1 <- unlist(.x)
  v2 <- unlist(.y)
  
  # Extract real part of the eigenvectors
  v1_real <- Re(v1)
  v2_real <- Re(v2)
  
  # Flip sign of the entire vector if any real part is negative
  #if (any(v1_real < 0)) v1 <- -v1
  #if (any(v2_real < 0)) v2 <- -v2
  
  # Optional: Normalize by sum or max value
  # v1 <- v1 / sum(v1)
  # v2 <- v2 / sum(v2)
  
  # Outer product of eigenvectors
  mat <- v1 %*% t(v2)
  mat / sum(mat)  # Normalize so the sum equals 1
})



# Example: Visualize the first matrix (for example)
m <- all_results_long2$outer_products[[1]]

df_heat <- as.data.frame(Re(m)) %>%  # <--- take the real part
  mutate(Row = row_number()) %>%
  pivot_longer(-Row, names_to = "Col", values_to = "Value")


ggplot(df_heat, aes(x = Col, y = factor(Row), fill = Value)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Heatmap of Outer Product Matrix", x = "Column", y = "Row") +
  theme_minimal()

##only ones that are happening are P21, P32, P33
# Visualize a specific outer product matrix
# Visualize a specific outer product matrix
m <- all_results_long2$outer_products[[1]]

df_heat <- as.data.frame(Re(m)) %>%
  mutate(Row = row_number()) %>%
  pivot_longer(-Row, names_to = "Col", values_to = "Value") %>%
  mutate(
    Col = as.numeric(gsub("V", "", Col))
  ) %>%
  filter((Row == 1 & Col == 2) |
           (Row == 2 & Col == 3) |
           (Row == 3 & Col == 3)) %>%
  mutate(
    Label = case_when(
      Row == 1 & Col == 2 ~ "P21",
      Row == 2 & Col == 3 ~ "P32",
      Row == 3 & Col == 3 ~ "P33"
    )
  )

ggplot(df_heat, aes(x = Col, y = Row, fill = Value)) +
  geom_tile(color = "white") +
  geom_label(aes(label = Label), fill = "white", color = "black", size = 5, label.size = 0) +
  scale_fill_viridis_c() +
  scale_y_reverse(breaks = 1:3, labels = c("Mañahak", "Dagge", "Hiteng kahlao")) +
  scale_x_continuous(
    breaks = 1:3,
    labels = c("Mañahak", "Dagge", "Hiteng kahlao"),
    limits = c(1, 3.5)  # Force inclusion of Stage 1
  ) +
  labs(
    x = "Destination", y = "Source"
  ) +
  theme_minimal()




plots <- map2(
  all_results_long2$outer_products,
  all_results_long2$scenario,
  ~ {
    m <- .x
    scenario_name <- .y
    
    df_heat <- as.data.frame(Re(m)) %>%
      mutate(Row = row_number()) %>%
      pivot_longer(-Row, names_to = "Col", values_to = "Value") %>%
      mutate(
        Col = as.numeric(gsub("V", "", Col))
      ) %>%
      filter((Row == 1 & Col == 2) |
               (Row == 2 & Col == 3) |
               (Row == 3 & Col == 3)) %>%
      mutate(
        Label = case_when(
          Row == 1 & Col == 2 ~ "P21",
          Row == 2 & Col == 3 ~ "P32",
          Row == 3 & Col == 3 ~ "P33"
        )
      )
    
    ggplot(df_heat, aes(x = Col, y = Row, fill = Value)) +
      geom_tile(color = "white") +
      geom_label(aes(label = Label), fill = "white", color = "black", size = 5, label.size = 0) +
      scale_fill_viridis_c(limits = c(0.09, 0.21)) +  # Set the same color scale limits
      scale_y_reverse(breaks = 1:3, labels = c("Mañahak", "Dagge", "Hiteng kahlao")) +
      scale_x_continuous(
        breaks = 1:3,
        labels = c("Mañahak", "Dagge", "Hiteng kahlao"),
        limits = c(1, 3.5)
      ) +
      labs(
        title = parse(text = scenario_name),
        # Use the dynamic scenario name with bquote
        x = "Destination", y = "Source"
      ) +
      theme_minimal()
  }
)



# Create a single legend using a dummy plot (use same limits and color scale)
legend_plot <- ggplot(df_heat, aes(x = Col, y = Row, fill = Value)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0.09, 0.21)) +  # Match the scale limits to the others
  theme_minimal() +
  theme(legend.position = "bottom")  # Add legend at the bottom

# Combine the plots with a common legend using patchwork
combined_plot <- wrap_plots(plots[1:nrow(all_results_long2)]) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# View the combined plot with the common legend
combined_plot


#ggsave("~/Desktop/transition_plot.png", combined_plot, width=18, height=8, bg="transparent")




# Split the plots by restocking = 0 and restocking = 10

library(tidyverse)
library(patchwork)

# Create helper function for plotting
make_heatmap_plot <- function(m, scenario_name) {
  df_heat <- as.data.frame(Re(m)) %>%
    mutate(Row = row_number()) %>%
    pivot_longer(-Row, names_to = "Col", values_to = "Value") %>%
    mutate(
      Col = as.numeric(gsub("V", "", Col))
    ) %>%
    filter((Row == 1 & Col == 2) |
             (Row == 2 & Col == 3) |
             (Row == 3 & Col == 3)) %>%
    mutate(
      Label = case_when(
        Row == 1 & Col == 2 ~ "P21",
        Row == 2 & Col == 3 ~ "P32",
        Row == 3 & Col == 3 ~ "P33"
      )
    )
  
  ggplot(df_heat, aes(x = Col, y = Row, fill = Value)) +
    geom_tile(color = "white") +
    geom_label(aes(label = Label), fill = "white", color = "black", size = 5, label.size = 0) +
    scale_fill_viridis_c(limits = c(0.09, 0.205)) +
    scale_y_reverse(breaks = 1:3, labels = c("Mañahak", "Dagge", "Hiteng kahlao")) +
    scale_x_continuous(
      breaks = 1:3,
      labels = c("Mañahak", "Dagge", "Hiteng kahlao"),
      limits = c(1, 3.5)
    ) +
    labs(
      title = parse(text = scenario_name),
      x = "Destination", y = "Source"
    ) +
    theme_minimal()
}

# Filter and plot restocking = 0
df_0 <- all_results_long2 %>% filter(Restocking_Percent == 0)
plots_0 <- map2(df_0$outer_products, df_0$scenario, make_heatmap_plot)

# Filter and plot restocking = 10
df_10 <- all_results_long2 %>% filter(Restocking_Percent == 10)
plots_10 <- map2(df_10$outer_products, df_10$scenario, make_heatmap_plot)




# Create custom titles using plotmath-compatible text
title_0 <- ggplot() + 
  annotate("text", x = 1, y = 1, 
           label = "Restocking~('% B'[0]*'')*' = 0'", 
           parse = TRUE, size = 5) +
  theme_void() +
  theme(plot.margin = margin(0, 0, -10, 0))  # reduce spacing below

title_10 <- ggplot() + 
  annotate("text", x = 1, y = 1, 
           label = "Restocking~('% B'[0]*'')*' = 10'", 
           parse = TRUE, size = 5) +
  theme_void() +
  theme(plot.margin = margin(0, 0, -10, 0))


col_0 <- title_0 / wrap_plots(plots_0, ncol = 1) + plot_layout(heights = c(0.08, 1))
col_10 <- title_10 / wrap_plots(plots_10, ncol = 1) + plot_layout(heights = c(0.08, 1))

combined_plot <- (col_0 | col_10) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


combined_plot

# Save to file
#ggsave("~/Desktop/combined_plot.png", combined_plot, width = 10, height = 12, bg = "transparent")




#Pairwise plot showing values at 0 and 10% restocking

extract_transitions <- function(m, scenario, restocking) {
  df <- as.data.frame(Re(m)) %>%
    mutate(Row = row_number()) %>%
    pivot_longer(-Row, names_to = "Col", values_to = "Value") %>%
    mutate(Col = as.numeric(gsub("V", "", Col))) %>%
    filter((Row == 1 & Col == 2) | (Row == 2 & Col == 3) | (Row == 3 & Col == 3)) %>%
    mutate(Transition = case_when(
      Row == 1 & Col == 2 ~ "P21",
      Row == 2 & Col == 3 ~ "P32",
      Row == 3 & Col == 3 ~ "P33"
    )) %>%
    select(Transition, Value) %>%
    mutate(Scenario = scenario, Restocking = restocking)
  return(df)
}

# Apply to all rows
df_transitions_all <- pmap_dfr(
  list(m = all_results_long2$outer_products,
       scenario = all_results_long2$scenario,
       restocking = all_results_long2$Restocking_Percent),
  extract_transitions
)


label_map <- c(
  "f[H]==0*','~f[M]==0" = "A*':'~f[H]==0*','~f[M]==0",
  "f[H]==0*','~f[M]==1" = "B*':'~f[H]==0*','~f[M]==1",
  "f[H]==1*','~f[M]==0" = "C*':'~f[H]==1*','~f[M]==0",
  "f[H]==1*','~f[M]==1" = "D*':'~f[H]==1*','~f[M]==1"
)

transition_plot <- ggplot(df_transitions_all, aes(x = Restocking, y = Value, color = Transition, group = Transition)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~Scenario, labeller = as_labeller(label_map, label_parsed)) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = expression("Restocking (% " * B[0] * ")"),
    y = "Transition sensitivity value"
  ) +
  theme_minimal()


ggsave("~/Desktop/Fig3_transition_plot.png", transition_plot, width = 6, height = 6, bg = "transparent")









