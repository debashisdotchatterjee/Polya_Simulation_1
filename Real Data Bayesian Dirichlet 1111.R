# Install and load necessary packages
required_packages <- c("vegan", "MCMCpack", "ggplot2", "dplyr", "tidyr", "gridExtra", "viridis")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if(length(new_packages)) install.packages(new_packages)

library(vegan)     # For ecological data and analysis
library(MCMCpack)  # For Dirichlet distribution sampling
library(ggplot2)   # For plotting
library(dplyr)     # For data manipulation
library(tidyr)     # For reshaping data
library(gridExtra) # For arranging multiple plots
library(viridis)   # For color scales

# Set seed for reproducibility
set.seed(123)

# Create a folder to save the results
results_folder <- "Bayesian_Polya_Urn_Real_Data_Results"
if(!dir.exists(results_folder)) {
  dir.create(results_folder)
}

# Load the dune dataset
data(dune)

# Sum the abundances of species across all sites to get a categorical distribution
species_abundance <- colSums(dune)

# Normalize the counts to get initial proportions
initial_proportions <- species_abundance / sum(species_abundance)

# Define parameters for the Bayesian Polya Urn process
k <- length(initial_proportions)  # Number of species
species_names <- names(initial_proportions)  # Species names
true_theta <- initial_proportions  # Use the observed proportions as the true values
beta <- 1  # Reinforcement factor
lambda <- 2  # Concentration parameter for Dirichlet prior
mu <- rep(1/k, k)  # Base measure for Dirichlet prior (non-informative)
t_max <- 10000  # Fixed number of iterations

# Assign fixed colors to each species
species_colors <- viridis::viridis_pal(option = "D")(k)
names(species_colors) <- species_names

# Function to simulate the Bayesian Polya Urn process
simulate_bayesian_polya_urn <- function(k, true_theta, beta, lambda, mu, t_max) {
  # Initialize counts based on the true proportions
  N <- round(true_theta * 100)  # Starting counts based on the true proportions
  theta <- N / sum(N)  # Calculate starting proportions
  
  # Initialize storage for posterior samples
  posterior_samples <- matrix(NA, nrow = t_max, ncol = k)
  colnames(posterior_samples) <- species_names
  
  # Run MCMC to sample from the posterior
  for (t in 1:t_max) {
    # Draw a species based on current proportions
    species <- sample(species_names, size = 1, prob = theta)
    N[species] <- N[species] + beta  # Reinforce the drawn species
    
    # Update posterior distribution using Dirichlet
    alpha <- lambda * mu + N
    theta <- as.numeric(rdirichlet(1, alpha))
    names(theta) <- species_names
    
    # Store the posterior sample
    posterior_samples[t, ] <- theta
  }
  
  return(posterior_samples)
}

# Simulate the Bayesian Polya Urn process with fixed t_max
simulation <- simulate_bayesian_polya_urn(k, true_theta, beta, lambda, mu, t_max)

# Function to plot posterior densities with true parameter values, with max 3 species per plot
plot_posterior_density <- function(simulation, title, species_subset) {
  theta_df <- as.data.frame(simulation)
  theta_df <- theta_df[, species_subset, drop = FALSE]  # Select only the specified species
  theta_df <- pivot_longer(theta_df, cols = everything(), names_to = "Species", values_to = "Proportion")
  
  # Add True Proportion
  theta_df <- theta_df %>%
    mutate(True_Proportion = true_theta[Species])
  
  # Generate the density plot
  density_plot <- ggplot(theta_df, aes(x = Proportion, fill = Species, color = Species)) +
    geom_density(alpha = 0.6, linewidth = 1) +
    geom_vline(aes(xintercept = True_Proportion), linetype = "dashed", size = 1.2) +
    scale_fill_manual(values = species_colors) +
    scale_color_manual(values = species_colors) +
    labs(title = paste(title),
         x = "Proportion", y = "Density") +
    theme_minimal() +
    theme(legend.position = "right")
  
  # Save the plot
  file_name <- paste0("Posterior_Density_", title, "_Species_", paste(species_subset, collapse = "_"), ".png")
  ggsave(filename = file.path(results_folder, file_name), plot = density_plot, width = 8, height = 6)
  
  # Print the plot
  print(density_plot)
}

# Function to plot convergence with max 3 species per plot
plot_convergence_faceted_rowwise <- function(simulation, species_subset) {
  simulation_df <- as.data.frame(simulation[, species_subset, drop = FALSE])
  simulation_df$Iteration <- 1:t_max
  
  # Pivot longer for plotting
  convergence_df <- pivot_longer(simulation_df, cols = species_subset, names_to = "Species", values_to = "Proportion")
  
  # Add True Proportion
  convergence_df <- convergence_df %>%
    mutate(True_Proportion = true_theta[Species])
  
  # Generate the convergence plot
  convergence_plot <- ggplot(convergence_df, aes(x = Iteration, y = Proportion, color = Species)) +
    geom_line(alpha = 0.7) +
    geom_hline(aes(yintercept = True_Proportion, color = Species), linetype = "dashed", size = 1.2) +
    scale_color_manual(values = species_colors) +
    labs(title = paste("Convergence of Posterior Means"),
         x = "Iteration", y = "Posterior Mean") +
    facet_wrap(~ Species, scales = "free_y", ncol = 1) +  # Row-wise layout with one column
    theme_minimal() +
    theme(legend.position = "right")
  
  # Save the plot
  file_name <- paste0("Posterior_Convergence_Species_", paste(species_subset, collapse = "_"), ".png")
  ggsave(filename = file.path(results_folder, file_name), plot = convergence_plot, width = 8, height = 12)
  
  # Print the plot
  print(convergence_plot)
}

# Function to evaluate and save final summary statistics and model fit for all species
save_final_summary <- function(simulation) {
  final_posterior_mean <- colMeans(simulation)
  
  # Create a data frame with formatted numbers
  summary_table <- data.frame(
    Species = species_names,
    True_Proportion = format(true_theta, digits = 4, nsmall = 4),
    Posterior_Mean = format(final_posterior_mean, digits = 4, nsmall = 4),
    Bias = format(final_posterior_mean - true_theta, digits = 4, nsmall = 4)
  )
  
  # Save summary table to CSV
  file_name <- "Final_Theta_Summary.csv"
  write.csv(summary_table, file.path(results_folder, file_name), row.names = FALSE)
  
  # Print summary table
  print(summary_table)
  
  # Model Fit: Mean Squared Error (MSE)
  mse <- mean((final_posterior_mean - true_theta)^2)
  cat("Mean Squared Error (MSE) of the model fit: ", format(mse, digits = 4, nsmall = 4), "\n")
}

# Define subsets of species to plot (max 3 per plot)
species_subsets <- split(species_names, ceiling(seq_along(species_names)/3))

# Plot and save the posterior densities for each subset of species
for (subset in species_subsets) {
  plot_posterior_density(simulation, "Species_Abundance", species_subset = subset)
}

# Plot and save the row-wise faceted convergence of posterior means for each subset of species
for (subset in species_subsets) {
  plot_convergence_faceted_rowwise(simulation, species_subset = subset)
}

# Save and print final summary statistics and model fit for all species
save_final_summary(simulation)

#######################

# Install and load necessary packages
required_packages <- c("vegan", "MCMCpack", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

if(length(new_packages)) install.packages(new_packages)

library(vegan)     # For ecological data and analysis
library(MCMCpack)  # For Dirichlet distribution sampling
library(ggplot2)   # For plotting

# Set seed for reproducibility
set.seed(123)

# Create a folder to save the results
results_folder <- "Bayesian_Polya_Urn_Results"
if(!dir.exists(results_folder)) {
  dir.create(results_folder)
}

# Load the dune dataset
data(dune)

# Sum the abundances of species across all sites to get a categorical distribution
species_abundance <- colSums(dune)

# Normalize the counts to get initial proportions
true_theta <- species_abundance / sum(species_abundance)

# Define parameters for the Bayesian Polya Urn process
k <- length(true_theta)  # Number of species
species_names <- names(true_theta)  # Species names
beta <- 1  # Reinforcement factor
lambda <- 2  # Concentration parameter for Dirichlet prior
mu <- rep(1/k, k)  # Base measure for Dirichlet prior (non-informative)
t_max <- 10000  # Fixed number of iterations

# Function to simulate the Bayesian Polya Urn process
simulate_bayesian_polya_urn <- function(k, true_theta, beta, lambda, mu, t_max) {
  # Initialize counts based on the true proportions
  N <- round(true_theta * 100)  # Starting counts based on the true proportions
  theta <- N / sum(N)  # Calculate starting proportions
  
  # Initialize storage for posterior samples
  posterior_samples <- matrix(NA, nrow = t_max, ncol = k)
  colnames(posterior_samples) <- species_names
  
  # Run MCMC to sample from the posterior
  for (t in 1:t_max) {
    # Draw a species based on current proportions
    species <- sample(species_names, size = 1, prob = theta)
    N[species] <- N[species] + beta  # Reinforce the drawn species
    
    # Update posterior distribution using Dirichlet
    alpha <- lambda * mu + N
    theta <- as.numeric(rdirichlet(1, alpha))
    names(theta) <- species_names
    
    # Store the posterior sample
    posterior_samples[t, ] <- theta
  }
  
  return(posterior_samples)
}

# Simulate the Bayesian Polya Urn process with fixed t_max
simulation <- simulate_bayesian_polya_urn(k, true_theta, beta, lambda, mu, t_max)

# Calculate the posterior mean for each species
posterior_mean <- colMeans(simulation)

# Create a data frame for plotting and summarizing, formatted without scientific notation
plot_data <- data.frame(
  Species = species_names,
  True_Proportion = formatC(true_theta, format = "f", digits = 6),
  Posterior_Mean = formatC(posterior_mean, format = "f", digits = 6),
  Bias = formatC(posterior_mean - true_theta, format = "f", digits = 6)
)

# Save summary table to CSV
summary_table_file <- file.path(results_folder, "Posterior_Mean_Summary.csv")
write.csv(plot_data, summary_table_file, row.names = FALSE)

# Print the summary table
print(plot_data)

# Plot Posterior Mean vs Actual Proportion with 45-degree line
posterior_plot <- ggplot(plot_data, aes(x = as.numeric(True_Proportion), y = as.numeric(Posterior_Mean), label = Species)) +
  geom_point(size = 3, color = "blue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +
  labs(title = "Posterior Mean vs Actual Proportion for All Species",
       x = "Actual Proportion",
       y = "Posterior Mean") +
  theme_minimal() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# Save the plot
posterior_plot_file <- file.path(results_folder, "Posterior_Mean_vs_Actual_Proportion.png")
ggsave(filename = posterior_plot_file, plot = posterior_plot, width = 10, height = 8)

# Print the plot
print(posterior_plot)

