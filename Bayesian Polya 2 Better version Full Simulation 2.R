# Load necessary libraries
library(MCMCpack)  # For Dirichlet distribution sampling
library(ggplot2)   # For plotting
library(dplyr)     # For data manipulation
library(tidyr)     # For reshaping data

# Set seed for reproducibility
set.seed(123)

# Create a folder to save the results
results_folder <- "Bayesian_Polya_Urn_Results"
dir.create(results_folder, showWarnings = FALSE)

# Function to simulate the Bayesian Polya Urn process
simulate_bayesian_polya_urn <- function(k, true_theta, beta, lambda, mu, t_max) {
  # Initialize counts based on the true proportions
  N <- round(true_theta * 100)  # Starting counts based on the true proportions
  theta <- N / sum(N)  # Calculate starting proportions
  
  # Initialize storage for posterior samples
  posterior_samples <- matrix(NA, nrow = t_max, ncol = k)
  
  # Run MCMC to sample from the posterior
  for (t in 1:t_max) {
    # Draw a color based on current proportions
    color <- sample(1:k, size = 1, prob = theta)
    N[color] <- N[color] + beta  # Reinforce the drawn color
    
    # Update posterior distribution using Dirichlet
    alpha <- lambda * mu + N
    theta <- rdirichlet(1, alpha)
    
    # Store the posterior sample
    posterior_samples[t, ] <- theta
  }
  
  return(list(posterior_samples = posterior_samples, true_theta = true_theta))
}

# Parameters
k <- 3  # Number of colors
true_theta <- c(0.3, 0.5, 0.2)  # True initial proportions for the simulation
beta <- 1  # Reinforcement factor
lambda <- 2  # Concentration parameter for Dirichlet prior
mu <- c(1/3, 1/3, 1/3)  # Base measure for Dirichlet prior
t_max <- 10000  # Number of iterations (increased to improve MCMC convergence)

# Simulate the Bayesian Polya Urn process
simulation <- simulate_bayesian_polya_urn(k, true_theta, beta, lambda, mu, t_max)

# Function to plot posterior densities with true parameter values
plot_posterior_density <- function(simulation, title) {
  theta_df <- as.data.frame(simulation$posterior_samples)
  theta_df <- pivot_longer(theta_df, cols = everything(), names_to = "Color", values_to = "Proportion")
  
  # Define color mapping for the true proportions
  true_color_mapping <- c("red", "green", "blue")
  
  density_plot <- ggplot(theta_df, aes(x = Proportion, fill = Color)) +
    geom_density(alpha = 0.6) +
    geom_vline(aes(xintercept = simulation$true_theta[1]), color = true_color_mapping[1], linetype = "solid", size = 1.5) +
    geom_vline(aes(xintercept = simulation$true_theta[2]), color = true_color_mapping[2], linetype = "solid", size = 1.5) +
    geom_vline(aes(xintercept = simulation$true_theta[3]), color = true_color_mapping[3], linetype = "solid", size = 1.5) +
    scale_fill_manual(values = true_color_mapping) +
    labs(title = title, x = "Proportion", y = "Density") +
    theme_minimal() +
    theme(legend.position = "top")
  
  # Save the plot
  file_name <- paste0("Posterior_Density_", title, ".png")
  ggsave(file.path(results_folder, file_name), plot = density_plot, width = 8, height = 6)
  
  # Print the plot
  print(density_plot)
}

# Plot and save the posterior densities
plot_posterior_density(simulation, "Posterior Density of Proportions with True Values")

# Function to evaluate and save final summary statistics
save_final_summary <- function(simulation) {
  final_posterior_mean <- colMeans(simulation$posterior_samples)
  summary_table <- data.frame(
    Color = factor(1:k),
    True_Proportion = simulation$true_theta,
    Posterior_Mean = final_posterior_mean,
    Bias = final_posterior_mean - simulation$true_theta
  )
  
  # Save summary table to CSV
  file_name <- "Final_Theta_Summary.csv"
  write.csv(summary_table, file.path(results_folder, file_name), row.names = FALSE)
  
  # Print summary table
  print(summary_table)
}

# Save and print final summary statistics
save_final_summary(simulation)
