# Load necessary libraries
library(MCMCpack)
library(ggplot2)
library(dplyr)
library(tidyr)

# Set seed for reproducibility
set.seed(123)

# Create a folder to save the results
results_folder <- "Bayesian_Polya_Urn_Results"
dir.create(results_folder, showWarnings = FALSE)

# Function to simulate the Bayesian Polya Urn process
simulate_bayesian_polya_urn <- function(k, alpha_0, beta, lambda, mu, t_max, known_initial = FALSE) {
  # Initialize counts and proportions
  if (known_initial) {
    N <- alpha_0
    theta <- N / sum(N)
  } else {
    theta <- rdirichlet(1, lambda * mu)
    N <- theta * alpha_0
  }
  
  # Store initial proportions for comparison
  initial_theta <- theta
  
  # Initialize storage for proportions and predictive distributions
  theta_t <- matrix(NA, nrow = t_max, ncol = k)
  theta_t[1, ] <- theta
  posterior_predictive_t <- matrix(NA, nrow = t_max, ncol = k)
  posterior_predictive_t[1, ] <- theta
  
  # Simulate the process
  for (t in 2:t_max) {
    # Draw a color based on current proportions
    color <- sample(1:k, size = 1, prob = theta)
    N[color] <- N[color] + beta
    
    # Update posterior distribution
    alpha <- lambda * mu + N
    theta <- rdirichlet(1, alpha)
    
    # Store the updated proportions and posterior predictive
    theta_t[t, ] <- theta
    posterior_predictive_t[t, ] <- alpha / sum(alpha)
  }
  
  return(list(theta_t = theta_t, N = N, posterior_predictive_t = posterior_predictive_t, initial_theta = initial_theta))
}

# Parameters
k <- 3  # Number of colors
alpha_0_known <- c(30, 50, 20)  # Known initial counts for the scenario where initial proportions are known
beta <- 1  # Reinforcement factor
lambda <- 2  # Concentration parameter for Dirichlet prior
mu <- c(1/3, 1/3, 1/3)  # Base measure for Dirichlet prior
t_max <- 100  # Number of iterations

# Simulate the process with known initial proportions
simulation_known <- simulate_bayesian_polya_urn(k, alpha_0_known, beta, lambda, mu, t_max, known_initial = TRUE)

# Simulate the process with unknown initial proportions
simulation_unknown <- simulate_bayesian_polya_urn(k, sum(alpha_0_known), beta, lambda, mu, t_max, known_initial = FALSE)

# Function to plot posterior densities with true initial parameters
plot_posterior_density <- function(simulation, title, initial_theta, known_initial = FALSE) {
  theta_df <- as.data.frame(simulation$theta_t)
  theta_df$Iteration <- 1:t_max
  theta_df <- pivot_longer(theta_df, cols = starts_with("V"), names_to = "Color", values_to = "Proportion")
  
  density_plot <- ggplot(theta_df, aes(x = Proportion, fill = Color)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = initial_theta, linetype = "dashed", color = "red", size = 1.5) +
    scale_fill_discrete(name = "Color", labels = paste("Color", 1:k)) +
    labs(title = title, x = "Proportion", y = "Density") +
    theme_minimal() +
    theme(legend.position = "top")
  
  # Save the plot
  file_name <- if (known_initial) {
    "Posterior_Density_Known.png"
  } else {
    "Posterior_Density_Unknown.png"
  }
  ggsave(file.path(results_folder, file_name), plot = density_plot, width = 8, height = 6)
  
  # Print the plot
  print(density_plot)
}

# Plot and save the posterior densities for known initial proportions
plot_posterior_density(simulation_known, "Posterior Density with Known Initial Proportions", simulation_known$initial_theta, known_initial = TRUE)

# Plot and save the posterior densities for unknown initial proportions
plot_posterior_density(simulation_unknown, "Posterior Density with Unknown Initial Proportions", simulation_unknown$initial_theta, known_initial = FALSE)

# Function to plot the evolution of posterior predictive distributions
plot_posterior_predictive <- function(simulation, title) {
  posterior_predictive_df <- as.data.frame(simulation$posterior_predictive_t)
  posterior_predictive_df$Iteration <- 1:t_max
  posterior_predictive_df <- pivot_longer(posterior_predictive_df, cols = starts_with("V"), names_to = "Color", values_to = "Probability")
  
  predictive_plot <- ggplot(posterior_predictive_df, aes(x = Iteration, y = Probability, color = Color)) +
    geom_line(size = 1) +
    scale_color_discrete(name = "Color", labels = paste("Color", 1:k)) +
    labs(title = title, x = "Iteration", y = "Probability") +
    theme_minimal() +
    theme(legend.position = "top")
  
  # Save the plot
  file_name <- paste0(title, ".png")
  ggsave(file.path(results_folder, file_name), plot = predictive_plot, width = 8, height = 6)
  
  # Print the plot
  print(predictive_plot)
}

# Plot and save the evolution of posterior predictive distributions
plot_posterior_predictive(simulation_known, "Posterior Predictive Evolution (Known Initial Proportions)")
plot_posterior_predictive(simulation_unknown, "Posterior Predictive Evolution (Unknown Initial Proportions)")

# Function to save final summary statistics
save_final_summary <- function(simulation, known_initial = FALSE) {
  final_theta <- simulation$theta_t[t_max, ]
  summary_table <- data.frame(
    Color = factor(1:k),
    Initial_Proportion = simulation$initial_theta,
    Final_Proportion = final_theta,
    Difference = final_theta - simulation$initial_theta
  )
  
  # Save summary table to CSV
  file_name <- if (known_initial) {
    "Final_Theta_Summary_Known.csv"
  } else {
    "Final_Theta_Summary_Unknown.csv"
  }
  write.csv(summary_table, file.path(results_folder, file_name), row.names = FALSE)
  
  # Print summary table
  print(summary_table)
}

# Save and print final summary statistics
save_final_summary(simulation_known, known_initial = TRUE)
save_final_summary(simulation_unknown, known_initial = FALSE)
