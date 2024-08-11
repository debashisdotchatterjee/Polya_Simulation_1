# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to simulate Pólya urn process
polya_urn <- function(a, b, c, n_draws) {
  # Initialize the urn
  urn <- c(a, b)
  draws <- numeric(n_draws)
  
  for (i in 1:n_draws) {
    # Draw a ball based on current urn proportions
    prob <- urn / sum(urn)
    draw <- sample(c(1, 2), size = 1, prob = prob)
    
    # Record the draw (1 = A, 2 = B)
    draws[i] <- draw
    
    # Reinforce the urn by adding c balls of the drawn color
    urn[draw] <- urn[draw] + c
  }
  
  return(draws)
}

# Simulation parameters
a <- 10    # Initial number of balls of color A
b <- 10    # Initial number of balls of color B
c <- 1     # Number of additional balls added each time
n_draws <- 1000 # Number of draws

# Run the simulation
set.seed(123) # Set seed for reproducibility
draws <- polya_urn(a, b, c, n_draws)

# Calculate the proportion of color A over time
prop_A <- cumsum(draws == 1) / 1:n_draws

# Plot the evolution of the proportion of color A
df <- data.frame(Draw = 1:n_draws, Prop_A = prop_A)
ggplot(df, aes(x = Draw, y = Prop_A)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = a / (a + b), linetype = "dashed", color = "red") +
  labs(title = "Proportion of Color A over Time",
       x = "Number of Draws",
       y = "Proportion of Color A") +
  theme_minimal()

# Simulate multiple urns to study the distribution of the final proportion
n_simulations <- 1000
final_props <- numeric(n_simulations)

for (i in 1:n_simulations) {
  draws <- polya_urn(a, b, c, n_draws)
  final_props[i] <- sum(draws == 1) / n_draws
}

# Compare the distribution of the final proportion with the Beta distribution
beta_shape1 <- a / c
beta_shape2 <- b / c
beta_samples <- rbeta(n_simulations, beta_shape1, beta_shape2)

# Plot the histogram of final proportions with the Beta distribution overlay
df_sim <- data.frame(Proportion = final_props)
ggplot(df_sim, aes(x = Proportion)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  stat_function(fun = dbeta, args = list(shape1 = beta_shape1, shape2 = beta_shape2), color = "red", size = 1) +
  labs(title = "Final Proportion Distribution vs Beta Distribution",
       x = "Final Proportion of Color A",
       y = "Density") +
  theme_minimal()

# Display summary statistics in a table
summary_df <- data.frame(
  Statistic = c("Mean", "Variance"),
  Simulation = c(mean(final_props), var(final_props)),
  Theoretical = c(beta_shape1 / (beta_shape1 + beta_shape2),
                  (beta_shape1 * beta_shape2) / ((beta_shape1 + beta_shape2)^2 * (beta_shape1 + beta_shape2 + 1)))
)

# Print the summary table
print(summary_df)

# Additional plots for visualizing the convergence behavior
n_steps <- seq(10, n_draws, by = 10)
proportions_over_time <- sapply(n_steps, function(n) {
  draws <- polya_urn(a, b, c, n)
  sum(draws == 1) / n
})

df_convergence <- data.frame(
  Steps = n_steps,
  Proportion = proportions_over_time
)

ggplot(df_convergence, aes(x = Steps, y = Proportion)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = a / (a + b), linetype = "dashed", color = "red") +
  labs(title = "Convergence of Proportion of Color A",
       x = "Number of Draws",
       y = "Proportion of Color A") +
  theme_minimal()

# Plot a comparison of theoretical vs simulated Beta distribution
ggplot(data.frame(x = beta_samples), aes(x)) +
  geom_density(fill = "lightblue", alpha = 0.7, color = "black") +
  geom_density(data = df_sim, aes(x = Proportion), color = "red", size = 1) +
  labs(title = "Theoretical vs Simulated Distribution",
       x = "Proportion of Color A",
       y = "Density") +
  theme_minimal()

