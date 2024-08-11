# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

# Function to simulate Pólya urn process
polya_urn <- function(a, b, c, n_draws) {
  urn <- c(a, b)
  draws <- numeric(n_draws)

  for (i in 1:n_draws) {
    prob <- urn / sum(urn)
    draw <- sample(c(1, 2), size = 1, prob = prob)
    draws[i] <- draw
    urn[draw] <- urn[draw] + c
  }

  return(draws)
}

# Simulation parameters
a <- 10
b <- 10
c <- 1
n_draws <- 1000

# Run the simulation
set.seed(123)
draws <- polya_urn(a, b, c, n_draws)

# Calculate the proportion of color A and B over time
prop_A <- cumsum(draws == 1) / 1:n_draws
prop_B <- cumsum(draws == 2) / 1:n_draws

################

beta_shape1 <- a / c
beta_shape2 <- b / c
######################
# Create a data frame for plotting
df <- data.frame(Draw = 1:n_draws, Prop_A = prop_A, Prop_B = prop_B)

# Plot the evolution of the proportion of color A and B over time
plot_A <- ggplot(df, aes(x = Draw, y = Prop_A)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = a / (a + b), linetype = "dashed", color = "red") +
  labs(title = "Proportion of Color A over Time",
       x = "Number of Draws",
       y = "Proportion of Color A") +
  theme_minimal() +
  annotate("text", x = n_draws * 0.7, y = 0.1, label = "Blue Line: Proportion of Color A\nRed Dashed Line: Theoretical Proportion", color = "black", size = 3, hjust = 0)

plot_B <- ggplot(df, aes(x = Draw, y = Prop_B)) +
  geom_line(color = "green") +
  geom_hline(yintercept = b / (a + b), linetype = "dashed", color = "purple") +
  labs(title = "Proportion of Color B over Time",
       x = "Number of Draws",
       y = "Proportion of Color B") +
  theme_minimal() +
  annotate("text", x = n_draws * 0.7, y = 0.1, label = "Green Line: Proportion of Color B\nPurple Dashed Line: Theoretical Proportion", color = "black", size = 3, hjust = 0)

# Save the plots
ggsave("proportion_A_over_time.png", plot_A, width = 8, height = 6)
ggsave("proportion_B_over_time.png", plot_B, width = 8, height = 6)

# Display the plots
grid.arrange(plot_A, plot_B, ncol = 2)

# Simulate multiple urns to study the distribution of the final proportion
n_simulations <- 1000
final_props_A <- numeric(n_simulations)
final_props_B <- numeric(n_simulations)

for (i in 1:n_simulations) {
  draws <- polya_urn(a, b, c, n_draws)
  final_props_A[i] <- sum(draws == 1) / n_draws
  final_props_B[i] <- sum(draws == 2) / n_draws
}

# Compare the distribution of the final proportion with the Beta distribution
beta_shape1 <- a / c
beta_shape2 <- b / c
beta_samples_A <- rbeta(n_simulations, beta_shape1, beta_shape2)
beta_samples_B <- rbeta(n_simulations, beta_shape2, beta_shape1)

# Plot the histogram of final proportions with the Beta distribution overlay for A
df_sim_A <- data.frame(Proportion = final_props_A)
plot_final_A <- ggplot(df_sim_A, aes(x = Proportion)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
  stat_function(fun = dbeta, args = list(shape1 = beta_shape1, shape2 = beta_shape2), color = "red", size = 1) +
  labs(title = "Final Proportion Distribution of Color A vs Beta Distribution",
       x = "Final Proportion of Color A",
       y = "Density") +
  theme_minimal() +
  annotate("text", x = 0.8, y = 2.5, label = "Light Blue: Simulated Proportion\nRed Line: Beta Distribution", color = "black", size = 3, hjust = 0)

# Plot the histogram of final proportions with the Beta distribution overlay for B
df_sim_B <- data.frame(Proportion = final_props_B)
plot_final_B <- ggplot(df_sim_B, aes(x = Proportion)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "lightgreen", color = "black", alpha = 0.7) +
  stat_function(fun = dbeta, args = list(shape1 = beta_shape2, shape2 = beta_shape1), color = "purple", size = 1) +
  labs(title = "Final Proportion Distribution of Color B vs Beta Distribution",
       x = "Final Proportion of Color B",
       y = "Density") +
  theme_minimal() +
  annotate("text", x = 0.8, y = 2.5, label = "Light Green: Simulated Proportion\nPurple Line: Beta Distribution", color = "black", size = 3, hjust = 0)

# Save the final proportion distribution plots
ggsave("final_distribution_A.png", plot_final_A, width = 8, height = 6)
ggsave("final_distribution_B.png", plot_final_B, width = 8, height = 6)

# Display the final distribution plots
grid.arrange(plot_final_A, plot_final_B, ncol = 2)

# Display summary statistics in a table
summary_df <- data.frame(
  Statistic = c("Mean", "Variance"),
  Simulation_A = c(mean(final_props_A), var(final_props_A)),
  Theoretical_A = c(beta_shape1 / (beta_shape1 + beta_shape2),
                    (beta_shape1 * beta_shape2) / ((beta_shape1 + beta_shape2)^2 * (beta_shape1 + beta_shape2 + 1))),
  Simulation_B = c(mean(final_props_B), var(final_props_B)),
  Theoretical_B = c(beta_shape2 / (beta_shape2 + beta_shape1),
                    (beta_shape2 * beta_shape1) / ((beta_shape2 + beta_shape1)^2 * (beta_shape2 + beta_shape1 + 1)))
)

# Print and save the summary table
print(summary_df)
write.csv(summary_df, "summary_statistics.csv", row.names = FALSE)

# Additional plots for visualizing the convergence behavior for both A and B
n_steps <- seq(10, n_draws, by = 10)
proportions_over_time_A <- sapply(n_steps, function(n) {
  draws <- polya_urn(a, b, c, n)
  sum(draws == 1) / n
})

proportions_over_time_B <- sapply(n_steps, function(n) {
  draws <- polya_urn(a, b, c, n)
  sum(draws == 2) / n
})

df_convergence_A <- data.frame(Steps = n_steps, Proportion = proportions_over_time_A)
df_convergence_B <- data.frame(Steps = n_steps, Proportion = proportions_over_time_B)

plot_convergence_A <- ggplot(df_convergence_A, aes(x = Steps, y = Proportion)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = a / (a + b), linetype = "dashed", color = "red") +
  labs(title = "Convergence of Proportion of Color A",
       x = "Number of Draws",
       y = "Proportion of Color A") +
  theme_minimal() +
  annotate("text", x = n_draws * 0.7, y = 0.1, label = "Blue Line: Proportion of Color A\nRed Dashed Line: Theoretical Proportion", color = "black", size = 3, hjust = 0)

plot_convergence_B <- ggplot(df_convergence_B, aes(x = Steps, y = Proportion)) +
  geom_line(color = "green") +
  geom_hline(yintercept = b / (a + b), linetype = "dashed", color = "purple") +
  labs(title = "Convergence of Proportion of Color B",
       x = "Number of Draws",
       y = "Proportion of Color B") +
  theme_minimal() +
  annotate("text", x = n_draws * 0.7, y = 0.1, label = "Green Line: Proportion of Color B\nPurple Dashed Line: Theoretical Proportion", color = "black", size = 3, hjust = 0)

# Save the convergence plots
ggsave("convergence_A.png", plot_convergence_A, width = 8, height = 6)
ggsave("convergence_B.png", plot_convergence_B, width = 8, height = 6)
xtable(summary_df)
# Display the convergence plots
grid.arrange(plot_convergence_A, plot_convergence_B, ncol = 2)
