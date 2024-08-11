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
n_draws <- 10000

# Run the simulation
set.seed(123)
draws <- polya_urn(a, b, c, n_draws)

# Calculate the expected long-term proportion (theta_A) for color A
theta_A <- (-b/c + sqrt((b/c)^2 + 4*a/c)) / 2

# Calculate the running proportion of color A and B over time
prop_A <- cumsum(draws == 1) / 1:n_draws
prop_B <- cumsum(draws == 2) / 1:n_draws

# Create a data frame for plotting
df <- data.frame(Draw = 1:n_draws, Prop_A = prop_A, Prop_B = prop_B)

# Plot the evolution of the running proportion of color A and B over time
plot_A <- ggplot(df, aes(x = Draw, y = Prop_A)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = theta_A, linetype = "dashed", color = "red") +
  labs(title = "Running Proportion of Color A over Time",
       x = "Number of Draws",
       y = "Proportion of Color A") +
  theme_minimal() +
  annotate("text", x = n_draws * 0.7, y = 0.1, label = "Blue Line: Running Proportion of Color A\nRed Dashed Line: Expected Proportion (Theta_A)", color = "black", size = 3, hjust = 0)

plot_B <- ggplot(df, aes(x = Draw, y = Prop_B)) +
  geom_line(color = "green") +
  geom_hline(yintercept = 1 - theta_A, linetype = "dashed", color = "purple") +
  labs(title = "Running Proportion of Color B over Time",
       x = "Number of Draws",
       y = "Proportion of Color B") +
  theme_minimal() +
  annotate("text", x = n_draws * 0.7, y = 0.1, label = "Green Line: Running Proportion of Color B\nPurple Dashed Line: Expected Proportion (Theta_B)", color = "black", size = 3, hjust = 0)

# Save the plots
ggsave("running_proportion_A_over_time.png", plot_A, width = 8, height = 6)
ggsave("running_proportion_B_over_time.png", plot_B, width = 8, height = 6)

# Display the plots
grid.arrange(plot_A, plot_B, ncol = 2)
