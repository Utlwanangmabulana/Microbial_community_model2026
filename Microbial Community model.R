# ============================================================
# Microbial community model: additive vs Monod growth functions
# Based on Jansma et al. (2026)
# ============================================================

library(deSolve)

# ---- Model function ----
model <- function(t, state, parameters) {
  bacterium   <- state[1]
  metabolites <- state[2:8]
  
  death_rate               <- parameters$death_rate
  carrying_capacity        <- parameters$carrying_capacity
  conversion_coefficient   <- parameters$conversion_coefficient
  interaction_coefficients <- parameters$interaction_coefficients
  half_saturation_constant <- parameters$half_saturation_constant
  growth_function_type     <- parameters$growth_function_type
  supply_rate              <- parameters$supply_rate
  metabolite_decay         <- parameters$metabolite_decay
  
  growth_additive <- conversion_coefficient *
    (t(pmax(interaction_coefficients, 0)) %*% metabolites)
  
  growth_monod <- conversion_coefficient *
    (t(pmax(interaction_coefficients, 0)) %*%
       (metabolites / (1 + metabolites / half_saturation_constant)))
  
  if (growth_function_type == "additive") {
    growth_rate <- growth_additive
  } else {
    growth_rate <- growth_monod
  }
  
  change_in_bacterium <- bacterium *
    (growth_rate - death_rate) *
    (carrying_capacity - bacterium)
  
  change_in_metabolites <- supply_rate -
    metabolite_decay * metabolites -
    metabolites * (interaction_coefficients %*% bacterium)
  
  return(list(c(change_in_bacterium, change_in_metabolites)))
}

# ---- Interaction coefficients ----
set.seed(42)
number_of_bacteria <- 1
number_of_produced <- 2
number_of_consumed <- 5

coefficients_produced <- matrix(runif(number_of_produced * number_of_bacteria, -1, 0),
                                nrow = number_of_produced, ncol = number_of_bacteria)
coefficients_consumed <- matrix(runif(number_of_consumed * number_of_bacteria, 0, 1),
                                nrow = number_of_consumed, ncol = number_of_bacteria)
interaction_coefficient_matrix <- rbind(coefficients_produced, coefficients_consumed)

# ---- Parameters ----
parameters <- list(
  death_rate               = 0.1,
  carrying_capacity        = 14,
  conversion_coefficient   = 0.1,
  half_saturation_constant = 0.5,
  interaction_coefficients = interaction_coefficient_matrix,
  growth_function_type     = "additive",
  supply_rate              = 0.1,
  metabolite_decay         = 0.1
)

# ---- Time sequence ----
time_sequence <- seq(0, 40, by = 0.1)

# ---- Single run (for the first plot) ----
set.seed(1)
initial_state <- c(
  bacterium            = runif(1, 0.05, 0.2),
  metabolites_produced = runif(2, exp(-3), exp(-2)),
  metabolites_consumed = runif(5, 0.5, 1)
)

parameters$growth_function_type <- "additive"
output_additive <- ode(y = initial_state, times = time_sequence, func = model, parms = parameters)

parameters$growth_function_type <- "monod"
output_monod <- ode(y = initial_state, times = time_sequence, func = model, parms = parameters)

# Plot single run comparison
par(mfrow = c(1, 2))
plot(output_additive[, "time"], output_additive[, "bacterium"], type = "l", lwd = 2, col = "blue",
     xlab = "Time", ylab = "Bacterial abundance", main = "Additive growth", ylim = c(0, 1.5))
plot(output_monod[, "time"], output_monod[, "bacterium"], type = "l", lwd = 2, col = "red",
     xlab = "Time", ylab = "Bacterial abundance", main = "Monod growth", ylim = c(0, 1.5))
par(mfrow = c(1, 1))

# ---- 50 Replicates & Boxplot ----
results_additive <- numeric(50)
results_monod    <- numeric(50)

for (i in 1:50) {
  set.seed(i + 100)
  initial_i <- c(
    bacterium            = runif(1, 0.05, 0.2),
    metabolites_produced = runif(2, exp(-3), exp(-2)),
    metabolites_consumed = runif(5, 0.5, 1)
  )
  
  parameters$growth_function_type <- "additive"
  out_add <- ode(y = initial_i, times = time_sequence, func = model, parms = parameters)
  results_additive[i] <- out_add[nrow(out_add), "bacterium"]
  
  parameters$growth_function_type <- "monod"
  out_mon <- ode(y = initial_i, times = time_sequence, func = model, parms = parameters)
  results_monod[i] <- out_mon[nrow(out_mon), "bacterium"]
}

# Print Summary
cat("--- ADDITIVE MODEL ---\n")
cat("Mean:", mean(results_additive), "| SD:", sd(results_additive), "\n\n")
cat("--- MONOD MODEL ---\n")
cat("Mean:", mean(results_monod), "| SD:", sd(results_monod), "\n\n")

# Draw Boxplot
plot_data <- data.frame(
  Growth_Function = rep(c("Additive", "Monod"), each = 50),
  Equilibrium_Abundance = c(results_additive, results_monod)
)

boxplot(Equilibrium_Abundance ~ Growth_Function, data = plot_data,
        col = c("lightblue", "lightcoral"),
        main = "Equilibrium Abundance after 50 Replicates",
        ylab = "Final Bacterial Abundance",
        xlab = "Growth Function")

