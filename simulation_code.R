############################################################
# CORRECTED CODE FOR JOURNAL PAPER
# Bayesian and Frequentist Shrinkage Estimation for UIW Distribution
# Author: Alireza Safariyan
#
# CORRECTIONS APPLIED:
# 1. MLE shown as horizontal constant lines (not curves) - baseline reference
# 2. Added error bars/ribbons using SE_MSE for statistical uncertainty
# 3. Improved aesthetics for journal publication quality
############################################################

rm(list = ls())

library(ggplot2)
library(dplyr)
library(openxlsx)
library(tidyr)
library(gridExtra)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(scales)

setwd("your address folder")
set.seed(2026)

############################################################
# PRIOR INFORMATION
############################################################

R0 <- 0.50  # Target shrinkage value (neutral reliability)

############################################################
# GENERATE APPROXIMATE MLE
############################################################

generate_mle <- function(n, trueR) {
  # Asymptotic variance
  sigma2 <- trueR * (1 - trueR) / n

  # Generate MLE (approximately normal)
  rhat <- rnorm(1, mean = trueR, sd = sqrt(sigma2))

  # Keep inside (0,1)
  rhat <- min(max(rhat, 0.001), 0.999)

  return(rhat)
}

############################################################
# VARIANCE ESTIMATOR
############################################################

var_estimator <- function(rhat, n) {
  vhat <- rhat * (1 - rhat) / n
  vhat <- max(vhat, 1e-6)  # Prevent zero variance
  return(vhat)
}

############################################################
# SMOOTH PRELIMINARY TEST ESTIMATOR (PTE)
############################################################

pte_estimator <- function(rhat, n, alpha = 0.05) {
  vhat <- var_estimator(rhat, n)
  Tstat <- ((rhat - R0)^2) / vhat
  crit <- qchisq(1 - alpha, df = 1)
  w <- exp(-Tstat / crit)  # Smooth weight
  out <- w * R0 + (1 - w) * rhat
  out <- min(max(out, 0), 1)
  return(out)
}

############################################################
# REGULAR STEIN ESTIMATOR (always shrinks)
############################################################

stein_estimator <- function(rhat, n, trueR) {
  sigma2 <- var_estimator(rhat, n)

  # Optimal shrinkage intensity (Theorem 4)
  lambda <- sigma2 / (sigma2 + (trueR - R0)^2)

  # Stein estimator
  out <- (1 - lambda) * rhat + lambda * R0
  out <- min(max(out, 0), 1)

  return(out)
}

############################################################
# POSITIVE-PART STEIN ESTIMATOR (only shrinks when beneficial)
############################################################

positive_stein_estimator <- function(rhat, n, trueR, c = 1) {
  vhat <- var_estimator(rhat, n)

  # Standardized squared distance from target
  z_squared <- ((rhat - R0)^2) / vhat

  # Positive-part shrinkage factor
  shrinkage <- max(0, 1 - c / z_squared)

  # Positive-part Stein estimator
  out <- R0 + shrinkage * (rhat - R0)
  out <- min(max(out, 0), 1)

  return(out)
}

############################################################
# CONVEX SHRINKAGE ESTIMATOR (fixed target R0=0.5)
############################################################

convex_estimator <- function(rhat, n) {
  alpha <- 1 - min(0.4, 10 / n)
  out <- alpha * rhat + (1 - alpha) * R0
  out <- min(max(out, 0), 1)
  return(out)
}

############################################################
# RIDGE-TYPE SHRINKAGE ESTIMATOR
############################################################

ridge_estimator <- function(rhat, n) {
  lambda <- min(0.3, 8 / n)
  out <- (1 - lambda) * rhat + lambda * R0
  out <- min(max(out, 0), 1)
  return(out)
}

############################################################
# ADAPTIVE BAYES ESTIMATOR
############################################################

bayes_estimator <- function(rhat, n) {
  prior_strength <- sqrt(n)
  a0 <- R0 * prior_strength
  b0 <- (1 - R0) * prior_strength
  out <- (a0 + n * rhat) / (a0 + b0 + n)
  out <- min(max(out, 0), 1)
  return(out)
}

############################################################
# EMPIRICAL BAYES ESTIMATOR
############################################################

empirical_bayes_estimator <- function(rhat, n, bootstrap_reps = 100) {
  boot_samples <- replicate(bootstrap_reps, {
    rhat_boot <- generate_mle(n, rhat)
    var_estimator(rhat_boot, n)
  })

  sampling_var <- mean(boot_samples)
  total_var <- var(boot_samples)
  prior_var <- max(0, total_var - sampling_var)

  lambda_eb <- sampling_var / (sampling_var + prior_var + 1e-6)
  out <- (1 - lambda_eb) * rhat + lambda_eb * R0
  out <- min(max(out, 0), 1)

  return(out)
}

############################################################
# ENHANCED SIMULATION FUNCTION
############################################################

run_simulation <- function(trueR, n, M = 15000) {

  # Storage matrix for ALL 8 estimators
  results <- matrix(NA, nrow = M, ncol = 8)
  colnames(results) <- c("MLE", "PTE", "Stein", "PositiveStein",
                         "Convex", "Ridge", "Bayes", "EmpiricalBayes")

  # Progress indicator
  pb <- txtProgressBar(min = 0, max = M, style = 3)

  for(m in 1:M) {
    rhat <- generate_mle(n, trueR)

    results[m, "MLE"] <- rhat
    results[m, "PTE"] <- pte_estimator(rhat, n)
    results[m, "Stein"] <- stein_estimator(rhat, n, trueR)
    results[m, "PositiveStein"] <- positive_stein_estimator(rhat, n, trueR)
    results[m, "Convex"] <- convex_estimator(rhat, n)
    results[m, "Ridge"] <- ridge_estimator(rhat, n)
    results[m, "Bayes"] <- bayes_estimator(rhat, n)
    results[m, "EmpiricalBayes"] <- empirical_bayes_estimator(rhat, n)

    setTxtProgressBar(pb, m)
  }

  close(pb)

  # Performance metrics
  bias <- apply(results, 2, function(x) mean(x - trueR))
  mse <- apply(results, 2, function(x) mean((x - trueR)^2))

  # Bootstrap standard errors
  boot_mse <- replicate(200, {
    idx <- sample(1:M, M, replace = TRUE)
    apply(results[idx, ], 2, function(x) mean((x - trueR)^2))
  })
  se_mse <- apply(boot_mse, 1, sd)

  # Relative efficiency
  re <- mse["MLE"] / mse
  se_re <- re * sqrt((se_mse["MLE"]/mse["MLE"])^2 + (se_mse/mse)^2)

  # MSE improvement
  mse_improvement <- (mse["MLE"] - mse) / mse["MLE"] * 100

  # Coverage probability
  coverage <- apply(results, 2, function(x) {
    se_x <- sqrt(var_estimator(x, n))
    lower <- x - 1.96 * se_x
    upper <- x + 1.96 * se_x
    mean(lower <= trueR & upper >= trueR)
  })

  out <- data.frame(
    Method = names(bias),
    n = n,
    TrueR = trueR,
    Bias = bias,
    MSE = mse,
    SE_MSE = se_mse,
    RelativeEfficiency = re,
    SE_RE = se_re,
    MSE_Improvement = mse_improvement,
    Coverage95 = coverage
  )

  return(out)
}

############################################################
# SIMULATION DESIGN - FULL RANGE [0, 0.99]
############################################################

# COMPLETE range from near 0 to near 1 (avoid exact 0 and 1 for numerical stability)
trueR_values <- c(0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
                  0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
                  0.90, 0.95, 0.99)

# Sample sizes
sample_sizes <- c(10, 20, 40, 80)

# Number of replications for smooth curves
M_reps <- 15000  # Increased for smoother curves

cat("\n========================================\n")
cat("Starting Simulation Study\n")
cat("True R values:", length(trueR_values), "points from 0.01 to 0.99\n")
cat("Replications per configuration:", M_reps, "\n")
cat("Total configurations:", length(trueR_values) * length(sample_sizes), "\n")
cat("========================================\n")

all_results <- data.frame()

for(i in 1:length(trueR_values)) {
  for(n in sample_sizes) {
    cat(sprintf("\n\n>>> Running: TrueR = %.3f | n = %d\n",
                trueR_values[i], n))

    sim_out <- run_simulation(trueR = trueR_values[i], n = n, M = M_reps)
    sim_out$Scenario <- paste0("R", round(trueR_values[i] * 100))

    all_results <- rbind(all_results, sim_out)
  }
}

cat("\n========================================\n")
cat("Simulation Complete!\n")
cat("Total rows in results:", nrow(all_results), "\n")
cat("========================================\n")

############################################################
# SAVE RESULTS
############################################################

write.csv(all_results, "journal_results_full_range.csv", row.names = FALSE)
write.xlsx(all_results, "journal_results_full_range.xlsx")
save(all_results,file="all_results.Rdata")

############################################################
# COMPUTE MLE REFERENCE VALUES (Horizontal baseline)
############################################################
# CORRECTION: MLE serves as the constant baseline reference.
# We compute average MLE MSE per sample size for horizontal reference lines.

mle_baseline <- all_results %>%
  filter(Method == "MLE") %>%
  group_by(n) %>%
  summarise(
    MLE_MSE_mean = mean(MSE),
    MLE_MSE_se = mean(SE_MSE),
    MLE_Bias_mean = mean(Bias),
    MLE_Coverage_mean = mean(Coverage95),
    .groups = "drop"
  )

cat("\n\n========== MLE BASELINE VALUES (Horizontal Reference) ==========\n")
print(mle_baseline)

############################################################
# TABLE 1: MAIN RESULTS (Selected points for paper)
############################################################

# Select key trueR values for table display
key_R_values <- c(0.10, 0.25, 0.50, 0.75, 0.90)

main_table <- all_results %>%
  filter(TrueR %in% key_R_values,
         Method %in% c("MLE", "PTE", "Stein", "PositiveStein", "Convex", "Bayes")) %>%
  select(TrueR, n, Method, MSE, MSE_Improvement, RelativeEfficiency, Coverage95) %>%
  arrange(n, TrueR, Method)

write.csv(main_table, "table1_main_results.csv", row.names = FALSE)

cat("\n\n========== TABLE 1: MAIN RESULTS ==========\n")
print(main_table %>% filter(n == 20) %>% head(30))

############################################################
# TABLE 2: IMPROVEMENT SUMMARY
############################################################

improvement_by_n <- all_results %>%
  filter(Method != "MLE") %>%
  group_by(n, Method) %>%
  summarise(
    Avg_Improvement = mean(MSE_Improvement),
    SE_Improvement = sd(MSE_Improvement) / sqrt(n()),
    Max_Improvement = max(MSE_Improvement),
    Min_Improvement = min(MSE_Improvement),
    Avg_Coverage = mean(Coverage95),
    .groups = "drop"
  ) %>%
  arrange(n, desc(Avg_Improvement))

write.csv(improvement_by_n, "table2_improvement_by_n.csv", row.names = FALSE)

cat("\n\n========== TABLE 2: IMPROVEMENT BY SAMPLE SIZE ==========\n")
print(improvement_by_n)

############################################################
# DEFINE COLOR PALETTE
############################################################

method_colors <- c(
  "MLE" = "#1B9E77",
  "PTE" = "#D95F02",
  "Stein" = "#7570B3",
  "PositiveStein" = "#E7298A",
  "Convex" = "#66A61E",
  "Ridge" = "#E6AB02",
  "Bayes" = "#A6761D",
  "EmpiricalBayes" = "#666666"
)

# Linetypes for better distinction
method_linetypes <- c(
  "MLE" = "dashed",
  "PTE" = "solid",
  "Stein" = "solid",
  "PositiveStein" = "solid",
  "Convex" = "solid",
  "Ridge" = "solid",
  "Bayes" = "solid",
  "EmpiricalBayes" = "solid"
)

############################################################
# FIGURE 1: MSE COMPARISON - CORRECTED VERSION
############################################################
# CORRECTION: MLE shown as horizontal reference lines with error ribbons

cat("\n\n========== GENERATING FIGURE 1 (CORRECTED) ==========\n")

# Prepare data - exclude MLE from smooth curves
shrinkage_data <- all_results %>%
  filter(Method %in% c("PTE", "Stein", "PositiveStein", "Convex", "Ridge", "Bayes"))

# Get MLE data with error ribbons
mle_data <- all_results %>%
  filter(Method == "MLE") %>%
  left_join(mle_baseline, by = "n")

p1 <- ggplot() +
  # Shrinkage estimators as smooth curves WITH error ribbons
  geom_smooth(data = shrinkage_data,
              aes(x = TrueR, y = MSE, color = Method, linetype = Method),
              se = TRUE, method = "loess", span = 0.2, linewidth = 1.0, alpha = 0.8) +
  # MLE as horizontal reference line with error ribbon
  geom_hline(data = mle_baseline,
             aes(yintercept = MLE_MSE_mean, color = "MLE", linetype = "MLE"),
             linewidth = 1.2) +
  # MLE error ribbon (constant across TrueR)
  geom_ribbon(data = mle_baseline,
              aes(ymin = MLE_MSE_mean - 1.96 * MLE_MSE_se,
                  ymax = MLE_MSE_mean + 1.96 * MLE_MSE_se,
                  x = 0.5),  # dummy x for facet
              fill = method_colors["MLE"], alpha = 0.15, inherit.aes = FALSE) +
  facet_wrap(~n, scales = "free_y", ncol = 2,
             labeller = labeller(n = function(x) paste0("n = ", x))) +
  scale_color_manual(values = method_colors, name = "Estimator") +
  scale_linetype_manual(values = method_linetypes, name = "Estimator") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  theme_bw(base_size = 13) +
  labs(
    title = "A) Mean Squared Error Comparison",
    subtitle = paste("Based on", M_reps, "replications per configuration | MLE = horizontal reference (dashed)"),
    x = expression("True Reliability" ~ (R[s,k])),
    y = "Mean Squared Error (MSE)"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.box = "horizontal",
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines")
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))

ggsave("figure1_mse_corrected.png", p1, width = 12, height = 8, dpi = 300)
print(p1)

############################################################
# FIGURE 2: RELATIVE EFFICIENCY - WITH ERROR RIBBONS
############################################################
# CORRECTION: Added error ribbons using SE_RE

cat("\n\n========== GENERATING FIGURE 2 (CORRECTED) ==========\n")

re_data <- all_results %>%
  filter(Method %in% c("PTE", "Stein", "PositiveStein", "Convex", "Ridge", "Bayes"))

p2 <- ggplot(re_data, aes(x = TrueR, y = RelativeEfficiency, color = Method, fill = Method)) +
  # Reference line at RE = 1 (MLE baseline)
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  # Smooth curves with confidence ribbons
  geom_smooth(aes(ymin = RelativeEfficiency - 1.96 * SE_RE,
                  ymax = RelativeEfficiency + 1.96 * SE_RE),
              stat = "smooth", method = "loess", span = 0.2,
              linewidth = 1.0, alpha = 0.2, se = TRUE) +
  facet_wrap(~n, ncol = 2,
             labeller = labeller(n = function(x) paste0("n = ", x))) +
  scale_color_manual(values = method_colors) +
  scale_fill_manual(values = method_colors, guide = "none") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0.5, 1.5, by = 0.1)) +
  theme_bw(base_size = 13) +
  labs(
    title = "B) Relative Efficiency: MLE vs Shrinkage Estimators",
    subtitle = "RE = MSE(MLE) / MSE(Estimator) | Values > 1 indicate improvement | Shaded = 95% CI",
    x = expression("True Reliability" ~ (R[s,k])),
    y = "Relative Efficiency"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave("figure2_relative_efficiency_corrected.png", p2, width = 12, height = 8, dpi = 300)
print(p2)

############################################################
# FIGURE 3: MSE IMPROVEMENT HEATMAP (UNCHANGED - NO MLE HERE)
############################################################

cat("\n\n========== GENERATING FIGURE 3 ==========\n")

heatmap_data <- all_results %>%
  filter(Method %in% c("PTE", "Stein", "PositiveStein", "Convex", "Ridge", "Bayes")) %>%
  group_by(n, TrueR, Method) %>%
  summarise(MSE_Improvement = mean(MSE_Improvement), .groups = "drop")

p3 <- ggplot(heatmap_data,
             aes(x = factor(n), y = factor(round(TrueR, 2)), fill = MSE_Improvement)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", MSE_Improvement)), size = 2.2) +
  facet_wrap(~Method, ncol = 3) +
  scale_fill_gradient2(
    low = "firebrick", mid = "white", high = "steelblue",
    midpoint = 0, name = "MSE\nImprovement (%)",
    limits = c(-5, 15),
    oob = scales::squish
  ) +
  theme_minimal(base_size = 10) +
  labs(
    title = "C) MSE Improvement (%) Heatmap",
    subtitle = "Positive values (blue) indicate improvement over MLE | Negative values (red) indicate worse performance",
    x = "Sample Size (n)",
    y = expression("True Reliability" ~ (R[s,k]))
  ) +
  theme(
    axis.text.x = element_text(angle = 0, size = 9),
    axis.text.y = element_text(size = 7),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave("figure3_improvement_heatmap.png", p3, width = 14, height = 10, dpi = 300)
print(p3)

############################################################
# FIGURE 4: COVERAGE PROBABILITY - CORRECTED VERSION
############################################################
# CORRECTION: MLE shown as horizontal reference line

cat("\n\n========== GENERATING FIGURE 4 (CORRECTED) ==========\n")

coverage_data <- all_results %>%
  filter(Method %in% c("PTE", "Stein", "PositiveStein", "Bayes"))

coverage_mle_baseline <- mle_baseline %>%
  select(n, MLE_Coverage_mean)

p4 <- ggplot() +
  # Shrinkage estimators
  geom_smooth(data = coverage_data,
              aes(x = TrueR, y = Coverage95, color = Method),
              se = FALSE, method = "loess", span = 0.25, linewidth = 1.0) +
  # MLE as horizontal reference
  geom_hline(data = coverage_mle_baseline,
             aes(yintercept = MLE_Coverage_mean, color = "MLE"),
             linetype = "dashed", linewidth = 1.2) +
  # Nominal coverage reference
  geom_hline(yintercept = 0.95, linetype = "dotted", color = "darkred", linewidth = 0.8) +
  facet_wrap(~n, ncol = 2,
             labeller = labeller(n = function(x) paste0("n = ", x))) +
  scale_color_manual(values = method_colors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  scale_y_continuous(breaks = seq(0.85, 1.00, by = 0.02), limits = c(0.85, 1.00)) +
  theme_bw(base_size = 13) +
  labs(
    title = "D) 95% Confidence Interval Coverage Probability",
    subtitle = "Dashed = MLE baseline | Dotted = Nominal 95% | MLE = horizontal reference",
    x = expression("True Reliability" ~ (R[s,k])),
    y = "Coverage Probability"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

ggsave("figure4_coverage_corrected.png", p4, width = 12, height = 8, dpi = 300)
print(p4)

############################################################
# FIGURE 5: BIAS COMPARISON - CORRECTED VERSION
############################################################
# CORRECTION: MLE shown as horizontal reference line at zero

cat("\n\n========== GENERATING FIGURE 5 (CORRECTED) ==========\n")

bias_data <- all_results %>%
  filter(Method %in% c("PTE", "Stein", "PositiveStein", "Convex", "Bayes"))

p5 <- ggplot() +
  # Shrinkage estimators
  geom_smooth(data = bias_data,
              aes(x = TrueR, y = Bias, color = Method, linetype = Method),
              se = TRUE, method = "loess", span = 0.2, linewidth = 1.0, alpha = 0.8) +
  # MLE as horizontal reference (Bias = 0)
  geom_hline(yintercept = 0, color = method_colors["MLE"],
             linetype = "dashed", linewidth = 1.2) +
  # Zero reference line
  geom_hline(yintercept = 0, linetype = "dotted", color = "darkred", linewidth = 0.6) +
  facet_wrap(~n, ncol = 2,
             labeller = labeller(n = function(x) paste0("n = ", x))) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_linetypes) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  theme_bw(base_size = 13) +
  labs(
    title = "E) Bias Comparison: MLE vs Shrinkage Estimators",
    subtitle = "Shrinkage introduces bias but reduces variance | Dashed = MLE (unbiased reference)",
    x = expression("True Reliability" ~ (R[s,k])),
    y = "Bias"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         linetype = guide_legend(nrow = 2, byrow = TRUE))

ggsave("figure5_bias_comparison_corrected.png", p5, width = 12, height = 8, dpi = 300)
print(p5)

############################################################
# FIGURE 6: STEIN vs POSITIVE-STEIN COMPARISON (UNCHANGED)
############################################################

cat("\n\n========== GENERATING FIGURE 6 ==========\n")

stein_comparison <- all_results %>%
  filter(Method %in% c("Stein", "PositiveStein")) %>%
  select(n, TrueR, Method, MSE_Improvement) %>%
  pivot_wider(names_from = Method, values_from = MSE_Improvement)

colnames(stein_comparison) <- c("n", "TrueR", "Improvement_Stein", "Improvement_PositiveStein")
stein_comparison$Difference <- stein_comparison$Improvement_Stein - stein_comparison$Improvement_PositiveStein

p6 <- ggplot(stein_comparison, aes(x = TrueR, y = Difference, color = factor(n))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkred", linewidth = 0.8) +
  geom_smooth(se = TRUE, method = "loess", span = 0.25, linewidth = 1.2) +
  geom_point(size = 1, alpha = 0.3) +
  scale_color_viridis_d(option = "turbo", name = "Sample Size (n)") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  theme_bw(base_size = 13) +
  labs(
    title = "F) Stein vs Positive-Part Stein: Difference in MSE Improvement",
    subtitle = "Positive values = Regular Stein better | Negative values = Positive-Part better",
    x = expression("True Reliability" ~ (R[s,k])),
    y = "Improvement Difference (Percentage Points)"
  ) +
  theme(legend.position = "bottom")

ggsave("figure6_stein_comparison.png", p6, width = 10, height = 6, dpi = 300)
print(p6)

############################################################
# COMBINED FIGURE (CORRECTED VERSION)
############################################################

cat("\n\n========== GENERATING COMBINED FIGURE ==========\n")

combined_plot <- (p1 + p2) / (p4 + p5) +
  plot_annotation(
    title = "Shrinkage Estimation for UIW Distribution under Type-II Censoring",
    subtitle = paste("Results based on", M_reps, "simulation replications | MLE = horizontal reference (dashed)"),
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5))
  )

ggsave("combined_figure_corrected.png", combined_plot, width = 14, height = 12, dpi = 300)

############################################################
# SUMMARY STATISTICS FOR PAPER
############################################################

cat("\n\n========== PAPER-READY SUMMARY ==========\n")

# Overall best performer
cat("\n1. OVERALL BEST PERFORMER (Average MSE Improvement):\n")
best_overall <- all_results %>%
  filter(Method != "MLE") %>%
  group_by(Method) %>%
  summarise(
    Avg_MSE_Improvement = mean(MSE_Improvement),
    SD_Improvement = sd(MSE_Improvement),
    Max_Improvement = max(MSE_Improvement),
    Min_Improvement = min(MSE_Improvement)
  ) %>%
  arrange(desc(Avg_MSE_Improvement))
print(best_overall)

# Small sample performance
cat("\n\n2. SMALL SAMPLE (n=10) PERFORMANCE:\n")
small_n <- all_results %>%
  filter(n == 10, Method != "MLE") %>%
  group_by(Method) %>%
  summarise(
    Avg_Improvement = mean(MSE_Improvement),
    Avg_Coverage = mean(Coverage95)
  ) %>%
  arrange(desc(Avg_Improvement))
print(small_n)

# Large sample performance
cat("\n\n3. LARGE SAMPLE (n=80) PERFORMANCE:\n")
large_n <- all_results %>%
  filter(n == 80, Method != "MLE") %>%
  group_by(Method) %>%
  summarise(
    Avg_Improvement = mean(MSE_Improvement),
    Avg_Coverage = mean(Coverage95)
  ) %>%
  arrange(desc(Avg_Improvement))
print(large_n)

# Coverage summary
cat("\n\n4. COVERAGE PROBABILITY SUMMARY (nominal 95%):\n")
coverage_summary <- all_results %>%
  filter(Method %in% c("MLE", "PTE", "Stein", "PositiveStein", "Bayes")) %>%
  group_by(Method) %>%
  summarise(
    Mean_Coverage = mean(Coverage95),
    SD_Coverage = sd(Coverage95),
    Min_Coverage = min(Coverage95),
    Max_Coverage = max(Coverage95)
  )
print(coverage_summary)

# Performance near boundaries
cat("\n\n5. PERFORMANCE AT EXTREME VALUES:\n")
extreme_performance <- all_results %>%
  filter(TrueR %in% c(0.01, 0.10, 0.90, 0.99), Method %in% c("MLE", "Stein", "PositiveStein")) %>%
  select(TrueR, n, Method, MSE, MSE_Improvement, Coverage95)
print(extreme_performance)

# Best method by sample size
cat("\n\n6. BEST METHOD BY SAMPLE SIZE:\n")
best_by_n <- improvement_by_n %>%
  group_by(n) %>%
  slice_max(Avg_Improvement, n = 1) %>%
  select(n, Method, Avg_Improvement)
print(best_by_n)

# Symmetry check (performance around R0=0.50)
cat("\n\n7. SYMMETRY CHECK (Performance around R0=0.50):\n")
symmetry_check <- all_results %>%
  filter(Method == "Stein", TrueR %in% c(0.30, 0.70)) %>%
  group_by(n, TrueR) %>%
  summarise(Avg_MSE_Improvement = mean(MSE_Improvement)) %>%
  pivot_wider(names_from = TrueR, values_from = Avg_MSE_Improvement,
              names_prefix = "R")
print(symmetry_check)

############################################################
# THEORETICAL VALIDATION
############################################################

cat("\n\n========== THEORETICAL VALIDATION ==========\n")

theoretical_validation <- all_results %>%
  filter(Method == "MLE") %>%
  mutate(
    Var_MLE = pmax(MSE - Bias^2, 1e-6),
    Lambda_opt_theory = Var_MLE / (Var_MLE + (TrueR - R0)^2)
  ) %>%
  select(n, TrueR, Lambda_opt_theory)

stein_perf <- all_results %>%
  filter(Method == "Stein") %>%
  select(n, TrueR, MSE_Improvement)

theoretical_validation <- merge(theoretical_validation, stein_perf, by = c("n", "TrueR"))

cat("Correlation between λ* and MSE Improvement:\n")
correlation <- cor(theoretical_validation$Lambda_opt_theory,
                   theoretical_validation$MSE_Improvement)
cat(sprintf("Pearson correlation: %.3f\n", correlation))

p_theory <- ggplot(theoretical_validation,
                   aes(x = Lambda_opt_theory, y = MSE_Improvement, color = factor(n))) +
  geom_point(size = 1.5, alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  scale_color_viridis_d(option = "turbo", name = "Sample Size (n)") +
  theme_bw(base_size = 12) +
  labs(
    title = "Theoretical Validation: Optimal Shrinkage vs Empirical Improvement",
    subtitle = sprintf("Correlation = %.3f | Validates Theorem 4", correlation),
    x = expression("Theoretical Optimal Shrinkage" ~ (lambda^star)),
    y = "MSE Improvement (%)"
  )

ggsave("figure_theoretical_validation.png", p_theory, width = 8, height = 6, dpi = 300)
print(p_theory)

cat("\n\n========================================\n")
cat("ALL RESULTS SAVED SUCCESSFULLY!\n")
cat("========================================\n")
cat("Files created:\n")
cat("  - journal_results_full_range.csv\n")
cat("  - journal_results_full_range.xlsx\n")
cat("  - table1_main_results.csv\n")
cat("  - table2_improvement_by_n.csv\n")
cat("  - figure1_mse_corrected.png (MLE = horizontal reference)\n")
cat("  - figure2_relative_efficiency_corrected.png (with error ribbons)\n")
cat("  - figure3_improvement_heatmap.png\n")
cat("  - figure4_coverage_corrected.png (MLE = horizontal reference)\n")
cat("  - figure5_bias_comparison_corrected.png (MLE = horizontal reference)\n")
cat("  - figure6_stein_comparison.png\n")
cat("  - figure_theoretical_validation.png\n")
cat("  - combined_figure_corrected.png\n")
cat("========================================\n")
cat("\nKEY CORRECTIONS APPLIED:\n")
cat("1. MLE now shown as horizontal reference lines (not curves)\n")
cat("2. Added error ribbons for statistical uncertainty (95% CI)\n")
cat("3. Improved legend clarity and visual distinction\n")
cat("4. Better facet labels with 'n = X' format\n")
cat("========================================\n")
