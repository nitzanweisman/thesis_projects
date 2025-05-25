# --- Multi-Model Averaging for Sorghum Circadian Plasticity ---

# Author: Nitzan Weisman
# Date: 25/5/25
# Description: Uses multi-model averaging (MMA) based on AICc to identify
#              environmental predictors of circadian trait plasticity (?? Amplitude,
#              ?? RAE, ?? RF) in sorghum ecotypes responding to temperature change.
#              This script assumes the input data for modeling is pre-cleaned and
#              contains NO MISSING VALUES (NAs) in the relevant columns.
# Data Source: Adapted from Master's Thesis, The Hebrew University of Jerusalem (March 2024)
# Input File: C:/Users/nitza/Documents/thesis_project_for_github/Sorghum_initial_screening_2.csv

# --- 1. Load Required Packages ---
# Load packages for data handling, statistics, and plotting.
library(tidyverse) # For data manipulation (dplyr, readr) and plotting (ggplot2)
library(car)       # For VIF (Variance Inflation Factor) calculation
library(MuMIn)     # For Multi-Model Inference (dredge, model.avg)
library(ggpubr)    # For 'ggplot2' based publication-ready plots (ggscatter, ggarrange)
library(gridExtra) # For arranging multiple grid-based plots
library(showtext)  # For custom font support in plots (optional)

# --- Optional: Custom Font Setup ---
# font_add(family = "Calibri", regular = "Calibri.ttf") # Ensure font file is accessible
# showtext_auto() # Enable custom fonts for plotting

# --- 2. Load and Prepare Data ---
# Load the primary dataset from a CSV file using an absolute path.
sorghum_data_raw <- readr::read_csv("C:/Users/nitza/Documents/thesis_project_for_github/Sorghum_initial_screening_2.csv")

# Adjust plasticity traits (delta values).
# The original CSV _delta columns represent (Trait at 35°C - Trait at 22°C).
# For this analysis, delta is defined as (Trait at 22°C - Trait at 35°C).
# This requires flipping the sign of the existing delta columns.
sorghum_plasticity_adjusted <- sorghum_data_raw %>%
  dplyr::mutate( # Modifies existing columns or creates new ones.
    Amplitude_delta = -Amplitude_delta,
    Period_delta = -Period_delta,
    RAE_delta = -RAE_delta,
    Rfd_avg_delta = -Rfd_avg_delta,
    Rhythmic.fraction_delta = -Rhythmic.fraction_delta
  )

# Select relevant columns for modeling:
# - Response variables (Y): The adjusted delta values for circadian traits.
# - Predictor variables (X): The 13 key environmental factors.
# - Also keeping 'Origin' and 'Climate.col' for plotting/identification later.
model_data <- sorghum_plasticity_adjusted %>%
  dplyr::select(
    # Response Variables (Plasticity)
    Amplitude_delta, RAE_delta, Rhythmic.fraction_delta, Period_delta, Rfd_avg_delta,
    # Predictor Variables (Environmental Factors)
    Ann.Mean.Tmp, t_season_delta, pet_annual, pet_season_delta,
    vpd_season_max, vpd_season_delta, Ann.Prc, prec_season_delta,
    PAR_season_max, Depth.to.Water.Table, soil_m, topsoil_pH,
    gs.length,
    # Columns for identification/grouping in plots
    Origin, Climate.col
  )

# Initial NA Check for modeling columns.
# This ensures the data is clean before proceeding with modeling.
# The dredge() function (used later) requires na.action = na.fail,
# which will stop the script if NAs are present in variables used by the global model.
modeling_columns <- c("Amplitude_delta", "RAE_delta", "Rhythmic.fraction_delta", "Period_delta", "Rfd_avg_delta",
                      "Ann.Mean.Tmp", "t_season_delta", "pet_annual", "pet_season_delta",
                      "vpd_season_max", "vpd_season_delta", "Ann.Prc", "prec_season_delta",
                      "PAR_season_max", "Depth.to.Water.Table", "soil_m", "topsoil_pH", "gs.length")
if (any(is.na(model_data[, modeling_columns]))) {
  warning("NAs detected in modeling columns! Clean data or implement NA handling for dredge to work.")
} else {
  cat("No NAs found in the selected modeling columns. Proceeding with analysis.\n")
}
# Set R's global NA action. 'na.fail' is crucial for MuMIn::dredge,
# ensuring all sub-models are compared on the same dataset.
options(na.action = "na.fail")


# --- 3. Modeling ?? Amplitude (Change in Amplitude: 22°C - 35°C) ---
cat("\n--- Modeling Delta Amplitude ---\n")
# Define the formula for the initial ("zero" or "global") linear model.
# This model includes all selected environmental parameters as potential predictors for ?? Amplitude.
formula_amp_zero <- Amplitude_delta ~ Ann.Mean.Tmp + t_season_delta + pet_annual +
  pet_season_delta + vpd_season_max + vpd_season_delta + Ann.Prc +
  prec_season_delta + PAR_season_max + Depth.to.Water.Table +
  soil_m + topsoil_pH + gs.length
# Build the initial linear model using lm().
zero_model_amp <- lm(formula_amp_zero, data = model_data)

# Iteratively reduce multicollinearity among predictors using Variance Inflation Factor (VIF).
# Multicollinearity can make model coefficients unstable and hard to interpret.
# A VIF > 4-5 (common threshold) suggests a predictor is highly correlated with others.
# The process removes the predictor with the highest VIF above the threshold, then recalculates VIFs.
cat("\nPerforming VIF reduction for Delta Amplitude model...\n")
# print(vif(zero_model_amp)) # Uncomment to see VIFs at each step
zero_model_amp <- update(zero_model_amp, ~ . - vpd_season_delta)
zero_model_amp <- update(zero_model_amp, ~ . - pet_annual)
zero_model_amp <- update(zero_model_amp, ~ . - vpd_season_max)
zero_model_amp <- update(zero_model_amp, ~ . - Depth.to.Water.Table)
zero_model_amp <- update(zero_model_amp, ~ . - prec_season_delta)
zero_model_amp <- update(zero_model_amp, ~ . - pet_season_delta)
cat("Final VIF after reduction (Delta Amplitude):\n"); print(vif(zero_model_amp)) # Show final VIFs.

# Perform multi-model selection using MuMIn::dredge.
# 'dredge' creates all possible simpler models from the VIF-reduced global model.
# Models are ranked by AICc (Akaike Information Criterion, corrected), balancing fit and complexity.
dredge_amp <- MuMIn::dredge(zero_model_amp)
# Select models with substantial empirical support (??AICc < 2 from the best model).
supported_amp <- subset(dredge_amp, delta < 2)
cat("\nSupported models (Delta Amplitude, ??AICc < 2):\n"); print(supported_amp)

# Perform multi-model averaging on the set of supported models.
# This averages parameter estimates (coefficients, SEs) weighted by each model's Akaike weight.
# It provides more robust estimates than relying on a single "best" model.
avg_model_amp <- MuMIn::model.avg(supported_amp, fit = TRUE)
cat("\n--- Averaged Model Summary (Delta Amplitude) ---\n")
# The summary shows averaged coefficients, their significance, and variable "importance"
# (sum of Akaike weights for models including that predictor).
print(summary(avg_model_amp))

# Assess the predictive performance of the final averaged model.
# This is done by correlating the model's predictions with the actual observed values.
pred_amp <- predict(avg_model_amp, newdata = model_data, type = "response", backtransform = FALSE)
obs_amp <- model_data$Amplitude_delta
if(length(na.omit(pred_amp)) != length(na.omit(obs_amp))) { # Sanity check for length
  warning("Prediction/observation length mismatch for Amplitude delta.")
}
# Fit a simple linear model of Predictions ~ Observations.
# A high R-squared from this indicates good predictive power of the averaged model.
fit_assessment_amp <- summary(lm(pred_amp ~ obs_amp, data=data.frame(pred_amp, obs_amp)))
cat("\n--- Fit Assessment R-squared (Delta Amplitude):", round(fit_assessment_amp$adj.r.squared, 4), "---\n")
cat("--- Fit Assessment p-value (Delta Amplitude):", format.pval(coefficients(fit_assessment_amp)[2,4]), "---\n")


# --- 4. Modeling ?? RAE (Change in Relative Amplitude Error: 22°C - 35°C) ---
# The same modeling process is repeated for the ?? RAE response variable.
cat("\n--- Modeling Delta RAE ---\n")
formula_rae_zero <- RAE_delta ~ Ann.Mean.Tmp + t_season_delta + pet_annual +
  pet_season_delta + vpd_season_max + vpd_season_delta + Ann.Prc +
  prec_season_delta + PAR_season_max + Depth.to.Water.Table +
  soil_m + topsoil_pH + gs.length
zero_model_rae <- lm(formula_rae_zero, data = model_data)
cat("\nPerforming VIF reduction for Delta RAE model...\n")
zero_model_rae <- update(zero_model_rae, ~ . - vpd_season_delta)
zero_model_rae <- update(zero_model_rae, ~ . - pet_annual)
zero_model_rae <- update(zero_model_rae, ~ . - vpd_season_max)
zero_model_rae <- update(zero_model_rae, ~ . - pet_season_delta)
zero_model_rae <- update(zero_model_rae, ~ . - Depth.to.Water.Table)
zero_model_rae <- update(zero_model_rae, ~ . - prec_season_delta)
cat("Final VIF after reduction (Delta RAE):\n"); print(vif(zero_model_rae))
dredge_rae <- MuMIn::dredge(zero_model_rae)
supported_rae <- subset(dredge_rae, delta < 2)
avg_model_rae <- MuMIn::model.avg(supported_rae, fit = TRUE)
cat("\n--- Averaged Model Summary (Delta RAE) ---\n")
print(summary(avg_model_rae))
pred_rae <- predict(avg_model_rae, newdata = model_data, type = "response", backtransform = FALSE)
obs_rae <- model_data$RAE_delta
if(length(na.omit(pred_rae)) != length(na.omit(obs_rae))) {
  warning("Prediction/observation length mismatch for RAE delta.")
}
fit_assessment_rae <- summary(lm(pred_rae ~ obs_rae, data=data.frame(pred_rae, obs_rae)))
cat("\n--- Fit Assessment R-squared (Delta RAE):", round(fit_assessment_rae$adj.r.squared, 4), "---\n")
cat("--- Fit Assessment p-value (Delta RAE):", format.pval(coefficients(fit_assessment_rae)[2,4]), "---\n")


# --- 5. Modeling ?? Rhythmic Fraction (Change in RF: 22°C - 35°C) ---
# The same modeling process is applied to the ?? Rhythmic Fraction response variable.
cat("\n--- Modeling Delta Rhythmic Fraction ---\n")
formula_rf_zero <- Rhythmic.fraction_delta ~ Ann.Mean.Tmp + t_season_delta + pet_annual +
  pet_season_delta + vpd_season_max + vpd_season_delta + Ann.Prc +
  prec_season_delta + PAR_season_max + Depth.to.Water.Table +
  soil_m + topsoil_pH + gs.length
zero_model_rf <- lm(formula_rf_zero, data = model_data)
cat("\nPerforming VIF reduction for Delta RF model...\n")
zero_model_rf <- update(zero_model_rf, ~ . - vpd_season_delta)
zero_model_rf <- update(zero_model_rf, ~ . - pet_annual)
zero_model_rf <- update(zero_model_rf, ~ . - vpd_season_max)
zero_model_rf <- update(zero_model_rf, ~ . - pet_season_delta)
zero_model_rf <- update(zero_model_rf, ~ . - Depth.to.Water.Table)
zero_model_rf <- update(zero_model_rf, ~ . - prec_season_delta)
cat("Final VIF after reduction (Delta RF):\n"); print(vif(zero_model_rf))
dredge_rf <- MuMIn::dredge(zero_model_rf)
supported_rf <- subset(dredge_rf, delta < 2)
avg_model_rf <- MuMIn::model.avg(supported_rf, fit = TRUE)
cat("\n--- Averaged Model Summary (Delta RF) ---\n")
print(summary(avg_model_rf))
pred_rf <- predict(avg_model_rf, newdata = model_data, type = "response", backtransform = FALSE)
obs_rf <- model_data$Rhythmic.fraction_delta
if(length(na.omit(pred_rf)) != length(na.omit(obs_rf))) {
  warning("Prediction/observation length mismatch for RF delta.")
}
fit_assessment_rf <- summary(lm(pred_rf ~ obs_rf, data=data.frame(pred_rf, obs_rf)))
cat("\n--- Fit Assessment R-squared (Delta RF):", round(fit_assessment_rf$adj.r.squared, 4), "---\n")
cat("--- Fit Assessment p-value (Delta RF):", format.pval(coefficients(fit_assessment_rf)[2,4]), "---\n")

# --- 6. Visualization ---
# This section visualizes key results:
# 1. Model Fit: Scatter plot of Observed vs. Predicted values for each delta trait.
# 2. Predictor Effects: Scatter plots of significant environmental predictors vs. each delta trait.

# --- 6.1 Visualization for Delta Amplitude ---
cat("\nGenerating plots for Delta Amplitude...\n")
# Prepare data for plotting observed vs. predicted values.
matched_indices_amp <- which(!is.na(model_data$Amplitude_delta) & !is.na(pred_amp))
a_d_prediction_plot_data <- data.frame(
  Observed = model_data$Amplitude_delta[matched_indices_amp],
  Predicted = pred_amp[matched_indices_amp],
  Climate.col = model_data$Climate.col[matched_indices_amp]
)
# Scatter plot of Observed vs. Predicted ?? Amplitude.
plot_obs_pred_amp <- ggplot(data=a_d_prediction_plot_data, aes(x=Observed, y=Predicted)) +
  geom_point(aes(color=Climate.col), size=3, show.legend = FALSE) + # Points colored by climate
  scale_color_identity() + # Use color values directly from 'Climate.col'
  stat_smooth(method=lm ,formula = y  ~ x, color="black", se=FALSE) + # Add linear regression line
  labs(x = "Observed ?? Amplitude [a.u] (22°C - 35°C)", y = "Predicted ?? Amplitude [a.u]",
       title= "Model Fit: Delta Amplitude Plasticity") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     text = element_text(size = 11.5, family = "Calibri"),
                     axis.text = element_text(size = 10, colour = 'black'))
print(plot_obs_pred_amp)

# Plot relationships between significant environmental predictors (identified from MMA) and ?? Amplitude.
# These plots show the individual effect of each key predictor.
plot_amp_prec <- ggpubr::ggscatter(model_data, x = "Ann.Prc", y = "Amplitude_delta",
                                   xlab = "Ann Precipitation\n[mm/Y]", ylab = "?? Amplitude [a.u] (22°C - 35°C)") +
  stat_smooth(method = lm, formula = y ~ x) + theme_bw() + theme(panel.grid = element_blank())
plot_amp_gsl <- ggpubr::ggscatter(model_data, x = "gs.length", y = "Amplitude_delta",
                                  xlab = "GSL\n[months]", ylab = "") + # Blank ylab for cleaner combined plot
  stat_smooth(method = lm, formula = y ~ x) + theme_bw() + theme(panel.grid = element_blank())
plot_amp_par <- ggpubr::ggscatter(model_data, x = "PAR_season_max", y = "Amplitude_delta",
                                  xlab = "Max PAR\n[W/m^2]", ylab = "") + # Blank ylab
  stat_smooth(method = lm, formula = y ~ x) + theme_bw() + theme(panel.grid = element_blank())
# Arrange the three predictor plots into a single figure.
combined_amp_predictors <- ggpubr::ggarrange(plot_amp_prec, plot_amp_gsl, plot_amp_par,
                                             labels = c("A", "B", "C"), ncol = 3, nrow = 1)
print(combined_amp_predictors)

# --- 6.2 Visualization for Delta RAE ---
cat("\nGenerating plots for Delta RAE...\n")
matched_indices_rae <- which(!is.na(model_data$RAE_delta) & !is.na(pred_rae))
rae_prediction_plot_data <- data.frame(
  Observed = model_data$RAE_delta[matched_indices_rae],
  Predicted = pred_rae[matched_indices_rae],
  Climate.col = model_data$Climate.col[matched_indices_rae]
)
plot_obs_pred_rae <- ggplot(data=rae_prediction_plot_data, aes(x=Observed, y=Predicted)) +
  geom_point(aes(color=Climate.col), size=3, show.legend = FALSE) +
  scale_color_identity() +
  stat_smooth(method=lm ,formula = y  ~ x, color="black", se=FALSE) +
  labs(x = "Observed ?? RAE [a.u] (22°C - 35°C)", y = "Predicted ?? RAE [a.u]",
       title= "Model Fit: Delta RAE Plasticity") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     text = element_text(size = 11.5, family = "Calibri"),
                     axis.text = element_text(size = 10, colour = 'black'))
print(plot_obs_pred_rae)

# Plot relationships between significant environmental predictors and ?? RAE.
plot_rae_par <- ggpubr::ggscatter(model_data, x = "PAR_season_max", y = "RAE_delta",
                                  xlab = "Max PAR\n[W/m^2]", ylab = "?? RAE [a.u] (22°C - 35°C)") +
  stat_smooth(method = lm, formula = y ~ x) + theme_bw() + theme(panel.grid = element_blank())
plot_rae_prec <- ggpubr::ggscatter(model_data, x = "Ann.Prc", y = "RAE_delta",
                                   xlab = "Ann Precipitation\n[mm/Y]", ylab = "") +
  stat_smooth(method = lm, formula = y ~ x) + theme_bw() + theme(panel.grid = element_blank())
plot_rae_gsl <- ggpubr::ggscatter(model_data, x = "gs.length", y = "RAE_delta",
                                  xlab = "GSL\n[months]", ylab = "") +
  stat_smooth(method = lm, formula = y ~ x) + theme_bw() + theme(panel.grid = element_blank())
combined_rae_predictors <- ggpubr::ggarrange(plot_rae_par, plot_rae_prec, plot_rae_gsl,
                                             labels = c("D", "E", "F"), ncol = 3, nrow = 1)
print(combined_rae_predictors)

# --- 6.3 Visualization for Delta Rhythmic Fraction ---
cat("\nGenerating plots for Delta Rhythmic Fraction...\n")
matched_indices_rf <- which(!is.na(model_data$Rhythmic.fraction_delta) & !is.na(pred_rf))
rf_prediction_plot_data <- data.frame(
  Observed = model_data$Rhythmic.fraction_delta[matched_indices_rf],
  Predicted = pred_rf[matched_indices_rf],
  Climate.col = model_data$Climate.col[matched_indices_rf]
)
plot_obs_pred_rf <- ggplot(data=rf_prediction_plot_data, aes(x=Observed, y=Predicted)) +
  geom_point(aes(color=Climate.col), size=3, show.legend = FALSE) +
  scale_color_identity() +
  stat_smooth(method=lm ,formula = y  ~ x, color="black", se=FALSE) +
  labs(x = "Observed ?? RF [a.u] (22°C - 35°C)", y = "Predicted ?? RF [a.u]",
       title= "Model Fit: Delta Rhythmic Fraction Plasticity") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     text = element_text(size = 11.5, family = "Calibri"),
                     axis.text = element_text(size = 10, colour = 'black'))
print(plot_obs_pred_rf)

# Plot relationship between the significant environmental predictor and ?? RF.
# For ?? RF, Max PAR was the primary significant predictor.
plot_rf_par <- ggpubr::ggscatter(model_data, x = "PAR_season_max", y = "Rhythmic.fraction_delta",
                                 xlab = "Max PAR\n[W/m^2]", ylab = "?? Rhythmic Fraction [a.u] (22°C - 35°C)") +
  stat_smooth(method = lm, formula = y ~ x) + theme_bw() + theme(panel.grid = element_blank())
print(plot_rf_par)


# --- 7. Reset global na.action ---
options(na.action = "na.omit")