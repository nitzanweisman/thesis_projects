# --- Multi-Model Averaging for Sorghum Circadian Plasticity ---

# Author: Nitzan Weisman
# Date: [Date you finalize the script]
# Description: Uses multi-model averaging (MMA) based on AICc to identify
#              environmental predictors of circadian trait plasticity (?? Amplitude,
#              ?? RAE, ?? RF) in sorghum ecotypes responding to temperature change.
# Data Source: Adapted from Master's Thesis, The Hebrew University of Jerusalem (March 2024)
# Input File: data/Sorghum_initial_screening_2.csv

# --- 1. Load Required Packages ---
# install.packages(c('tidyverse', 'car', 'MuMIn', 'ggpubr', 'gridExtra', 'showtext'))
library(tidyverse)
library(car)
library(MuMIn)
library(ggpubr)
library(gridExtra)
library(showtext)

# --- Optional: Custom Font Setup ---
# font_add(family = "Calibri", regular = "Calibri.ttf")
# showtext_auto()

# --- 2. Load and Prepare Data ---
sorghum_data <- readr::read_csv("data/Sorghum_initial_screening_2.csv")

# Calculate plasticity traits (?? = 22°C - 35°C)
# Important: The script section flips the sign of deltas.
# This assumes your input CSV might have delta calculated as (35C - 22C).
# Verify this against your data and adjust if your '..._delta' columns already represent (22C - 35C).
sorghum_plasticity <- sorghum_data %>%
  dplyr::mutate(
    Amplitude_delta = -Amplitude_delta,
    Period_delta = -Period_delta,
    RAE_delta = -RAE_delta,
    Rfd_avg_delta = -Rfd_avg_delta,
    Rhythmic.fraction_delta = -Rhythmic.fraction_delta
  )

# Select response variables and the 13 environmental predictors
model_data_full <- sorghum_plasticity %>%
  dplyr::select(
    # Response Variables
    Amplitude_delta, RAE_delta, Rhythmic.fraction_delta, Period_delta, Rfd_avg_delta,
    # Predictor Variables
    Ann.Mean.Tmp, t_season_delta, pet_annual, pet_season_delta,
    vpd_season_max, vpd_season_delta, Ann.Prc, prec_season_delta,
    PAR_season_max, Depth.to.Water.Table, soil_m, topsoil_pH,
    gs.length,
    # Keep 'Origin' and 'Climate.col' for plotting if needed later, but exclude from models
    Origin, Climate.col
  )

# Models require complete cases for predictors and the specific response variable.
# This will be handled per model after formula definition.


# --- 3. Modeling ?? Amplitude ---

cat("\n--- Modeling Delta Amplitude ---\n")

# Initial full model
# Ensure na.action = na.fail for model building steps before dredge
zero_model_amp <- lm(
  Amplitude_delta ~ Ann.Mean.Tmp + t_season_delta + pet_annual +
    pet_season_delta + vpd_season_max + vpd_season_delta + Ann.Prc +
    prec_season_delta + PAR_season_max + Depth.to.Water.Table +
    soil_m + topsoil_pH + gs.length,
  data = model_data_full, # Use the full data here, lm handles NAs by casewise deletion by default if not na.fail
  na.action = na.fail # Critical for dredge
)
cat("Initial model summary (Delta Amplitude):\n"); print(summary(zero_model_amp))

# Iterative VIF reduction
cat("\nPerforming VIF reduction for Delta Amplitude model...\n")
print(vif(zero_model_amp))
zero_model_amp <- update(zero_model_amp, ~ . - vpd_season_delta); print(vif(zero_model_amp))
zero_model_amp <- update(zero_model_amp, ~ . - pet_annual); print(vif(zero_model_amp))
zero_model_amp <- update(zero_model_amp, ~ . - vpd_season_max); print(vif(zero_model_amp))
zero_model_amp <- update(zero_model_amp, ~ . - Depth.to.Water.Table); print(vif(zero_model_amp))
zero_model_amp <- update(zero_model_amp, ~ . - prec_season_delta); print(vif(zero_model_amp))
zero_model_amp <- update(zero_model_amp, ~ . - pet_season_delta); print(vif(zero_model_amp))
cat("Final VIF after reduction (Delta Amplitude):\n"); print(vif(zero_model_amp))
cat("Reduced model summary (Delta Amplitude):\n"); print(summary(zero_model_amp))

# Model selection and averaging
options(na.action = "na.fail") # Required for dredge
dredge_amp <- MuMIn::dredge(zero_model_amp)
options(na.action = "na.omit") # Reset to default

supported_amp <- subset(dredge_amp, delta < 2)
cat("\nSupported models (Delta Amplitude):\n"); print(supported_amp)

avg_model_amp <- MuMIn::model.avg(supported_amp, fit = TRUE)
cat("\n--- Averaged Model Summary (Delta Amplitude) ---\n")
print(summary(avg_model_amp))

# Assess model fit (Observed vs Predicted)
# Need to use data with complete cases for THIS specific model's variables
# Get complete cases for the variables in the final averaged model for prediction
final_amp_vars <- all.vars(formula(avg_model_amp$formula[[1]])) # Get vars from one of the component models
model_data_amp_complete <- model_data_full[complete.cases(model_data_full[, final_amp_vars]), ]

pred_amp <- predict(avg_model_amp, newdata = model_data_amp_complete, type = "response", backtransform = FALSE)
obs_amp <- model_data_amp_complete$Amplitude_delta
fit_assessment_amp <- summary(lm(pred_amp ~ obs_amp))
cat("\n--- Fit Assessment R-squared (Delta Amplitude):", round(fit_assessment_amp$adj.r.squared, 4), "---\n")
cat("--- Fit Assessment p-value (Delta Amplitude):", format.pval(fit_assessment_amp$coefficients[2, 4]), "---\n")


# --- 4. Modeling ?? RAE ---

cat("\n--- Modeling Delta RAE ---\n")
zero_model_rae <- lm(
  RAE_delta ~ Ann.Mean.Tmp + t_season_delta + pet_annual +
    pet_season_delta + vpd_season_max + vpd_season_delta + Ann.Prc +
    prec_season_delta + PAR_season_max + Depth.to.Water.Table +
    soil_m + topsoil_pH + gs.length,
  data = model_data_full,
  na.action = na.fail
)
cat("Initial model summary (Delta RAE):\n"); print(summary(zero_model_rae))

cat("\nPerforming VIF reduction for Delta RAE model...\n")
print(vif(zero_model_rae))
zero_model_rae <- update(zero_model_rae, ~ . - vpd_season_delta); print(vif(zero_model_rae))
zero_model_rae <- update(zero_model_rae, ~ . - pet_annual); print(vif(zero_model_rae))
zero_model_rae <- update(zero_model_rae, ~ . - vpd_season_max); print(vif(zero_model_rae))
zero_model_rae <- update(zero_model_rae, ~ . - pet_season_delta); print(vif(zero_model_rae))
zero_model_rae <- update(zero_model_rae, ~ . - Depth.to.Water.Table); print(vif(zero_model_rae))
zero_model_rae <- update(zero_model_rae, ~ . - prec_season_delta); print(vif(zero_model_rae))
cat("Final VIF after reduction (Delta RAE):\n"); print(vif(zero_model_rae))
cat("Reduced model summary (Delta RAE):\n"); print(summary(zero_model_rae))

options(na.action = "na.fail")
dredge_rae <- MuMIn::dredge(zero_model_rae)
options(na.action = "na.omit")

supported_rae <- subset(dredge_rae, delta < 2)
cat("\nSupported models (Delta RAE):\n"); print(supported_rae)

avg_model_rae <- MuMIn::model.avg(supported_rae, fit = TRUE)
cat("\n--- Averaged Model Summary (Delta RAE) ---\n")
print(summary(avg_model_rae))

# Assess model fit
final_rae_vars <- all.vars(formula(avg_model_rae$formula[[1]]))
model_data_rae_complete <- model_data_full[complete.cases(model_data_full[, final_rae_vars]), ]
pred_rae <- predict(avg_model_rae, newdata = model_data_rae_complete, type = "response", backtransform = FALSE)
obs_rae <- model_data_rae_complete$RAE_delta
fit_assessment_rae <- summary(lm(pred_rae ~ obs_rae))
cat("\n--- Fit Assessment R-squared (Delta RAE):", round(fit_assessment_rae$adj.r.squared, 4), "---\n")
cat("--- Fit Assessment p-value (Delta RAE):", format.pval(fit_assessment_rae$coefficients[2, 4]), "---\n")


# --- 5. Modeling ?? Rhythmic Fraction ---

cat("\n--- Modeling Delta Rhythmic Fraction ---\n")
zero_model_rf <- lm(
  Rhythmic.fraction_delta ~ Ann.Mean.Tmp + t_season_delta + pet_annual +
    pet_season_delta + vpd_season_max + vpd_season_delta + Ann.Prc +
    prec_season_delta + PAR_season_max + Depth.to.Water.Table +
    soil_m + topsoil_pH + gs.length,
  data = model_data_full,
  na.action = na.fail
)
cat("Initial model summary (Delta RF):\n"); print(summary(zero_model_rf))

cat("\nPerforming VIF reduction for Delta RF model...\n")
print(vif(zero_model_rf))
zero_model_rf <- update(zero_model_rf, ~ . - vpd_season_delta); print(vif(zero_model_rf))
zero_model_rf <- update(zero_model_rf, ~ . - pet_annual); print(vif(zero_model_rf))
zero_model_rf <- update(zero_model_rf, ~ . - vpd_season_max); print(vif(zero_model_rf))
zero_model_rf <- update(zero_model_rf, ~ . - pet_season_delta); print(vif(zero_model_rf))
zero_model_rf <- update(zero_model_rf, ~ . - Depth.to.Water.Table); print(vif(zero_model_rf))
zero_model_rf <- update(zero_model_rf, ~ . - prec_season_delta); print(vif(zero_model_rf))
cat("Final VIF after reduction (Delta RF):\n"); print(vif(zero_model_rf))
cat("Reduced model summary (Delta RF):\n"); print(summary(zero_model_rf))

options(na.action = "na.fail")
dredge_rf <- MuMIn::dredge(zero_model_rf)
options(na.action = "na.omit")

supported_rf <- subset(dredge_rf, delta < 2)
cat("\nSupported models (Delta RF):\n"); print(supported_rf)

avg_model_rf <- MuMIn::model.avg(supported_rf, fit = TRUE)
cat("\n--- Averaged Model Summary (Delta RF) ---\n")
print(summary(avg_model_rf))

# Assess model fit
final_rf_vars <- all.vars(formula(avg_model_rf$formula[[1]]))
model_data_rf_complete <- model_data_full[complete.cases(model_data_full[, final_rf_vars]), ]
pred_rf <- predict(avg_model_rf, newdata = model_data_rf_complete, type = "response", backtransform = FALSE)
obs_rf <- model_data_rf_complete$Rhythmic.fraction_delta
fit_assessment_rf <- summary(lm(pred_rf ~ obs_rf))
cat("\n--- Fit Assessment R-squared (Delta RF):", round(fit_assessment_rf$adj.r.squared, 4), "---\n")
cat("--- Fit Assessment p-value (Delta RF):", format.pval(fit_assessment_rf$coefficients[2, 4]), "---\n")


# --- 6. Visualization (Example for Delta Amplitude) ---
# Create data frame for plotting Observed vs Predicted for Amplitude, using complete cases for this model
a_d_prediction_plot_data <- data.frame(
  Observed = obs_amp,
  Predicted = pred_amp,
  Climate.col = model_data_amp_complete$Climate.col # Use climate info from the data subset used for prediction
)

plot_obs_pred_amp <- ggplot(data=a_d_prediction_plot_data, aes(x=Observed, y=Predicted)) +
  geom_point(aes(color=Climate.col), size=3, show.legend = FALSE) + # color by climate, hide direct legend
  scale_color_identity() + # Use colors directly from Climate.col
  stat_smooth(method=lm ,formula = y  ~ x, color="black", se=FALSE) +
  labs(x = "Observed ?? Amplitude [a.u]", y = "Predicted ?? Amplitude [a.u]",
       title= "Model Fit: Delta Amplitude Plasticity") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    text = element_text(size = 11.5, family = "Calibri"), # Use "sans" if needed
    axis.text = element_text(size = 10, colour = 'black')
  )
print(plot_obs_pred_amp)

# Example plots for significant predictors of Delta Amplitude
# (Assuming 'Ann.Prc', 'gs.length', 'PAR_season_max' are generally significant predictors from your full results)
plot_amp_prec <- ggpubr::ggscatter(model_data_full, x = "Ann.Prc", y = "Amplitude_delta",
                                   xlab = "Ann Precipitation\n[mm/Y]", ylab = "?? Amplitude [a.u]") +
  stat_smooth(method = lm, formula = y ~ x) + theme_bw() + theme(panel.grid = element_blank())

plot_amp_gsl <- ggpubr::ggscatter(model_data_full, x = "gs.length", y = "Amplitude_delta",
                                  xlab = "GSL\n[months]", ylab = "") +
  stat_smooth(method = lm, formula = y ~ x) + theme_bw() + theme(panel.grid = element_blank())

plot_amp_par <- ggpubr::ggscatter(model_data_full, x = "PAR_season_max", y = "Amplitude_delta",
                                  xlab = "Max PAR\n[W/m^2]", ylab = "") +
  stat_smooth(method = lm, formula = y ~ x) + theme_bw() + theme(panel.grid = element_blank())

combined_amp_predictors <- ggpubr::ggarrange(plot_amp_prec, plot_amp_gsl, plot_amp_par,
                                             labels = c("A", "B", "C"),
                                             ncol = 3, nrow = 1)
print(combined_amp_predictors)

# Plot RAE vs Amplitude_delta from original script
r2_rae_amp <- format(summary(lm(RAE_delta ~ Amplitude_delta, data = model_data_full))$r.squared, digits = 3)
plot_rae_vs_amp <- ggpubr::ggscatter(model_data_full, x = "Amplitude_delta", y = "RAE_delta",
                                     xlab = "?? Amplitude [a.u]", ylab = "?? RAE [a.u]") +
  stat_smooth(method = lm, formula = y ~ x) +
  annotate("text", x = min(model_data_full$Amplitude_delta, na.rm=T)*0.9, y = max(model_data_full$RAE_delta, na.rm=T)*0.9, label = paste("R-squared =", r2_rae_amp),
           hjust = 0, vjust = 1, size = 4, color = "black") +
  theme_bw() +
  theme(panel.grid = element_blank(), text = element_text(size = 12, family = "Calibri"))
print(plot_rae_vs_amp)


# --- 7. Save Outputs (Optional) ---
# if (!dir.exists("plots")) dir.create("plots")
# if (!dir.exists("output")) dir.create("output")
# ggsave("plots/Obs_vs_Pred_AmplitudeDelta_MMA.png", plot = plot_obs_pred_amp, width = 7, height = 5, dpi = 300)
# ggsave("plots/Predictor_Plots_AmplitudeDelta_MMA.png", plot = combined_amp_predictors, width = 10, height = 4, dpi = 300)
# ggsave("plots/RAE_vs_AmplitudeDelta_MMA.png", plot = plot_rae_vs_amp, width = 6, height = 5, dpi = 300)
# capture.output(summary(avg_model_amp), file = "output/MMA_Summary_AmplitudeDelta.txt")
# capture.output(summary(avg_model_rae), file = "output/MMA_Summary_RAEDelta.txt")
# capture.output(summary(avg_model_rf), file = "output/MMA_Summary_RFDelta.txt")
# saveRDS(avg_model_amp, file = "output/avg_model_amp_delta.rds") # etc.