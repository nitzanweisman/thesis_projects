# --- PCA of Sorghum Ecotype Environmental Factors ---

# Author: Nitzan Weisman
# Date: [Date you finalize the script]
# Description: Explores the environmental landscape of sorghum ecotype origins.
#              Performs Principal Component Analysis (PCA) on selected environmental
#              variables to visualize major gradients and relationships.
#              Includes initial screening for multicollinearity among environmental parameters.
# Data Source: Adapted from Master's Thesis, The Hebrew University of Jerusalem (March 2024)
# Input File: data/Sorghum_initial_screening_2.csv

# --- 1. Load Required Packages ---
# install.packages(c('tidyverse', 'vegan', 'factoextra', 'ggrepel', 'GGally', 'showtext'))
library(tidyverse) # Includes readr, dplyr, ggplot2
library(vegan)     # For decostand() standardization
library(factoextra) # For PCA visualization (fviz_pca_biplot)
library(ggrepel)   # For non-overlapping text labels in plots
library(GGally)    # For ggpairs() correlation matrix plots
library(showtext)  # For custom font support if needed

# --- Optional: Custom Font Setup ---
# font_add(family = "Calibri", regular = "Calibri.ttf")
# showtext_auto()

# --- 2. Load and Prepare Data ---

# Load the main dataset
sorghum_data <- readr::read_csv("data/Sorghum_initial_screening_2.csv")

# --- 3. Initial Exploration of All Environmental Parameters (for Collinearity Assessment) ---
# This section demonstrates the process of examining relationships between the original
# set of environmental parameters to inform the selection of a reduced set for PCA.

cat("\n--- 3.1. Collinearity Exploration: Temperature-related Parameters ---\n")
temp_params <- sorghum_data %>%
  dplyr::select(Alt, tmax_1, tmax_8, tmin_1, tmin_8, tmean_1, tmean_8,
                t_season_max, t_season_range, t_season_delta, Ann.Mean.Tmp)
print(round(cor(temp_params, use = "pairwise.complete.obs"), 2)) # Correlation matrix
# print(ggpairs(temp_params, title = "Temperature Parameter Correlations")) # Visual matrix

cat("\n--- 3.2. Collinearity Exploration: Precipitation-related Parameters ---\n")
prec_params <- sorghum_data %>%
  dplyr::select(prec_1, prec_8, prec_season_delta, prec_season_max, Ann.Prc)
print(round(cor(prec_params, use = "pairwise.complete.obs"), 2))
# print(ggpairs(prec_params, title = "Precipitation Parameter Correlations"))

cat("\n--- 3.3. Collinearity Exploration: PET and VPD Parameters ---\n")
pet_vpd_params <- sorghum_data %>%
  dplyr::select(pet_1, pet_8, pet_season_max, pet_season_delta, pet_annual,
                vpd_1, vpd_8, vpd_season_max, vpd_season_delta)
print(round(cor(pet_vpd_params, use = "pairwise.complete.obs"), 2))
# print(ggpairs(pet_vpd_params, title = "PET & VPD Parameter Correlations"))

cat("\n--- 3.4. Collinearity Exploration: PAR Parameters ---\n")
par_params <- sorghum_data %>%
  dplyr::select(PAR_season_max, PAR_season_delta, `PAROct-Dec`, `PARJul-Sep`) # Note backticks for special characters
print(round(cor(par_params, use = "pairwise.complete.obs"), 2))
# print(ggpairs(par_params, title = "PAR Parameter Correlations"))

cat("\n--- 3.5. Collinearity Exploration: Other Parameters ---\n")
other_params <- sorghum_data %>%
  dplyr::select(Depth.to.Water.Table, soil_m, topsoil_pH, gs.length)
print(round(cor(other_params, use = "pairwise.complete.obs"), 2))
# print(ggpairs(other_params, title = "Other Parameter Correlations"))

# Rationale for Parameter Selection (from original script's intent):
# Based on the above (and more detailed external examination if done for thesis),
# a subset of 13 parameters was chosen to represent key environmental axes while
# reducing redundancy for a clearer PCA. This selection aims to balance
# comprehensive environmental characterization with interpretability.

# --- 4. Select Final Set of Environmental Parameters for PCA ---
env_parameters_selected <- sorghum_data %>%
  dplyr::select(
    Ann.Mean.Tmp, t_season_delta, pet_annual, pet_season_delta,
    vpd_season_max, vpd_season_delta, Ann.Prc, prec_season_delta,
    PAR_season_max, Depth.to.Water.Table, soil_m, topsoil_pH,
    gs.length
  )

# Rename columns for clearer labels in the PCA plot
env_parameters_renamed <- env_parameters_selected %>%
  dplyr::rename(
    "Ann Mean Temp" = Ann.Mean.Tmp,
    "Temp Seasonal Delta" = t_season_delta,
    "Ann PET" = pet_annual,
    "PET Seasonal Delta" = pet_season_delta,
    "VPD Max" = vpd_season_max,
    "VPD Seasonal Delta" = vpd_season_delta,
    "Ann Precipitation" = Ann.Prc,
    "Precip Seasonal Delta" = prec_season_delta,
    "PAR Max" = PAR_season_max,
    "Water Table Depth" = Depth.to.Water.Table,
    "Soil Depth" = soil_m,
    "Topsoil pH" = topsoil_pH,
    "Growing Season Length" = gs.length
  )

# --- 5. Standardize Data for PCA ---
# Standardize variables (mean=0, sd=1)
env_data_scaled <- vegan::decostand(env_parameters_renamed, 'standardize', margin = 2)

# Check for NAs introduced by standardization
if(any(is.na(env_data_scaled))) {
  warning("NAs introduced during standardization. Check variables with zero variance.")
  env_data_scaled <- na.omit(env_data_scaled) # Example: remove rows with NAs
  # This might cause issues if it misaligns with climate_groups
}

# --- 6. Perform PCA ---
# Perform PCA using the standardized environmental data
pca_results <- stats::prcomp(env_data_scaled, scale. = FALSE, center = TRUE)

# Display summary to see variance explained by each PC
print(summary(pca_results))

# --- 7. Visualize PCA Results ---
# Define climate groups for coloring points
# Ensure row alignment if na.omit was used on env_data_scaled
if (nrow(env_data_scaled) != nrow(sorghum_data)) {
  warning("Number of rows changed after NA removal in scaled data. Grouping might be misaligned.")
  # Consider filtering sorghum_data to match rows in env_data_scaled if this happens
  # For simplicity, assuming no rows were removed, or they are handled appropriately.
  # If rows were removed by na.omit(env_data_scaled), you'd need to filter sorghum_data accordingly:
  # sorghum_data_filtered <- sorghum_data[complete.cases(env_parameters_renamed), ]
  # climate_groups <- as.factor(sorghum_data_filtered$Climate)
} else {
  climate_groups <- as.factor(sorghum_data$Climate)
}


# Create the PCA biplot using factoextra
pca_biplot <- fviz_pca_biplot(
  pca_results,
  geom = c("point"),
  col.ind = climate_groups,
  palette = c("#0078ff", "#53b1fc", "#faa800", "#96ff96", "#64C864"),
  col.var = "black",
  addEllipses = TRUE,
  ellipse.type = "convex",
  repel = TRUE,
  label = "var",
  title = "PCA of Selected Sorghum Environmental Factors"
)

# Customize plot appearance
pca_biplot_customized <- pca_biplot +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 11.5, family = "Calibri"), # Use "sans" if Calibri not available
    axis.text = element_text(size = 10, colour = 'black'),
    legend.text = element_text(size = 10, colour = 'black'),
    legend.title = element_text(size=11)
  ) +
  labs(title = NULL) +
  scale_fill_manual(
    name = "Climate classification:",
    values = c("#0078ff", "#53b1fc", "#faa800", "#96ff96", "#64C864"),
    labels = c("Am | Tropical", "Aw | Tropical", "BSh | Hot Semi-Arid", "Cwa | Temperate", "Cwb | Temperate")
  ) +
  scale_color_manual(
    name = "Climate classification:",
    values = c("#0078ff", "#53b1fc", "#faa800", "#96ff96", "#64C864"),
    labels = c("Am | Tropical", "Aw | Tropical", "BSh | Hot Semi-Arid", "Cwa | Temperate", "Cwb | Temperate")
  ) +
  scale_shape_discrete(
    name = "Climate classification:",
    labels = c("Am | Tropical", "Aw | Tropical", "BSh | Hot Semi-Arid", "Cwa | Temperate", "Cwb | Temperate")
  )

# Display the final plot
print(pca_biplot_customized)

# --- 8. Save Output (Optional) ---
# if (!dir.exists("plots")) dir.create("plots")
# if (!dir.exists("output")) dir.create("output")
# ggsave("plots/PCA_Biplot_Sorghum_Environment_Selected.png", plot = pca_biplot_customized, width = 8, height = 6, dpi = 300)
# saveRDS(pca_results, file = "output/pca_results_sorghum_env_selected.rds")