# --- PCA of Sorghum Ecotype Environmental Factors ---

# Author: Nitzan Weisman
# Date: 25/5/25
# Description: This script explores the environmental parameters associated with
#              sorghum ecotype origin sites. It includes steps for assessing and
#              reducing co-linearity and redundancy among numerous environmental
#              variables. Finally, Principal Component Analysis (PCA) is performed
#              on a selected set of environmental variables to visualize the major
#              axes of environmental variation and the relationships between these variables.
# Data Source: Adapted from Master's Thesis, The Hebrew University of Jerusalem (March 2024)
# Input File: C:/Users/nitza/Documents/thesis_project_for_github/Sorghum_initial_screening_2.csv

# --- 1. Load Required Packages ---
# Load packages for data handling, statistics, and plotting.
# install.packages(c('tidyverse', 'vegan', 'factoextra', 'ggrepel', 'GGally', 'showtext'))
library(tidyverse)
library(vegan)
library(factoextra)
library(ggrepel)
library(GGally)
library(showtext)

# --- Optional: Custom Font Setup ---
# font_add(family = "Calibri", regular = "Calibri.ttf")
# showtext_auto()

# --- 2. Load and Prepare Data ---
# Load the main dataset containing sorghum ecotype information and associated environmental parameters.
sorghum_data_raw <- readr::read_csv("C:/Users/nitza/Documents/thesis_project_for_github/Sorghum_initial_screening_2.csv")

# --- 3. Initial Exploration & Collinearity Assessment of Environmental Parameters ---
# The goal here is to understand the relationships between the full set of
# environmental parameters to identify highly correlated (collinear) and potentially
# redundant variables. This informs the selection of a more concise set for PCA.

cat("\n--- 3.1. Collinearity Exploration: Temperature-related Parameters ---\n")
# Select and examine temperature-related variables.
temp_params <- sorghum_data_raw %>%
  dplyr::select(Alt, tmax_1, tmax_8, tmin_1, tmin_8, tmean_1, tmean_8,
                t_season_max, t_season_range, t_season_delta, Ann.Mean.Tmp)
# Print correlation matrix to identify highly correlated pairs.
print(round(cor(temp_params, use = "pairwise.complete.obs"), 2))
# Visual correlation matrix for a graphical overview.
print(ggpairs(temp_params, title = "Temperature Parameter Correlations"))

cat("\n--- 3.2. Collinearity Exploration: Precipitation-related Parameters ---\n")
prec_params <- sorghum_data_raw %>%
  dplyr::select(prec_1, prec_8, prec_season_delta, prec_season_max, Ann.Prc)
print(round(cor(prec_params, use = "pairwise.complete.obs"), 2))
print(ggpairs(prec_params, title = "Precipitation Parameter Correlations"))

cat("\n--- 3.3. Collinearity Exploration: PET and VPD Parameters ---\n")
pet_vpd_params <- sorghum_data_raw %>%
  dplyr::select(pet_1, pet_8, pet_season_max, pet_season_delta, pet_annual,
                vpd_1, vpd_8, vpd_season_max, vpd_season_delta)
print(round(cor(pet_vpd_params, use = "pairwise.complete.obs"), 2))
print(ggpairs(pet_vpd_params, title = "PET & VPD Parameter Correlations"))

cat("\n--- 3.4. Collinearity Exploration: PAR Parameters ---\n")
par_params <- sorghum_data_raw %>%
  dplyr::select(PAR_season_max, PAR_season_delta, `PAROct-Dec`, `PARJul-Sep`)
print(round(cor(par_params, use = "pairwise.complete.obs"), 2))
print(ggpairs(par_params, title = "PAR Parameter Correlations"))

cat("\n--- 3.5. Collinearity Exploration: Other Parameters ---\n")
other_params <- sorghum_data_raw %>%
  dplyr::select(Depth.to.Water.Table, soil_m, topsoil_pH, gs.length)
print(round(cor(other_params, use = "pairwise.complete.obs"), 2))
print(ggpairs(other_params, title = "Other Parameter Correlations"))

# Rationale for Parameter Selection:
# Following the collinearity assessment, a subset of 13 parameters is chosen.
# This selection aims to retain key environmental information while reducing
# redundancy, leading to a more interpretable PCA.

# --- 4. Select Final Set of Environmental Parameters for PCA ---
# These 13 parameters are selected to represent diverse environmental aspects
# after considering their inter-correlations.
env_parameters_selected <- sorghum_data_raw %>%
  dplyr::select(
    Ann.Mean.Tmp, t_season_delta, pet_annual, pet_season_delta,
    vpd_season_max, vpd_season_delta, Ann.Prc, prec_season_delta,
    PAR_season_max, Depth.to.Water.Table, soil_m, topsoil_pH,
    gs.length
  )

# Rename columns for more descriptive and reader-friendly labels in the PCA plot.
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
# Standardization (e.g., to mean 0, sd 1) is crucial before PCA when variables
# are measured on different scales or have vastly different variances.
# This ensures all variables contribute more equally to the analysis.
env_data_scaled <- vegan::decostand(env_parameters_renamed, 'standardize', margin = 2)

# Check for NAs after standardization, which can occur if a column had zero variance.
if(any(is.na(env_data_scaled))) {
  warning("NAs introduced during standardization. This might indicate columns with zero variance. Review data.")
  # PCA cannot be performed on columns with no variance or all NAs.
}

# --- 6. Perform PCA ---
# Perform Principal Component Analysis on the standardized set of selected environmental variables.
# 'scale. = FALSE' because data was already standardized.
# 'center = TRUE' ensures data is centered before rotation (default for prcomp).
pca_results <- stats::prcomp(env_data_scaled, scale. = FALSE, center = TRUE)

# Display a summary of the PCA results.
# This shows the proportion of variance explained by each Principal Component (PC),
# which helps in deciding how many PCs are important to interpret.
print(summary(pca_results))

# --- 7. Visualize PCA Results ---
# Define climate groups from the original data for coloring points in the PCA plot.
climate_groups <- as.factor(sorghum_data_raw$Climate)

# Create the PCA biplot using factoextra's fviz_pca_biplot.
# A biplot visualizes:
#   - Principal component scores: Position of each sorghum ecotype origin site in the new PC space.
#   - Principal component loadings: Arrows representing the original environmental variables,
#     showing their contribution and correlation with the PCs.
pca_biplot <- fviz_pca_biplot(
  pca_results,
  geom = c("point"),                # Show sites as points.
  col.ind = climate_groups,         # Color points by 'Climate' group.
  palette = c("#0078ff", "#53b1fc", "#faa800", "#96ff96", "#64C864"), # Custom colors.
  col.var = "black",                # Variable arrows in black.
  addEllipses = TRUE,               # Draw ellipses around climate groups.
  ellipse.type = "convex",          # Ellipse type (convex hull).
  repel = TRUE,                     # Avoid overlapping text labels for variables.
  label = "var",                    # Display labels for variables (arrows) only.
  title = "PCA of Selected Sorghum Environmental Factors"
)

# Customize plot appearance for better readability and presentation.
pca_biplot_customized <- pca_biplot +
  theme_bw() + # Use a black and white theme.
  theme(
    panel.grid.major = element_blank(), # Remove major grid lines.
    panel.grid.minor = element_blank(), # Remove minor grid lines.
    text = element_text(size = 11.5, family = "Calibri"), # Set font style (ensure font is available).
    axis.text = element_text(size = 10, colour = 'black'),
    legend.text = element_text(size = 10, colour = 'black'),
    legend.title = element_text(size=11)
  ) +
  labs(title = NULL) + # Remove default title from fviz_pca_biplot if a custom one is preferred.
  # Manually define legend for consistency.
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

# Display the final customized PCA plot.
# This plot will show how ecotype origins cluster based on the selected environmental variables
# and how those variables relate to the principal axes of variation.
print(pca_biplot_customized)