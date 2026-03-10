# =============================================================================
# Script 04: Sensitivity Analysis and Residual Diagnostics
# =============================================================================
# Manuscript: "Flood Exposure and Spatial Patterns of Leptospirosis in
#              Kelantan, Malaysia"
# Journal:    GeoHealth (American Geophysical Union)
# Authors:    Che Norhalila Che Mohamed, Anis Kausar Ghazali, et al.
# -----------------------------------------------------------------------------
# PURPOSE:
#   - Residual diagnostics for non-spatial and final spatial GAMLSS models:
#     histogram, Q-Q plot, GAMLSS worm plot
#   - Moran's I test on model residuals to assess residual spatial
#     autocorrelation (pre- and post-spatial adjustment)
#   - Marginal effect plots for all four selected covariates
#   - Absolute prediction error maps (NBI spatial vs. Poisson spatial)
# REPRODUCES: Supplementary diagnostics in the manuscript
# -----------------------------------------------------------------------------
# INPUT:
#   data/cases_by_mukim_with_climate.rds
#   data/model_spatial_gamlss_final.rds  -- saved from Script 03
#   data/shape_gamlss_output.rds         -- shapefile with fitted values
#   shapefiles/kelantan.shp
# OUTPUT:
#   outputs/figures/FigS1_Residuals_NoSpatial.tiff
#   outputs/figures/FigS2_Residuals_Spatial.tiff
#   outputs/figures/FigS3_Marginal_Effects.tiff
#   outputs/figures/FigS4_Prediction_Error_Map.tiff
# =============================================================================

library(tidyverse)
library(sf)
library(spdep)
library(gamlss)
library(gamlss.add)
library(gamlss.dist)
library(ggplot2)
library(viridis)
library(patchwork)
library(here)

figs_dir <- here("outputs", "figures")
dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# STEP 1: Load data and models
# =============================================================================
cases            <- readRDS(here("data", "cases_by_mukim_with_climate.rds"))
model_nb_spatial <- readRDS(here("data", "model_spatial_gamlss_final.rds"))
shape_mukim      <- readRDS(here("data", "shape_gamlss_output.rds"))

offset_val <- log(cases$Mean_Population)

# Re-fit non-spatial model for pre-spatial diagnostics
model_nosmooth <- gamlss(
  total_cases ~ male_ratio + melayu_ratio + flood1km_ratio + Temp_C,
  family = NBI, data = cases, offset = offset_val
)

# Spatial weights (queen contiguity)
nb_q   <- poly2nb(shape_mukim, queen = TRUE)
listw  <- nb2listw(nb_q, style = "W", zero.policy = TRUE)

# =============================================================================
# STEP 2: Residual diagnostics — Non-spatial model
# =============================================================================
resid_nosmooth <- residuals(model_nosmooth)

# Moran's I on non-spatial residuals
cat("--- Moran's I on non-spatial model residuals ---\n")
moran_nosmooth <- moran.test(resid_nosmooth, listw, zero.policy = TRUE)
print(moran_nosmooth)

# Diagnostic plots
tiff(file.path(figs_dir, "FigS1_Residuals_NoSpatial.tiff"),
     width = 10, height = 4, units = "in", res = 600, compression = "lzw")
par(mfrow = c(1, 3))
hist(resid_nosmooth, col = "lightblue", main = "Pearson Residuals\n(Non-spatial NBI)",
     xlab = "Residuals", cex.main = 1)
qqnorm(resid_nosmooth, main = "Q-Q Plot\n(Non-spatial NBI)", cex.main = 1)
qqline(resid_nosmooth, col = "red", lwd = 2)
wp(model_nosmooth, ylim.all = 2.5, main = "Worm Plot\n(Non-spatial NBI)")
dev.off()

# =============================================================================
# STEP 3: Residual diagnostics — Final spatial model
# =============================================================================
resid_spatial <- residuals(model_nb_spatial)

# Moran's I on spatial model residuals
cat("--- Moran's I on final spatial model residuals ---\n")
moran_spatial <- moran.test(resid_spatial, listw, zero.policy = TRUE)
print(moran_spatial)

# Diagnostic plots
tiff(file.path(figs_dir, "FigS2_Residuals_Spatial.tiff"),
     width = 10, height = 4, units = "in", res = 600, compression = "lzw")
par(mfrow = c(1, 3))
hist(resid_spatial, col = "lightblue",
     main = "Pearson Residuals\n(NBI + Spatial Random)",
     xlab = "Residuals", cex.main = 1)
qqnorm(resid_spatial, main = "Q-Q Plot\n(NBI + Spatial Random)", cex.main = 1)
qqline(resid_spatial, col = "red", lwd = 2)
wp(model_nb_spatial, ylim.all = 2.5, main = "Worm Plot\n(NBI + Spatial Random)")
dev.off()

# =============================================================================
# STEP 4: Marginal effect plots for each covariate (from non-spatial model)
# =============================================================================
# Marginal effects computed from model_nosmooth, holding other variables
# constant at their sample means, to aid interpretability of fixed effects

newdata_base <- list(
  male_ratio    = mean(cases$male_ratio),
  melayu_ratio  = mean(cases$melayu_ratio),
  flood1km_ratio= mean(cases$flood1km_ratio),
  Temp_C        = mean(cases$Temp_C)
)

# flood1km_ratio effect
df_flood <- data.frame(
  flood1km_ratio = seq(min(cases$flood1km_ratio), max(cases$flood1km_ratio), length.out = 100),
  male_ratio     = newdata_base$male_ratio,
  melayu_ratio   = newdata_base$melayu_ratio,
  Temp_C         = newdata_base$Temp_C
)
df_flood$fit <- predict(model_nosmooth, newdata = df_flood, type = "response")

p1 <- ggplot(df_flood, aes(x = flood1km_ratio, y = fit)) +
  geom_line(color = "#1f78b4", linewidth = 1) +
  labs(title = "Flood Exposure Within 1 km", x = "Flood Exposure Proportion",
       y = "Predicted Cases (mu)") +
  theme_minimal(base_size = 11)

# male_ratio effect
df_male <- data.frame(
  male_ratio     = seq(min(cases$male_ratio), max(cases$male_ratio), length.out = 100),
  flood1km_ratio = newdata_base$flood1km_ratio,
  melayu_ratio   = newdata_base$melayu_ratio,
  Temp_C         = newdata_base$Temp_C
)
df_male$fit <- predict(model_nosmooth, newdata = df_male, type = "response")

p2 <- ggplot(df_male, aes(x = male_ratio, y = fit)) +
  geom_line(color = "#e31a1c", linewidth = 1) +
  labs(title = "Male Population Ratio", x = "Male Ratio",
       y = "Predicted Cases (mu)") +
  theme_minimal(base_size = 11)

# melayu_ratio effect
df_melayu <- data.frame(
  melayu_ratio   = seq(min(cases$melayu_ratio), max(cases$melayu_ratio), length.out = 100),
  flood1km_ratio = newdata_base$flood1km_ratio,
  male_ratio     = newdata_base$male_ratio,
  Temp_C         = newdata_base$Temp_C
)
df_melayu$fit <- predict(model_nosmooth, newdata = df_melayu, type = "response")

p3 <- ggplot(df_melayu, aes(x = melayu_ratio, y = fit)) +
  geom_line(color = "#33a02c", linewidth = 1) +
  labs(title = "Malay Ethnicity Ratio", x = "Malay Ratio",
       y = "Predicted Cases (mu)") +
  theme_minimal(base_size = 11)

# Temp_C effect
df_temp <- data.frame(
  Temp_C         = seq(min(cases$Temp_C), max(cases$Temp_C), length.out = 100),
  flood1km_ratio = newdata_base$flood1km_ratio,
  male_ratio     = newdata_base$male_ratio,
  melayu_ratio   = newdata_base$melayu_ratio
)
df_temp$fit <- predict(model_nosmooth, newdata = df_temp, type = "response")

p4 <- ggplot(df_temp, aes(x = Temp_C, y = fit)) +
  geom_line(color = "#ff7f00", linewidth = 1) +
  labs(title = "Average Temperature", x = "Temperature (degrees C)",
       y = "Predicted Cases (mu)") +
  theme_minimal(base_size = 11)

p_marginal <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title    = "Marginal Effect Plots: Non-spatial NBI GAMLSS Model",
    subtitle = "Predicted mu as a function of each covariate, others held at mean"
  )

ggsave(file.path(figs_dir, "FigS3_Marginal_Effects.tiff"),
       p_marginal, width = 10, height = 8, units = "in", dpi = 600,
       compression = "lzw")

# =============================================================================
# STEP 5: Absolute prediction error maps
# =============================================================================
# Re-fit Poisson spatial model for comparison
model_pois_spatial <- gamlss(
  total_cases ~ male_ratio + melayu_ratio + flood1km_ratio + Temp_C + random(mukim_id),
  family = PO, data = cases, offset = offset_val
)

cases$abs_error_nb   <- abs(cases$total_cases - fitted(model_nb_spatial))
cases$abs_error_pois <- abs(cases$total_cases - fitted(model_pois_spatial))

shape_mukim <- shape_mukim %>%
  left_join(
    cases %>% select(MUKIM, abs_error_nb, abs_error_pois),
    by = "MUKIM"
  )

p_err_nb <- ggplot(shape_mukim) +
  geom_sf(aes(fill = abs_error_nb), color = "black", linewidth = 0.2) +
  scale_fill_viridis(name = "Absolute\nError", direction = -1) +
  labs(title = "NBI + Spatial Random Effect") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

p_err_pois <- ggplot(shape_mukim) +
  geom_sf(aes(fill = abs_error_pois), color = "black", linewidth = 0.2) +
  scale_fill_viridis(name = "Absolute\nError", direction = -1) +
  labs(title = "Poisson + Spatial Random Effect") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

p_err <- p_err_nb + p_err_pois +
  plot_annotation(
    title    = "Absolute Prediction Error Maps by Subdistrict",
    subtitle = "Lower values indicate better model fit"
  )

ggsave(file.path(figs_dir, "FigS4_Prediction_Error_Map.tiff"),
       p_err, width = 14, height = 6, units = "in", dpi = 600,
       compression = "lzw")

cat("\nScript 04 complete.\n")
cat("Diagnostic figures saved to:", figs_dir, "\n")

# Summary: Moran's I results
cat("\n--- Summary: Spatial autocorrelation in residuals ---\n")
cat("Non-spatial NBI model: I =",
    round(moran_nosmooth$estimate["Moran I statistic"], 3),
    "  p =", round(moran_nosmooth$p.value, 3), "\n")
cat("Spatial NBI model:     I =",
    round(moran_spatial$estimate["Moran I statistic"], 3),
    "  p =", round(moran_spatial$p.value, 3), "\n")
