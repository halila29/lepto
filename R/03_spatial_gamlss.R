# =============================================================================
# Script 03: Spatial GAMLSS Modeling
# =============================================================================
# Manuscript: "Flood Exposure and Spatial Patterns of Leptospirosis in
#              Kelantan, Malaysia"
# Journal:    GeoHealth (American Geophysical Union)
# Authors:    Che Norhalila Che Mohamed, Anis Kausar Ghazali, et al.
# -----------------------------------------------------------------------------
# PURPOSE:
#   - Fit eight candidate GAMLSS models (Poisson, NBI, ZINBI x
#     non-spatial / spatial random effect / penalized B-spline)
#   - Select final model based on AIC, RMSE, and overdispersion diagnostics
#   - Estimate subdistrict-level relative risk (RR) and exceedance
#     probability P(RR > 1)
#   - Produce risk maps: raw RR, quintile classification, exceedance probability
# FINAL MODEL: NBI GAMLSS with spatial random effect random(mukim_id)
# REPRODUCES:  Table 2, Table 3, Figure 7(a-c) in the manuscript
# -----------------------------------------------------------------------------
# INPUT:
#   data/cases_by_mukim_with_climate.rds  -- aggregated analysis dataset
#   shapefiles/kelantan.shp               -- subdistrict boundaries
# OUTPUT:
#   outputs/tables/model_comparison.csv
#   outputs/tables/final_model_coefficients.csv
#   outputs/tables/rr_exceedance_by_mukim.csv
#   outputs/figures/Figure7a_RR_Map.tiff
#   outputs/figures/Figure7b_RR_Quintile_Map.tiff
#   outputs/figures/Figure7c_Exceedance_Map.tiff
#   data/model_spatial_gamlss_final.rds   -- saved final model object
# =============================================================================

library(tidyverse)
library(sf)
library(spdep)
library(gamlss)
library(gamlss.add)
library(gamlss.dist)
library(car)
library(ggplot2)
library(viridis)
library(knitr)
library(here)

# Output directories
tables_dir <- here("outputs", "tables")
figs_dir   <- here("outputs", "figures")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figs_dir,   showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# STEP 1: Load data
# =============================================================================
cases       <- readRDS(here("data", "cases_by_mukim_with_climate.rds"))
shape_mukim <- st_read(here("shapefiles", "kelantan.shp"))

# Assign numeric mukim ID for spatial random effect
shape_mukim <- shape_mukim %>%
  mutate(mukim_id = as.numeric(as.factor(MUKIM)))

cases <- left_join(
  cases,
  shape_mukim %>% st_drop_geometry() %>% select(MUKIM, mukim_id),
  by = "MUKIM"
) %>%
  mutate(mukim_id = as.factor(MUKIM))   # character factor for random()

# Spatial centroid coordinates for penalized spline model
coords        <- st_centroid(shape_mukim) %>% st_coordinates()
cases$x_coord <- coords[, 1]
cases$y_coord <- coords[, 2]

# =============================================================================
# STEP 2: Overdispersion check
# =============================================================================
cat("--- Overdispersion check ---\n")
cat("Mean total cases:", round(mean(cases$total_cases), 2), "\n")
cat("Variance total cases:", round(var(cases$total_cases), 2), "\n")
cat("Variance/Mean ratio:", round(var(cases$total_cases) / mean(cases$total_cases), 2), "\n")
hist(cases$total_cases, breaks = 15, col = "skyblue",
     main = "Distribution of Total Cases by Subdistrict",
     xlab = "Total Cases (2011-2024)")

# =============================================================================
# STEP 3: Covariate selection (stepwise) and multicollinearity check
# =============================================================================
# Full model with all candidate covariates
full_nb <- gamlss(
  total_cases ~ median_age + male_ratio + melayu_ratio + warganegara_ratio +
    flood500m_ratio + flood1km_ratio + flood2km_ratio +
    flood30d_ratio + flood60d_ratio +
    Avg_Rainfall + Temp_C + RH + Dew_C,
  family = NBI,
  data   = cases,
  offset = log(cases$Mean_Population)
)

# Stepwise selection using BIC penalty (k = log n)
step_nb <- stepGAIC(full_nb, direction = "both", k = log(nrow(cases)))
summary(step_nb)

# Variance Inflation Factor check
vif_model <- lm(total_cases ~ male_ratio + melayu_ratio + flood1km_ratio + Temp_C,
                data = cases)
cat("\n--- VIF check ---\n")
print(vif(vif_model))

# =============================================================================
# STEP 4: Fit eight candidate models
# =============================================================================
# Population offset: log(mean annual population per subdistrict)
offset_val <- log(cases$Mean_Population)

# -- Non-spatial models --
model_pois  <- gamlss(total_cases ~ male_ratio + melayu_ratio + flood1km_ratio + Temp_C,
                      family = PO,    data = cases, offset = offset_val)

model_nb    <- gamlss(total_cases ~ male_ratio + melayu_ratio + flood1km_ratio + Temp_C,
                      family = NBI,   data = cases, offset = offset_val)

model_zinb  <- gamlss(total_cases ~ male_ratio + melayu_ratio + flood1km_ratio + Temp_C,
                      family = ZINBI, data = cases, offset = offset_val)

# -- Spatial random effect models (random intercept for mukim) --
model_pois_spatial <- gamlss(
  total_cases ~ male_ratio + melayu_ratio + flood1km_ratio + Temp_C + random(mukim_id),
  family = PO,    data = cases, offset = offset_val)

model_nb_spatial <- gamlss(
  total_cases ~ male_ratio + melayu_ratio + flood1km_ratio + Temp_C + random(mukim_id),
  family = NBI,   data = cases, offset = offset_val)

model_zinb_spatial <- gamlss(
  total_cases ~ male_ratio + melayu_ratio + flood1km_ratio + Temp_C + random(mukim_id),
  family = ZINBI, data = cases, offset = offset_val)

# -- Non-spatial NBI (for residual Moran's I diagnostic) --
model_nosmooth <- gamlss(
  total_cases ~ male_ratio + melayu_ratio + flood1km_ratio + Temp_C,
  family = NBI, data = cases, offset = offset_val)

# -- Penalized B-spline spatial model (sensitivity check) --
model_pb <- gamlss(
  total_cases ~ male_ratio + melayu_ratio + flood1km_ratio + Temp_C +
    pb(x_coord) + pb(y_coord),
  family = NBI, data = cases, offset = offset_val)

# =============================================================================
# STEP 5: Model comparison (Table 2 in manuscript)
# =============================================================================
eval_model <- function(model, obs, name) {
  preds <- fitted(model)
  data.frame(
    Model    = name,
    RMSE     = sqrt(mean((obs - preds)^2)),
    MAE      = mean(abs(obs - preds)),
    Pseudo_R2 = 1 - sum((obs - preds)^2) / sum((obs - mean(obs))^2),
    AIC      = AIC(model),
    GAIC     = GAIC(model, k = log(length(obs)))
  )
}

obs <- cases$total_cases

results <- bind_rows(
  eval_model(model_pois,         obs, "Model 1: Poisson"),
  eval_model(model_nb,           obs, "Model 2: NBI"),
  eval_model(model_zinb,         obs, "Model 3: ZINBI"),
  eval_model(model_pois_spatial, obs, "Model 4: Poisson + Spatial Random"),
  eval_model(model_nb_spatial,   obs, "Model 5: NBI + Spatial Random [FINAL]"),
  eval_model(model_zinb_spatial, obs, "Model 6: ZINBI + Spatial Random"),
  eval_model(model_nosmooth,     obs, "Model 7: NBI (No Spatial)"),
  eval_model(model_pb,           obs, "Model 8: NBI + pb(x,y)")
) %>%
  arrange(AIC) %>%
  mutate(Delta_AIC = round(AIC - min(AIC), 2))

print(knitr::kable(results, digits = 3,
      caption = "Table 2. Model comparison for eight candidate GAMLSS models"))
write_csv(results, file.path(tables_dir, "model_comparison.csv"))

# =============================================================================
# STEP 6: Final model — NBI + spatial random effect (Model 5)
# =============================================================================
# Model 5 selected: lowest RMSE (3.06), lowest MAE (2.20),
# pseudo-R2 = 0.999, theoretically appropriate for overdispersed count data
# (variance = 16,456; mean = 130.6; variance/mean ratio >> 1)

cat("\n--- Final Model Summary (Model 5: NBI + Spatial Random) ---\n")
summary(model_nb_spatial)

# Extract and save fixed-effect coefficients
coef_tbl <- data.frame(
  Term        = names(coef(model_nb_spatial)),
  Estimate    = coef(model_nb_spatial),
  Std_Error   = coef(model_nb_spatial, what = "se"),
  row.names   = NULL
) %>%
  mutate(
    IRR       = exp(Estimate),
    IRR_lower = exp(Estimate - 1.96 * Std_Error),
    IRR_upper = exp(Estimate + 1.96 * Std_Error)
  )

write_csv(coef_tbl, file.path(tables_dir, "final_model_coefficients.csv"))

# =============================================================================
# STEP 7: Subdistrict relative risk and exceedance probability P(RR > 1)
# =============================================================================
# P(RR > 1): exceedance probability that a subdistrict's fitted risk
# exceeds 1.0 (the null/average risk), computed from the linear predictor
# and its standard error via the normal approximation.
# This probabilistic measure avoids frequentist confidence interval framing
# and is consistent with the Bayesian spatial epidemiology literature
# (Moraga, 2023).

pred_link <- predict(model_nb_spatial, type = "link",   se.fit = TRUE)
pred_resp <- predict(model_nb_spatial, type = "response")

# Normalize RR relative to population offset
cases$fitted_mu  <- pred_resp
cases$fitted_rr  <- pred_resp / exp(offset_val)   # RR relative to null

# Exceedance probability P(RR > 1): P(linear predictor > 0)
# Note: mu = exp(eta + offset); RR > 1 iff eta > -offset
cases$eta    <- pred_link$fit
cases$eta_se <- pred_link$se.fit

# P(RR > 1) = P(eta > -log(offset)) via normal CDF
log_offset      <- offset_val
cases$exc_prob  <- pnorm((cases$eta - (-log_offset)) / cases$eta_se,
                          lower.tail = TRUE)

# Table 3: RR and P(RR > 1) per subdistrict (sorted by exceedance probability)
rr_tbl <- cases %>%
  select(MUKIM, total_cases, Mean_Population, fitted_mu, fitted_rr, exc_prob) %>%
  arrange(desc(exc_prob)) %>%
  mutate(
    fitted_rr  = round(fitted_rr, 3),
    exc_prob   = round(exc_prob, 3)
  )

write_csv(rr_tbl, file.path(tables_dir, "rr_exceedance_by_mukim.csv"))
cat("\nTop 10 high-risk subdistricts (Table 3):\n")
print(head(rr_tbl, 10))

# =============================================================================
# STEP 8: Risk maps (Figure 7a-c in manuscript)
# =============================================================================
# Join fitted values to shapefile
shape_mukim <- shape_mukim %>%
  left_join(
    cases %>% select(MUKIM, fitted_mu, fitted_rr, exc_prob),
    by = "MUKIM"
  )

# Figure 7a: Raw fitted RR map
p_rr <- ggplot(shape_mukim) +
  geom_sf(aes(fill = fitted_rr), color = "black", linewidth = 0.2) +
  scale_fill_viridis(name = "Relative Risk", direction = -1) +
  labs(title = "Fitted Relative Risk of Leptospirosis by Subdistrict, Kelantan (2011-2024)",
       subtitle = "Spatial GAMLSS: NBI + spatial random effect") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

ggsave(file.path(figs_dir, "Figure7a_RR_Map.tiff"),
       p_rr, width = 8, height = 7, units = "in", dpi = 600,
       compression = "lzw")

# Figure 7b: Quintile classification map
shape_mukim$rr_quintile <- cut(
  shape_mukim$fitted_rr,
  breaks         = quantile(shape_mukim$fitted_rr, probs = seq(0, 1, 0.2), na.rm = TRUE),
  labels         = c("Q1 (Lowest)", "Q2", "Q3", "Q4", "Q5 (Highest)"),
  include.lowest = TRUE
)

p_quintile <- ggplot(shape_mukim) +
  geom_sf(aes(fill = rr_quintile), color = "black", linewidth = 0.2) +
  scale_fill_viridis_d(name = "RR Quintile", direction = -1) +
  labs(title = "Relative Risk Quintile Classification by Subdistrict, Kelantan",
       subtitle = "Q5 = Highest 20% of predicted leptospirosis risk") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

ggsave(file.path(figs_dir, "Figure7b_RR_Quintile_Map.tiff"),
       p_quintile, width = 8, height = 7, units = "in", dpi = 600,
       compression = "lzw")

# Figure 7c: Exceedance probability map P(RR > 1)
p_exc <- ggplot(shape_mukim) +
  geom_sf(aes(fill = exc_prob), color = "black", linewidth = 0.2) +
  scale_fill_viridis(name = "P(RR > 1)", limits = c(0, 1), direction = -1) +
  labs(title = "Exceedance Probability P(RR > 1) by Subdistrict, Kelantan",
       subtitle = "Probability that a subdistrict's leptospirosis risk exceeds the null (RR = 1)") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

ggsave(file.path(figs_dir, "Figure7c_Exceedance_Probability_Map.tiff"),
       p_exc, width = 8, height = 7, units = "in", dpi = 600,
       compression = "lzw")

# =============================================================================
# STEP 9: Save final model object
# =============================================================================
saveRDS(model_nb_spatial, here("data", "model_spatial_gamlss_final.rds"))
saveRDS(shape_mukim,      here("data", "shape_gamlss_output.rds"))

cat("\nScript 03 complete.\n")
cat("Final model saved: data/model_spatial_gamlss_final.rds\n")
cat("Risk tables saved to:", tables_dir, "\n")
cat("Risk maps saved to:", figs_dir, "\n")
