# =============================================================================
# Script 02: Descriptive Spatial Analysis
# =============================================================================
# Manuscript: "Flood Exposure and Spatial Patterns of Leptospirosis in
#              Kelantan, Malaysia"
# Journal:    GeoHealth (American Geophysical Union)
# Authors:    Che Norhalila Che Mohamed, Anis Kausar Ghazali, et al.
# -----------------------------------------------------------------------------
# PURPOSE:
#   - Global Moran's I for leptospirosis incidence and flood frequency
#     by year (2011-2024), using 999-permutation Monte Carlo tests
#   - Local Indicators of Spatial Association (LISA): High-High, Low-Low,
#     High-Low, Low-High cluster maps
#   - Getis-Ord Gi* hotspot and coldspot maps (95% and 99% thresholds)
#   - Bivariate Moran's I: leptospirosis incidence vs. flood frequency
# REPRODUCES: Table 1, Figure 4, Figure 5, Figure 6 in the manuscript
# -----------------------------------------------------------------------------
# INPUT:
#   data/cases_by_mukim_with_climate.rds  -- from Script 01 + climate merge
#   data/LAPORAN_BANJIR.rds               -- flood event records
#   shapefiles/kelantan.shp               -- subdistrict boundaries
# OUTPUT:
#   outputs/tables/  -- Moran's I tables, LISA class counts, Gi* summaries
#   outputs/figures/ -- LISA facet map (TIFF), Gi* facet map (TIFF),
#                       Bivariate Moran time series (TIFF)
# =============================================================================

library(tidyverse)
library(sf)
library(spdep)
library(ggplot2)
library(viridis)
library(knitr)
library(kableExtra)
library(here)
library(purrr)

# Output directories
tables_dir <- here("outputs", "tables")
figs_dir   <- here("outputs", "figures")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figs_dir,   showWarnings = FALSE, recursive = TRUE)

set.seed(123)

# =============================================================================
# STEP 1: Load and harmonize data
# =============================================================================
case_count   <- readRDS(here("data", "Leptos.rds"))
mukim_sf     <- st_read(here("shapefiles", "kelantan.shp"))
flood_events <- readRDS(here("data", "LAPORAN_BANJIR.rds"))

# Standardize identifiers and CRS (UTM Zone 47N)
mukim_sf <- mukim_sf %>%
  mutate(
    DAERAH = toupper(trimws(DAERAH)),
    MUKIM  = toupper(trimws(MUKIM))
  ) %>%
  st_transform(crs = 32647)

case_count <- case_count %>%
  mutate(
    Tahun  = as.integer(as.character(Tahun)),
    DAERAH = toupper(trimws(DAERAH)),
    MUKIM  = toupper(trimws(MUKIM))
  ) %>%
  filter(between(Tahun, 2011, 2024))

flood_events <- flood_events %>%
  rename(Tahun = TAHUN, DAERAH = DAERAH.y) %>%
  mutate(
    Tahun  = as.integer(as.character(Tahun)),
    DAERAH = toupper(trimws(DAERAH)),
    MUKIM  = toupper(trimws(MUKIM))
  ) %>%
  filter(between(Tahun, 2011, 2024)) %>%
  st_transform(crs = 32647)

# =============================================================================
# STEP 2: Aggregate to subdistrict-year level
# =============================================================================
# Leptospirosis incidence rate per 100,000 population
# NOTE: 2024 is a partial year (January-April); it is retained but NOT
# annualized in the incidence computation.
cases_by_mukimyear <- case_count %>%
  st_drop_geometry() %>%
  group_by(DAERAH, MUKIM, Tahun) %>%
  summarise(
    total_cases = n(),
    pop = suppressWarnings(max(as.numeric(JUM_JANTIN), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    pop            = if_else(is.finite(pop) & pop > 0, pop, NA_real_),
    incidence_rate = if_else(!is.na(pop) & pop > 0,
                             (total_cases / pop) * 1e5, NA_real_)
  )

# Flood event counts per subdistrict-year
flood_by_mukimyear <- flood_events %>%
  st_drop_geometry() %>%
  group_by(DAERAH, MUKIM, Tahun) %>%
  summarise(flood_events = n(), .groups = "drop")

# =============================================================================
# STEP 3: Merge with shapefile and build spatial weights
# =============================================================================
mukim_year_sf <- mukim_sf %>%
  right_join(
    bind_rows(
      cases_by_mukimyear %>% select(DAERAH, MUKIM, Tahun),
      flood_by_mukimyear %>% select(DAERAH, MUKIM, Tahun)
    ) %>% distinct(),
    by = c("DAERAH", "MUKIM")
  ) %>%
  left_join(cases_by_mukimyear, by = c("DAERAH", "MUKIM", "Tahun")) %>%
  left_join(flood_by_mukimyear, by = c("DAERAH", "MUKIM", "Tahun")) %>%
  arrange(Tahun, DAERAH, MUKIM) %>%
  mutate(
    Tahun          = as.integer(Tahun),
    incidence_rate = if_else(is.finite(incidence_rate), incidence_rate, NA_real_),
    flood_events   = if_else(is.finite(flood_events),   flood_events,   NA_real_)
  )

# Build stable queen contiguity weights ONCE on static subdistrict geometry
# This ensures consistent spatial weights across all 14 years
mukim_base <- mukim_sf %>%
  select(DAERAH, MUKIM, geometry) %>%
  arrange(DAERAH, MUKIM)

nb_queen <- poly2nb(mukim_base, queen = TRUE)
lw_queen <- nb2listw(nb_queen, style = "W", zero.policy = TRUE)

# Helper: align year-specific attributes to stable mukim_base ordering
subset_year_align <- function(year_val) {
  dat_attr <- mukim_year_sf %>%
    filter(Tahun == year_val) %>%
    st_drop_geometry() %>%
    select(DAERAH, MUKIM, Tahun, total_cases, pop, incidence_rate, flood_events)

  out <- mukim_base %>%
    left_join(dat_attr, by = c("DAERAH", "MUKIM"))

  out$incidence_rate[is.na(out$incidence_rate)] <- 0
  out$flood_events[is.na(out$flood_events)]     <- 0
  return(out)
}

years <- mukim_year_sf %>%
  st_drop_geometry() %>%
  distinct(Tahun) %>%
  arrange(Tahun) %>%
  pull(Tahun)

# =============================================================================
# STEP 4: Global Moran's I by year (leptospirosis incidence and flood events)
# =============================================================================
# One-sided permutation p-value with add-one correction (Phipson & Smyth, 2010)
perm_pval_greater <- function(mc_obj) {
  obs <- unname(mc_obj$statistic)
  res <- as.numeric(mc_obj$res)
  (sum(res >= obs) + 1) / (length(res) + 1)
}

moran_incidence_tbl <- map_dfr(years, function(yr) {
  dat <- subset_year_align(yr)
  mc  <- moran.mc(dat$incidence_rate, listw = lw_queen, nsim = 999,
                  zero.policy = TRUE, alternative = "greater")
  tibble(Year = yr, Moran_I = as.numeric(mc$statistic),
         P_value = perm_pval_greater(mc), N = length(dat$incidence_rate))
}) %>% arrange(Year)

moran_flood_tbl <- map_dfr(years, function(yr) {
  dat <- subset_year_align(yr)
  mc  <- moran.mc(dat$flood_events, listw = lw_queen, nsim = 999,
                  zero.policy = TRUE, alternative = "greater")
  tibble(Year = yr, Moran_I = as.numeric(mc$statistic),
         P_value = perm_pval_greater(mc), N = length(dat$flood_events))
}) %>% arrange(Year)

write_csv(moran_incidence_tbl, file.path(tables_dir, "moran_incidence_by_year.csv"))
write_csv(moran_flood_tbl,     file.path(tables_dir, "moran_flood_by_year.csv"))

# Combined Moran's I table (Table 1 in manuscript)
moran_combined <- moran_incidence_tbl %>%
  rename(MoranI_Incidence = Moran_I, P_Incidence = P_value) %>%
  left_join(
    moran_flood_tbl %>% rename(MoranI_Flood = Moran_I, P_Flood = P_value),
    by = c("Year", "N")
  ) %>%
  mutate(across(c(MoranI_Incidence, P_Incidence, MoranI_Flood, P_Flood),
                ~round(., 3)))

write_csv(moran_combined, file.path(tables_dir, "moran_combined_table.csv"))

# =============================================================================
# STEP 5: LISA - Local Indicators of Spatial Association
# =============================================================================
compute_lisa_year <- function(yr, listw = lw_queen) {
  dat       <- subset_year_align(yr)
  dat$Tahun <- as.integer(yr)
  x         <- dat$incidence_rate

  lm_obj  <- localmoran(x, listw = listw, zero.policy = TRUE)
  dat$lisa_i <- lm_obj[, "Ii"]
  dat$lisa_z <- (lm_obj[, "Ii"] - lm_obj[, "E.Ii"]) / sqrt(lm_obj[, "Var.Ii"])
  dat$lisa_p <- 2 * pnorm(-abs(dat$lisa_z))

  zx  <- scale(x)[, 1]
  zwx <- scale(lag.listw(listw, x, zero.policy = TRUE))[, 1]

  quad <- ifelse(zx >= 0 & zwx >= 0, "High-High",
          ifelse(zx <  0 & zwx <  0, "Low-Low",
          ifelse(zx >= 0 & zwx <  0, "High-Low", "Low-High")))

  dat$quadrant <- factor(
    ifelse(dat$lisa_p <= 0.05, quad, "Not Significant"),
    levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not Significant")
  )
  dat
}

lisa_all <- lapply(years, compute_lisa_year) %>% bind_rows()
lisa_all$Tahun <- factor(lisa_all$Tahun, levels = 2011:2024)

# Save LISA tables
lisa_named_long <- lisa_all %>%
  st_drop_geometry() %>%
  select(DAERAH, MUKIM, Tahun, lisa_i, lisa_z, lisa_p, quadrant)

lisa_counts_wide <- lisa_named_long %>%
  count(Tahun, quadrant, name = "n") %>%
  complete(Tahun, quadrant, fill = list(n = 0)) %>%
  mutate(Tahun = as.integer(as.character(Tahun))) %>%
  arrange(Tahun, quadrant) %>%
  pivot_wider(names_from = quadrant, values_from = n)

lisa_hotspot_recurring <- lisa_named_long %>%
  filter(quadrant == "High-High") %>%
  count(MUKIM, sort = TRUE, name = "Years_as_Hotspot")

write_csv(lisa_counts_wide,       file.path(tables_dir, "lisa_class_counts_by_year.csv"))
write_csv(lisa_named_long,        file.path(tables_dir, "lisa_local_moran_all.csv"))
write_csv(lisa_hotspot_recurring, file.path(tables_dir, "lisa_recurring_hotspots.csv"))

# LISA facet map (Figure 4 in manuscript)
pal_lisa <- c("High-High"       = "#d7191c",
              "Low-Low"         = "#2c7bb6",
              "High-Low"        = "#fdae61",
              "Low-High"        = "#abdda4",
              "Not Significant" = "grey80")

p_lisa <- ggplot(lisa_all) +
  geom_sf(aes(fill = quadrant), color = "black", linewidth = 0.15) +
  scale_fill_manual(values = pal_lisa, drop = FALSE, name = "Cluster Type") +
  facet_wrap(~ Tahun, ncol = 5) +
  labs(
    title    = "LISA Cluster Map of Leptospirosis Incidence, Kelantan (2011-2024)",
    subtitle = "High-High, Low-Low, and spatial outliers (p <= 0.05)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text   = element_text(size = 11, face = "bold"),
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = "right"
  )

ggsave(file.path(figs_dir, "Figure4_LISA_Cluster_2011_2024.tiff"),
       p_lisa, width = 16, height = 10, units = "in", dpi = 600,
       compression = "lzw")

# =============================================================================
# STEP 6: Getis-Ord Gi* hotspot and coldspot analysis
# =============================================================================
compute_gi_year <- function(yr, listw = lw_queen) {
  dat       <- subset_year_align(yr)
  dat$Tahun <- as.integer(yr)
  x         <- dat$incidence_rate
  giz       <- as.numeric(localG(x, listw = listw, zero.policy = TRUE))
  dat$GiZ   <- giz

  dat$Gi_class <- factor(
    case_when(
      giz >=  2.58 ~ "Hotspot 99%",
      giz >=  1.96 ~ "Hotspot 95%",
      giz <= -2.58 ~ "Coldspot 99%",
      giz <= -1.96 ~ "Coldspot 95%",
      TRUE         ~ "Not significant"
    ),
    levels = c("Hotspot 99%", "Hotspot 95%",
               "Coldspot 99%", "Coldspot 95%", "Not significant")
  )
  dat
}

gi_all <- lapply(years, compute_gi_year) %>% bind_rows()
gi_all$Tahun <- factor(gi_all$Tahun, levels = 2011:2024)

gi_named_long <- gi_all %>%
  st_drop_geometry() %>%
  select(DAERAH, MUKIM, Tahun, GiZ, Gi_class)

gi_counts_wide <- gi_named_long %>%
  count(Tahun, Gi_class, name = "n") %>%
  complete(Tahun, Gi_class, fill = list(n = 0)) %>%
  mutate(Tahun = as.integer(as.character(Tahun))) %>%
  arrange(Tahun, Gi_class) %>%
  pivot_wider(names_from = Gi_class, values_from = n)

gi_hotspot_recurring <- gi_named_long %>%
  filter(Gi_class %in% c("Hotspot 99%", "Hotspot 95%")) %>%
  count(MUKIM, sort = TRUE, name = "Years_as_Hotspot")

write_csv(gi_counts_wide,       file.path(tables_dir, "gi_star_class_counts_by_year.csv"))
write_csv(gi_named_long,        file.path(tables_dir, "gi_star_named_long.csv"))
write_csv(gi_hotspot_recurring, file.path(tables_dir, "gi_star_recurring_hotspots.csv"))

# Gi* facet map (Figure 5 in manuscript)
pal_gi <- c("Hotspot 99%"    = "#B2182B",
            "Hotspot 95%"    = "#EF8A62",
            "Coldspot 99%"   = "#2166AC",
            "Coldspot 95%"   = "#67A9CF",
            "Not significant"= "#DDDDDD")

p_gi <- ggplot(gi_all) +
  geom_sf(aes(fill = Gi_class), color = "black", linewidth = 0.15) +
  scale_fill_manual(values = pal_gi, drop = FALSE, name = "Gi* class") +
  facet_wrap(~ Tahun, ncol = 5) +
  labs(
    title    = "Getis-Ord Gi* Hotspot/Coldspot Map of Leptospirosis Incidence, Kelantan (2011-2024)",
    subtitle = "Z-scores classified at 95% and 99% confidence thresholds"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text   = element_text(size = 11, face = "bold"),
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.ticks   = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    legend.position = "right"
  )

ggsave(file.path(figs_dir, "Figure5_GiStar_Cluster_2011_2024.tiff"),
       p_gi, width = 16, height = 10, units = "in", dpi = 600,
       compression = "lzw")

# =============================================================================
# STEP 7: Bivariate Moran's I (leptospirosis incidence vs. flood frequency)
# =============================================================================
# Bivariate Moran's I measures spatial co-clustering of two variables:
# how much a subdistrict's leptospirosis incidence covaries with its
# neighbors' flood event frequency.

bivar_moran_once <- function(x, y, listw) {
  zx  <- scale(x)[, 1]
  zy  <- scale(y)[, 1]
  wzy <- lag.listw(listw, zy, zero.policy = TRUE)
  sum(zx * wzy, na.rm = TRUE) / sum(zx * zx, na.rm = TRUE)
}

bivar_moran_mc <- function(x, y, listw, nsim = 999) {
  obs  <- bivar_moran_once(x, y, listw)
  sims <- replicate(nsim, bivar_moran_once(x, sample(y, replace = FALSE), listw))
  pval <- (sum(sims >= obs) + 1) / (nsim + 1)   # add-one correction
  list(stat = obs, sims = sims, p = pval)
}

bivar_tbl <- map_dfr(years, function(yr) {
  dat <- subset_year_align(yr)
  x   <- dat$incidence_rate
  y   <- dat$flood_events
  mc  <- bivar_moran_mc(x, y, lw_queen, nsim = 999)
  wy  <- lag.listw(lw_queen, y, zero.policy = TRUE)
  tibble(
    Year        = yr,
    Bivariate_I = as.numeric(mc$stat),
    P_value     = as.numeric(mc$p),
    Cor_x_Wy    = suppressWarnings(cor(x, wy, use = "complete.obs")),
    N_mukim     = length(x)
  )
}) %>%
  arrange(Year) %>%
  mutate(Signif = case_when(
    P_value < 0.001 ~ "***",
    P_value < 0.01  ~ "**",
    P_value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

write_csv(bivar_tbl, file.path(tables_dir, "bivariate_moran_by_year.csv"))

# Bivariate Moran time-series plot (Figure 6 in manuscript)
p_bivar <- ggplot(bivar_tbl, aes(x = Year, y = Bivariate_I)) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_line(linewidth = 0.8, color = "steelblue") +
  geom_point(size = 2.5, color = "steelblue") +
  scale_x_continuous(breaks = 2011:2024) +
  labs(
    title    = "Bivariate Moran's I: Leptospirosis Incidence vs. Flood Frequency, Kelantan (2011-2024)",
    subtitle = "Positive values indicate spatial co-clustering; asterisks denote p < 0.05 (*), 0.01 (**), 0.001 (***)",
    y        = "Bivariate Moran's I",
    x        = "Year"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(figs_dir, "Figure6_Bivariate_Moran_TimeSeries.tiff"),
       p_bivar, width = 10, height = 5, units = "in", dpi = 600,
       compression = "lzw")

cat("\nScript 02 complete.\n")
cat("Tables saved to:", tables_dir, "\n")
cat("Figures saved to:", figs_dir, "\n")
