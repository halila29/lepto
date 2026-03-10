# =============================================================================
# Script 01: Data Preparation and Aggregation
# =============================================================================
# Manuscript: "Flood Exposure and Spatial Patterns of Leptospirosis in
#              Kelantan, Malaysia"
# Journal:    GeoHealth (American Geophysical Union)
# Authors:    Che Norhalila Che Mohamed, Anis Kausar Ghazali, et al.
# Ethical approvals: USM/JEPeM/KK/23110900 | NMRR-23-03438-R4P
# -----------------------------------------------------------------------------
# PURPOSE:
#   1. Import and harmonize leptospirosis case records, flood reports,
#      subdistrict shapefile, population, and climate data
#   2. Aggregate to subdistrict (mukim) level
#   3. Compute annual incidence rates and flood exposure proportions
#   4. Export analysis-ready dataset: cases_by_mukim_with_climate.rds
# -----------------------------------------------------------------------------
# INPUT FILES (place in project root as described in README):
#   data/Leptos.rds              -- individual leptospirosis case records
#   data/LAPORAN_BANJIR.rds      -- flood event records
#   data/mean_population_by_mukim.csv -- mean annual population by mukim
#   shapefiles/kelantan.shp      -- Kelantan subdistrict boundaries
# OUTPUT:
#   data/cases_by_mukim_with_climate.rds  -- aggregated analysis dataset
# =============================================================================

library(tidyverse)
library(sf)
library(here)
library(moments)
library(stringr)
library(writexl)

# -----------------------------------------------------------------------------
# 1.1  Import subdistrict shapefile
# -----------------------------------------------------------------------------
shp_mukim <- st_read(here("shapefiles", "kelantan.shp"))

# Standardize identifiers
shp_mukim <- shp_mukim %>%
  mutate(
    DAERAH = toupper(trimws(DAERAH)),
    MUKIM  = toupper(trimws(MUKIM))
  )

# Summary table of subdistricts
senarai_mukim <- shp_mukim %>%
  st_drop_geometry() %>%
  select(DAERAH, MUKIM, LELAKI, PEREMPUAN, JUM_JANTIN) %>%
  arrange(DAERAH, MUKIM)

write.csv(senarai_mukim,
          here("data", "Senarai_Mukim_Mengikut_Daerah.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 1.2  Import case records and population
# -----------------------------------------------------------------------------
lepto    <- readRDS(here("data", "Leptos.rds"))
pop_mean <- read.csv(here("data", "mean_population_by_mukim.csv"))

glimpse(lepto)
summary(lepto)

# -----------------------------------------------------------------------------
# 1.3  Outlier and age distribution check
# -----------------------------------------------------------------------------
skew <- skewness(lepto$Umur, na.rm = TRUE)
kurt <- kurtosis(lepto$Umur, na.rm = TRUE)
cat("Skewness:", round(skew, 3), "\nKurtosis:", round(kurt, 3), "\n")

# Kolmogorov-Smirnov normality test (appropriate for n >= 5000)
ks.test(lepto$Umur, "pnorm",
        mean(lepto$Umur, na.rm = TRUE),
        sd(lepto$Umur, na.rm = TRUE))

# -----------------------------------------------------------------------------
# 1.4  Convert categorical variables to factors
# -----------------------------------------------------------------------------
lepto <- lepto %>%
  mutate(
    Jantina              = as.factor(Jantina),
    Keturunan            = as.factor(Keturunan),
    Kewarganegaraan      = as.factor(Kewarganegaraan),
    Status.Diagnosis     = as.factor(Status.Diagnosis),
    Status.Pesakit       = as.factor(Status.Pesakit),
    Jenis.Rawatan        = as.factor(Jenis.Rawatan),
    month                = as.factor(month),
    Tahun                = as.factor(Tahun),
    in_flood_500m        = as.factor(in_flood_500m),
    in_flood_1km         = as.factor(in_flood_1km),
    in_flood_2km         = as.factor(in_flood_2km),
    related_to_flood_30d = as.factor(related_to_flood_30d),
    related_to_flood_60d = as.factor(related_to_flood_60d)
  )

# -----------------------------------------------------------------------------
# 1.5  Re-categorize categorical variables
# -----------------------------------------------------------------------------
# Ethnicity: Malay vs. Other
lepto <- lepto %>%
  mutate(
    Keturunan_recat = ifelse(Keturunan == "Melayu", "Melayu", "Lain-lain"),
    Keturunan_recat = factor(Keturunan_recat,
                             levels = c("Melayu", "Lain-lain"))
  )

# Treatment type: consolidate minor categories
lepto <- lepto %>%
  mutate(
    Rawatan_recat = case_when(
      Jenis.Rawatan == "Jabatan Kecemasan & Kemalangan" ~ "Kecemasan",
      Jenis.Rawatan == "Wad Perubatan"                  ~ "Wad Perubatan",
      Jenis.Rawatan == "Jabatan Pesakit Luar"           ~ "Pesakit Luar",
      TRUE                                               ~ "Lain-lain"
    ),
    Rawatan_recat = factor(Rawatan_recat,
                           levels = c("Kecemasan", "Wad Perubatan",
                                      "Pesakit Luar", "Lain-lain"))
  )

# -----------------------------------------------------------------------------
# 1.6  Aggregate to subdistrict (mukim) level
# -----------------------------------------------------------------------------
cases_by_mukim <- lepto %>%
  group_by(MUKIM) %>%
  summarise(
    total_cases         = n(),
    median_age          = median(Umur, na.rm = TRUE),
    male_ratio          = mean(Jantina == "Lelaki", na.rm = TRUE),
    melayu_ratio        = mean(Keturunan == "Melayu", na.rm = TRUE),
    warganegara_ratio   = mean(Kewarganegaraan == "Warganegara", na.rm = TRUE),
    probable_ratio      = mean(Status.Diagnosis == "Probable", na.rm = TRUE),
    hidup_ratio         = mean(Status.Pesakit == "Hidup", na.rm = TRUE),
    kecemasan_ratio     = mean(Rawatan_recat == "Kecemasan", na.rm = TRUE),
    wadperubatan_ratio  = mean(Rawatan_recat == "Wad Perubatan", na.rm = TRUE),
    pesakitluar_ratio   = mean(Rawatan_recat == "Pesakit Luar", na.rm = TRUE),
    rawatanlain_ratio   = mean(Rawatan_recat == "Lain-lain", na.rm = TRUE),
    flood500m_ratio     = mean(in_flood_500m == "TRUE", na.rm = TRUE),
    flood1km_ratio      = mean(in_flood_1km == "TRUE", na.rm = TRUE),
    flood2km_ratio      = mean(in_flood_2km == "TRUE", na.rm = TRUE),
    flood30d_ratio      = mean(related_to_flood_30d == "TRUE", na.rm = TRUE),
    flood60d_ratio      = mean(related_to_flood_60d == "TRUE", na.rm = TRUE)
  )

# -----------------------------------------------------------------------------
# 1.7  Merge population and compute log offset
# -----------------------------------------------------------------------------
clean_mukim <- function(x) x %>% str_trim() %>% str_squish() %>% str_to_lower()

cases_by_mukim <- cases_by_mukim %>%
  mutate(MUKIM_key = clean_mukim(MUKIM))

pop_mean_clean <- pop_mean %>%
  mutate(mukim_key = clean_mukim(mukim)) %>%
  group_by(mukim_key) %>%
  summarise(Mean_Population = mean(Mean_Population, na.rm = TRUE), .groups = "drop")

cases_by_mukim <- cases_by_mukim %>%
  select(-any_of("Mean_Population")) %>%
  left_join(pop_mean_clean, by = c("MUKIM_key" = "mukim_key")) %>%
  mutate(offset_log_pop = log(Mean_Population)) %>%
  select(-MUKIM_key)

stopifnot("Mean_Population" %in% names(cases_by_mukim))

# -----------------------------------------------------------------------------
# 1.8  Compute annual incidence rate
# -----------------------------------------------------------------------------
# Study period: January 2011 – April 2024 = 13 years and 4 months
n_years <- 13 + 4 / 12

cases_by_mukim <- cases_by_mukim %>%
  mutate(
    annual_incidence_rate     = (total_cases / Mean_Population / n_years) * 100000,
    cumulative_incidence_rate = (total_cases / Mean_Population) * 100000
  )

summary(cases_by_mukim$annual_incidence_rate)
summary(cases_by_mukim$cumulative_incidence_rate)

# Cross-tabulation: cases by year, district, subdistrict
lepto_df <- lepto %>%
  st_drop_geometry() %>%
  filter(!is.na(Tahun)) %>%
  mutate(
    Tahun  = as.integer(as.character(Tahun)),
    DAERAH = str_to_title(DAERAH),
    MUKIM  = str_to_title(MUKIM)
  )

table_lepto <- lepto_df %>%
  group_by(DAERAH, MUKIM, Tahun) %>%
  summarise(Kes = n(), .groups = "drop") %>%
  pivot_wider(names_from = Tahun, values_from = Kes, values_fill = 0) %>%
  arrange(DAERAH, MUKIM)

# -----------------------------------------------------------------------------
# 1.9  Export analysis-ready dataset
# -----------------------------------------------------------------------------
if (!dir.exists(here("data")))    dir.create(here("data"),    recursive = TRUE)
if (!dir.exists(here("outputs"))) dir.create(here("outputs"), recursive = TRUE)

saveRDS(cases_by_mukim, here("data", "cases_by_mukim.rds"))
write.csv(cases_by_mukim, here("data", "DATA_ALL.csv"), row.names = FALSE)
write_xlsx(table_lepto, here("outputs", "CrossTab_Kes_Lepto_ByTahun_Daerah_Mukim.xlsx"))

cat("\nScript 01 complete. Output: data/cases_by_mukim.rds\n")
cat("Note: Climate variables (Temp_C) must be merged externally from ERA5.\n")
cat("Final analysis dataset: data/cases_by_mukim_with_climate.rds\n")
