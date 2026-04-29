#Title : Monitoring the colonisation dynamics of Rugulopteryx okamurae on algal communities in the Northwestern Mediterranean sea
#Authors : Marie Borriglione, Thierry Thibaut, Bastien Thouroude, Aurélie Blanfuné,
#           Frédéric Zuberer, Dorian Guillemain, Delphine Thibault, Charles-François Boudouresque, 
#           Sandrine Ruitton
# Journal : Aquatic conservation
# Date    : 2026
#
## =============================================================================
# Analysis of Rugulopteryx okamurae biomass across depth and time
# =============================================================================
# Description:
#   This script analyses the dry weight biomass of the invasive brown alga
#   Rugulopteryx okamurae across two depths (5m and 15m) and nine sampling
#   dates over one year at two sites in Marseille (France).
#
#   Statistical approach:
#   - Generalized Linear Mixed Model (GLMM) with Gamma distribution and
#     log-link function (glmmTMB), appropriate for continuous positive,
#     right-skewed biomass data with no zeros.
#   - Site included as a random effect (random representation of the
#     broader study area, not a fixed treatment of interest).
#   - Model selection via AIC comparing nested random-effect structures.
#   - Significance of fixed effects tested by likelihood ratio tests (LRT),
#     comparing nested models (robust to large Z-statistics from strong
#     seasonal effects).
#   - Post-hoc pairwise comparisons (Tukey-adjusted) via emmeans to
#     decompose the significant Depth x Time interaction.
#
# Data structure:
#   Each row = one quadrat measured once (independent quadrats per date).
#   Columns used: MS_Rug (dry weight, g/25cm²), SITE, PROFONDEUR (depth: 5 or 15m),
#                 DATE2 (sampling date), Quadrat (unique quadrat ID)
# =============================================================================


# --- 0. Packages --------------------------------------------------------------

library(glmmTMB)   # GLMM fitting
library(DHARMa)    # Residual diagnostics via simulation
library(emmeans)   # Estimated marginal means and post-hoc tests
library(ggplot2)   # Visualization
library(dplyr)     # Data manipulation
library(tidyr)
# --- 1. Data preparation ------------------------------------------------------

# Convert date to factor and remove the March 2022 sampling date
# (excluded due to insufficient sampling effort at this date)
raw <- readxl::read_xlsx("Biomasse_rugulopteryx.xlsx")
biom_rugu <- raw |>
  mutate(date = as.factor(DATE2)) |>
  filter(date != "2022-03-03")

# Response variable: MS_Rug * 25 converts g/25cm² to g/m²
# Check data structure
str(biom_rugu)
table(biom_rugu$SITE, biom_rugu$PROFONDEUR, biom_rugu$DATE2)


# --- 2. Model fitting ---------------------------------------------------------

# --- 2a. Full model: site + quadrat as random effects
# (1 | SITE) : among-site variation
# (1 | SITE:Quadrat) : quadrat nested within site
# Retained for comparison only — see AIC below
mod_glmm <- glmmTMB(
  (MS_Rug * 25) ~ as.factor(PROFONDEUR) * as.factor(DATE2) +
    (1 | SITE) +
    (1 | SITE:Quadrat),
  data   = biom_rugu,
  family = Gamma(link = "log")
)

# --- 2b. Simplified model: site only as random effect
# Justified because each quadrat is measured only once (no repeated measures
# at the quadrat level) → (1 | SITE:Quadrat) variance ≈ 0
mod_final <- glmmTMB(
  (MS_Rug * 25) ~ as.factor(PROFONDEUR) * as.factor(DATE2) +
    (1 | SITE),
  data   = biom_rugu,
  family = Gamma(link = "log")
)

# --- 2c. Model selection: compare random effect structures
AIC(mod_glmm, mod_final)
# ΔAIC < 2 → mod_final (simpler) is preferred


# --- 3. Model diagnostics (DHARMa) -------------------------------------------

# Simulate residuals (n = 1000 for stable results)
residus_glmm <- simulateResiduals(fittedModel = mod_final, n = 1000)

# Overall diagnostic plot (QQ + residuals vs. fitted)
plot(residus_glmm)

# Individual tests
testDispersion(residus_glmm)    # Should be non-significant (p > 0.05)
testZeroInflation(residus_glmm) # Should be non-significant (no zeros expected)
testOutliers(residus_glmm)      # Flag any extreme residuals

# Identify outlying observations if present
outliers(residus_glmm)

# Check for convergence issues and large Z-statistics
# Note: large Z-statistics on seasonal effects reflect strong biological signal,
# not model failure. LRT-based inference (see below) is robust to this.
diagnose(mod_final)


# --- 4. Significance of fixed effects (Likelihood Ratio Tests) ---------------
#
# Because the Depth x Time interaction is present, main effects cannot be
# suppressed independently (principle of marginality). We therefore compare
# nested models to test each term.
#
# LRT is preferred over Wald tests here (robust to large Z-statistics).

# Additive model (no interaction) — used as baseline for main effect tests
mod_add <- glmmTMB(
  (MS_Rug * 25) ~ as.factor(PROFONDEUR) + as.factor(DATE2) +
    (1 | SITE),
  data = biom_rugu, family = Gamma(link = "log")
)

# Model without depth
mod_nodepth <- glmmTMB(
  (MS_Rug * 25) ~ as.factor(DATE2) + (1 | SITE),
  data = biom_rugu, family = Gamma(link = "log")
)

# Model without time
mod_notime <- glmmTMB(
  (MS_Rug * 25) ~ as.factor(PROFONDEUR) + (1 | SITE),
  data = biom_rugu, family = Gamma(link = "log")
)

# LRT: Depth x Time interaction (full vs. additive)
anova(mod_add, mod_final)

# LRT: effect of Depth (additive vs. no depth)
anova(mod_nodepth, mod_add)

# LRT: effect of Time (additive vs. no time)
anova(mod_notime, mod_add)

# Summary table of LRT results (reported in Table 1):
#   Effect        | Df | chi²   | p
#   Depth         |  1 | 94.63  | < 0.0001
#   Time          |  7 | 134.64 | < 0.0001
#   Depth x Time  |  7 | 34.85  | < 0.0001


# --- 5. Post-hoc pairwise comparisons ----------------------------------------
#
# The significant Depth x Time interaction means that depth effects depend on
# date (and vice versa). Post-hoc tests decompose this interaction.

# --- 5a. Effect of TIME within each depth
# (How does biomass change across dates, separately at 5m and 15m?)
emm_time <- emmeans(mod_final, ~ DATE2 | PROFONDEUR)
pairs(emm_time, adjust = "tukey")
# Results on the log scale (model scale)

# --- 5b. Effect of DEPTH at each date
# (Is 5m different from 15m, and does this vary across dates?)
emm_depth <- emmeans(mod_final, ~ PROFONDEUR | DATE2)
pairs(emm_depth, adjust = "tukey")
writexl::write_xlsx(as.data.frame(pairs(emm_depth, adjust = "tukey")), "posthoc_depth_by_date.xlsx")
# Results on the log scale (model scale)


# --- 6. Figure: estimated marginal means -------------------------------------
#Back-transformed to the response scale (g/m²) with 95% confidence intervals

 #emm_fig <- as.data.frame(
 #emmeans(mod_final, ~ DATE2 * PROFONDEUR, type = "response")
 

 #ggplot(emm_fig,
        #aes(x     = as.Date(DATE2),
            #y     = response,
            #color = as.factor(PROFONDEUR),
            #group = as.factor(PROFONDEUR))) +
   #geom_point(size = 3) +
   #geom_line() +
   #geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 5) +
   #scale_color_manual(
     #values = c("5" = "#2196F3", "15" = "#FF5722"),
     #name   = "Depth (m)"
   #) +
   #labs(
     #x = "Date",
     #y = expression("Estimated biomass (g m"^{-2}*")")
   #) +
   #theme_bw(base_size = 13) +
   #theme(axis.text.x = element_text(hjust = 0.5, vjust = 0.5))

# =============================================================================
# PART 2: BIOMASS VISUALIZATION ON RAW DATA
# =============================================================================

# --- 6. Compute summary statistics -------------------------------------------

r_biomasse_site <- biom_rugu %>%
  group_by(SITE, PROFONDEUR, DATE2) %>%
  summarise(
    mean_Rug_MS = mean(MS_Rug, na.rm = TRUE),
    se = sd(MS_Rug, na.rm = TRUE) / sqrt(sum(!is.na(MS_Rug))),
    .groups = "drop"
  ) %>%
  mutate(PROFONDEUR = factor(PROFONDEUR, levels = c(5, 15), labels = c("5 m", "15 m")))

# Temperature data in long format
r_temp_long <- biom_rugu %>%
  select(Date_temp, temp_5, temp_15) %>%
  distinct() %>%
  pivot_longer(cols = starts_with("temp_"), names_to = "PROFONDEUR_TEMP", values_to = "temp") %>%
  mutate(PROFONDEUR_TEMP = recode(PROFONDEUR_TEMP, "temp_5" = "5 m", "temp_15" = "15 m"))


# --- 7. Biomass × Temperature plot function ----------------------------------

plot_biomasse_temp_site <- function(site_name) {
  
  dat_bio <- r_biomasse_site %>% filter(SITE == site_name)
  
  # Scaling factor for dual y-axis
  range_bio  <- range(dat_bio$mean_Rug_MS, na.rm = TRUE)
  range_temp <- range(r_temp_long$temp, na.rm = TRUE)
  coef <- diff(range_bio) / diff(range_temp)
  if (!is.finite(coef) || coef == 0) coef <- 1
  coef <- coef * 0.65
  
  # Axis breaks
  
  bio_max <- max(dat_bio$mean_Rug_MS + dat_bio$se, na.rm = TRUE)
  bio_breaks <- seq(0, ceiling(bio_max / 100) * 100, by = 100)
  temp_breaks <- seq(floor(min(r_temp_long$temp, na.rm = TRUE) / 5) * 5,
                     ceiling(max(r_temp_long$temp, na.rm = TRUE) / 5) * 5, by = 5)
  
  ggplot() +
    geom_hline(yintercept = 0, linewidth = 0.6, color = "black") +
    # Biomass
    geom_line(data = dat_bio,
              aes(x = DATE2, y = mean_Rug_MS, color = PROFONDEUR, group = PROFONDEUR),
              linewidth = 0.7) +
    geom_point(data = dat_bio, aes(x = DATE2, y = mean_Rug_MS, color = PROFONDEUR), size = 2) +
    geom_errorbar(data = dat_bio,
                  aes(x = DATE2, ymin = mean_Rug_MS - se, ymax = mean_Rug_MS + se, color = PROFONDEUR),
                  width = 3 * 24 * 3600) +
    # Temperature (secondary axis)
    geom_line(data = r_temp_long,
              aes(x = Date_temp, y = temp * coef, linetype = PROFONDEUR_TEMP),
              linewidth = 0.7, color = "grey20") +
    scale_y_continuous(
      name = expression("Dry biomass (g m"^-2*")"),
      breaks = bio_breaks,
      sec.axis = sec_axis(~ . / coef, name = "Temperature (°C)", breaks = temp_breaks)
    ) +
    scale_color_manual(values = c("5 m" = "#1E90FF", "15 m" = "#00008B"), name = "Depth (biomass)") +
    scale_linetype_manual(values = c("5 m" = "dashed", "15 m" = "solid"), name = "Depth (temp.)") +
    labs(title = paste("Biomass & temperature –", site_name), x = "") +
    theme_minimal() +
    theme(legend.position = "top", legend.title = element_text(face = "bold"), panel.grid.minor = element_blank())
}

# Generate plots
p_callelongue <- plot_biomasse_temp_site("Callelongue")
p_maire       <- plot_biomasse_temp_site("Maire")

p_callelongue
p_maire


# =============================================================================
# PART 3: COMMUNITY COMPOSITION ANALYSIS
# =============================================================================

library(vegan)
library(tidyr)
library(forcats)
library(readxl)

data <- read_excel("data_quadrat-colonisation.xlsx", col_names = TRUE)

# Set chronological order of sampling months
data$mois <- fct_relevel(data$mois,
                         "juillet_1", "septembre", "novembre", "janvier",
                         "avril", "mai", "juin", "juillet_2"
)


# --- 8. Stacked barplot of relative abundances -------------------------------

# Reshape to long format
data_bar <- data %>%
  pivot_longer(cols = -c(site, quadrat, mois), names_to = "espece", values_to = "abondance")

# Compute mean abundance per site × month × species
data_bar2 <- data_bar %>%
  group_by(site, mois, espece) %>%
  summarise(mean_abondance = mean(abondance), .groups = "drop")

# Identify and pool rare species (<1% in all combinations)
rare_species <- data_bar2 %>%
  group_by(espece) %>%
  summarise(always_rare = all(mean_abondance < 1)) %>%
  filter(always_rare) %>%
  pull(espece)

data_bar2 <- data_bar2 %>%
  mutate(espece = if_else(espece %in% rare_species, "Others", espece)) %>%
  group_by(site, mois, espece) %>%
  summarise(mean_abondance = sum(mean_abondance), .groups = "drop")

# Order species by total abundance
ordre_especes <- data_bar2 %>%
  filter(espece != "Others") %>%
  group_by(espece) %>%
  summarise(total = sum(mean_abondance)) %>%
  arrange(total) %>%
  pull(espece)
ordre_especes <- c("Others", ordre_especes)

# Color palette
colors <- c(
  "rock" = "black", "Dictyota spp." = "#D9C5B6FF", "asparagopsis_armata" = "#A56A3EFF",
  "Coralline" = "#CFB267FF", "Turf" = "#693829FF", "rugulopteryx_okamurae" = "#5480B5FF",
  "lithophyllum_incrustans" = "#345084FF"
)
final_colors <- c("Others" = "gray50", colors)

# Plot
ggplot(data_bar2, aes(x = mois, y = mean_abondance, fill = espece)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ site) +
  scale_fill_manual(values = final_colors) +
  labs(x = "", y = "Mean relative abundance (%)", fill = "") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Export summary table
tableau_stats <- data_bar %>%
  group_by(site, mois, espece) %>%
  summarise(mean = mean(abondance), sd = sd(abondance), .groups = "drop")
writexl::write_xlsx(tableau_stats, "tableau_abondances_stats.xlsx")


# --- 9. NMDS ordination ------------------------------------------------------

data_nmds <- data.frame(data)
data_nmds$id <- paste0(data_nmds$site, "-", data_nmds$quadrat, "-", data_nmds$mois)
rownames(data_nmds) <- data_nmds$id

meta_data <- data_nmds %>% select(id, site, quadrat, mois)
mat_nmds  <- data_nmds %>% select(-site, -mois, -quadrat, -id)

res_nmds <- metaMDS(mat_nmds, k = 2, distance = "bray")

# Species vectors
species_fit <- envfit(res_nmds, mat_nmds)
sp_vectors <- as.data.frame(species_fit$vectors$arrows)
sp_vectors$r <- species_fit$vectors$r

# Site scores
scores_sites <- as.data.frame(scores(res_nmds, display = "sites"))
scores_sites$id <- rownames(scores_sites)
scores_sites <- left_join(scores_sites, meta_data, by = "id")

# Select species to display
species_to_show <- c(1, 2, 3, 4, 5, 9, 10)

# NMDS plot
ggplot(scores_sites, aes(x = NMDS1, y = NMDS2, shape = site, color = mois)) +
  geom_point(size = 2.8, alpha = 0.9) +
  scale_color_manual(values = c(
    "juillet_1" = "#752305FF", "septembre" = "#AC420AFF", "novembre" = "#CE8A37FF",
    "janvier" = "#F8CA7CFF", "avril" = "#76A1CDFF", "mai" = "#4E79A5FF",
    "juin" = "#1B3A6BFF", "juillet_2" = "#0D1D3CFF"
  )) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_segment(data = sp_vectors[species_to_show, ],
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", inherit.aes = FALSE) +
  geom_text(data = sp_vectors[species_to_show, ],
            aes(x = NMDS1, y = NMDS2, label = rownames(sp_vectors)[species_to_show]),
            inherit.aes = FALSE, hjust = 1, vjust = 1, size = 3) +
  labs(x = "NMDS1", y = "NMDS2", title = "NMDS (Bray-Curtis)", color = "Month", shape = "Site") +
  theme_bw()


# --- 10. PERMANOVA and PERMDISP ----------------------------------------------

# Prepare data
meta_perm <- data %>%
  mutate(
    site = factor(site),
    quadrat = factor(quadrat),
    mois = factor(mois),
    quadrat_id = interaction(site, quadrat, drop = TRUE)
  ) %>%
  select(site, quadrat, mois, quadrat_id)

mat_commu <- data %>%
  select(-site, -quadrat, -mois) %>%
  mutate(across(everything(), as.numeric))

# Remove empty rows
rows_ok <- rowSums(mat_commu) > 0
mat_commu <- mat_commu[rows_ok, ]
meta_perm <- meta_perm[rows_ok, ]

dist_bray <- vegdist(mat_commu, method = "bray")

# PERMANOVA: test effect of Site × Month (permutations within quadrats)
res_permanova <- adonis2(
  dist_bray ~ site * mois,
  data = meta_perm,
  permutations = 999,
  strata = meta_perm$quadrat_id,
  by = "term"
)
print(res_permanova)

# PERMDISP: test homogeneity of dispersion among months
disp_month <- betadisper(dist_bray, meta_perm$mois)
permutest(disp_month, permutations = 999)
TukeyHSD(disp_month)

