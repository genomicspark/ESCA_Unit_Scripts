# ==============================================================================
# Figure 2 - Updated R Script (v3)
#
# Key changes from v2:
#   - CSV loaded with check.names = TRUE (default) to match original script
#   - All column references use dot-sanitised names (Days.from.exp.start, etc.)
#   - Removed make.names() rename-back loop from LMM for-loop (no longer needed)
#   - Formula string now identical to original: paste(marker, "~ Time * ...")
#   - K-means cluster colours changed from Set3 (contains yellow) to a
#     high-contrast palette with no yellow
#
# All other fixes from v2 retained:
#   - Deduplication of duplicate Tag+Time rows
#   - Post-inflection heatmap shows row labels
#   - Boxplots show exactly 2 boxes per facet (AL vs CR)
#   - Permutation tests as robust alternative for small n
#   - Cohen's d effect sizes on post-hoc contrasts
#   - K-means clustering panel
#   - FDR correction within test families
#   - Corrected 3-step post-hoc contrast (Time x Diet x AgeGroup)
#
# Required packages:
#   install.packages(c("lme4","lmerTest","emmeans","tidyr","dplyr","ggplot2",
#                      "ggrepel","patchwork","pheatmap","vegan","RColorBrewer",
#                      "scales","cowplot","coin"))
# ==============================================================================

library(lme4)
library(lmerTest)
library(emmeans)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(pheatmap)
library(vegan)
library(RColorBrewer)
library(scales)
library(cowplot)
library(coin)

# ── Font size control ──────────────────────────────────────────────────────────
BASE_FONT <- 19   # 50% increase from original 13 for final publication version

# All figure, table, and supplementary outputs are saved here.
out_dir <- "Figure2_R_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Colour palettes ────────────────────────────────────────────────────────────
# Four-group colours for PCA: Age x Diet combination
# Pre-inflection (Early MA): skyblue = AL, darkblue = CR
# Post-inflection (Late MA): orange  = AL, red      = CR
group_colors <- c(
  "Early middle age_Ad lib" = "#87CEEB",   # skyblue
  "Early middle age_CR"     = "#00008B",   # darkblue
  "Late middle age_Ad lib"  = "#FFA500",   # orange
  "Late middle age_CR"      = "#CC0000"    # red
)

# Kept for heatmap annotation and boxplots
age_colors  <- c("Early middle age" = "#87CEEB",
                 "Late middle age"  = "#FFA500",
                 "Old"              = "#4575b4")
diet_colors <- c("Ad lib" = "#555555", "CR" = "#2166ac")

# High-contrast K-means cluster colours — no yellow, all clearly distinguishable
kmeans_colors <- c("1" = "#e41a1c",   # red
                   "2" = "#377eb8",   # blue
                   "3" = "#4daf4a",   # green
                   "4" = "#984ea3")   # purple

# ── Marker display name mapping ────────────────────────────────────────────────
# Maps dot-sanitised R column names (check.names=TRUE) back to clean display
# names matching the original figure exactly
# Marker display names — keys are the EXACT R check.names=TRUE column names
# (verified by running names(df)[14:49] in R directly)
marker_labels <- c(
  "T.cells.FoL"                  = "T cells FoL",
  "CD4..of.T.Cells"              = "CD4+ of T Cells",
  "CD4..T.cells.FoL"             = "CD4+ T cells FoL",
  "CD25..of.CD4.T.cells"         = "CD25+ of CD4 T cells",
  "CD3..CD4..CD25..FoL"          = "CD3+/CD4+/CD25+ FoL",
  "CD8..of.T.cells"              = "CD8+ of T cells",
  "CD3..CD8..FoL"                = "CD3+/CD8+ FoL",
  "B.cells.FoL"                  = "B cells FoL",
  "NK.cells.FoL"                 = "NK cells FoL",
  "Monocytes.FoL"                = "Monocytes FoL",
  "Neutrophils.FoL"              = "Neutrophils FoL",
  "Stimulated.monocytes.FoP"     = "Stimulated monocytes FoP",
  "Stimulated.monocytes.FoL"     = "Stimulated monocytes FoL",
  "Stimulated.regular.monocytes" = "Stimulated/regular monocytes",
  "Myeloid.lymphoid.ratio..Flow." = "Myeloid/lymphoid ratio (Flow)",
  "WBC..k.ul."                   = "WBC (k/ul)",
  "NE..k.ul."                    = "NE (k/ul)",
  "LY..k.ul."                    = "LY (k/ul)",
  "MO..k.ul."                    = "MO (k/ul)",
  "EO..k.ul."                    = "EO (k/ul)",
  "BA..k.ul."                    = "BA (k/ul)",
  "NE...."                       = "NE (%)",
  "LY...."                       = "LY (%)",
  "MO...."                       = "MO (%)",
  "EO...."                       = "EO (%)",
  "BA...."                       = "BA (%)",
  "RBC..M.ul."                   = "RBC (M/ul)",
  "Hb..g.dl."                    = "Hb (g/dl)",
  "HCT...."                      = "HCT (%)",
  "MCV..fl."                     = "MCV (fl)",
  "MCH..pg."                     = "MCH (pg)",
  "MCHC..g.dl."                  = "MCHC (g/dl)",
  "RDW...."                      = "RDW (%)",
  "PLT..K.ul."                   = "PLT (K/ul)",
  "MPV..fl."                     = "MPV (fl)",
  "Myeloid.Lymphoid.ratio..CBC." = "Myeloid/Lymphoid ratio (CBC)"
)

# Helper: convert dot-names to display names (falls back to original if not found)
# Uses direct vector lookup + unname() to preserve parentheses and special characters.
# ifelse() was avoided because it strips the names attribute of named vectors,
# causing the dot-names to be returned instead of the display names.
to_display <- function(x) {
  result           <- as.character(marker_labels[x])
  result[is.na(result)] <- x[is.na(result)]   # fallback for unmapped names
  unname(result)
}

# ==============================================================================
# 1. LOAD AND PREPARE DATA
# ==============================================================================
# check.names = TRUE (default): converts spaces/special chars to dots
# This matches the original script so column references are consistent
df <- read.csv("Complete data draft2.csv")

# Explicit reference ordering makes the sign of the direct interaction contrast
# unambiguous: Post − Baseline, CR − Ad libitum, and Late MA − Early MA.
df$Time      <- factor(ifelse(df$Days.from.exp.start <= 0, "Baseline", "Post"),
                       levels = c("Baseline", "Post"))
df$Age_Group <- factor(df$Age.group,
                       levels = c("Early middle age", "Late middle age", "Old"))
df$Diet      <- factor(df$Diet, levels = c("Ad lib", "CR"))
df$Tag       <- as.factor(df$Tag.)

# Marker columns — dot-sanitised names from check.names = TRUE
markers <- names(df)[14:49]

# Ensure all marker columns are numeric
for (m in markers) df[[m]] <- suppressWarnings(as.numeric(df[[m]]))

# ── Deduplication ──────────────────────────────────────────────────────────────
# Remove only exact same-day duplicates (same Tag + same Day, one row all-NA).
# Do NOT deduplicate by Tag+Time: each rat has multiple legitimate post-
# intervention measurement dates that are correctly averaged in Section 2.
# Calculate the number of available marker measurements with base R first.
# This avoids a dplyr data-masking error from rowSums(across(...)).
df$.n_valid <- rowSums(!is.na(as.matrix(df[, markers, drop = FALSE])))
df <- df %>%
  group_by(Tag, Days.from.exp.start) %>%
  slice_max(order_by = .n_valid, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-.n_valid)

cat("After deduplication:", nrow(df), "rows (original: 347)\n")

# ── Exclude Old age group ───────────────────────────────────────────────────────
# Analysis restricted to Early middle age (Pre-inflection) and
# Late middle age (Post-inflection) only.
df <- df %>% filter(Age_Group != "Old")
df$Age_Group <- droplevels(df$Age_Group)
cat("After excluding Old:", nrow(df), "rows\n")

# ==============================================================================
# 2. COMPUTE PER-RAT DELTAS (Day 105 − Baseline)
# ==============================================================================
# Original figure legend: "delta between study start (baseline) and study end (day 105)"
# Use only the last available measurement per rat as "study end" (day 105 or nearest)

baseline_df <- df %>%
  filter(Time == "Baseline") %>%
  group_by(Tag, Age_Group, Diet) %>%
  summarise(across(all_of(markers), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

post_df <- df %>%
  filter(Time == "Post") %>%
  group_by(Tag, Age_Group, Diet) %>%
  slice_max(order_by = Days.from.exp.start, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(Tag, Age_Group, Diet) %>%
  summarise(across(all_of(markers), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

cat("Post-intervention: using last measurement day per rat (study end)\n")
cat("Last days per rat:\n")
print(df %>% filter(Time == "Post") %>%
        group_by(Tag, Age_Group, Diet) %>%
        summarise(Last_Day = max(Days.from.exp.start), .groups = "drop"))

delta_df <- baseline_df %>%
  inner_join(post_df, by = c("Tag", "Age_Group", "Diet"), suffix = c("_base", "_post"))

for (m in markers) {
  delta_df[[m]] <- delta_df[[paste0(m, "_post")]] - delta_df[[paste0(m, "_base")]]
}
delta_df <- delta_df %>% select(Tag, Age_Group, Diet, all_of(markers))
delta_df <- delta_df[complete.cases(delta_df[, markers]), ]

cat("Delta data:", nrow(delta_df), "rats\n")

# ==============================================================================
# 3. PRIMARY LMM: DIRECT TIME × AGE GROUP × DIET INTERACTION CONTRAST
# ==============================================================================
# For every marker, the following model is fitted:
#   marker ~ Time * Age_Group * Diet + (1 | Tag)
#
# The primary inferential test is the direct difference-in-differences contrast:
#   [(Post - Baseline)_CR - (Post - Baseline)_AL]_Late MA
#   - [(Post - Baseline)_CR - (Post - Baseline)_AL]_Early MA
#
# This is the direct Time × Age_Group × Diet interaction. A positive estimate
# means that the CR-versus-AL change over time is larger in Late MA than in
# Early MA. The 95% CI and unadjusted p-value are obtained from the same LMM;
# Benjamini-Hochberg FDR correction is then performed across all 36 markers.

cat("Running direct Time x Age_Group x Diet interaction analysis...\n")
contrast_definition <- paste0(
  "[(Post-Baseline)_CR - (Post-Baseline)_AL]_Late MA - ",
  "[(Post-Baseline)_CR - (Post-Baseline)_AL]_Early MA"
)
direct_3way_results <- data.frame()

for (marker in markers) {
  formula_str <- paste(marker, "~ Time * Age_Group * Diet + (1 | Tag)")

  tryCatch({
    model <- lmer(as.formula(formula_str), data = df, REML = TRUE)

    # The factor levels are explicitly ordered above. revpairwise yields:
    # Post - Baseline, CR - Ad lib, and Late MA - Early MA.
    emm_full <- emmeans(
      model, ~ Time * Diet * Age_Group,
      lmer.df = "satterthwaite"
    )
    direct_contrast <- contrast(
      emm_full,
      interaction = list(
        Time      = "revpairwise",
        Diet      = "revpairwise",
        Age_Group = "revpairwise"
      )
    )
    direct_summary <- as.data.frame(summary(
      direct_contrast,
      infer  = c(TRUE, TRUE),
      level  = 0.95,
      adjust = "none"
    ))

    # There is exactly one contrast because Time, Diet, and Age_Group each
    # have two levels. The standard error, CI, and p-value are all LMM-derived.
    direct_3way_results <- rbind(direct_3way_results, data.frame(
      Marker   = marker,
      Contrast = contrast_definition,
      estimate = direct_summary$estimate[1],
      SE       = direct_summary$SE[1],
      df       = direct_summary$df[1],
      CI_lower = direct_summary$lower.CL[1],
      CI_upper = direct_summary$upper.CL[1],
      p_value  = direct_summary$p.value[1],
      stringsAsFactors = FALSE
    ))

  }, error = function(e) {
    cat("  LMM error for", marker, ":", conditionMessage(e), "\n")
    direct_3way_results <<- rbind(direct_3way_results, data.frame(
      Marker = marker, Contrast = contrast_definition,
      estimate = NA, SE = NA, df = NA, CI_lower = NA, CI_upper = NA,
      p_value = NA, stringsAsFactors = FALSE
    ))
  })
}

# One BH-FDR family of 36 pre-specified direct interaction tests.
direct_3way_results <- direct_3way_results %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH", n = length(markers)),
    Significance = case_when(
      !is.na(q_value) & q_value < 0.001 ~ "***",
      !is.na(q_value) & q_value < 0.01  ~ "**",
      !is.na(q_value) & q_value < 0.05  ~ "*",
      TRUE                              ~ "ns"
    )
  )

# Compatibility aliases used in the plotting and workbook sections below.
sig_df          <- direct_3way_results
updated_results <- direct_3way_results

write.csv(
  direct_3way_results,
  file.path(out_dir, "Supplementary_Table_Direct_ThreeWay_Interaction.csv"),
  row.names = FALSE
)
cat("Direct interaction table saved to Supplementary_Table_Direct_ThreeWay_Interaction.csv\n\n")

# ==============================================================================
# 3A. HISTORICAL WITHIN-AGE-GROUP ANALYSES (DESCRIPTIVE / EXPLORATORY)
# ==============================================================================
# These results reproduce the two earlier visual summaries:
#   (1) Welch t-test on per-rat delta values, within each age group;
#   (2) within-age-group LMM Time × Diet tests, using raw repeated measures.
# They are retained for transparent presentation of the complete analysis history.
# IMPORTANT: Neither result is a test of differential response between EMA and
# LMA. Only the direct 3-way contrast in Section 3 tests that primary question.

cat("Running historical within-age-group t-test and LMM summaries...\n")

run_delta_ttest <- function(age_value, column_label) {
  result_list <- lapply(markers, function(marker) {
    sub <- delta_df %>% filter(Age_Group == age_value)
    al_vals <- sub %>% filter(Diet == "Ad lib") %>% pull(all_of(marker))
    cr_vals <- sub %>% filter(Diet == "CR") %>% pull(all_of(marker))

    tryCatch({
      tt <- t.test(cr_vals, al_vals, var.equal = FALSE)
      data.frame(
        Marker = marker, Analysis = column_label, Method = "Welch t-test on delta",
        estimate = unname(tt$estimate[1] - tt$estimate[2]),
        SE = NA_real_, df = unname(tt$parameter),
        CI_lower = unname(tt$conf.int[1]), CI_upper = unname(tt$conf.int[2]),
        p_value = tt$p.value, stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        Marker = marker, Analysis = column_label, Method = "Welch t-test on delta",
        estimate = NA_real_, SE = NA_real_, df = NA_real_,
        CI_lower = NA_real_, CI_upper = NA_real_, p_value = NA_real_,
        stringsAsFactors = FALSE
      )
    })
  })

  bind_rows(result_list) %>%
    mutate(
      q_value = p.adjust(p_value, method = "BH", n = length(markers)),
      Significance = case_when(
        !is.na(q_value) & q_value < 0.001 ~ "***",
        !is.na(q_value) & q_value < 0.01  ~ "**",
        !is.na(q_value) & q_value < 0.05  ~ "*",
        TRUE                              ~ "ns"
      )
    )
}

run_within_age_lmm <- function(age_value, column_label) {
  result_list <- lapply(markers, function(marker) {
    sub <- df %>% filter(Age_Group == age_value) %>% droplevels()
    formula_str <- paste(marker, "~ Time * Diet + (1 | Tag)")

    tryCatch({
      model <- lmer(as.formula(formula_str), data = sub, REML = TRUE)
      av <- anova(model, ddf = "Satterthwaite")
      p_td <- av["Time:Diet", "Pr(>F)"]
      data.frame(
        Marker = marker, Analysis = column_label,
        Method = "Within-age LMM Time × Diet",
        p_value = as.numeric(p_td), stringsAsFactors = FALSE
      )
    }, error = function(e) {
      data.frame(
        Marker = marker, Analysis = column_label,
        Method = "Within-age LMM Time × Diet",
        p_value = NA_real_, stringsAsFactors = FALSE
      )
    })
  })

  bind_rows(result_list) %>%
    mutate(
      q_value = p.adjust(p_value, method = "BH", n = length(markers)),
      Significance = case_when(
        !is.na(q_value) & q_value < 0.001 ~ "***",
        !is.na(q_value) & q_value < 0.01  ~ "**",
        !is.na(q_value) & q_value < 0.05  ~ "*",
        TRUE                              ~ "ns"
      )
    )
}

# Historical columns corresponding to the initial submission and first revision.
ttest_ema <- run_delta_ttest("Early middle age", "EMA: delta t-test")
ttest_lma <- run_delta_ttest("Late middle age",  "LMA: delta t-test")
lmm_ema   <- run_within_age_lmm("Early middle age", "EMA: within-age LMM")
lmm_lma   <- run_within_age_lmm("Late middle age",  "LMA: within-age LMM")

historical_within_age_results <- bind_rows(ttest_ema, ttest_lma, lmm_ema, lmm_lma)
write.csv(
  historical_within_age_results,
  file.path(out_dir, "Supplementary_Historical_WithinAge_Results.csv"),
  row.names = FALSE
)
cat("Historical within-age-group results saved to Supplementary_Historical_WithinAge_Results.csv\n\n")

# ==============================================================================
# 3B. EXPLORATORY PERMUTATION SENSITIVITY ANALYSIS ON DELTAS
# ==============================================================================
# This optional sensitivity analysis is reported separately and is not used for
# primary inference, Figure 3 significance labels, or claims of age-dependent CR effects.
cat("Running exploratory permutation sensitivity analysis on deltas...\n")
perm_results <- data.frame()

for (marker in markers) {
  sub_df <- delta_df %>%
    filter(Age_Group %in% c("Early middle age", "Late middle age")) %>%
    filter(!is.na(.data[[marker]])) %>%
    mutate(Age_Group = droplevels(Age_Group),
           Diet      = droplevels(Diet))

  if (nrow(sub_df) >= 8) {
    tryCatch({
      it     <- independence_test(
                  as.formula(paste0("`", marker, "` ~ Diet | Age_Group")),
                  data     = sub_df,
                  teststat = "quadratic"
                )
      p_perm <- as.numeric(pvalue(it))

      perm_results <- rbind(perm_results, data.frame(
        Marker     = marker,
        Comparison = "Permutation: Diet effect diff (Early vs Late MA)",
        p_value    = p_perm,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      perm_results <<- rbind(perm_results, data.frame(
        Marker = marker, Comparison = "Permutation Error",
        p_value = NA, stringsAsFactors = FALSE
      ))
    })
  }
}

if (nrow(perm_results) > 0) {
  perm_results$q_value_FDR <- p.adjust(perm_results$p_value, method = "fdr")
  write.csv(perm_results, file.path(out_dir, "Permutation_Test_Results.csv"), row.names = FALSE)
  cat("Exploratory permutation results saved to Permutation_Test_Results.csv\n\n")
}

# ==============================================================================
# 4. ANALYTIC CLARIFICATION
# ==============================================================================
# The direct Time × Age_Group × Diet contrast defined in Section 3 is the sole
# primary inferential result used for Figure 3, its boxplot labels, and the
# supplementary direct-interaction table. Separate within-age-group p-values
# are not used to claim an age-dependent treatment response.
#
# The permutation analysis above is retained only as an exploratory sensitivity
# analysis on delta values and is reported separately from the primary LMM.
# ==============================================================================

# ==============================================================================
# 7. DELTA-PCA + PERMANOVA
# ==============================================================================
cat("Running Delta-PCA...\n")

delta_matrix <- as.matrix(delta_df[, markers])
pca_res      <- prcomp(delta_matrix, center = TRUE, scale. = TRUE)
var_exp      <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

pca_plot_df <- data.frame(
  PC1       = pca_res$x[, 1],
  PC2       = pca_res$x[, 2],
  Age_Group = delta_df$Age_Group,
  Diet      = delta_df$Diet,
  Tag       = delta_df$Tag
)
# Create combined Age x Diet group label for 4-colour PCA
pca_plot_df$Group <- paste(pca_plot_df$Age_Group, pca_plot_df$Diet, sep = "_")

set.seed(123)
dist_mat      <- dist(scale(delta_matrix), method = "euclidean")
permanova_res <- adonis2(dist_mat ~ Age_Group * Diet, data = delta_df, permutations = 999)
p_permanova   <- permanova_res$`Pr(>F)`[3]
cat("PERMANOVA Age:Diet interaction p =", round(p_permanova, 3), "\n")

# ==============================================================================
# 7B. K-MEANS CLUSTERING
# ==============================================================================
cat("Running K-means clustering...\n")

kmeans_df  <- df[complete.cases(df[, markers]), ]
kmeans_mat <- scale(as.matrix(kmeans_df[, markers]))

set.seed(42)
k_res             <- kmeans(kmeans_mat, centers = 3, nstart = 25)
kmeans_df$Cluster <- as.factor(k_res$cluster)

k_pca             <- prcomp(kmeans_mat, center = FALSE, scale. = FALSE)
kmeans_df$PC1     <- k_pca$x[, 1]
kmeans_df$PC2     <- k_pca$x[, 2]

kmeans_df$Inflection <- factor(
  ifelse(kmeans_df$Age_Group == "Early middle age", "Pre-inflection", "Post-inflection"),
  levels = c("Pre-inflection", "Post-inflection")
)

cluster_summary <- kmeans_df %>%
  group_by(Inflection, Diet, Cluster) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Inflection, Diet) %>%
  mutate(Prop = Count / sum(Count))

p_k_bar <- ggplot(cluster_summary, aes(x = Diet, y = Prop, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  facet_wrap(~ Inflection) +
  scale_fill_manual(values = kmeans_colors) +
  labs(y = "Proportion", x = NULL, title = "K-means Cluster Distribution") +
  theme_bw(base_size = BASE_FONT) +
  theme(legend.position    = "right",
        legend.key.size    = unit(0.3, "cm"),
        strip.text         = element_text(size = BASE_FONT - 2))

p_k_pca <- ggplot(kmeans_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(aes(shape = Time), size = 2, alpha = 0.8) +
  stat_ellipse(type = "norm", level = 0.9, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = kmeans_colors) +
  labs(title = "K-means on PCA (raw values)") +
  theme_bw(base_size = BASE_FONT) +
  theme(legend.key.size = unit(0.3, "cm"))

panel_kmeans <- plot_grid(p_k_bar, p_k_pca, ncol = 1, rel_heights = c(1, 1.5))

# ==============================================================================
# 8. BUILD FIGURE 2 PANELS
# ==============================================================================
cat("Building Figure 2 panels...\n")

# ── Panel a: Delta-PCA ────────────────────────────────────────────────────────
# 4 colours: Early MA AL (skyblue), Early MA CR (darkblue),
#            Late MA AL (orange),   Late MA CR (red)
# Each group gets its own ellipse; legend shows combined Age+Diet label
group_labels <- c(
  "Early middle age_Ad lib" = "Early MA  Al",
  "Early middle age_CR"     = "Early MA  CR",
  "Late middle age_Ad lib"  = "Late MA  Al",
  "Late middle age_CR"      = "Late MA  CR"
)

panel_a <- ggplot(pca_plot_df, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3, alpha = 0.85, stroke = 0.3) +
  stat_ellipse(aes(group = Group),
               type = "norm", level = 0.9,
               linetype = "dashed", linewidth = 0.6) +
  scale_color_manual(values = group_colors,
                     labels = group_labels,
                     name   = "Group") +
  labs(
    title    = "Delta-PCA\n(Post \u2212 Baseline)",
    subtitle = paste0("PERMANOVA Age\u00d7Diet p = ", round(p_permanova, 3)),
    x        = paste0("PC1 (", var_exp[1], "%)"),
    y        = paste0("PC2 (", var_exp[2], "%)") 
  ) +
  theme_bw(base_size = BASE_FONT) +
  theme(
    plot.title      = element_text(face = "bold", size = BASE_FONT + 1),
    plot.subtitle   = element_text(size = BASE_FONT - 1, color = "grey40"),
    legend.key.size = unit(0.3, "cm"),
    legend.text     = element_text(size = BASE_FONT - 2)
  ) +
  annotate("text", x = -Inf, y = -Inf, label = "\u03b1 = 0.9 ellipses",
           hjust = -0.1, vjust = -0.5, size = BASE_FONT / 3, color = "grey50")

# ── Panel b: Hierarchical clustering heatmaps ─────────────────────────────────
# Input: DELTA values (Post - Baseline), normalised per marker to [-1, 1]
# Row labels shown on BOTH heatmaps

pre_mat <- delta_df %>%
  filter(Age_Group == "Early middle age") %>%
  select(Tag, Diet, all_of(markers))

post_mat <- delta_df %>%
  filter(Age_Group == "Late middle age") %>%
  select(Tag, Diet, all_of(markers))

make_pheatmap <- function(group_df, title_str, show_rownames = TRUE) {
  mat     <- t(as.matrix(group_df[, markers]))
  row_max <- apply(abs(mat), 1, max, na.rm = TRUE)
  row_max[row_max == 0] <- 1
  mat_norm           <- mat / row_max
  colnames(mat_norm) <- as.character(group_df$Tag)
  # Apply clean display names to row labels
  rownames(mat_norm) <- to_display(rownames(mat_norm))

  ann_col    <- data.frame(Diet = group_df$Diet,
                           row.names = as.character(group_df$Tag))
  ann_colors <- list(Diet = diet_colors)
  # Original legend: blue (decreased) to grey (no change) to red (increased)
  bwr_colors <- colorRampPalette(c("#2166ac", "#808080", "#d73027"))(100)

  pheatmap(mat_norm,
           color             = bwr_colors,
           breaks            = seq(-1, 1, length.out = 101),
           cluster_rows      = TRUE,
           cluster_cols      = TRUE,
           clustering_method = "ward.D2",
           annotation_col    = ann_col,
           annotation_colors = ann_colors,
           show_colnames     = FALSE,
           show_rownames     = show_rownames,
           fontsize_row      = BASE_FONT - 1,
           fontsize          = BASE_FONT,
           angle_col         = 0,
           main              = title_str,
           border_color      = NA,
           silent            = TRUE)
}

ph_pre  <- make_pheatmap(pre_mat,  "Pre-inflection",  show_rownames = TRUE)
ph_post <- make_pheatmap(post_mat, "Post-inflection", show_rownames = TRUE)

# ── Panel c: Primary direct-interaction significance heatmap ──────────────────
# This one-column heatmap remains the main Figure 3a panel because it displays
# the pre-specified test that answers the central biological question.
max_val <- 3
heat_colors <- colorRampPalette(c("#f0f0f0",   # white/light grey = not significant
                                   "#fee08b",   # yellow = q approximately 0.05
                                   "#f46d43",   # orange = q approximately 0.01
                                   "#a50026"))(100)  # dark red = q below 0.001

primary_sig_mat <- direct_3way_results %>%
  transmute(Marker, neg_log_q = -log10(pmax(q_value, 1e-10)))
primary_sig_mat <- primary_sig_mat[match(markers, primary_sig_mat$Marker), ]
primary_sig_mat <- primary_sig_mat[!is.na(primary_sig_mat$Marker), ]

primary_mat_capped <- as.matrix(primary_sig_mat$neg_log_q)
primary_mat_capped <- pmin(primary_mat_capped, max_val)
rownames(primary_mat_capped) <- to_display(primary_sig_mat$Marker)
colnames(primary_mat_capped) <- "Primary direct\nTime×Age×Diet"

panel_c_pheatmap <- pheatmap(
  primary_mat_capped,
  color         = heat_colors,
  breaks        = seq(0, max_val, length.out = 101),
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row  = BASE_FONT,
  fontsize_col  = BASE_FONT - 2,
  fontsize      = BASE_FONT,
  main          = "Direct Time × Age × Diet interaction\n−log10(FDR q-value)",
  border_color  = "white",
  legend_breaks = c(0, -log10(0.05), -log10(0.01), 3),
  legend_labels = c("ns", "q=0.05", "q=0.01", "q≤0.001"),
  angle_col     = 45,
  silent        = TRUE)

# ── Supplementary five-column heatmap: full analysis history ──────────────────
# Columns 1–4 reproduce historical within-age-group analyses. Column 5 is the
# PRIMARY direct interaction. The first four columns are descriptive/exploratory
# and must not be interpreted as tests of a difference between EMA and LMA.
heatmap_column_order <- c(
  "EMA: delta t-test",
  "LMA: delta t-test",
  "EMA: within-age LMM",
  "LMA: within-age LMM",
  "Primary: direct 3-way interaction"
)
heatmap_column_labels <- c(
  "EMA: delta t-test"                 = "EMA\nΔ t-test",
  "LMA: delta t-test"                 = "LMA\nΔ t-test",
  "EMA: within-age LMM"               = "EMA\nTime×Diet LMM",
  "LMA: within-age LMM"               = "LMA\nTime×Diet LMM",
  "Primary: direct 3-way interaction" = "Primary\nTime×Age×Diet"
)

five_column_results <- bind_rows(
  ttest_ema %>% transmute(Marker, Analysis, q_value),
  ttest_lma %>% transmute(Marker, Analysis, q_value),
  lmm_ema   %>% transmute(Marker, Analysis, q_value),
  lmm_lma   %>% transmute(Marker, Analysis, q_value),
  direct_3way_results %>% transmute(
    Marker, Analysis = "Primary: direct 3-way interaction", q_value
  )
) %>%
  mutate(Analysis = factor(Analysis, levels = heatmap_column_order))

five_column_mat <- five_column_results %>%
  pivot_wider(names_from = Analysis, values_from = q_value) %>%
  select(Marker, all_of(heatmap_column_order))
five_column_mat <- five_column_mat[match(markers, five_column_mat$Marker), ]
five_column_mat <- five_column_mat[!is.na(five_column_mat$Marker), ]

five_column_capped <- as.matrix(five_column_mat[, heatmap_column_order, drop = FALSE])
storage.mode(five_column_capped) <- "numeric"
five_column_capped <- -log10(pmax(five_column_capped, 1e-10))
five_column_capped <- pmin(five_column_capped, max_val)
rownames(five_column_capped) <- to_display(five_column_mat$Marker)
colnames(five_column_capped) <- unname(heatmap_column_labels[heatmap_column_order])

panel_c_5col_pheatmap <- pheatmap(
  five_column_capped,
  color         = heat_colors,
  breaks        = seq(0, max_val, length.out = 101),
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row  = BASE_FONT - 1,
  fontsize_col  = BASE_FONT - 4,
  fontsize      = BASE_FONT,
  main          = "Historical within-age analyses and primary direct interaction\n−log10(FDR q-value)",
  border_color  = "white",
  na_col        = "#D9D9D9",
  legend_breaks = c(0, -log10(0.05), -log10(0.01), 3),
  legend_labels = c("ns", "q=0.05", "q=0.01", "q≤0.001"),
  angle_col     = 45,
  silent        = TRUE)

# ── Panel d: Delta boxplots for top 8 markers ─────────────────────────────────
# Top 8 markers by interaction contrast q-value
top8 <- sig_df %>%
  arrange(q_value) %>%
  slice_head(n = 8) %>%
  pull(Marker)

cat("Top 8 markers for boxplots:", paste(top8, collapse = ", "), "\n")

# Significance annotations from the direct Time × Age_Group × Diet contrast
sig_annot <- sig_df %>%
  mutate(stars = case_when(
    q_value < 0.001 ~ "***",
    q_value < 0.01  ~ "**",
    q_value < 0.05  ~ "*",
    TRUE            ~ "ns"
  ))

# ── make_boxplot: compact design with individual points + boxplot ──────────────
# Layout: x = Inflection group (Pre / Post), fill = Diet (AL / CR)
# This is more space-efficient than faceting — both groups on one x-axis
make_boxplot <- function(marker_name) {
  # Build plot data: one value per rat per group
  plot_df <- delta_df %>%
    filter(!is.na(.data[[marker_name]])) %>%
    mutate(
      Inflection = ifelse(Age_Group == "Early middle age",
                          "Pre", "Post"),
      Inflection = factor(Inflection, levels = c("Pre", "Post")),
      value      = .data[[marker_name]],
      # x-position: group by Inflection x Diet for dodging
      Group      = interaction(Inflection, Diet, sep = "\n")
    )

  # Significance labels from Interaction Contrast
  int_q <- sig_annot %>% filter(Marker == marker_name) %>% pull(stars)
  int_q <- if (length(int_q) == 0) "ns" else int_q[1]

  # Clean display label for y-axis
  display_name <- to_display(marker_name)

  ggplot(plot_df, aes(x = Inflection, y = value, fill = Diet, color = Diet)) +
    geom_boxplot(alpha = 0.55, outlier.shape = NA,
                 position = position_dodge(0.7), width = 0.55, linewidth = 0.35) +
    geom_jitter(aes(group = Diet),
                position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
                size = 1.2, alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.35) +
    # Single significance bracket for the interaction
    labs(
      title = paste0("Direct 3-way interaction: ", int_q),
      y     = paste0("\u0394 ", display_name),
      x     = NULL
    ) +
    scale_fill_manual(values  = diet_colors) +
    scale_color_manual(values = diet_colors) +
    theme_bw(base_size = BASE_FONT) +
    theme(
      legend.position  = "none",
      strip.text       = element_text(size = BASE_FONT - 2),
      axis.text.x      = element_text(size = BASE_FONT - 1),
      axis.title.y     = element_text(size = BASE_FONT - 2),
      plot.margin      = unit(c(4, 4, 2, 4), "pt"),
      plot.title       = element_text(size = BASE_FONT - 2, hjust = 0.5,
                                      color = ifelse(int_q == "ns", "grey50", "darkred"),
                                      fontface = "bold")
    )
}

# Build 8 plots and arrange as 2 columns x 4 rows
bp_plots <- lapply(top8, make_boxplot)

# Add a shared legend from the first plot
legend_plot <- ggplot(
  data.frame(Diet = c("Ad lib", "CR"), x = 1:2, y = 1:2),
  aes(x = x, y = y, fill = Diet, color = Diet)
) +
  geom_boxplot() +
  scale_fill_manual(values = diet_colors, name = "Diet") +
  scale_color_manual(values = diet_colors, name = "Diet") +
  theme_bw(base_size = BASE_FONT) +
  theme(legend.position = "bottom",
        legend.key.size = unit(0.3, "cm"),
        legend.text     = element_text(size = BASE_FONT - 1))
shared_legend <- cowplot::get_legend(legend_plot)

# Left column: top 1-4, Right column: top 5-8
left_col  <- wrap_plots(bp_plots[1:4], ncol = 1)
right_col <- wrap_plots(bp_plots[5:8], ncol = 1)

panel_d <- plot_grid(
  plot_grid(left_col, right_col, ncol = 2, labels = c("Top 1-4", "Top 5-8"),
            label_size = BASE_FONT - 1, label_colour = "grey40"),
  shared_legend,
  ncol = 1, rel_heights = c(1, 0.05)
)

# ==============================================================================
# 9. ASSEMBLE AND SAVE FIGURE AS SINGLE-PAGE PDF
# ==============================================================================
# Output directory was created near the beginning of this script.

wrap_pheatmap <- function(grob) ggdraw() + draw_grob(grob)

panel_b_pre  <- wrap_pheatmap(ph_pre$gtable)
panel_b_post <- wrap_pheatmap(ph_post$gtable)
panel_c_gg      <- wrap_pheatmap(panel_c_pheatmap$gtable)
panel_c_5col_gg <- wrap_pheatmap(panel_c_5col_pheatmap$gtable)

# ── FIGURE 2 ──────────────────────────────────────────────────────────────────
# Top row:    Panel a (Delta-PCA)  |  Panel b (K-means)
# Bottom row: Panel c (Heatmaps: Pre-inflection + Post-inflection side by side)

fig2_top <- plot_grid(
  panel_a,
  panel_kmeans,
  ncol = 2, rel_widths = c(1, 1),
  labels = c("a", "b"), label_size = 14
)

fig2_bottom <- plot_grid(
  panel_b_pre, panel_b_post,
  ncol = 2, rel_widths = c(1, 1),
  labels = c("c", ""), label_size = 14
)

figure2 <- plot_grid(
  fig2_top,
  fig2_bottom,
  nrow = 2, rel_heights = c(1, 1.4)
)

ggsave(
  filename = file.path(out_dir, "Figure2.pdf"),
  plot     = figure2,
  width    = 22, height = 26, limitsize = FALSE,
  device   = "pdf"
)
cat("Figure 2 saved:", file.path(out_dir, "Figure2.pdf"), "\n")

# ── FIGURE 3 ──────────────────────────────────────────────────────────────────
# Panel a (top):    Significance heatmap (formerly panel d)
# Panel b (bottom): Boxplots top 8 markers (formerly panel e)

figure3 <- plot_grid(
  panel_c_gg,
  panel_d,
  nrow = 2, rel_heights = c(1, 1.6),
  labels = c("a", "b"), label_size = 14
)

ggsave(
  filename = file.path(out_dir, "Figure3.pdf"),
  plot     = figure3,
  width    = 18, height = 30, limitsize = FALSE,
  device   = "pdf"
)
cat("Figure 3 saved:", file.path(out_dir, "Figure3.pdf"), "\n")

# Individual panel PDFs
ggsave(file.path(out_dir, "Fig2a_delta_pca.pdf"),       panel_a,      width = 11,  height = 9,  device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Fig2b_kmeans.pdf"),          panel_kmeans, width = 11,  height = 13, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Fig2c_heatmap_pre.pdf"),     panel_b_pre,  width = 11,  height = 16, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Fig2c_heatmap_post.pdf"),    panel_b_post, width = 11,  height = 16, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Fig3a_significance.interaction.pdf"),    panel_c_gg,      width = 14, height = 16, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Fig3a_significance.pdf"),
       panel_c_5col_gg, width = 16, height = 16, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Fig3b_boxplots.pdf"),        panel_d,      width = 14,  height = 20, device = "pdf", limitsize = FALSE)

cat("\nAll files saved to:", out_dir, "\n")

# ==============================================================================
# 10. CLEAN MARKER LABELS IN ALL CSV OUTPUTS + SUPPLEMENTARY EXCEL
# ==============================================================================
cat("Generating supplementary Excel file...\n")

if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
library(openxlsx)

# Helper: apply display names to Marker column of any data frame
clean_markers <- function(df) {
  if ("Marker" %in% colnames(df)) df$Marker <- to_display(df$Marker)
  df
}

# ── Sheet 1: Primary direct Time × Age Group × Diet interaction ───────────────
# This is the complete numerical table requested by the reviewer. Estimate, SE,
# and the unadjusted 95% CI are derived from the LMM; q-value is BH-FDR adjusted
# across the 36 pre-specified direct interaction tests.
direct_interaction_sheet <- direct_3way_results %>%
  transmute(
    Marker          = to_display(Marker),
    Contrast        = Contrast,
    Estimate        = round(estimate, 4),
    `Standard error` = round(SE, 4),
    `Degrees of freedom` = round(df, 2),
    `95% CI lower`  = round(CI_lower, 4),
    `95% CI upper`  = round(CI_upper, 4),
    `Unadjusted p-value` = round(p_value, 4),
    `FDR q-value`   = round(q_value, 4),
    Significance    = Significance
  )

# ── Sheet 2: Historical within-age-group results (descriptive / exploratory) ──
# These are included for transparency and correspond to the first four columns
# of Figure 3a. They are not evidence for a differential EMA-versus-LMA effect.
historical_within_age_sheet <- historical_within_age_results %>%
  mutate(
    Marker = to_display(Marker),
    Interpretation = "Within-age-group result only; not a test of EMA versus LMA difference.",
    across(where(is.numeric), \(x) round(x, 4))
  ) %>%
  select(Marker, Analysis, Method, estimate, SE, df, CI_lower, CI_upper,
         p_value, q_value, Significance, Interpretation)

# ── Sheet 3: Delta values per rat (matching heatmap b) ────────────────────────
delta_sheet <- delta_df %>%
  select(Tag, Age_Group, Diet, all_of(markers)) %>%
  mutate(
    Inflection = ifelse(Age_Group == "Early middle age",
                        "Pre-inflection", "Post-inflection")
  ) %>%
  select(Tag, Age_Group, Inflection, Diet, all_of(markers))
# Rename marker columns to display names
colnames(delta_sheet)[5:ncol(delta_sheet)] <- to_display(markers)

# ── Sheet 4: Top 8 markers summary (matching panel d boxplots) ────────────────
top8_summary <- delta_df %>%
  select(Tag, Age_Group, Diet, all_of(top8)) %>%
  mutate(
    Inflection = ifelse(Age_Group == "Early middle age",
                        "Pre-inflection", "Post-inflection")
  ) %>%
  select(Tag, Age_Group, Inflection, Diet, all_of(top8))
colnames(top8_summary)[5:ncol(top8_summary)] <- to_display(top8)

# ── Sheet 5: Full direct 3-way LMM results (clean labels) ─────────────────────
lmm_full_sheet <- clean_markers(updated_results) %>%
  mutate(across(where(is.numeric), \(x) round(x, 4)))

# ── Sheet 6: Exploratory permutation test results (clean labels) ──────────────
perm_sheet <- clean_markers(perm_results) %>%
  mutate(across(where(is.numeric), \(x) round(x, 4)))

# ── Build workbook with formatting ────────────────────────────────────────────
wb <- createWorkbook()

# Styles
header_style <- createStyle(fontColour = "#FFFFFF", fgFill = "#2c3e50",
                             halign = "CENTER", textDecoration = "Bold",
                             border = "Bottom", borderColour = "#FFFFFF")
sig_style    <- createStyle(fgFill = "#c0392b", fontColour = "#FFFFFF",
                             textDecoration = "Bold", halign = "CENTER")
border_style <- createStyle(fgFill = "#f9e4e4", halign = "CENTER")
ns_style     <- createStyle(fgFill = "#f0f0f0", halign = "CENTER")

add_sheet_formatted <- function(wb, sheet_name, data, sig_cols = NULL) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data, headerStyle = header_style)
  setColWidths(wb, sheet_name, cols = 1:ncol(data), widths = "auto")
  # Freeze top row
  freezePane(wb, sheet_name, firstRow = TRUE)
  # Colour significant cells
  if (!is.null(sig_cols)) {
    for (col_name in sig_cols) {
      col_idx <- which(colnames(data) == col_name)
      if (length(col_idx) == 0) next
      for (row_idx in seq_len(nrow(data))) {
        val <- data[[col_name]][row_idx]
        if (!is.na(val) && val != "ns") {
          addStyle(wb, sheet_name, sig_style,
                   rows = row_idx + 1, cols = col_idx, stack = TRUE)
        } else {
          addStyle(wb, sheet_name, ns_style,
                   rows = row_idx + 1, cols = col_idx, stack = TRUE)
        }
      }
    }
  }
}

add_sheet_formatted(wb, "Direct_3Way_LMM", direct_interaction_sheet,
                    sig_cols = c("Significance"))
add_sheet_formatted(wb, "Historical_WithinAge", historical_within_age_sheet,
                    sig_cols = c("Significance"))
add_sheet_formatted(wb, "Delta_Values",           delta_sheet)
add_sheet_formatted(wb, "Top8_Markers_Boxplots",  top8_summary)
add_sheet_formatted(wb, "LMM_Full_Results",       lmm_full_sheet)
add_sheet_formatted(wb, "Permutation_Tests",      perm_sheet)

# Save Excel
excel_path <- file.path(out_dir, "Supplementary_Statistics.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)
cat("Supplementary Excel saved:", excel_path, "\n")

cat("\nCSV outputs:\n")
cat("  Supplementary_Table_Direct_ThreeWay_Interaction.csv - primary direct LMM contrast\n")
cat("  Supplementary_Historical_WithinAge_Results.csv       - descriptive historical t-test/LMM results\n")
cat("  Supplementary_Figure3_FiveColumn_Heatmap.pdf           - historical t-test/LMM columns plus primary direct interaction\n")
cat("  Permutation_Test_Results.csv                          - exploratory sensitivity analysis\n")
cat("  Supplementary_Statistics.xlsx                          - formatted supplementary table\n")
