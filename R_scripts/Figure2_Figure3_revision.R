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
#                      "scales","cowplot","coin","cluster"))
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
# cluster is a recommended R package used for silhouette-width calculations.
library(cluster)

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

# Reproducible K-means settings. k = 3 is evaluated using the elbow curve,
# average silhouette width, and bootstrap Jaccard stability (Supplementary Fig.).
KMEANS_K              <- 3
KMEANS_NSTART         <- 25
KMEANS_ITER_MAX       <- 100
KMEANS_SEED           <- 42
KMEANS_STABILITY_SEED <- 2026
KMEANS_BOOTSTRAPS     <- 1000

# High-contrast K-means cluster colours — no yellow, all clearly distinguishable
kmeans_colors <- c("1" = "#e41a1c",   # red
                   "2" = "#377eb8",   # blue
                   "3" = "#4daf4a")   # green

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
# 7. SUPPLEMENTARY DELTA-PCA + PERMANOVA
# ==============================================================================
# Each row represents one rat. For each marker, delta is defined as day 105 minus
# baseline. The 36 delta variables are z-standardized before PCA, K-means, and
# Euclidean-distance calculations. Baseline and follow-up are therefore NOT
# treated as separate observations in any multivariate or clustering analysis.
cat("Running Delta-PCA...\n")

delta_matrix <- as.matrix(delta_df[, markers])
delta_scaled <- scale(delta_matrix, center = TRUE, scale = TRUE)

# PCA calculates all possible components. K-means is NOT applied to PCA scores;
# PC1 and PC2 are used only to visualize clusters in two dimensions.
pca_res <- prcomp(delta_scaled, center = FALSE, scale. = FALSE)
var_exp <- 100 * pca_res$sdev^2 / sum(pca_res$sdev^2)
var_exp_display <- round(var_exp, 1)
n_pcs_calculated <- length(pca_res$sdev)

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
dist_mat      <- dist(delta_scaled, method = "euclidean")
permanova_res <- adonis2(dist_mat ~ Age_Group * Diet, data = delta_df, permutations = 999)
p_permanova   <- permanova_res$`Pr(>F)`[3]
pca_variance_table <- data.frame(
  Principal_component = paste0("PC", seq_along(var_exp)),
  Variance_explained_percent = round(var_exp, 3),
  Cumulative_variance_explained_percent = round(cumsum(var_exp), 3)
)
write.csv(pca_variance_table,
          file.path(out_dir, "Supplementary_PCA_VarianceExplained.csv"),
          row.names = FALSE)

cat("PCA components calculated:", n_pcs_calculated,
    "| PC1:", round(var_exp[1], 1), "% | PC2:", round(var_exp[2], 1), "%\n")
cat("PERMANOVA Age:Diet interaction p =", round(p_permanova, 3), "\n")

# ==============================================================================
# 7B. K-MEANS CLUSTERING ON THE STANDARDIZED DELTA MATRIX
# ==============================================================================
# K-means is applied directly to the 36-variable standardized delta matrix, NOT
# to PC scores. This produces one clustering observation per rat (n = nrow(delta_df)).
# The Hartigan-Wong algorithm, 25 random starts, 100 maximum iterations, and a
# fixed seed are specified for reproducibility.
cat("Running delta-based K-means clustering...\n")

kmeans_df  <- delta_df
kmeans_mat <- delta_scaled

set.seed(KMEANS_SEED)
k_res <- kmeans(
  kmeans_mat,
  centers = KMEANS_K,
  nstart = KMEANS_NSTART,
  iter.max = KMEANS_ITER_MAX,
  algorithm = "Hartigan-Wong"
)
kmeans_df$Cluster <- factor(k_res$cluster, levels = seq_len(KMEANS_K))
kmeans_df$PC1 <- pca_res$x[, 1]
kmeans_df$PC2 <- pca_res$x[, 2]
kmeans_df$Group <- paste(kmeans_df$Age_Group, kmeans_df$Diet, sep = "_")
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
  scale_fill_manual(values = kmeans_colors, drop = FALSE) +
  labs(y = "Proportion of rats", x = NULL,
       title = paste0("K-means cluster distribution (k = ", KMEANS_K, ")")) +
  theme_bw(base_size = BASE_FONT) +
  theme(legend.position = "right",
        legend.key.size = unit(0.3, "cm"),
        strip.text = element_text(size = BASE_FONT - 2))

p_k_pca <- ggplot(kmeans_df, aes(x = PC1, y = PC2, color = Cluster)) +
  geom_point(aes(shape = Group), size = 2.8, alpha = 0.85) +
  stat_ellipse(type = "norm", level = 0.9, linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = kmeans_colors, drop = FALSE) +
  labs(
    title = "K-means clusters projected on Delta-PCA",
    subtitle = "K-means fitted to the standardized 36-marker delta matrix; PC1/PC2 shown for visualization",
    x = paste0("PC1 (", var_exp_display[1], "%)"),
    y = paste0("PC2 (", var_exp_display[2], "%)"),
    shape = "Experimental group"
  ) +
  theme_bw(base_size = BASE_FONT) +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = BASE_FONT - 3))

panel_kmeans <- plot_grid(p_k_bar, p_k_pca, ncol = 1, rel_heights = c(1, 1.5))

# ==============================================================================
# 7C. K-MEANS QUALITY AND BOOTSTRAP STABILITY ASSESSMENT
# ==============================================================================
# The elbow curve describes within-cluster compactness across k = 1–8. The
# silhouette width describes separation of each rat from the nearest cluster.
# Stability is assessed with subject-level non-parametric bootstrap resampling.
# In each resample, K-means is re-run and each bootstrap cluster is matched to
# the original cluster using the largest Jaccard similarity of rat membership.

cat("Assessing K-means quality and bootstrap stability...\n")

set.seed(KMEANS_SEED)
k_grid <- 1:8
wcss <- vapply(k_grid, function(k) {
  kmeans(kmeans_mat, centers = k, nstart = KMEANS_NSTART,
         iter.max = KMEANS_ITER_MAX, algorithm = "Hartigan-Wong")$tot.withinss
}, numeric(1))

elbow_df <- data.frame(k = k_grid, WCSS = wcss)

silhouette_obj <- cluster::silhouette(as.integer(kmeans_df$Cluster), dist(kmeans_mat))
silhouette_df <- data.frame(
  Tag = as.character(kmeans_df$Tag),
  Cluster = factor(silhouette_obj[, "cluster"], levels = seq_len(KMEANS_K)),
  Silhouette = silhouette_obj[, "sil_width"]
) %>%
  arrange(Cluster, Silhouette) %>%
  mutate(PlotOrder = row_number())

average_silhouette <- mean(silhouette_df$Silhouette, na.rm = TRUE)
cluster_silhouette <- silhouette_df %>%
  group_by(Cluster) %>%
  summarise(
    n = n(),
    mean_silhouette = mean(Silhouette, na.rm = TRUE),
    .groups = "drop"
  )

# Bootstrap Jaccard stability without treating repeated time points as separate
# observations. B = 1000 resamples provides stable empirical estimates while
# retaining the same K-means settings as the final clustering solution.
bootstrap_kmeans_jaccard <- function(data_mat, reference_clusters, k, B,
                                     seed, nstart, iter.max) {
  set.seed(seed)
  n <- nrow(data_mat)
  reference_clusters <- as.integer(reference_clusters)
  jaccard_values <- matrix(NA_real_, nrow = B, ncol = k)

  for (b in seq_len(B)) {
    sampled_index <- sample.int(n, size = n, replace = TRUE)
    boot_fit <- tryCatch(
      kmeans(data_mat[sampled_index, , drop = FALSE], centers = k,
             nstart = nstart, iter.max = iter.max, algorithm = "Hartigan-Wong"),
      error = function(e) NULL
    )
    if (is.null(boot_fit)) next

    for (reference_cluster in seq_len(k)) {
      original_members <- which(reference_clusters == reference_cluster)
      boot_jaccard <- vapply(seq_len(k), function(bootstrap_cluster) {
        bootstrap_members <- unique(sampled_index[boot_fit$cluster == bootstrap_cluster])
        union_members <- union(original_members, bootstrap_members)
        if (length(union_members) == 0) return(NA_real_)
        length(intersect(original_members, bootstrap_members)) / length(union_members)
      }, numeric(1))
      jaccard_values[b, reference_cluster] <- max(boot_jaccard, na.rm = TRUE)
    }
  }

  data.frame(
    Cluster = factor(seq_len(k), levels = seq_len(k)),
    Mean_Jaccard = colMeans(jaccard_values, na.rm = TRUE),
    Median_Jaccard = apply(jaccard_values, 2, median, na.rm = TRUE),
    Lower_2.5pct = apply(jaccard_values, 2, quantile, probs = 0.025, na.rm = TRUE),
    Upper_97.5pct = apply(jaccard_values, 2, quantile, probs = 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

stability_df <- bootstrap_kmeans_jaccard(
  data_mat = kmeans_mat,
  reference_clusters = k_res$cluster,
  k = KMEANS_K,
  B = KMEANS_BOOTSTRAPS,
  seed = KMEANS_STABILITY_SEED,
  nstart = KMEANS_NSTART,
  iter.max = KMEANS_ITER_MAX
) %>%
  mutate(
    Stability_interpretation = case_when(
      Mean_Jaccard >= 0.75 ~ "Stable",
      Mean_Jaccard >= 0.60 ~ "Moderate stability",
      TRUE                 ~ "Low stability"
    )
  )

explained_between_ss <- 100 * k_res$betweenss / k_res$totss

# Value is intentionally stored as character because this summary table contains
# both numeric metrics and the text algorithm name ("Hartigan-Wong").
metric_row <- function(metric, value, details) {
  data.frame(
    Metric = metric,
    Value = as.character(value),
    Details = details,
    stringsAsFactors = FALSE
  )
}

kmeans_metrics <- bind_rows(
  metric_row("Input observations", nrow(kmeans_mat), "One standardized delta profile per rat"),
  metric_row("Input variables", ncol(kmeans_mat), "36 peripheral blood parameters"),
  metric_row("K-means algorithm", "Hartigan-Wong", "Base R kmeans implementation"),
  metric_row("Chosen number of clusters", KMEANS_K, "Evaluated using elbow, silhouette, and stability metrics"),
  metric_row("Random seed", KMEANS_SEED, "Fixed before final K-means fit"),
  metric_row("Random starts", KMEANS_NSTART, "Independent K-means initializations"),
  metric_row("Maximum iterations", KMEANS_ITER_MAX, "Per random start"),
  metric_row("Between/total SS (%)", round(explained_between_ss, 2), "Variance explained by clustering"),
  metric_row("Average silhouette width", round(average_silhouette, 3), "Higher values indicate better separation"),
  metric_row("Bootstrap resamples", KMEANS_BOOTSTRAPS, "Subject-level non-parametric bootstrap"),
  stability_df %>% transmute(
    Metric = paste0("Mean Jaccard similarity: cluster ", Cluster),
    Value = as.character(round(Mean_Jaccard, 3)),
    Details = Stability_interpretation
  )
)

write.csv(kmeans_metrics,
          file.path(out_dir, "Supplementary_KmeansMetrics.csv"), row.names = FALSE)
write.csv(cluster_silhouette,
          file.path(out_dir, "Supplementary_KmeansSilhouette_ByCluster.csv"), row.names = FALSE)
write.csv(stability_df,
          file.path(out_dir, "Supplementary_KmeansBootstrapStability.csv"), row.names = FALSE)
write.csv(kmeans_df %>% select(Tag, Age_Group, Diet, Cluster, PC1, PC2),
          file.path(out_dir, "Supplementary_KmeansClusterAssignments.csv"), row.names = FALSE)

# Make a dedicated plotting table with simple names. This avoids non-standard
# evaluation issues with long column names inside ggplot's aes() call.
scree_plot_df <- data.frame(
  PC = seq_along(var_exp),
  Variance = as.numeric(var_exp),
  PC_label = paste0("PC", seq_along(var_exp))
)

p_scree <- ggplot(scree_plot_df, aes(x = PC, y = Variance)) +
  geom_line(linewidth = 0.65, colour = "grey40") +
  geom_point(size = 2.4, colour = "black") +
  scale_x_continuous(breaks = scree_plot_df$PC,
                     labels = scree_plot_df$PC_label) +
  labs(title = "Delta-PCA scree plot", x = "Principal component",
       y = "Variance explained (%)") +
  theme_bw(base_size = BASE_FONT) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

p_elbow <- ggplot(elbow_df, aes(x = k, y = WCSS)) +
  geom_line(linewidth = 0.65, colour = "grey40") +
  geom_point(size = 2.6, colour = "black") +
  scale_x_continuous(breaks = k_grid) +
  labs(title = "K-means elbow plot", x = "Number of clusters (k)",
       y = "Within-cluster sum of squares") +
  theme_bw(base_size = BASE_FONT) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_silhouette <- ggplot(silhouette_df,
                       aes(x = PlotOrder, y = Silhouette, fill = Cluster)) +
  geom_col(width = 0.9) +
  geom_hline(yintercept = average_silhouette, linetype = "dashed",
             colour = "#D73027", linewidth = 0.55) +
  scale_fill_manual(values = kmeans_colors, drop = FALSE) +
  labs(title = paste0("Silhouette plot (k = ", KMEANS_K, ")"),
       subtitle = paste0("Average silhouette width = ", round(average_silhouette, 3)),
       x = "Rats, sorted by cluster and silhouette width", y = "Silhouette width",
       fill = "Cluster") +
  theme_bw(base_size = BASE_FONT) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_stability <- ggplot(stability_df,
                      aes(x = Cluster, y = Mean_Jaccard, fill = Cluster)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = Lower_2.5pct, ymax = Upper_97.5pct), width = 0.16) +
  geom_hline(yintercept = 0.75, linetype = "dashed", colour = "#D73027", linewidth = 0.55) +
  annotate("text", x = 0.6, y = 0.77, label = "Stable ≥ 0.75", hjust = 0,
           size = BASE_FONT / 3.5, colour = "#D73027") +
  scale_fill_manual(values = kmeans_colors, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Bootstrap cluster stability", subtitle = paste0(KMEANS_BOOTSTRAPS, " resamples"),
       x = "K-means cluster", y = "Mean Jaccard similarity") +
  theme_bw(base_size = BASE_FONT) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

quality_plot <- (p_elbow + p_silhouette) + plot_annotation(
  title = "K-means Cluster Quality Metrics",
  theme = theme(plot.title = element_text(face = "bold", size = BASE_FONT + 1, hjust = 0.5))
)
stability_plot <- p_stability + plot_annotation(
  title = "K-means Bootstrap Stability",
  subtitle = "Cluster membership stability under subject-level resampling",
  theme = theme(
    plot.title = element_text(face = "bold", size = BASE_FONT + 1, hjust = 0.5),
    plot.subtitle = element_text(size = BASE_FONT - 2, hjust = 0.5)
  )
)

ggsave(file.path(out_dir, "Supplementary_PCA_ScreePlot.pdf"), p_scree,
       width = 11, height = 7.5, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Supplementary_KmeansQuality.pdf"), quality_plot,
       width = 18, height = 8.5, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Supplementary_KmeansStability.pdf"), stability_plot,
       width = 9, height = 7.5, device = "pdf", limitsize = FALSE)

cat("Delta-based K-means quality, stability, and cluster assignment files saved.\n")

# ==============================================================================
# 7D. PARALLEL RAW-MEASUREMENT K-MEANS (COMPARISON / SENSITIVITY ANALYSIS)
# ==============================================================================
# This section intentionally reproduces clustering of the full standardized raw
# measurement table so it can be compared with delta-based K-means. Each row is
# one observation at baseline or follow-up; the same rat therefore contributes
# repeated observations. This raw analysis is supplementary/descriptive only and
# must NOT be used as the primary intervention-response clustering result.

cat("Running parallel raw-measurement K-means comparison...\n")

raw_kmeans_df <- df[complete.cases(df[, markers]), ]
raw_kmeans_mat <- scale(as.matrix(raw_kmeans_df[, markers]), center = TRUE, scale = TRUE)
raw_pca_res <- prcomp(raw_kmeans_mat, center = FALSE, scale. = FALSE)
raw_var_exp <- 100 * raw_pca_res$sdev^2 / sum(raw_pca_res$sdev^2)

set.seed(KMEANS_SEED)
raw_k_res <- kmeans(
  raw_kmeans_mat,
  centers = KMEANS_K,
  nstart = KMEANS_NSTART,
  iter.max = KMEANS_ITER_MAX,
  algorithm = "Hartigan-Wong"
)
raw_kmeans_df$Cluster <- factor(raw_k_res$cluster, levels = seq_len(KMEANS_K))
raw_kmeans_df$PC1 <- raw_pca_res$x[, 1]
raw_kmeans_df$PC2 <- raw_pca_res$x[, 2]
raw_kmeans_df$Inflection <- factor(
  ifelse(raw_kmeans_df$Age_Group == "Early middle age", "Pre-inflection", "Post-inflection"),
  levels = c("Pre-inflection", "Post-inflection")
)

# Main Figure 2A: raw-measurement PCA and K-means presentation.
# PCA is calculated from the raw standardized measurement matrix only for the
# two-dimensional display; K-means itself is fit directly to all 36 variables.
raw_cluster_summary <- raw_kmeans_df %>%
  group_by(Inflection, Diet, Cluster) %>%
  summarise(Count = n(), .groups = "drop") %>%
  group_by(Inflection, Diet) %>%
  mutate(Proportion = Count / sum(Count))

panel_kmeans_raw_bar <- ggplot(raw_cluster_summary,
                               aes(x = Diet, y = Proportion, fill = Cluster)) +
  geom_col(width = 0.6) +
  facet_wrap(~ Inflection) +
  scale_fill_manual(values = kmeans_colors, drop = FALSE, name = "Cluster") +
  labs(title = "K-means Cluster Distribution", y = "Proportion", x = NULL) +
  theme_bw(base_size = BASE_FONT) +
  theme(
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    strip.text = element_text(size = BASE_FONT - 2)
  )

panel_raw_pca_main <- ggplot(raw_kmeans_df,
                             aes(x = PC1, y = PC2, colour = Cluster, shape = Time)) +
  geom_point(size = 2.7, alpha = 0.85) +
  stat_ellipse(aes(group = Cluster), type = "norm", level = 0.9,
               linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = kmeans_colors, drop = FALSE, name = "Cluster") +
  scale_shape_manual(values = c("Baseline" = 16, "Post" = 17), name = "Time") +
  labs(
    title = "PCA of raw peripheral blood measurements",
    subtitle = "Raw measurement matrix; clusters shown for descriptive visualization",
    x = paste0("PC1 (", round(raw_var_exp[1], 1), "% variance)"),
    y = paste0("PC2 (", round(raw_var_exp[2], 1), "% variance)")
  ) +
  theme_bw(base_size = BASE_FONT) +
  theme(
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = BASE_FONT - 3)
  )

# Figure 2A combines the raw-data PCA (left) and K-means distribution (right).
panel_figure2a <- plot_grid(
  panel_raw_pca_main, panel_kmeans_raw_bar,
  ncol = 2, rel_widths = c(1.15, 0.85)
)

# Quality: elbow curve and silhouette width for raw observations.
set.seed(KMEANS_SEED)
raw_wcss <- vapply(k_grid, function(k) {
  kmeans(raw_kmeans_mat, centers = k, nstart = KMEANS_NSTART,
         iter.max = KMEANS_ITER_MAX, algorithm = "Hartigan-Wong")$tot.withinss
}, numeric(1))
raw_elbow_df <- data.frame(k = k_grid, WCSS = raw_wcss)

raw_silhouette_obj <- cluster::silhouette(as.integer(raw_kmeans_df$Cluster),
                                          dist(raw_kmeans_mat))
raw_silhouette_df <- data.frame(
  Tag = as.character(raw_kmeans_df$Tag),
  Time = as.character(raw_kmeans_df$Time),
  Cluster = factor(raw_silhouette_obj[, "cluster"], levels = seq_len(KMEANS_K)),
  Silhouette = raw_silhouette_obj[, "sil_width"]
) %>%
  arrange(Cluster, Silhouette) %>%
  mutate(PlotOrder = row_number())
raw_average_silhouette <- mean(raw_silhouette_df$Silhouette, na.rm = TRUE)
raw_cluster_silhouette <- raw_silhouette_df %>%
  group_by(Cluster) %>%
  summarise(n = n(), mean_silhouette = mean(Silhouette, na.rm = TRUE), .groups = "drop")

# Stability: resample animals as blocks, retaining all measurements belonging
# to a selected rat. This prevents individual baseline and follow-up rows from
# being resampled as if they were independent animals.
bootstrap_raw_subject_jaccard <- function(data_mat, subject_id, reference_clusters,
                                          k, B, seed, nstart, iter.max) {
  set.seed(seed)
  subject_id <- as.character(subject_id)
  subject_levels <- unique(subject_id)
  reference_clusters <- as.integer(reference_clusters)
  jaccard_values <- matrix(NA_real_, nrow = B, ncol = k)

  for (b in seq_len(B)) {
    sampled_subjects <- sample(subject_levels, size = length(subject_levels), replace = TRUE)
    sampled_index <- unlist(lapply(sampled_subjects, function(id) which(subject_id == id)),
                            use.names = FALSE)
    boot_fit <- tryCatch(
      kmeans(data_mat[sampled_index, , drop = FALSE], centers = k,
             nstart = nstart, iter.max = iter.max, algorithm = "Hartigan-Wong"),
      error = function(e) NULL
    )
    if (is.null(boot_fit)) next

    for (reference_cluster in seq_len(k)) {
      original_members <- which(reference_clusters == reference_cluster)
      boot_jaccard <- vapply(seq_len(k), function(bootstrap_cluster) {
        bootstrap_members <- unique(sampled_index[boot_fit$cluster == bootstrap_cluster])
        union_members <- union(original_members, bootstrap_members)
        if (length(union_members) == 0) return(NA_real_)
        length(intersect(original_members, bootstrap_members)) / length(union_members)
      }, numeric(1))
      jaccard_values[b, reference_cluster] <- max(boot_jaccard, na.rm = TRUE)
    }
  }

  data.frame(
    Cluster = factor(seq_len(k), levels = seq_len(k)),
    Mean_Jaccard = colMeans(jaccard_values, na.rm = TRUE),
    Median_Jaccard = apply(jaccard_values, 2, median, na.rm = TRUE),
    Lower_2.5pct = apply(jaccard_values, 2, quantile, probs = 0.025, na.rm = TRUE),
    Upper_97.5pct = apply(jaccard_values, 2, quantile, probs = 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

raw_stability_df <- bootstrap_raw_subject_jaccard(
  data_mat = raw_kmeans_mat,
  subject_id = raw_kmeans_df$Tag,
  reference_clusters = raw_k_res$cluster,
  k = KMEANS_K,
  B = KMEANS_BOOTSTRAPS,
  seed = KMEANS_STABILITY_SEED,
  nstart = KMEANS_NSTART,
  iter.max = KMEANS_ITER_MAX
) %>%
  mutate(
    Stability_interpretation = case_when(
      Mean_Jaccard >= 0.75 ~ "Stable",
      Mean_Jaccard >= 0.60 ~ "Moderate stability",
      TRUE                 ~ "Low stability"
    )
  )
raw_explained_between_ss <- 100 * raw_k_res$betweenss / raw_k_res$totss

# Comparable raw-measurement summary tables.
raw_kmeans_metrics <- bind_rows(
  metric_row("Input observations", nrow(raw_kmeans_mat), "Baseline and follow-up rows; repeated observations per rat"),
  metric_row("Input variables", ncol(raw_kmeans_mat), "36 peripheral blood parameters"),
  metric_row("K-means algorithm", "Hartigan-Wong", "Base R kmeans implementation"),
  metric_row("Chosen number of clusters", KMEANS_K, "Same k as delta analysis for side-by-side comparison"),
  metric_row("Random seed", KMEANS_SEED, "Fixed before final K-means fit"),
  metric_row("Random starts", KMEANS_NSTART, "Independent K-means initializations"),
  metric_row("Maximum iterations", KMEANS_ITER_MAX, "Per random start"),
  metric_row("Between/total SS (%)", round(raw_explained_between_ss, 2), "Variance explained by clustering"),
  metric_row("Average silhouette width", round(raw_average_silhouette, 3), "Higher values indicate better separation"),
  metric_row("Bootstrap resamples", KMEANS_BOOTSTRAPS, "Subject-level non-parametric bootstrap; all measurements within a rat retained together"),
  raw_stability_df %>% transmute(
    Metric = paste0("Mean Jaccard similarity: cluster ", Cluster),
    Value = as.character(round(Mean_Jaccard, 3)),
    Details = Stability_interpretation
  )
)

raw_pca_variance_table <- data.frame(
  Principal_component = paste0("PC", seq_along(raw_var_exp)),
  Variance_explained_percent = round(raw_var_exp, 3),
  Cumulative_variance_explained_percent = round(cumsum(raw_var_exp), 3)
)

write.csv(raw_kmeans_metrics,
          file.path(out_dir, "Supplementary_KmeansRaw_Metrics.csv"), row.names = FALSE)
write.csv(raw_cluster_silhouette,
          file.path(out_dir, "Supplementary_KmeansRaw_Silhouette_ByCluster.csv"), row.names = FALSE)
write.csv(raw_stability_df,
          file.path(out_dir, "Supplementary_KmeansRaw_BootstrapStability.csv"), row.names = FALSE)
write.csv(raw_kmeans_df %>% select(Tag, Time, Age_Group, Diet, Cluster, PC1, PC2),
          file.path(out_dir, "Supplementary_KmeansRaw_ClusterAssignments.csv"), row.names = FALSE)
write.csv(raw_pca_variance_table,
          file.path(out_dir, "Supplementary_KmeansRaw_PCA_VarianceExplained.csv"), row.names = FALSE)

# Raw PCA projection of the raw K-means assignments.
p_raw_kmeans_pca <- ggplot(raw_kmeans_df, aes(x = PC1, y = PC2, colour = Cluster)) +
  geom_point(aes(shape = Time), size = 2.7, alpha = 0.85) +
  stat_ellipse(type = "norm", level = 0.9, linetype = "dashed", linewidth = 0.5) +
  scale_colour_manual(values = kmeans_colors, drop = FALSE) +
  labs(
    title = "Raw-measurement K-means projected on PCA",
    subtitle = "Supplementary comparison only: baseline and follow-up rows are separate observations",
    x = paste0("PC1 (", round(raw_var_exp[1], 1), "%)"),
    y = paste0("PC2 (", round(raw_var_exp[2], 1), "%)")
  ) +
  theme_bw(base_size = BASE_FONT) +
  theme(legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = BASE_FONT - 3))

# Raw-data PCA scree plot. This reports the variance explained by every
# component; PC1 and PC2 values are also printed on the Figure 2A PCA axes.
raw_scree_plot_df <- data.frame(
  PC = seq_along(raw_var_exp),
  Variance = as.numeric(raw_var_exp),
  PC_label = paste0("PC", seq_along(raw_var_exp))
)

p_raw_scree <- ggplot(raw_scree_plot_df, aes(x = PC, y = Variance)) +
  geom_line(linewidth = 0.65, colour = "grey40") +
  geom_point(size = 2.4, colour = "black") +
  scale_x_continuous(breaks = raw_scree_plot_df$PC,
                     labels = raw_scree_plot_df$PC_label) +
  labs(
    title = "Raw-data PCA scree plot",
    subtitle = paste0("PC1 = ", round(raw_var_exp[1], 1),
                      "% variance; PC2 = ", round(raw_var_exp[2], 1), "% variance"),
    x = "Principal component",
    y = "Variance explained (%)"
  ) +
  theme_bw(base_size = BASE_FONT) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = BASE_FONT - 3),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(out_dir, "Supplementary_Fig2_RawPCA_ScreePlot.pdf"),
  plot = p_raw_scree,
  width = 12, height = 8, device = "pdf", limitsize = FALSE
)

p_raw_elbow <- ggplot(raw_elbow_df, aes(x = k, y = WCSS)) +
  geom_line(linewidth = 0.65, colour = "grey40") +
  geom_point(size = 2.6, colour = "black") +
  scale_x_continuous(breaks = k_grid) +
  labs(title = "Raw-measurement K-means elbow plot", x = "Number of clusters (k)",
       y = "Within-cluster sum of squares") +
  theme_bw(base_size = BASE_FONT) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

p_raw_silhouette <- ggplot(raw_silhouette_df,
                           aes(x = PlotOrder, y = Silhouette, fill = Cluster)) +
  geom_col(width = 0.9) +
  geom_hline(yintercept = raw_average_silhouette, linetype = "dashed",
             colour = "#D73027", linewidth = 0.55) +
  scale_fill_manual(values = kmeans_colors, drop = FALSE) +
  labs(title = paste0("Raw-measurement silhouette plot (k = ", KMEANS_K, ")"),
       subtitle = paste0("Average silhouette width = ", round(raw_average_silhouette, 3)),
       x = "Observations, sorted by cluster and silhouette width", y = "Silhouette width",
       fill = "Cluster") +
  theme_bw(base_size = BASE_FONT) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

p_raw_stability <- ggplot(raw_stability_df,
                          aes(x = Cluster, y = Mean_Jaccard, fill = Cluster)) +
  geom_col(width = 0.65) +
  geom_errorbar(aes(ymin = Lower_2.5pct, ymax = Upper_97.5pct), width = 0.16) +
  geom_hline(yintercept = 0.75, linetype = "dashed", colour = "#D73027", linewidth = 0.55) +
  annotate("text", x = 0.6, y = 0.77, label = "Stable ≥ 0.75", hjust = 0,
           size = BASE_FONT / 3.5, colour = "#D73027") +
  scale_fill_manual(values = kmeans_colors, drop = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Raw-measurement bootstrap stability",
       subtitle = paste0(KMEANS_BOOTSTRAPS, " subject-level resamples; all rows per rat retained together"),
       x = "K-means cluster", y = "Mean Jaccard similarity") +
  theme_bw(base_size = BASE_FONT) +
  theme(legend.position = "none", plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

raw_quality_plot <- (p_raw_elbow + p_raw_silhouette) + plot_annotation(
  title = "Raw-measurement K-means Cluster Quality Metrics",
  subtitle = "Supplementary comparison: repeated baseline and follow-up observations are included separately",
  theme = theme(
    plot.title = element_text(face = "bold", size = BASE_FONT + 1, hjust = 0.5),
    plot.subtitle = element_text(size = BASE_FONT - 2, hjust = 0.5)
  )
)

# Side-by-side numerical summary of the two input choices.
kmeans_input_comparison <- bind_rows(
  data.frame(
    Analysis = "Per-rat delta matrix (primary exploratory clustering)",
    Observation_unit = "One day-105 minus baseline profile per rat",
    N_observations = nrow(kmeans_mat),
    Average_silhouette = round(average_silhouette, 3),
    Between_total_SS_percent = round(explained_between_ss, 2),
    Mean_Jaccard_across_clusters = round(mean(stability_df$Mean_Jaccard, na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Analysis = "Raw measurement table (supplementary sensitivity comparison)",
    Observation_unit = "Baseline and follow-up rows; repeated observations per rat",
    N_observations = nrow(raw_kmeans_mat),
    Average_silhouette = round(raw_average_silhouette, 3),
    Between_total_SS_percent = round(raw_explained_between_ss, 2),
    Mean_Jaccard_across_clusters = round(mean(raw_stability_df$Mean_Jaccard, na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  )
)
write.csv(kmeans_input_comparison,
          file.path(out_dir, "Supplementary_KmeansRaw_vs_Delta_Comparison.csv"),
          row.names = FALSE)

comparison_long <- kmeans_input_comparison %>%
  select(Analysis, Average_silhouette, Between_total_SS_percent, Mean_Jaccard_across_clusters) %>%
  pivot_longer(-Analysis, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = recode(
    Metric,
    Average_silhouette = "Average silhouette width",
    Between_total_SS_percent = "Between/total SS (%)",
    Mean_Jaccard_across_clusters = "Mean Jaccard stability"
  ))

p_kmeans_comparison <- ggplot(comparison_long, aes(x = Analysis, y = Value, fill = Analysis)) +
  geom_col(width = 0.65) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(
    "Per-rat delta matrix (primary exploratory clustering)" = "#377eb8",
    "Raw measurement table (supplementary sensitivity comparison)" = "#999999"
  )) +
  labs(title = "K-means comparison: raw measurements versus per-rat delta profiles",
       subtitle = "The two analyses use different observation units and should not be treated as equivalent",
       x = NULL, y = NULL) +
  theme_bw(base_size = BASE_FONT) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = BASE_FONT - 2))

ggsave(file.path(out_dir, "Supplementary_KmeansRaw_PCAProjection.pdf"),
       p_raw_kmeans_pca, width = 11, height = 9, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Supplementary_KmeansRaw_Quality.pdf"),
       raw_quality_plot, width = 18, height = 8.5, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Supplementary_KmeansRaw_Stability.pdf"),
       p_raw_stability, width = 9, height = 7.5, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Supplementary_KmeansRaw_vs_Delta_Comparison.pdf"),
       p_kmeans_comparison, width = 18, height = 8, device = "pdf", limitsize = FALSE)

cat("Raw-measurement K-means comparison files saved.\n")

# ==============================================================================
# 7E. LEGACY REPRODUCTION: UNSCALED PER-RAT DELTA K-MEANS QUALITY FIGURE
# ==============================================================================
# This section reproduces the design of the submitted K-means quality figure.
# It clusters the one-row-per-rat delta matrix in ORIGINAL measurement units,
# not standardized values and not PCA scores. It is retained to reproduce the
# historical figure only; the standardized delta clustering above is preferable
# for comparing markers with different units.

cat("Reproducing legacy unscaled-delta K-means quality figure...\n")

LEGACY_KMEANS_K <- 3
LEGACY_KMEANS_SEED <- 42
LEGACY_KMEANS_NSTART <- 25
LEGACY_KMEANS_ITER_MAX <- 100
legacy_kmeans_colors <- c("1" = "red", "2" = "blue", "3" = "green")
legacy_kmeans_mat <- as.matrix(delta_df[, markers])

set.seed(LEGACY_KMEANS_SEED)
legacy_wcss <- vapply(1:8, function(k) {
  kmeans(legacy_kmeans_mat, centers = k, nstart = LEGACY_KMEANS_NSTART,
         iter.max = LEGACY_KMEANS_ITER_MAX, algorithm = "Hartigan-Wong")$tot.withinss
}, numeric(1))
legacy_elbow_df <- data.frame(k = 1:8, WCSS = legacy_wcss)

set.seed(LEGACY_KMEANS_SEED)
legacy_k_res <- kmeans(legacy_kmeans_mat, centers = LEGACY_KMEANS_K,
                       nstart = LEGACY_KMEANS_NSTART,
                       iter.max = LEGACY_KMEANS_ITER_MAX,
                       algorithm = "Hartigan-Wong")
legacy_silhouette <- cluster::silhouette(legacy_k_res$cluster, dist(legacy_kmeans_mat))
legacy_silhouette_df <- data.frame(
  Tag = as.character(delta_df$Tag),
  Cluster = factor(legacy_silhouette[, "cluster"], levels = seq_len(LEGACY_KMEANS_K)),
  Silhouette = legacy_silhouette[, "sil_width"]
) %>%
  arrange(Cluster, Silhouette) %>%
  mutate(PlotOrder = row_number())
legacy_average_silhouette <- mean(legacy_silhouette_df$Silhouette, na.rm = TRUE)
legacy_cluster_silhouette <- legacy_silhouette_df %>%
  group_by(Cluster) %>%
  summarise(n = n(), mean_silhouette = mean(Silhouette, na.rm = TRUE), .groups = "drop")
legacy_explained_between_ss <- 100 * legacy_k_res$betweenss / legacy_k_res$totss

legacy_metrics <- bind_rows(
  metric_row("Input observations", nrow(legacy_kmeans_mat), "One day-105 minus baseline profile per rat"),
  metric_row("Input variables", ncol(legacy_kmeans_mat), "36 peripheral blood parameters"),
  metric_row("Input scaling", "None", "Original measurement units retained"),
  metric_row("K-means input", "Per-rat delta matrix", "K-means was not applied to PCA scores"),
  metric_row("K-means algorithm", "Hartigan-Wong", "Base R kmeans implementation"),
  metric_row("Number of clusters", LEGACY_KMEANS_K, "Historical/reproduced figure"),
  metric_row("Random seed", LEGACY_KMEANS_SEED, "Fixed before final K-means fit"),
  metric_row("Random starts", LEGACY_KMEANS_NSTART, "Independent K-means initializations"),
  metric_row("Maximum iterations", LEGACY_KMEANS_ITER_MAX, "Per random start"),
  metric_row("Between/total SS (%)", round(legacy_explained_between_ss, 2), "Variance explained by clustering"),
  metric_row("Average silhouette width", round(legacy_average_silhouette, 3), "Higher values indicate better separation")
)

p_legacy_elbow <- ggplot(legacy_elbow_df, aes(x = k, y = WCSS)) +
  geom_line(linewidth = 0.65, colour = "grey40") +
  geom_point(size = 2.7, colour = "black") +
  scale_x_continuous(breaks = 1:8) +
  labs(title = "K-means Elbow Plot", x = "Number of clusters (k)",
       y = "Within-cluster sum of squares") +
  theme_bw(base_size = BASE_FONT) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_line(colour = "grey92", linewidth = 0.25))

p_legacy_silhouette <- ggplot(legacy_silhouette_df,
                              aes(x = PlotOrder, y = Silhouette, fill = Cluster)) +
  geom_col(width = 0.9) +
  geom_hline(yintercept = legacy_average_silhouette, linetype = "dashed",
             colour = "#D73027", linewidth = 0.55) +
  scale_fill_manual(values = legacy_kmeans_colors, drop = FALSE, name = "Cluster") +
  labs(title = paste0("Silhouette Plot (k = ", LEGACY_KMEANS_K, ")"),
       subtitle = paste0("Average silhouette width = ", round(legacy_average_silhouette, 3)),
       x = "Rats (sorted by cluster)", y = "Silhouette width") +
  theme_bw(base_size = BASE_FONT) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.minor = element_line(colour = "grey92", linewidth = 0.25))

legacy_quality_core <- p_legacy_elbow + p_legacy_silhouette + plot_annotation(
  title = "K-means Cluster Quality Metrics",
  theme = theme(plot.title = element_text(face = "bold", size = BASE_FONT, hjust = 0.5))
)
legacy_quality_plot <- cowplot::ggdraw() +
  cowplot::draw_label("Figure2. Supplementary Figure2", x = 0.01, y = 0.99,
                      hjust = 0, vjust = 1, size = BASE_FONT / 1.8) +
  cowplot::draw_plot(legacy_quality_core, x = 0, y = 0, width = 1, height = 0.95)

ggsave(file.path(out_dir, "Figure2_SupplementaryFigure2_LegacyReproduction.pdf"),
       legacy_quality_plot, width = 18, height = 10, device = "pdf", limitsize = FALSE)
write.csv(legacy_metrics,
          file.path(out_dir, "Supplementary_KmeansLegacyUnscaledDelta_Metrics.csv"), row.names = FALSE)
write.csv(legacy_cluster_silhouette,
          file.path(out_dir, "Supplementary_KmeansLegacyUnscaledDelta_Silhouette_ByCluster.csv"), row.names = FALSE)
write.csv(data.frame(Tag = as.character(delta_df$Tag), Age_Group = delta_df$Age_Group,
                     Diet = delta_df$Diet, Cluster = legacy_k_res$cluster,
                     stringsAsFactors = FALSE),
          file.path(out_dir, "Supplementary_KmeansLegacyUnscaledDelta_Assignments.csv"), row.names = FALSE)

cat("Legacy unscaled-delta K-means reproduction files saved.\n")

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
    x        = paste0("PC1 (", var_exp_display[1], "%)"),
    y        = paste0("PC2 (", var_exp_display[2], "%)") 
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

# ── Panel c: Significance Heatmap (Interaction Contrast) ──────────────────────
# 1 column: Early MA vs Late MA Interaction
# Colour: -log10(q) scale so gradient is spread across meaningful range:
#   q=1.0  → -log10 = 0.0  (grey, not significant)
#   q=0.05 → -log10 = 1.3  (light colour)
#   q=0.01 → -log10 = 2.0  (medium)
#   q=0.001→ -log10 = 3.0  (dark, highly significant)
# Colour scheme: white (not sig) → yellow → orange → red (highly sig)

# Build a single-column matrix
sig_mat <- sig_df %>%
  mutate(neg_log_q = -log10(pmax(q_value, 1e-10))) %>%
  select(Marker, neg_log_q)

# Fix marker order to match original figure
sig_mat <- sig_mat[match(markers, sig_mat$Marker), ]
sig_mat <- sig_mat[!is.na(sig_mat$Marker), ]

mat_capped <- as.matrix(sig_mat$neg_log_q)
rownames(mat_capped) <- to_display(sig_mat$Marker)
colnames(mat_capped) <- c("Direct Time × Age × Diet\ninteraction")

max_val <- 3
mat_capped <- pmin(mat_capped, max_val)

heat_colors <- colorRampPalette(c("#f0f0f0",   # white/light grey = not sig
                                   "#fee08b",   # yellow = p~0.05
                                   "#f46d43",   # orange = p~0.01
                                   "#a50026"))(100)  # dark red = p<0.001

panel_c_pheatmap <- pheatmap(mat_capped,
  color         = heat_colors,
  breaks        = seq(0, max_val, length.out = 101),
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row  = BASE_FONT,
  fontsize_col  = BASE_FONT,
  fontsize      = BASE_FONT,
  main          = "Direct Time × Age × Diet interaction\n-log10(q)",
  border_color  = "white",
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
panel_c_gg   <- wrap_pheatmap(panel_c_pheatmap$gtable)

# ── FIGURE 2: PRESERVED ORIGINAL STRUCTURE ────────────────────────────────────
# Panel a: raw-data PCA (left; PC1/PC2 variance labelled) and raw-data K-means
#          cluster distribution (right).
# Panel b: pre-inflection and post-inflection heatmaps.
# Delta-PCA is saved separately as a supplementary figure below.

fig2_panel_a <- panel_figure2a

fig2_panel_b <- plot_grid(
  panel_b_pre, panel_b_post,
  ncol = 2, rel_widths = c(1, 1)
)

figure2 <- plot_grid(
  fig2_panel_a,
  fig2_panel_b,
  nrow = 2, rel_heights = c(1, 1.2),
  labels = c("a", "b"), label_size = 14
)

ggsave(
  filename = file.path(out_dir, "Figure2.pdf"),
  plot     = figure2,
  width    = 22, height = 22, limitsize = FALSE,
  device   = "pdf"
)
# Save individual Figure 2A components for transparent review.
ggsave(
  filename = file.path(out_dir, "Fig2a_raw_PCA.pdf"),
  plot     = panel_raw_pca_main,
  width    = 12, height = 10, limitsize = FALSE,
  device   = "pdf"
)
ggsave(
  filename = file.path(out_dir, "Fig2a_cluster_distribution.pdf"),
  plot     = panel_kmeans_raw_bar,
  width    = 10, height = 8, limitsize = FALSE,
  device   = "pdf"
)
# All delta-based outputs are supplementary analyses.
ggsave(
  filename = file.path(out_dir, "Supplementary_Fig2_DeltaPCA.pdf"),
  plot     = panel_a,
  width    = 12, height = 10, limitsize = FALSE,
  device   = "pdf"
)
ggsave(
  filename = file.path(out_dir, "Supplementary_KmeansDelta_ClusterDistribution.pdf"),
  plot     = p_k_bar,
  width    = 11, height = 7, limitsize = FALSE,
  device   = "pdf"
)
ggsave(
  filename = file.path(out_dir, "Supplementary_KmeansDelta_PCAProjection.pdf"),
  plot     = p_k_pca,
  width    = 11, height = 9, limitsize = FALSE,
  device   = "pdf"
)
cat("Figure 2 saved with raw-data PCA + cluster distribution as panel a and heatmaps as panel b. Delta-PCA and delta K-means outputs saved as supplementary figures.\n")

# ── FIGURE 3: ORIGINAL WITHIN-AGE-GROUP PRESENTATION ─────────────────────────
# This reproduces the original Figure 3 structure: four historical within-age
# comparison columns and delta box-dot plots. It is retained as a descriptive
# main figure, while the direct Time × Age × Diet result below remains the
# supplementary formal comparison of differential CR response between ages.

cat("Generating original within-age-group Figure 3...\n")

within_age_delta_test <- function(marker, age_level) {
  dat <- delta_df %>%
    filter(Age_Group == age_level, !is.na(.data[[marker]])) %>%
    mutate(Diet = droplevels(Diet))

  test_out <- tryCatch(
    t.test(as.formula(paste0("`", marker, "` ~ Diet")), data = dat),
    error = function(e) NULL
  )

  data.frame(
    Marker = marker,
    Age_Group = age_level,
    Analysis = "Delta t-test",
    estimate = if (is.null(test_out)) NA_real_ else unname(test_out$estimate[2] - test_out$estimate[1]),
    p_value = if (is.null(test_out)) NA_real_ else test_out$p.value,
    stringsAsFactors = FALSE
  )
}

within_age_lmm_test <- function(marker, age_level) {
  dat <- df %>%
    filter(Age_Group == age_level, !is.na(.data[[marker]])) %>%
    mutate(Time = droplevels(Time), Diet = droplevels(Diet))

  test_out <- tryCatch({
    fit <- lmer(
      as.formula(paste0("`", marker, "` ~ Time * Diet + (1 | Tag)")),
      data = dat, REML = TRUE
    )
    emm <- emmeans(fit, ~ Time * Diet, lmer.df = "satterthwaite")
    as.data.frame(summary(
      contrast(emm, interaction = list(Time = "revpairwise", Diet = "revpairwise")),
      infer = c(TRUE, TRUE), adjust = "none"
    ))
  }, error = function(e) NULL)

  data.frame(
    Marker = marker,
    Age_Group = age_level,
    Analysis = "Time × Diet LMM",
    estimate = if (is.null(test_out)) NA_real_ else test_out$estimate[1],
    p_value = if (is.null(test_out)) NA_real_ else test_out$p.value[1],
    stringsAsFactors = FALSE
  )
}

age_levels_fig3 <- c("Early middle age", "Late middle age")
historical_figure3_results <- bind_rows(lapply(markers, function(marker) {
  bind_rows(lapply(age_levels_fig3, function(age_level) {
    bind_rows(
      within_age_delta_test(marker, age_level),
      within_age_lmm_test(marker, age_level)
    )
  }))
})) %>%
  group_by(Age_Group, Analysis) %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH", n = length(markers)),
    Significance = case_when(
      !is.na(q_value) & q_value < 0.001 ~ "***",
      !is.na(q_value) & q_value < 0.01 ~ "**",
      !is.na(q_value) & q_value < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  ungroup() %>%
  mutate(
    Age_Label = ifelse(Age_Group == "Early middle age", "Pre", "Post"),
    Column = paste(Age_Label, ifelse(Analysis == "Delta t-test", "t-test", "LMM"))
  )

write.csv(
  historical_figure3_results,
  file.path(out_dir, "Supplementary_Fig3_OriginalWithinAge_Results.csv"),
  row.names = FALSE
)

historical_columns <- c("Pre t-test", "Post t-test", "Pre LMM", "Post LMM")
historical_heatmap_df <- historical_figure3_results %>%
  select(Marker, Column, q_value) %>%
  tidyr::pivot_wider(names_from = Column, values_from = q_value) %>%
  select(Marker, all_of(historical_columns))

historical_heatmap_matrix <- as.matrix(historical_heatmap_df[, historical_columns])
rownames(historical_heatmap_matrix) <- to_display(historical_heatmap_df$Marker)
historical_heatmap_values <- pmin(-log10(historical_heatmap_matrix), 3)
historical_heatmap_values[is.infinite(historical_heatmap_values)] <- 3

original_fig3_heatmap <- pheatmap(
  historical_heatmap_values,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#f0f0f0", "#fee08b", "#f46d43", "#a50026"))(100),
  breaks = seq(0, 3, length.out = 101),
  main = "CR effect\n−log10(FDR q-value)",
  fontsize = BASE_FONT - 6,
  fontsize_row = BASE_FONT - 7,
  fontsize_col = BASE_FONT - 6,
  border_color = "white",
  legend_breaks = c(0, -log10(0.05), -log10(0.01), 3),
  legend_labels = c("ns", "q=0.05", "q=0.01", "q≤0.001"),
  silent = TRUE
)
panel_original_fig3_heatmap <- wrap_pheatmap(original_fig3_heatmap$gtable)

original_top8 <- historical_figure3_results %>%
  group_by(Marker) %>%
  summarise(min_q = min(q_value, na.rm = TRUE), .groups = "drop") %>%
  arrange(min_q) %>%
  slice_head(n = 8) %>%
  pull(Marker)

make_original_fig3_boxplot <- function(marker) {
  plot_df <- delta_df %>%
    filter(!is.na(.data[[marker]])) %>%
    mutate(
      Age_Label = factor(
        ifelse(Age_Group == "Early middle age", "Pre", "Post"),
        levels = c("Pre", "Post")
      ),
      Value = .data[[marker]]
    )

  annotation_df <- historical_figure3_results %>%
    filter(Marker == marker, Analysis == "Delta t-test") %>%
    transmute(
      Age_Label = factor(ifelse(Age_Group == "Early middle age", "Pre", "Post"),
                         levels = c("Pre", "Post")),
      label = Significance
    )

  y_max <- max(plot_df$Value, na.rm = TRUE)
  y_min <- min(plot_df$Value, na.rm = TRUE)
  annotation_df$y <- y_max + max(0.08 * (y_max - y_min), 0.05)

  ggplot(plot_df, aes(x = Age_Label, y = Value, fill = Diet, colour = Diet)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.35) +
    geom_boxplot(position = position_dodge(width = 0.72), width = 0.6,
                 alpha = 0.58, outlier.shape = NA, linewidth = 0.35) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.11, dodge.width = 0.72),
      size = 1.5, alpha = 0.85
    ) +
    geom_text(data = annotation_df, aes(x = Age_Label, y = y, label = label),
              inherit.aes = FALSE, fontface = "bold", size = BASE_FONT / 4.8) +
    scale_fill_manual(values = diet_colors) +
    scale_colour_manual(values = diet_colors) +
    labs(
      title = to_display(marker),
      subtitle = "Within-age delta t-test (exploratory)",
      x = NULL,
      y = paste0("Δ ", to_display(marker))
    ) +
    theme_bw(base_size = BASE_FONT - 4) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = BASE_FONT - 2),
      plot.subtitle = element_text(size = BASE_FONT - 6, colour = "grey45"),
      axis.text.x = element_text(face = "bold")
    )
}

original_boxplot_grid <- wrap_plots(lapply(original_top8, make_original_fig3_boxplot), ncol = 2)
original_legend_plot <- ggplot(
  data.frame(Diet = factor(c("Ad lib", "CR"), levels = c("Ad lib", "CR")), x = 1:2, y = 1:2),
  aes(x = x, y = y, fill = Diet, colour = Diet)
) +
  geom_point(size = 3) +
  scale_fill_manual(values = diet_colors, name = "Diet") +
  scale_colour_manual(values = diet_colors, name = "Diet") +
  theme_void() + theme(legend.position = "bottom")
original_shared_legend <- cowplot::get_legend(original_legend_plot)
original_boxplots_with_legend <- plot_grid(
  original_boxplot_grid, original_shared_legend,
  ncol = 1, rel_heights = c(1, 0.05)
)

figure3_original <- plot_grid(
  panel_original_fig3_heatmap,
  original_boxplots_with_legend,
  ncol = 2, rel_widths = c(1, 1.6),
  labels = c("a", "b"), label_size = 14
)

ggsave(
  filename = file.path(out_dir, "Figure3.pdf"),
  plot = figure3_original,
  width = 25, height = 20, limitsize = FALSE, device = "pdf"
)
ggsave(
  filename = file.path(out_dir, "Fig3a_original_within_age_heatmap.pdf"),
  plot = panel_original_fig3_heatmap,
  width = 13, height = 17, limitsize = FALSE, device = "pdf"
)
ggsave(
  filename = file.path(out_dir, "Fig3b_original_within_age_boxplots.pdf"),
  plot = original_boxplots_with_legend,
  width = 17, height = 17, limitsize = FALSE, device = "pdf"
)
cat("Original Figure 3 saved with within-age t-test/LMM heatmap and delta boxplots.\n")

# ── SUPPLEMENTARY FIGURE: DIRECT INTERACTION ANALYSIS ─────────────────────────
# The revised direct Time × Age Group × Diet analysis is reported as a
# supplementary exploratory analysis. The submitted original Figure 3 remains
# the main figure and is deliberately not overwritten by this script.

supplementary_figure3_direct_interaction <- plot_grid(
  panel_c_gg,
  panel_d,
  nrow = 2, rel_heights = c(1, 1.6),
  labels = c("a", "b"), label_size = 14
)

ggsave(
  filename = file.path(out_dir, "Supplementary_Figure3_DirectInteraction.pdf"),
  plot     = supplementary_figure3_direct_interaction,
  width    = 18, height = 30, limitsize = FALSE,
  device   = "pdf"
)
cat("Supplementary direct-interaction Figure 3 saved:",
    file.path(out_dir, "Supplementary_Figure3_DirectInteraction.pdf"), "\n")

# Individual supplementary/main-panel PDFs
ggsave(file.path(out_dir, "Fig2c_heatmap_pre.pdf"),     panel_b_pre,      width = 11, height = 16, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Fig2c_heatmap_post.pdf"),    panel_b_post,     width = 11, height = 16, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Supplementary_Fig3a_DirectInteraction_Heatmap.pdf"),
       panel_c_gg, width = 14, height = 16, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "Supplementary_Fig3b_DirectInteraction_Boxplots.pdf"),
       panel_d, width = 14, height = 20, device = "pdf", limitsize = FALSE)

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

# ── Sheet 2: PCA variance explained ──────────────────────────────────────────
pca_variance_sheet <- pca_variance_table

# ── Sheet 3: K-means configuration and quality metrics ────────────────────────
kmeans_metrics_sheet <- kmeans_metrics

# ── Sheet 4: Per-cluster silhouette widths ────────────────────────────────────
kmeans_silhouette_sheet <- cluster_silhouette %>%
  mutate(across(where(is.numeric), \(x) round(x, 4)))

# ── Sheet 5: Bootstrap stability metrics ──────────────────────────────────────
kmeans_stability_sheet <- stability_df %>%
  mutate(across(where(is.numeric), \(x) round(x, 4)))

# ── Sheet 6: K-means cluster assignments ──────────────────────────────────────
kmeans_assignments_sheet <- kmeans_df %>%
  select(Tag, Age_Group, Diet, Cluster, PC1, PC2) %>%
  mutate(PC1 = round(PC1, 4), PC2 = round(PC2, 4))

# ── Sheets 7–11: Raw-measurement K-means sensitivity comparison ──────────────
# These outputs retain baseline/follow-up rows separately and are therefore
# supplementary/descriptive. They are provided solely to compare input choices.
raw_kmeans_metrics_sheet <- raw_kmeans_metrics
raw_kmeans_silhouette_sheet <- raw_cluster_silhouette %>%
  mutate(across(where(is.numeric), \(x) round(x, 4)))
raw_kmeans_stability_sheet <- raw_stability_df %>%
  mutate(across(where(is.numeric), \(x) round(x, 4)))
raw_kmeans_assignments_sheet <- raw_kmeans_df %>%
  select(Tag, Time, Age_Group, Diet, Cluster, PC1, PC2) %>%
  mutate(PC1 = round(PC1, 4), PC2 = round(PC2, 4))
kmeans_comparison_sheet <- kmeans_input_comparison

# ── Sheets 12–14: Historical figure reproduction (unscaled per-rat deltas) ───
legacy_kmeans_metrics_sheet <- legacy_metrics
legacy_kmeans_silhouette_sheet <- legacy_cluster_silhouette %>%
  mutate(across(where(is.numeric), \(x) round(x, 4)))
legacy_kmeans_assignments_sheet <- data.frame(
  Tag = as.character(delta_df$Tag),
  Age_Group = delta_df$Age_Group,
  Diet = delta_df$Diet,
  Cluster = legacy_k_res$cluster,
  stringsAsFactors = FALSE
)

# ── Sheet 15: Delta values per rat (matching heatmap b) ───────────────────────
delta_sheet <- delta_df %>%
  select(Tag, Age_Group, Diet, all_of(markers)) %>%
  mutate(
    Inflection = ifelse(Age_Group == "Early middle age",
                        "Pre-inflection", "Post-inflection")
  ) %>%
  select(Tag, Age_Group, Inflection, Diet, all_of(markers))
# Rename marker columns to display names
colnames(delta_sheet)[5:ncol(delta_sheet)] <- to_display(markers)

# ── Sheet 16: Top 8 markers summary (matching panel d boxplots) ───────────────
top8_summary <- delta_df %>%
  select(Tag, Age_Group, Diet, all_of(top8)) %>%
  mutate(
    Inflection = ifelse(Age_Group == "Early middle age",
                        "Pre-inflection", "Post-inflection")
  ) %>%
  select(Tag, Age_Group, Inflection, Diet, all_of(top8))
colnames(top8_summary)[5:ncol(top8_summary)] <- to_display(top8)

# ── Sheet 17: Full LMM results (updated, clean labels) ────────────────────────
lmm_full_sheet <- clean_markers(updated_results) %>%
  mutate(across(where(is.numeric), \(x) round(x, 4)))

# ── Sheet 18: Permutation test results (clean labels) ─────────────────────────
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
add_sheet_formatted(wb, "PCA_Variance",        pca_variance_sheet)
add_sheet_formatted(wb, "Kmeans_Metrics_Delta", kmeans_metrics_sheet)
add_sheet_formatted(wb, "Kmeans_Silhouette_Delta", kmeans_silhouette_sheet)
add_sheet_formatted(wb, "Kmeans_Stability_Delta", kmeans_stability_sheet)
add_sheet_formatted(wb, "Kmeans_Assignments_Delta", kmeans_assignments_sheet)
add_sheet_formatted(wb, "Kmeans_Metrics_Raw", raw_kmeans_metrics_sheet)
add_sheet_formatted(wb, "Kmeans_Silhouette_Raw", raw_kmeans_silhouette_sheet)
add_sheet_formatted(wb, "Kmeans_Stability_Raw", raw_kmeans_stability_sheet)
add_sheet_formatted(wb, "Kmeans_Assignments_Raw", raw_kmeans_assignments_sheet)
add_sheet_formatted(wb, "Kmeans_Raw_vs_Delta", kmeans_comparison_sheet)
add_sheet_formatted(wb, "Kmeans_Legacy_Metrics", legacy_kmeans_metrics_sheet)
add_sheet_formatted(wb, "Kmeans_Legacy_Silhouette", legacy_kmeans_silhouette_sheet)
add_sheet_formatted(wb, "Kmeans_Legacy_Assignments", legacy_kmeans_assignments_sheet)
add_sheet_formatted(wb, "Delta_Values",        delta_sheet)
add_sheet_formatted(wb, "Top8_Markers_Boxplots", top8_summary)
add_sheet_formatted(wb, "LMM_Full_Results",    lmm_full_sheet)
add_sheet_formatted(wb, "Permutation_Tests",   perm_sheet)

# Save Excel
excel_path <- file.path(out_dir, "Supplementary_Statistics.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)
cat("Supplementary Excel saved:", excel_path, "\n")

cat("\nCSV outputs:\n")
cat("  Supplementary_Table_Direct_ThreeWay_Interaction.csv - primary direct LMM contrast\n")
cat("  Supplementary_PCA_VarianceExplained.csv            - variance explained by all PCs\n")
cat("  Supplementary_KmeansMetrics.csv                    - settings and quality metrics\n")
cat("  Supplementary_KmeansSilhouette_ByCluster.csv       - cluster separation results\n")
cat("  Supplementary_KmeansBootstrapStability.csv         - bootstrap Jaccard stability\n")
cat("  Supplementary_KmeansClusterAssignments.csv         - one cluster per rat\n")
cat("  Supplementary_PCA_ScreePlot.pdf                    - PCA variance plot\n")
cat("  Supplementary_KmeansQuality.pdf                    - elbow and silhouette plots\n")
cat("  Supplementary_KmeansStability.pdf                  - delta-based bootstrap stability plot\n")
cat("  Supplementary_KmeansRaw_PCAProjection.pdf          - raw-measurement K-means PCA projection\n")
cat("  Supplementary_KmeansRaw_Quality.pdf                - raw-measurement elbow and silhouette plots\n")
cat("  Supplementary_KmeansRaw_Stability.pdf              - raw-measurement bootstrap stability plot\n")
cat("  Supplementary_KmeansRaw_vs_Delta_Comparison.pdf    - side-by-side quality/stability summary\n")
cat("  Supplementary_KmeansRaw_vs_Delta_Comparison.csv    - numerical comparison of raw and delta inputs\n")
cat("  Figure2_SupplementaryFigure2_LegacyReproduction.pdf - historical unscaled-delta figure reproduction\n")
cat("  Supplementary_KmeansLegacyUnscaledDelta_Metrics.csv - historical figure settings and metrics\n")
cat("  Permutation_Test_Results.csv                       - exploratory sensitivity analysis\n")
cat("  Supplementary_Statistics.xlsx                      - formatted supplementary table\n")
