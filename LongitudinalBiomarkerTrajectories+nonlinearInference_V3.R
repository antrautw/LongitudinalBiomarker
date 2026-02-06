############################################################
# Longitudinal biomarker trajectories (publication-leaning)
# Updates in this version:
# - If VISIT_MODE="visit_number" and Perc_Dur exists, Visit := Perc_Dur (preferred time axis)
# - SLOPES (and slope labels in PDF) use the *continuous* Visit axis:
#     -> If Perc_Dur is present, slopes are per 1% disease progression (since Visit is 0–100)
# - Builds and writes a Feature×Disease slope table, and also merges slopes (wide) into results_df
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(mgcv)
  library(lme4)
})

# -----------------------------
# Paths (EDIT THESE)
# -----------------------------
traits_path <- "numericMeta_forBoxes.csv"  # or "Traits_v3_FirstVLast.csv"
expr_path   <- "3.Regressed(Batch)LongitudinalCSF_cleanDat-1909x143-ALS_TAMPOR+Regression_04042024_subset.csv"  # or "LongitudinalBiomarker_candidates.csv"

pdf_out     <- "Longitudinal_Biomarker_Trajectories_TargetALSSelected_ALLpoints_newslope.pdf"
res_out     <- "Longitudinal_Biomarker_model_results_TargetALSSelected_ALLpoints_newslope.csv"
slope_out   <- "Longitudinal_Biomarker_slopes_byDiseases_TargetALSSelected_ALLpoints_newslope.csv"

# -----------------------------
# Core options
# -----------------------------
VISIT_MODE <- "visit_number"    # "firstvlast" or "visit_number"

# If VISIT_MODE="visit_number" AND Perc_Dur exists, Visit becomes Perc_Dur (preferred)
# Perc_Dur scale auto-detected: 0-1 -> converted to 0-100

# Binning for Perc_Dur (continuous) so boxplots/tests remain interpretable
BIN_PERC_DUR   <- TRUE
BIN_WIDTH_PCT  <- 20            # 20% bins

# Abundance transform within subject (baseline = earliest timepoint per subject)
TRANSFORM_TO_CHANGE <- FALSE    # TRUE -> baseline-referenced change
ANCHOR_TO <- "zero"             # "zero" => % change baseline=0; "one" => fold change baseline=1

# Plot toggles
ADD_GAM_CURVES       <- TRUE
SHOW_BOXPLOTS        <- TRUE
SHOW_SUBJECT_LINES   <- TRUE
ADD_PER_VISIT_TESTS  <- TRUE
ADD_SLOPE_LABELS     <- TRUE

LEGEND_MODE <- "bottom_compact" # "bottom_compact", "right_compact", "none"

# Alpha controls
ALPHA_SUBJECT <- 0.25
ALPHA_BOX     <- 0.35
ALPHA_MEAN    <- 1.0
ALPHA_RIBBON  <- 0.12

# Smooth complexity + testing
k_max   <- 6
ci_mult <- 1.96
ASSUMPTION_ALPHA <- 0.05  # if Shapiro or Fligner p < this => use Kruskal–Wallis

# -----------------------------
# 1) Load data
# -----------------------------
traits_raw <- read.csv(traits_path, header = TRUE, row.names = 1, check.names = FALSE)
expr       <- read.csv(expr_path,   header = TRUE, row.names = 1, check.names = FALSE)

colnames(traits_raw) <- make.unique(colnames(traits_raw), sep = "__dup")
if ("sample_id" %in% colnames(traits_raw)) {
  colnames(traits_raw)[colnames(traits_raw) == "sample_id"] <- "sample_id__from_column"
}

traits <- traits_raw %>%
  tibble::rownames_to_column("sample_id") %>%
  mutate(sample_id = as.character(sample_id))

# -----------------------------
# 2) Align expression cols to traits sample_id
# -----------------------------
keep_samples <- intersect(traits$sample_id, colnames(expr))
if (length(keep_samples) == 0) stop("No overlap between traits rownames (sample_id) and expression colnames.")

expr <- expr[, keep_samples, drop = FALSE]

traits <- traits %>%
  filter(sample_id %in% colnames(expr)) %>%
  slice(match(colnames(expr), sample_id))
stopifnot(identical(traits$sample_id, colnames(expr)))

# -----------------------------
# 3) Required columns + Disease/subject/Visit
# -----------------------------
req_cols <- c("subject.id", "Disease")
missing_req <- setdiff(req_cols, names(traits))
if (length(missing_req) > 0) stop("Traits missing required columns: ", paste(missing_req, collapse = ", "))

traits$Disease    <- as.factor(traits$Disease)
traits$subject.id <- as.factor(traits$subject.id)

VISIT_MODE <- tolower(VISIT_MODE)
TIME_SOURCE <- NA_character_

if (VISIT_MODE == "firstvlast") {
  
  if (!("FirstVLast" %in% names(traits))) stop("VISIT_MODE='firstvlast' requires 'FirstVLast' in traits.")
  traits$Visit <- ifelse(traits$FirstVLast == "First", 1,
                         ifelse(traits$FirstVLast == "Last", 2, NA))
  if (any(is.na(traits$Visit))) stop("FirstVLast must be exactly 'First'/'Last'.")
  TIME_SOURCE <- "FirstVLast"
  
} else if (VISIT_MODE == "visit_number") {
  
  # Prefer Perc_Dur if present
  if ("Perc_Dur" %in% names(traits)) {
    traits$Visit <- suppressWarnings(as.numeric(traits$Perc_Dur))
    if (any(!is.finite(traits$Visit))) stop("Perc_Dur exists but could not be coerced to numeric cleanly (NAs introduced).")
    
    # Auto-detect: if looks like 0-1, convert to 0-100
    if (max(traits$Visit, na.rm = TRUE) <= 1.5) traits$Visit <- traits$Visit * 100
    TIME_SOURCE <- "Perc_Dur"
    
  } else if ("visit.number" %in% names(traits)) {
    traits$Visit <- suppressWarnings(as.numeric(traits$visit.number))
    if (any(!is.finite(traits$Visit))) stop("visit.number could not be coerced to numeric cleanly.")
    TIME_SOURCE <- "visit.number"
    
  } else if ("Visit" %in% names(traits)) {
    traits$Visit <- suppressWarnings(as.numeric(traits$Visit))
    if (any(!is.finite(traits$Visit))) stop("Visit could not be coerced to numeric cleanly.")
    TIME_SOURCE <- "Visit"
    
  } else {
    stop("VISIT_MODE='visit_number' requires Perc_Dur, visit.number, or Visit in traits.")
  }
  
} else {
  stop("VISIT_MODE must be either 'firstvlast' or 'visit_number'.")
}

traits$Visit <- as.numeric(traits$Visit)

# -----------------------------
# 3b) If using Perc_Dur, create bins for boxplots/tests (continuous axis stays Visit)
# -----------------------------
if (TIME_SOURCE == "Perc_Dur" && isTRUE(BIN_PERC_DUR)) {
  
  traits$Visit <- pmin(pmax(traits$Visit, 0), 100)
  
  breaks <- seq(0, 100, by = BIN_WIDTH_PCT)
  if (tail(breaks, 1) < 100) breaks <- c(breaks, 100)
  mids <- (breaks[-1] + breaks[-length(breaks)]) / 2
  
  idx <- findInterval(traits$Visit, vec = breaks, rightmost.closed = TRUE, all.inside = TRUE)
  traits$Visit_bin_mid   <- mids[idx]
  traits$Visit_bin_label <- paste0(breaks[idx], "-", breaks[idx + 1], "%")
  
} else {
  traits$Visit_bin_mid   <- traits$Visit
  traits$Visit_bin_label <- as.character(traits$Visit)
}

# -----------------------------
# 4) Long format expression
# -----------------------------
expr_long <- expr %>%
  rownames_to_column("Feature") %>%
  pivot_longer(cols = -Feature, names_to = "sample_id", values_to = "Abundance") %>%
  left_join(traits %>% select(sample_id, subject.id, Disease, Visit, Visit_bin_mid, Visit_bin_label), by = "sample_id") %>%
  filter(is.finite(Abundance))

# -----------------------------
# 4b) Optional baseline-referenced change within subject
# -----------------------------
ANCHOR_TO <- tolower(ANCHOR_TO)
if (isTRUE(TRANSFORM_TO_CHANGE)) {
  
  expr_long <- expr_long %>%
    group_by(Feature, subject.id) %>%
    mutate(
      baseline_visit = min(Visit, na.rm = TRUE),
      baseline_abund = mean(Abundance[Visit == baseline_visit], na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(Abundance = ifelse(!is.finite(baseline_abund) | baseline_abund == 0, NA_real_, Abundance / baseline_abund))
  
  if (ANCHOR_TO == "zero") {
    expr_long <- expr_long %>% mutate(Abundance = 100 * (Abundance - 1))
    ylab_txt <- "% change vs baseline"
  } else if (ANCHOR_TO == "one") {
    ylab_txt <- "Fold change vs baseline"
  } else {
    stop("ANCHOR_TO must be 'zero' or 'one'.")
  }
  
  expr_long <- expr_long %>% filter(is.finite(Abundance))
  
} else {
  ylab_txt <- "Abundance"
}

# Disease means for plotting (NOTE: uses continuous Visit, which is Perc_Dur if present)
mean_long <- expr_long %>%
  group_by(Feature, Disease, Visit) %>%
  summarize(mean_abund = mean(Abundance, na.rm = TRUE), .groups = "drop")

# -----------------------------
# 5) Global model helpers: GAMM with fallback to LMM
# -----------------------------
choose_k <- function(df_f, k_max = 6) {
  n_visits <- length(unique(df_f$Visit))
  min(k_max, max(3, n_visits))
}

fit_feature_lmm_test <- function(df_f) {
  red  <- try(lmer(Abundance ~ Disease + Visit + (1 | subject.id), data = df_f, REML = FALSE), silent = TRUE)
  full <- try(lmer(Abundance ~ Disease * Visit + (1 | subject.id), data = df_f, REML = FALSE), silent = TRUE)
  
  if (inherits(red, "try-error") || inherits(full, "try-error")) {
    return(list(ok = FALSE, p_global = NA_real_, note = "lmm_fit_failed", model_used = "LMM", k_used = NA_real_, fit_reml = NULL))
  }
  
  lrt <- try(anova(red, full), silent = TRUE)
  if (inherits(lrt, "try-error")) {
    return(list(ok = FALSE, p_global = NA_real_, note = "lmm_anova_failed", model_used = "LMM", k_used = NA_real_, fit_reml = NULL))
  }
  
  p <- suppressWarnings(as.numeric(lrt$`Pr(>Chisq)`[2]))
  list(ok = TRUE, p_global = p, note = "", model_used = "LMM", k_used = NA_real_, fit_reml = NULL)
}

fit_feature_gam_or_fallback <- function(df_f, k_max = 6) {
  n_visits <- length(unique(df_f$Visit))
  if (n_visits <= 2) return(fit_feature_lmm_test(df_f))
  
  k_use <- choose_k(df_f, k_max = k_max)
  
  full_form <- Abundance ~ Disease +
    s(Visit, by = Disease, k = k_use) +
    s(subject.id, bs = "re")
  
  red_form  <- Abundance ~ Disease +
    s(Visit, k = k_use) +
    s(subject.id, bs = "re")
  
  full_ml <- try(gam(full_form, data = df_f, method = "ML"), silent = TRUE)
  red_ml  <- try(gam(red_form,  data = df_f, method = "ML"), silent = TRUE)
  
  if (inherits(full_ml, "try-error") || inherits(red_ml, "try-error")) {
    out <- fit_feature_lmm_test(df_f)
    out$note <- paste0(out$note, "|fallback_from_GAM")
    return(out)
  }
  
  lrt <- try(anova(red_ml, full_ml, test = "Chisq"), silent = TRUE)
  if (inherits(lrt, "try-error")) {
    out <- fit_feature_lmm_test(df_f)
    out$note <- paste0(out$note, "|fallback_from_GAM_anova")
    return(out)
  }
  
  p_global <- suppressWarnings(as.numeric(lrt[2, "Pr(>Chi)"]))
  
  full_reml <- try(gam(full_form, data = df_f, method = "REML"), silent = TRUE)
  if (inherits(full_reml, "try-error")) full_reml <- NULL
  
  list(ok = TRUE, p_global = p_global, note = "", model_used = "GAM", k_used = k_use, fit_reml = full_reml)
}

predict_gam_curves <- function(fit_reml, df_f, ci_mult = 1.96, n_grid = 200) {
  if (is.null(fit_reml)) return(NULL)
  
  visit_min <- min(df_f$Visit, na.rm = TRUE)
  visit_max <- max(df_f$Visit, na.rm = TRUE)
  
  grid <- expand.grid(
    Visit = seq(visit_min, visit_max, length.out = n_grid),
    Disease = levels(droplevels(df_f$Disease))
  )
  grid$subject.id <- df_f$subject.id[1]
  
  pr <- try(predict(fit_reml, newdata = grid, se.fit = TRUE, exclude = "s(subject.id)"), silent = TRUE)
  if (inherits(pr, "try-error")) return(NULL)
  
  grid$fit <- as.numeric(pr$fit)
  grid$se  <- as.numeric(pr$se.fit)
  grid$lwr <- grid$fit - ci_mult * grid$se
  grid$upr <- grid$fit + ci_mult * grid$se
  
  as_tibble(grid)
}

# -----------------------------
# 6) Per-bin/per-visit group test labels: ANOVA with auto KW fallback + stars + p
# -----------------------------
p_to_stars <- function(p) {
  if (!is.finite(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  return("")
}
fmt_p <- function(p) {
  if (!is.finite(p)) return("NA")
  format.pval(p, digits = 3, eps = 1e-4)
}

safe_group_test <- function(d, assumption_alpha = 0.05) {
  d <- d %>% filter(is.finite(Abundance)) %>% mutate(Disease = droplevels(Disease))
  if (nlevels(d$Disease) < 2) return(list(test = NA_character_, p = NA_real_, note = "single_group"))
  
  ns <- table(d$Disease)
  if (any(ns < 3)) {
    kw <- try(kruskal.test(Abundance ~ Disease, data = d), silent = TRUE)
    if (inherits(kw, "try-error")) return(list(test = "KW", p = NA_real_, note = "kw_failed"))
    return(list(test = "KW", p = as.numeric(kw$p.value), note = "small_n"))
  }
  
  lm_fit <- try(lm(Abundance ~ Disease, data = d), silent = TRUE)
  if (inherits(lm_fit, "try-error")) {
    kw <- try(kruskal.test(Abundance ~ Disease, data = d), silent = TRUE)
    if (inherits(kw, "try-error")) return(list(test = "KW", p = NA_real_, note = "lm_failed_kw_failed"))
    return(list(test = "KW", p = as.numeric(kw$p.value), note = "lm_failed"))
  }
  
  sh_p <- try(shapiro.test(residuals(lm_fit))$p.value, silent = TRUE)
  if (inherits(sh_p, "try-error")) sh_p <- NA_real_
  
  fl_p <- try(fligner.test(Abundance ~ Disease, data = d)$p.value, silent = TRUE)
  if (inherits(fl_p, "try-error")) fl_p <- NA_real_
  
  shaky <- (is.finite(sh_p) && sh_p < assumption_alpha) || (is.finite(fl_p) && fl_p < assumption_alpha)
  
  if (shaky) {
    kw <- try(kruskal.test(Abundance ~ Disease, data = d), silent = TRUE)
    if (inherits(kw, "try-error")) return(list(test = "KW", p = NA_real_, note = "kw_failed"))
    return(list(test = "KW", p = as.numeric(kw$p.value), note = "assumptions_failed"))
  }
  
  aov_fit <- try(aov(Abundance ~ Disease, data = d), silent = TRUE)
  if (inherits(aov_fit, "try-error")) {
    kw <- try(kruskal.test(Abundance ~ Disease, data = d), silent = TRUE)
    if (inherits(kw, "try-error")) return(list(test = "KW", p = NA_real_, note = "aov_failed_kw_failed"))
    return(list(test = "KW", p = as.numeric(kw$p.value), note = "aov_failed"))
  }
  
  p <- try(summary(aov_fit)[[1]][["Pr(>F)"]][1], silent = TRUE)
  if (inherits(p, "try-error")) p <- NA_real_
  
  list(test = "ANOVA", p = as.numeric(p), note = "")
}

per_visit_labels <- function(df_f, assumption_alpha = 0.05) {
  out <- df_f %>%
    group_by(Visit) %>%
    group_modify(~{
      gt <- safe_group_test(.x, assumption_alpha = assumption_alpha)
      tibble(test = gt$test, p = gt$p, note = gt$note)
    }) %>%
    ungroup()
  
  y_max <- df_f %>%
    group_by(Visit) %>%
    summarize(y_pos = max(Abundance, na.rm = TRUE), .groups = "drop")
  
  y_range <- max(df_f$Abundance, na.rm = TRUE) - min(df_f$Abundance, na.rm = TRUE)
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  
  out %>%
    left_join(y_max, by = "Visit") %>%
    mutate(
      stars = vapply(p, p_to_stars, character(1)),
      p_txt = vapply(p, fmt_p, character(1)),
      label = ifelse(is.na(test) | test == "", "NA", paste0(stars, " ", test, " p=", p_txt)),
      y_pos = y_pos + 0.06 * y_range
    )
}

per_bin_labels <- function(df_f, assumption_alpha = 0.05) {
  out <- df_f %>%
    group_by(Visit_bin_mid) %>%
    group_modify(~{
      gt <- safe_group_test(.x, assumption_alpha = assumption_alpha)
      tibble(test = gt$test, p = gt$p, note = gt$note)
    }) %>%
    ungroup()
  
  y_max <- df_f %>%
    group_by(Visit_bin_mid) %>%
    summarize(y_pos = max(Abundance, na.rm = TRUE), .groups = "drop")
  
  y_range <- max(df_f$Abundance, na.rm = TRUE) - min(df_f$Abundance, na.rm = TRUE)
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  
  out %>%
    left_join(y_max, by = "Visit_bin_mid") %>%
    mutate(
      stars = vapply(p, p_to_stars, character(1)),
      p_txt = vapply(p, fmt_p, character(1)),
      label = ifelse(is.na(test) | test == "", "NA", paste0(stars, " ", test, " p=", p_txt)),
      y_pos = y_pos + 0.06 * y_range
    )
}

# -----------------------------
# 7) Slope/trend labels (Option 2): average GAM derivative per Disease
#     - If model_used == "GAM": compute mean derivative of fitted curve over Visit range
#     - Else fallback to per-disease endpoint slope (robust patch from before)
#     Notes:
#       * If Perc_Dur is present, Visit is 0–100 so "trend" is per 1% progression.
#       * This aligns the numeric trend with the trajectory shown in the figure.
# -----------------------------

# Fallback endpoint slope (per-disease endpoints; fixes missing C9, avoids global min/max)
slope_labels_endpoint <- function(mean_f) {
  
  first_df <- mean_f %>%
    group_by(Disease) %>%
    slice_min(order_by = Visit, n = 1, with_ties = FALSE) %>%
    transmute(Disease, first_visit = Visit, y1 = mean_abund)
  
  last_df <- mean_f %>%
    group_by(Disease) %>%
    slice_max(order_by = Visit, n = 1, with_ties = FALSE) %>%
    transmute(Disease, last_visit = Visit, y2 = mean_abund)
  
  first_df %>%
    inner_join(last_df, by = "Disease") %>%
    mutate(
      slope = ifelse(last_visit == first_visit, NA_real_, (y2 - y1) / (last_visit - first_visit)),
      x = (first_visit + last_visit) / 2,
      y = (y1 + y2) / 2,
      label = ifelse(is.finite(slope), paste0("slope=", signif(slope, 3)), "slope=NA")
    ) %>%
    ungroup()
}

# Option 2: mean derivative of GAM curve (per Disease)
gam_trend_slopes <- function(fit_reml, df_f, n_grid = 200) {
  
  if (is.null(fit_reml)) return(NULL)
  
  visit_min <- min(df_f$Visit, na.rm = TRUE)
  visit_max <- max(df_f$Visit, na.rm = TRUE)
  
  grid <- expand.grid(
    Visit = seq(visit_min, visit_max, length.out = n_grid),
    Disease = levels(droplevels(df_f$Disease))
  )
  # dummy subject so random effect isn't included; we exclude it explicitly too
  grid$subject.id <- df_f$subject.id[1]
  
  pr <- try(
    predict(fit_reml, newdata = grid, exclude = "s(subject.id)"),
    silent = TRUE
  )
  if (inherits(pr, "try-error")) return(NULL)
  
  grid$fit <- as.numeric(pr)
  
  # numerical derivative + mean derivative ("average trend")
  out <- as_tibble(grid) %>%
    group_by(Disease) %>%
    arrange(Visit) %>%
    mutate(
      dfit  = c(NA, diff(fit)),
      dtime = c(NA, diff(Visit)),
      deriv = dfit / dtime
    ) %>%
    summarize(
      slope = mean(deriv, na.rm = TRUE),
      first_visit = min(Visit, na.rm = TRUE),
      last_visit  = max(Visit, na.rm = TRUE),
      .groups = "drop"
    )
  
  out
}

# Helper to build label-ready slope_df from either GAM-derivative or endpoint fallback
build_trend_labels <- function(df_f, mean_f, fit_obj) {
  
  # If GAM available, use average derivative (Option 2)
  if (identical(fit_obj$model_used, "GAM") && !is.null(fit_obj$fit_reml)) {
    
    gsl <- gam_trend_slopes(fit_obj$fit_reml, df_f)
    if (!is.null(gsl) && nrow(gsl) > 0) {
      
      # place label at mid-time, near mid-mean (unobtrusive but visible)
      mean_mid <- mean(mean_f$mean_abund, na.rm = TRUE)
      
      return(
        gsl %>%
          mutate(
            x = (first_visit + last_visit) / 2,
            y = mean_mid,
            label = ifelse(is.finite(slope), paste0("trend=", signif(slope, 3)), "trend=NA"),
            slope_type = "GAM_mean_derivative"
          )
      )
    }
  }
  
  # Fallback: per-disease endpoint slope
  ep <- slope_labels_endpoint(mean_f)
  ep$slope_type <- "endpoint_fallback"
  ep
}

# -----------------------------
# 8) Legend theme helper
# -----------------------------
apply_legend_theme <- function(p, mode = "bottom_compact") {
  mode <- tolower(mode)
  if (mode == "none") return(p + theme(legend.position = "none"))
  
  if (mode == "right_compact") {
    return(
      p +
        theme(
          legend.position = "right",
          legend.title = element_text(size = 10),
          legend.text  = element_text(size = 9),
          legend.key.size = unit(0.35, "cm"),
          legend.spacing.y = unit(0.05, "cm")
        ) +
        guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))
    )
  }
  
  p +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.title = element_text(size = 10),
      legend.text  = element_text(size = 9),
      legend.key.size = unit(0.35, "cm"),
      legend.margin = margin(t = 2, r = 2, b = 2, l = 2)
    ) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(alpha = 1, linewidth = 1)))
}

# -----------------------------
# 9) Loop features: plot + stats + collect slopes (Feature×Disease)
#     PATCH: slope_df now built by build_trend_labels() (Option 2) and saved in slope_results
# -----------------------------
features <- unique(expr_long$Feature)
results  <- vector("list", length(features))
slope_results <- vector("list", length(features))

pdf(pdf_out, width = 7.2, height = 5.8)

for (ii in seq_along(features)) {
  
  f      <- features[ii]
  df_f   <- expr_long %>% filter(Feature == f)
  mean_f <- mean_long %>% filter(Feature == f)
  
  fit <- fit_feature_gam_or_fallback(df_f, k_max = k_max)
  
  results[[ii]] <- tibble(
    Feature = f,
    n_samples = nrow(df_f),
    n_subjects = length(unique(df_f$subject.id)),
    n_visits = length(unique(df_f$Visit)),
    n_disease_levels = nlevels(droplevels(df_f$Disease)),
    time_source = TIME_SOURCE,
    model_used = fit$model_used,
    k_used = fit$k_used,
    p_global = fit$p_global,
    note = fit$note
  )
  
  subtitle_txt <- if (is.finite(fit$p_global)) {
    paste0(
      fit$model_used, " global Disease×Time: p=", format.pval(fit$p_global, digits = 3, eps = 1e-4),
      if (fit$model_used == "GAM") paste0(" | k=", fit$k_used) else "",
      " | time=", TIME_SOURCE,
      if (TRANSFORM_TO_CHANGE) paste0(" | change=", ANCHOR_TO) else ""
    )
  } else {
    paste0(fit$model_used, " global Disease×Time: p=NA | time=", TIME_SOURCE)
  }
  
  # Optional GAM curves + ribbons
  curve_df <- NULL
  if (isTRUE(ADD_GAM_CURVES) && identical(fit$model_used, "GAM")) {
    curve_df <- predict_gam_curves(fit$fit_reml, df_f, ci_mult = ci_mult)
  }
  
  # Per-"timepoint" labels: visit vs bin (bin when Perc_Dur)
  lbl_df <- NULL
  if (isTRUE(ADD_PER_VISIT_TESTS)) {
    if (TIME_SOURCE == "Perc_Dur" && isTRUE(BIN_PERC_DUR)) {
      lbl_df <- per_bin_labels(df_f, assumption_alpha = ASSUMPTION_ALPHA)
    } else {
      lbl_df <- per_visit_labels(df_f, assumption_alpha = ASSUMPTION_ALPHA)
    }
  }
  
  # ---- NEW: trend slope labels (Option 2) + collection table
  slope_df <- NULL
  if (isTRUE(ADD_SLOPE_LABELS) && nrow(mean_f) > 0) {
    
    slope_df <- build_trend_labels(df_f, mean_f, fit)
    
    slope_results[[ii]] <- slope_df %>%
      transmute(
        Feature = f,
        Disease,
        time_source = TIME_SOURCE,
        slope_type,
        slope,
        first_visit,
        last_visit
      )
    
  } else {
    slope_results[[ii]] <- tibble(
      Feature = f, Disease = NA_character_, time_source = TIME_SOURCE,
      slope_type = NA_character_, slope = NA_real_,
      first_visit = NA_real_, last_visit = NA_real_
    )
  }
  
  dodge_w <- 0.6
  
  p <- ggplot(df_f, aes(x = Visit, y = Abundance)) +
    
    # GAM ribbon/line (explicit x mapping)
    { if (!is.null(curve_df)) geom_ribbon(
      data = curve_df,
      aes(x = Visit, ymin = lwr, ymax = upr, fill = Disease, group = Disease),
      alpha = ALPHA_RIBBON,
      inherit.aes = FALSE
    ) } +
    { if (!is.null(curve_df)) geom_line(
      data = curve_df,
      aes(x = Visit, y = fit, color = Disease, group = Disease),
      linewidth = 1.05,
      alpha = 1,
      inherit.aes = FALSE
    ) } +
    
    # Boxplots: binned midpoints when Perc_Dur; otherwise Visit_bin_mid==Visit
    { if (isTRUE(SHOW_BOXPLOTS)) geom_boxplot(
      data = df_f,
      aes(x = Visit_bin_mid, y = Abundance, group = interaction(Disease, Visit_bin_mid), color = Disease),
      position = position_dodge(width = dodge_w),
      alpha = ALPHA_BOX,
      outlier.alpha = 0.25,
      linewidth = 0.55,
      inherit.aes = FALSE
    ) } +
    
    # Subject trajectories
    { if (isTRUE(SHOW_SUBJECT_LINES)) geom_line(
      aes(group = subject.id, color = Disease),
      alpha = ALPHA_SUBJECT,
      linewidth = 0.45
    ) } +
    
    # Disease mean line/points
    geom_line(
      data = mean_f,
      aes(x = Visit, y = mean_abund, color = Disease, group = Disease),
      linewidth = 1.2,
      alpha = ALPHA_MEAN
    ) +
    geom_point(
      data = mean_f,
      aes(x = Visit, y = mean_abund, color = Disease),
      size = 3,
      alpha = ALPHA_MEAN
    ) +
    
    # Labels: use Visit_bin_mid when binned Perc_Dur
    { if (!is.null(lbl_df) && ("Visit_bin_mid" %in% names(lbl_df))) geom_text(
      data = lbl_df,
      aes(x = Visit_bin_mid, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 3.0,
      vjust = 0
    ) } +
    { if (!is.null(lbl_df) && ("Visit" %in% names(lbl_df))) geom_text(
      data = lbl_df,
      aes(x = Visit, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 3.0,
      vjust = 0
    ) } +
    
    # Trend label (Option 2): average GAM derivative if available, else endpoint slope
    { if (!is.null(slope_df)) geom_text(
      data = slope_df,
      aes(x = x, y = y, label = label, color = Disease),
      inherit.aes = FALSE,
      size = 2.8,
      alpha = 0.85,
      show.legend = FALSE
    ) } +
    
    labs(
      title = f,
      subtitle = subtitle_txt,
      x = ifelse(TIME_SOURCE == "Perc_Dur", "Disease progression (%)", "Visit"),
      y = ylab_txt
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 9.5)
    )
  
  if (TIME_SOURCE == "Perc_Dur") {
    p <- p + scale_x_continuous(breaks = seq(0, 100, by = 20), limits = c(0, 100))
  } else {
    p <- p + scale_x_continuous(breaks = sort(unique(df_f$Visit)))
  }
  
  p <- apply_legend_theme(p, LEGEND_MODE)
  print(p)
}

dev.off()


# -----------------------------
# 10) Save results table (global test p + FDR) and slopes table
# -----------------------------
results_df <- bind_rows(results) %>%
  mutate(q_global = p.adjust(p_global, method = "fdr")) %>%
  arrange(q_global, p_global)

slope_df_all <- bind_rows(slope_results) %>%
  filter(!is.na(Disease)) %>%
  arrange(Feature, Disease)

# Write slopes in long format (Feature×Disease)
write.csv(slope_df_all, slope_out, row.names = FALSE)

# Also merge slopes into results_df as wide columns (slope_<Disease>)
slope_wide <- slope_df_all %>%
  mutate(Disease = as.character(Disease)) %>%
  select(Feature, Disease, slope) %>%
  pivot_wider(names_from = Disease, values_from = slope, names_prefix = "slope_")

results_df <- results_df %>%
  left_join(slope_wide, by = "Feature")

write.csv(results_df, res_out, row.names = FALSE)

message(
  "Done.\nSaved PDF: ", pdf_out,
  "\nSaved results: ", res_out,
  "\nSaved slopes: ", slope_out,
  "\nVISIT_MODE = ", VISIT_MODE,
  "\nTIME_SOURCE used = ", TIME_SOURCE,
  if (TIME_SOURCE == "Perc_Dur") paste0("\nBIN_PERC_DUR = ", BIN_PERC_DUR, " (", BIN_WIDTH_PCT, "% bins)") else ""
)
