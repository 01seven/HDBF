# ==============================================================================
# Real-data GO analysis for Table 2
#
# Manuscript:
#   Power-Enhanced Cauchy Combination Tests for the Two-Sample
#   High-Dimensional Behrens-Fisher Problem
#
# Purpose:
#   Reproduce the real-data gene ontology (GO) analysis based on the ALL data
#   set and summarize the number of significant GO sets for the eight methods
#   reported in Table 2.
#
# Strategy:
#   Full-spectrum comparison across ZWL, CQ, SKK, CLX, CCTL, WCCT, EWCCT,
#   and EWCCTS.
#
# Tested with:
#   R 4.4.3 on Windows
#
# Main entry points:
#   source("HDBF_real_data_GO_analysis.R")
#   results <- run_real_data_analysis()
#
# Notes:
#   - This script requires both CRAN and Bioconductor packages.
#   - It writes WCCT_Advantage_Report.txt by default.
#   - To reduce CPU use, modify hdbf_real_data_config() after sourcing this
#     file, for example cfg <- hdbf_real_data_config(); cfg$workers <- 1.
# ==============================================================================

#### Package requirements ####

required_bioc <- c("ALL", "genefilter", "hgu95av2.db", "GO.db", "AnnotationDbi")
required_cran <- c(
  "foreach", "doParallel",
  "highDmean", "highmean",
  "matrixStats", "stats", "knitr"
)
required_packages <- c(required_bioc, required_cran)

# Stop early if a package required by the real-data analysis is not installed.
check_packages <- function(packages = required_packages) {
  missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0L) {
    stop(
      "Please install the following R packages before running this script: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

# Load packages that provide data objects, annotation databases, and parallel backends.
load_required_packages <- function() {
  check_packages()
  invisible(lapply(required_packages, library, character.only = TRUE))
}

# Return the real-data settings used for the Table 2 GO analysis.
hdbf_real_data_config <- function() {
  list(
    target_probe_count = 2391L,
    min_go_size = 10L,
    max_go_size = 3000L,
    fdr_level = 0.05,
    top_unique_sets = 15L,
    workers = max(1L, as.integer(parallel::detectCores()) - 1L, na.rm = TRUE),
    output_file = "WCCT_Advantage_Report.txt"
  )
}

#### Test statistics and p-value combinations ####

# Compute the Chen-Qin quadratic-form statistic using precomputed Gram matrices.
calculate_Tn_optimized <- function(H1, H2, K12, n1, n2) {
  term1 <- (sum(H1) - sum(diag(H1))) / (n1 * (n1 - 1))
  term2 <- (sum(H2) - sum(diag(H2))) / (n2 * (n2 - 1))
  term3 <- -2 * sum(K12) / (n1 * n2)
  term1 + term2 + term3
}

# Estimate tr(Sigma^2) from a Gram matrix for one sample.
est_sq_R_from_H <- function(H, m) {
  if (m < 4L) {
    return(NA_real_)
  }

  H_sq <- H^2
  r_sum <- rowSums(H)
  r_sum_sq <- rowSums(H_sq)
  HH <- H %*% H

  S_sq <- outer(r_sum_sq, r_sum_sq, "+") - 2 * HH
  S_val <- outer(r_sum, r_sum, "-")

  diag_H <- diag(H)
  Di <- matrix(diag_H, m, m)
  Dj <- t(Di)

  diff_i <- Di - H
  diff_j <- H - Dj

  S_sq_v <- S_sq - diff_i^2 - diff_j^2
  S_val_v <- S_val - diff_i - diff_j

  Inner <- 2 * ((m - 2L) * S_sq_v - S_val_v^2)
  diag(Inner) <- 0

  sum(Inner) / (4 * m * (m - 1) * (m - 2) * (m - 3))
}

# Estimate the cross trace tr(Sigma_1 Sigma_2) from two centered samples.
cross_trace_est <- function(X, Y) {
  m <- nrow(X)
  n <- nrow(Y)
  Xc <- t(t(X) - matrixStats::colMeans2(X))
  Yc <- t(t(Y) - matrixStats::colMeans2(Y))
  cross_terms <- tcrossprod(Yc, Xc)
  sum(cross_terms^2) / ((m - 1) * (n - 1))
}

# Compute the WCCT p-value from coordinate-wise Welch t-tests.
calc_wcct_R <- function(X1, X2) {
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  m1 <- matrixStats::colMeans2(X1)
  m2 <- matrixStats::colMeans2(X2)
  v1 <- matrixStats::colVars(X1)
  v2 <- matrixStats::colVars(X2)

  vn1 <- v1 / n1
  vn2 <- v2 / n2
  se <- sqrt(vn1 + vn2)
  se[se < 1e-12] <- 1e-12

  t_stat <- (m1 - m2) / se
  df <- (vn1 + vn2)^2 / ((vn1^2) / (n1 - 1) + (vn2^2) / (n2 - 1))

  p_vals <- 2 * pt(abs(t_stat), df = df, lower.tail = FALSE)
  tan_vals <- tan((0.5 - p_vals) * pi)
  tan_vals[!is.finite(tan_vals) | tan_vals > 1e10] <- 1e10

  stat <- mean(tan_vals)
  0.5 - atan(stat) / pi
}

# Combine two p-values using the Cauchy combination rule.
cauchy_comb <- function(p1, p2) {
  if (is.na(p1) || is.na(p2)) {
    return(NA_real_)
  }

  t1 <- tan((0.5 - p1) * pi)
  t2 <- tan((0.5 - p2) * pi)

  if (!is.finite(t1)) {
    t1 <- 1e10
  }
  if (!is.finite(t2)) {
    t2 <- 1e10
  }

  pcauchy(0.5 * t1 + 0.5 * t2, lower.tail = FALSE)
}

# Compute the EWCCT p-value using the finite enhancement term in the manuscript.
ewcct_pvalue <- function(p_wcct, Qn, n1, n2, d, threshold_coef = 2.1) {
  if (is.na(p_wcct) || is.na(Qn)) {
    return(NA_real_)
  }

  if (abs(Qn) <= sqrt(threshold_coef * log(d))) {
    return(p_wcct)
  }

  T_wcct <- tan((0.5 - p_wcct) * pi)
  T_ewcct <- T_wcct + n1 + n2
  0.5 - atan(T_ewcct) / pi
}

#### Data loading and GO mapping ####

# Load and filter the ALL expression data to match the Yu et al. probe count.
prepare_real_data <- function(config = hdbf_real_data_config()) {
  cat(">> [Step 1] Loading Data & Preprocessing (Strict Match to Yu et al.)...\n")
  data(ALL)

  bcell <- ALL[, ALL$BT %in% c("B", "B1", "B2", "B3", "B4")]
  mol.biol <- bcell$mol.biol
  subset_idx <- which(mol.biol %in% c("BCR/ABL", "NEG"))
  eset <- bcell[, subset_idx]
  mol.biol <- factor(mol.biol[subset_idx])

  cat(sprintf(
    "   Filtering probes to match paper count (Target: %s)...\n",
    format(config$target_probe_count, big.mark = ",", scientific = FALSE)
  ))
  all_iqrs <- apply(exprs(eset), 1, IQR)
  valid_probes <- names(sort(all_iqrs, decreasing = TRUE))[seq_len(config$target_probe_count)]
  X_filtered <- exprs(eset)[valid_probes, ]
  cat(
    "   Final Matrix Dimensions:",
    nrow(X_filtered),
    "probes x",
    ncol(X_filtered),
    "samples\n"
  )

  grp1_idx <- which(mol.biol == "BCR/ABL")
  grp2_idx <- which(mol.biol == "NEG")
  X1_full <- t(X_filtered[, grp1_idx])
  X2_full <- t(X_filtered[, grp2_idx])

  list(
    valid_probes = valid_probes,
    X1_full = X1_full,
    X2_full = X2_full
  )
}

# Build the probe-to-GO mapping for one ontology.
build_go_map <- function(ontology, valid_probes, config = hdbf_real_data_config()) {
  go_df <- AnnotationDbi::select(
    hgu95av2.db,
    keys = valid_probes,
    columns = c("GO", "ONTOLOGY"),
    keytype = "PROBEID"
  )
  go_df <- go_df[go_df$ONTOLOGY == ontology, ]
  go_df <- na.omit(go_df)
  go_list <- split(go_df$PROBEID, go_df$GO)
  valid_ids <- names(go_list)[vapply(go_list, function(x) {
    len <- length(unique(x))
    len >= config$min_go_size && len <= config$max_go_size
  }, logical(1))]
  go_list[valid_ids]
}

# Build all BP, CC, and MF GO mappings used in the analysis.
build_go_maps <- function(valid_probes, config = hdbf_real_data_config()) {
  cat(">> [Step 2] Building GO Mappings...\n")
  go_maps <- list(
    BP = build_go_map("BP", valid_probes, config),
    CC = build_go_map("CC", valid_probes, config),
    MF = build_go_map("MF", valid_probes, config)
  )
  cat(
    "   Mapped Sets -> BP:",
    length(go_maps$BP),
    "| CC:",
    length(go_maps$CC),
    "| MF:",
    length(go_maps$MF),
    "\n"
  )
  go_maps
}

#### GO-set testing ####

# Run all nine tests for one GO set.
run_tests_optimized <- function(go_id, X1, X2, gene_map) {
  probes <- unique(gene_map[[go_id]])
  valid <- intersect(probes, colnames(X1))
  if (length(valid) < 4L) {
    return(rep(NA_real_, 8L))
  }

  d1 <- X1[, valid, drop = FALSE]
  d2 <- X2[, valid, drop = FALSE]
  n1 <- nrow(d1)
  n2 <- nrow(d2)
  p <- ncol(d1)

  res <- c(
    ZWL = NA_real_,
    CQ = NA_real_,
    SKK = NA_real_,
    CLX = NA_real_,
    CCTL = NA_real_,
    WCCT = NA_real_,
    EWCCT = NA_real_,
    EWCCTS = NA_real_
  )

  try({
    res["ZWL"] <- highDmean::zwl_test(d1, d2)$pvalue
  }, silent = TRUE)
  try({
    res["CQ"] <- highmean::apval_Chen2010(d1, d2, eq.cov = FALSE)$pval
  }, silent = TRUE)
  try({
    res["SKK"] <- highDmean::SKK_test(d1, d2)$pvalue
  }, silent = TRUE)
  try({
    res["CLX"] <- highmean::apval_Cai2014(d1, d2, eq.cov = FALSE)$pval
  }, silent = TRUE)

  res["CCTL"] <- cauchy_comb(res["CLX"], res["CQ"])
  try({
    res["WCCT"] <- calc_wcct_R(d1, d2)
  }, silent = TRUE)

  try({
    H1 <- tcrossprod(d1)
    H2 <- tcrossprod(d2)
    K12 <- tcrossprod(d1, d2)
    Tn <- calculate_Tn_optimized(H1, H2, K12, n1, n2)
    tr1 <- est_sq_R_from_H(H1, n1)
    tr2 <- est_sq_R_from_H(H2, n2)
    tr12 <- cross_trace_est(d1, d2)
    var_est <- max(
      1e-12,
      (2 / (n1 * (n1 - 1))) * tr1 +
        (2 / (n2 * (n2 - 1))) * tr2 +
        (4 / (n1 * n2)) * tr12
    )
    Qn <- Tn / sqrt(var_est)

    res["EWCCT"] <- ewcct_pvalue(res["WCCT"], Qn, n1, n2, p)
  }, silent = TRUE)
  res["EWCCTS"] <- cauchy_comb(res["EWCCT"], res["SKK"])
  res
}

# Run the GO-set tests in parallel for all ontology categories.
run_go_analysis <- function(real_data, go_maps, config = hdbf_real_data_config()) {
  cat(">> [Step 3] Running Parallel Analysis...\n")

  cl <- parallel::makeCluster(config$workers)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  export_funcs <- c(
    "run_tests_optimized",
    "calculate_Tn_optimized",
    "est_sq_R_from_H",
    "cross_trace_est",
    "calc_wcct_R",
    "cauchy_comb"
  )

  all_results <- list()
  for (cat in names(go_maps)) {
    ids <- names(go_maps[[cat]])
    if (length(ids) > 0L) {
      res <- foreach::foreach(
        g = ids,
        .combine = rbind,
        .packages = c("highDmean", "highmean", "matrixStats"),
        .export = export_funcs
      ) %dopar% {
        run_tests_optimized(g, real_data$X1_full, real_data$X2_full, go_maps[[cat]])
      }
      if (is.vector(res)) {
        res <- t(res)
      }
      rownames(res) <- ids
      all_results[[cat]] <- res
    }
  }

  all_results
}

#### Summaries and output ####

# Count significant GO sets for each method after BH adjustment.
significant_counts <- function(all_results, go_maps, config = hdbf_real_data_config()) {
  sig_counts <- data.frame(Method = colnames(all_results$BP))
  for (cat in names(go_maps)) {
    if (!is.null(all_results[[cat]])) {
      sig_counts[[cat]] <- apply(
        all_results[[cat]],
        2,
        function(x) sum(p.adjust(x, "BH") < config$fdr_level, na.rm = TRUE)
      )
    }
  }
  sig_counts
}

# Write the WCCT advantage report used for the real-data summary.
write_advantage_report <- function(all_results, go_maps, config = hdbf_real_data_config()) {
  cat(sprintf(
    ">> [Step 4] Generating Advantage Report (%s)...\n",
    config$output_file
  ))

  sig_counts <- significant_counts(all_results, go_maps, config)

  sink(config$output_file)
  on.exit(sink(), add = TRUE)

  cat("================================================================\n")
  cat("COMPARATIVE ADVANTAGE REPORT: Eight-Method Analysis\n")
  cat("================================================================\n\n")

  cat("[1] Overall Performance (Total Significant Sets, FDR < 0.05)\n")
  print(sig_counts)

  cat("\n[2] Reliability-Oriented Interpretation\n")
  cat("EWCCTS is retained as the proposed combined method because it pairs EWCCT with the scale-invariant SKK component.\n")
  cat("Interpret raw discovery counts together with the simulation evidence on type I error control.\n")

  invisible(list(
    output_file = config$output_file,
    significant_counts = sig_counts
  ))
}

#### Public entry points ####

# Run the full real-data GO analysis for Table 2.
run_real_data_analysis <- function(config = hdbf_real_data_config(), write_report = TRUE) {
  load_required_packages()
  real_data <- prepare_real_data(config)
  go_maps <- build_go_maps(real_data$valid_probes, config)
  all_results <- run_go_analysis(real_data, go_maps, config)

  report <- NULL
  if (isTRUE(write_report)) {
    report <- write_advantage_report(all_results, go_maps, config)
    cat(">> Report saved to:", config$output_file, "\n")
  }

  list(
    real_data = real_data,
    go_maps = go_maps,
    all_results = all_results,
    report = report
  )
}

#### Main ####

# Run the real-data analysis when this file is executed directly.
main <- function() {
  config <- hdbf_real_data_config()
  cat("HDBF real-data GO analysis\n")
  cat(sprintf(
    "target_probe_count = %d, GO size range = [%d, %d]\n",
    config$target_probe_count,
    config$min_go_size,
    config$max_go_size
  ))
  cat(sprintf("workers = %d\n", config$workers))
  start_time <- proc.time()
  run_real_data_analysis(config, write_report = TRUE)
  elapsed <- (proc.time() - start_time)[3]
  cat(sprintf("Completed in %.2f seconds.\n", elapsed))
}

if (sys.nframe() == 0L) {
  main()
}
