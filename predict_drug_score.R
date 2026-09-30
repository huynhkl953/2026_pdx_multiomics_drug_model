#!/usr/bin/env Rscript
# =============================================================================
# 01_run_fig8_models.R   (checkpointed / resumable)  -- v2: 1 uM + 10 uM everywhere
# -----------------------------------------------------------------------------
# Keeps the manuscript pipeline (custom calcPhenotype: ComBat + removeLowVaryingGenes
# = 0.2 + per-drug ridge); changes only the EVALUATION to per-PDX x compound k-fold CV.
#   (1) per-PDX x compound predicted vs observed      -> 1 + 10 uM
#   (3) matched-feature (RNA n protein symbols)       -> 1 (NEW) + 10 uM
#   (2) tail-based per-sample imputation (optional)   -> 1 + 10 uM
# Checkpoints: inputs.rds | heldout_pred_<name>.rds/.csv | cache/<name>/fold_XX.rds
# Existing runs load from cache; only matched_*_1uM is new.
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(matrixStats)
  library(caret); library(sva); library(ridge)
})

# ---------------------------------- CONFIG -----------------------------------
base_dir <- "/lustre/home/huynhk4/BioinformaticsCore/Project 2"
paths <- list(
  fn_oncopredict = file.path(base_dir, "code", "function_oncoPredict.R"),
  protein   = file.path(base_dir, "data_version_2", "20251203_DetectedProteins_MaxSum_RepAvgImpute_RepAveraged.csv"),
  rna       = file.path(base_dir, "data_version_2", "25.11.20_HUMAN_U54_AllSamples_nf-JaxPDXNet_RSEM_Gene_TPM_U54manuscript.csv"),
  drug_1um  = file.path(base_dir, "data_version_2", "20251113_All_U54_1uM_inhibition.xlsx"),
  drug_10um = file.path(base_dir, "data_version_2", "20251113_All_U54_10uM_inhibition.xlsx"))
outdir    <- file.path(base_dir, "results", "results_v2", "Fig8_rerun_artifact")
cache_dir <- file.path(outdir, "cache"); dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

resume      <- TRUE
clear_cache <- FALSE   # TRUE only if a modelling parameter / input changed

protein_id_col <- "Protein_ID"; protein_symbol_col <- "Gene"; protein_meta_cols <- 2
rna_symbol_col <- "geneSym";    rna_sample_regex   <- "^VCU"

top_n_features <- 8000; cv_folds <- 9; 
min_nonNA_train <- 3;   min_frac_nonNA <- 0.60
sample_rename <- c("VCUBC003_052223" = "VCU-BC-003")
preimpute_protein_path <- NA   # PRE-imputation protein matrix (NAs kept) to enable point (2)

if (clear_cache) {
  unlink(cache_dir, recursive = TRUE); dir.create(cache_dir, recursive = TRUE)
  file.remove(list.files(outdir, pattern = "\\.rds$", full.names = TRUE))
}
source(paths$fn_oncopredict); stopifnot(exists("calcPhenotype"))

# ---------------------------------- HELPERS ----------------------------------
clean_samples <- function(x) {
  x <- gsub("\\.", "-", x); x <- sub("_.*", "", x)
  ix <- x %in% names(sample_rename); x[ix] <- sample_rename[x[ix]]; x
}
build_expr <- function(df, symbol_vec, sample_cols, top_n = NULL) {
  m <- as.matrix(df[, sample_cols, drop = FALSE]); storage.mode(m) <- "numeric"
  colnames(m) <- clean_samples(colnames(m)); sym <- as.character(symbol_vec)
  ord <- order(matrixStats::rowVars(m, na.rm = TRUE), decreasing = TRUE)
  m <- m[ord, , drop = FALSE]; sym <- sym[ord]
  keep <- !duplicated(sym) & !is.na(sym) & sym != ""
  m <- m[keep, , drop = FALSE]; rownames(m) <- sym[keep]
  if (!is.null(top_n) && top_n < nrow(m))
    m <- m[order(matrixStats::rowVars(m, na.rm = TRUE), decreasing = TRUE)[seq_len(top_n)], , drop = FALSE]
  m
}
load_drug <- function(path) {
  d <- as.data.frame(readxl::read_excel(path)); colnames(d)[2] <- "Target"
  map <- d[, c("Drug", "Target")]
  wide <- d %>% dplyr::select(-Target) %>%
    tidyr::pivot_longer(-Drug, names_to = "Sample", values_to = "Value") %>%
    tidyr::pivot_wider(names_from = Drug, values_from = Value) %>% as.data.frame()
  rownames(wide) <- clean_samples(wide$Sample); wide$Sample <- NULL
  list(mat = as.matrix(wide), map = map)
}
drug_keep <- function(mat, samples) {
  m <- mat[samples, , drop = FALSE]
  names(which(colMeans(!is.na(m)) >= min_frac_nonNA & matrixStats::colVars(m, na.rm = TRUE) > 0))
}

run_cv <- function(expr, ptype, k = 9,  batchCorrect = "eb",
                   powerTransform = FALSE, removeLowVaryingGenes = 0.2,
                   selection = 1, removeLowVaringGenesFrom = "homogenizeData",
                   percent = 80, min_test = 2, min_nonNA = 3, printOutput = FALSE,
                   fold_cache = NULL, resume = TRUE) {
  samples <- colnames(expr); n <- length(samples)
  ptype <- ptype[samples, , drop = FALSE]; stopifnot(all(rownames(ptype) == samples))
  k_use <- max(2, min(k, floor(n / min_test)))
  fold <- caret::createFolds(seq_len(n), k = k_use, list = FALSE)
  if (!is.null(fold_cache)) dir.create(fold_cache, recursive = TRUE, showWarnings = FALSE)
  preds <- list()
  for (f in sort(unique(fold))) {
    fc <- if (!is.null(fold_cache)) file.path(fold_cache, sprintf("fold_%02d.rds", f)) else NULL
    if (!is.null(fc) && resume && file.exists(fc)) {
      message("    fold ", f, "/", k_use, " (cached)"); preds[[as.character(f)]] <- readRDS(fc); next
    }
    message("    fold ", f, "/", k_use, " ...")
    te <- which(fold == f); tr <- setdiff(seq_len(n), te)
    tr_ptype <- ptype[tr, , drop = FALSE]
    ok <- vapply(seq_len(ncol(tr_ptype)), function(j) {
      x <- tr_ptype[, j]; sum(!is.na(x)) >= min_nonNA && stats::sd(x, na.rm = TRUE) > 0
    }, logical(1))
    p <- calcPhenotype(trainingExprData = as.matrix(expr[, tr, drop = FALSE]),
                       trainingPtype = as.matrix(tr_ptype[, ok, drop = FALSE]),
                       testExprData = as.matrix(expr[, te, drop = FALSE]),
                       batchCorrect = batchCorrect, powerTransformPhenotype = powerTransform,
                       removeLowVaryingGenes = removeLowVaryingGenes, minNumSamples = 10,
                       selection = selection, printOutput = printOutput, pcr = FALSE,
                       removeLowVaringGenesFrom = removeLowVaringGenesFrom,
                       report_pc = FALSE, cc = FALSE, percent = percent, rsq = FALSE, folder = FALSE)
    p <- as.data.frame(p, check.names = FALSE); p$.PDX <- rownames(p)
    if (!is.null(fc)) saveRDS(p, fc)
    preds[[as.character(f)]] <- p
  }
  out <- dplyr::bind_rows(preds); rownames(out) <- out$.PDX; out$.PDX <- NULL
  as.matrix(out[samples, , drop = FALSE])
}
cv_cached <- function(name, expr, ptype, ...) {
  rds <- file.path(outdir, paste0("heldout_pred_", name, ".rds"))
  if (resume && file.exists(rds)) { message("[skip] ", name, " (cached run)"); return(readRDS(rds)) }
  message("[run ] ", name, " ...")
  pred <- run_cv(expr, ptype, fold_cache = file.path(cache_dir, name), resume = resume, ...)
  saveRDS(pred, rds); write.csv(pred, file.path(outdir, paste0("heldout_pred_", name, ".csv")))
  pred
}
per_observation <- function(pred, drug_mat, modality, conc) {
  d <- intersect(colnames(pred), colnames(drug_mat)); s <- intersect(rownames(pred), rownames(drug_mat))
  df <- data.frame(PDX = rep(s, times = length(d)), Drug = rep(d, each = length(s)),
                   Predicted = as.vector(pred[s, d, drop = FALSE]),
                   Observed  = as.vector(drug_mat[s, d, drop = FALSE]),
                   Modality = modality, Conc = conc)
  df[stats::complete.cases(df), ]
}
per_drug_metrics <- function(pred, drug_mat, modality, conc, min_n = 3) {
  d <- intersect(colnames(pred), colnames(drug_mat)); s <- intersect(rownames(pred), rownames(drug_mat))
  P <- pred[s, d, drop = FALSE]; O <- drug_mat[s, d, drop = FALSE]
  do.call(rbind, lapply(d, function(k) {
    p <- P[, k]; o <- O[, k]; ok <- !is.na(p) & !is.na(o)
    if (sum(ok) < min_n) return(data.frame(Drug=k, n=sum(ok), r=NA, MSE=NA, RMSE=NA, Modality=modality, Conc=conc))
    e <- p[ok] - o[ok]
    data.frame(Drug=k, n=sum(ok), r=suppressWarnings(cor(p[ok], o[ok])),
               MSE=mean(e^2), RMSE=sqrt(mean(e^2)), Modality=modality, Conc=conc)
  }))
}

# ------------------------------ INPUTS (cached) ------------------------------
inputs_rds <- file.path(outdir, "inputs.rds")
if (resume && file.exists(inputs_rds)) {
  message("[skip] inputs (cached)"); list2env(readRDS(inputs_rds), envir = environment())
} else {
  protein_raw <- readr::read_csv(paths$protein, show_col_types = FALSE)
  rna_raw     <- readr::read_csv(paths$rna,     show_col_types = FALSE)
  protein_scol <- colnames(protein_raw)[(protein_meta_cols + 1):ncol(protein_raw)]
  protein_expr_full <- build_expr(protein_raw, protein_raw[[protein_symbol_col]], protein_scol)
  protein_expr      <- build_expr(protein_raw, protein_raw[[protein_symbol_col]], protein_scol, top_n_features)
  rna_scol <- grep(rna_sample_regex, colnames(rna_raw), value = TRUE)
  rna_expr_full <- build_expr(rna_raw, rna_raw[[rna_symbol_col]], rna_scol)
  rna_expr      <- build_expr(rna_raw, rna_raw[[rna_symbol_col]], rna_scol, top_n_features)
  d1 <- load_drug(paths$drug_1um); d10 <- load_drug(paths$drug_10um)
  drug_map <- unique(rbind(d1$map, d10$map))
  common_1um  <- Reduce(intersect, list(colnames(rna_expr), colnames(protein_expr), rownames(d1$mat)))
  common_10um <- Reduce(intersect, list(colnames(rna_expr), colnames(protein_expr), rownames(d10$mat)))
  keep10 <- intersect(drug_keep(d1$mat, common_10um), drug_keep(d10$mat, common_10um))
  keep1  <- intersect(drug_keep(d1$mat, common_1um),  drug_keep(d10$mat, common_1um))
  saveRDS(list(protein_expr=protein_expr, protein_expr_full=protein_expr_full,
               rna_expr=rna_expr, rna_expr_full=rna_expr_full, d1=d1, d10=d10,
               drug_map=drug_map, common_1um=common_1um, common_10um=common_10um,
               keep1=keep1, keep10=keep10), inputs_rds)
  write.csv(drug_map, file.path(outdir, "drug_class_map.csv"), row.names = FALSE)
}
message("Common samples 1uM: ", length(common_1um), " | 10uM: ", length(common_10um))

# ------------------------ MAIN CV (top-n) 1 + 10 uM --------------------------
pred_rna_10  <- cv_cached("RNA_10uM",     rna_expr[, common_10um],     d10$mat[common_10um, keep10], k=cv_folds, min_nonNA=min_nonNA_train)
pred_prot_10 <- cv_cached("protein_10uM", protein_expr[, common_10um], d10$mat[common_10um, keep10], k=cv_folds, min_nonNA=min_nonNA_train)
pred_rna_1   <- cv_cached("RNA_1uM",      rna_expr[, common_1um],      d1$mat[common_1um, keep1],    k=cv_folds, min_nonNA=min_nonNA_train)
pred_prot_1  <- cv_cached("protein_1uM",  protein_expr[, common_1um],  d1$mat[common_1um, keep1],    k=cv_folds, min_nonNA=min_nonNA_train)

# ------------- MATCHED FEATURES (Reviewer point 3) 1 (NEW) + 10 uM -----------
common_genes <- intersect(rownames(rna_expr_full), rownames(protein_expr_full))
message("Matched features (shared symbols): ", length(common_genes))
rna_m <- rna_expr_full[common_genes, ]; prot_m <- protein_expr_full[common_genes, ]
pred_rna_m_10  <- cv_cached("matched_RNA_10uM",     rna_m[, common_10um],  d10$mat[common_10um, keep10], k=cv_folds, removeLowVaryingGenes=0, min_nonNA=min_nonNA_train)
pred_prot_m_10 <- cv_cached("matched_protein_10uM", prot_m[, common_10um], d10$mat[common_10um, keep10], k=cv_folds, removeLowVaryingGenes=0, min_nonNA=min_nonNA_train)
pred_rna_m_1   <- cv_cached("matched_RNA_1uM",      rna_m[, common_1um],   d1$mat[common_1um, keep1],    k=cv_folds, removeLowVaryingGenes=0, min_nonNA=min_nonNA_train)
pred_prot_m_1  <- cv_cached("matched_protein_1uM",  prot_m[, common_1um],  d1$mat[common_1um, keep1],    k=cv_folds, removeLowVaryingGenes=0, min_nonNA=min_nonNA_train)

# ------------- OPTIONAL tail imputation (Reviewer point 2) 1 + 10 uM ---------
impute_tail <- function(mat, width = 0.3, downshift = 1.8) {
  apply(mat, 2, function(col) {
    miss <- is.na(col); if (!any(miss)) return(col)
    mu <- mean(col, na.rm = TRUE); sdv <- sd(col, na.rm = TRUE)
    col[miss] <- rnorm(sum(miss), mu - downshift * sdv, width * sdv); col
  })
}
pred_tail_10 <- pred_tail_1 <- NULL
if (!is.na(preimpute_protein_path) && file.exists(preimpute_protein_path)) {
  pre <- readr::read_csv(preimpute_protein_path, show_col_types = FALSE)
  pre_sym <- if (protein_symbol_col %in% colnames(pre)) pre[[protein_symbol_col]] else sub("^.*_", "", pre[[protein_id_col]])
  pre_mat <- build_expr(pre, pre_sym, colnames(pre)[(protein_meta_cols + 1):ncol(pre)])
  tail_mat <- impute_tail(pre_mat)
  tail_expr <- build_expr(as.data.frame(tail_mat, check.names = FALSE), rownames(tail_mat), colnames(tail_mat), top_n_features)
  ct10 <- intersect(colnames(tail_expr), common_10um); ct1 <- intersect(colnames(tail_expr), common_1um)
  pred_tail_10 <- cv_cached("protein_tail_10uM", tail_expr[, ct10], d10$mat[ct10, keep10], k=cv_folds, min_nonNA=min_nonNA_train)
  pred_tail_1  <- cv_cached("protein_tail_1uM",  tail_expr[, ct1],  d1$mat[ct1, keep1],    k=cv_folds, min_nonNA=min_nonNA_train)
} else message("[skip] tail-imputation sensitivity (set preimpute_protein_path).")

# ------------------------------ DERIVED TABLES -------------------------------
obs_long <- bind_rows(
  per_observation(pred_rna_10, d10$mat, "RNA", "10uM"), per_observation(pred_prot_10, d10$mat, "Protein", "10uM"),
  per_observation(pred_rna_1,  d1$mat,  "RNA", "1uM"),  per_observation(pred_prot_1,  d1$mat,  "Protein", "1uM"))
write.csv(obs_long, file.path(outdir, "per_observation_long.csv"), row.names = FALSE)
headline_r <- obs_long %>% group_by(Modality, Conc) %>% summarise(r = cor(Predicted, Observed), n = n(), .groups = "drop")
write.csv(headline_r, file.path(outdir, "headline_r.csv"), row.names = FALSE)
message("Per-observation r:"); print(as.data.frame(headline_r))

per_drug_all <- bind_rows(
  per_drug_metrics(pred_rna_10, d10$mat, "RNA", "10uM"), per_drug_metrics(pred_prot_10, d10$mat, "Protein", "10uM"),
  per_drug_metrics(pred_rna_1,  d1$mat,  "RNA", "1uM"),  per_drug_metrics(pred_prot_1,  d1$mat,  "Protein", "1uM"))
write.csv(per_drug_all, file.path(outdir, "per_drug_metrics.csv"), row.names = FALSE)

matched <- bind_rows(
  per_drug_metrics(pred_rna_m_10, d10$mat, "RNA", "10uM"), per_drug_metrics(pred_prot_m_10, d10$mat, "Protein", "10uM"),
  per_drug_metrics(pred_rna_m_1,  d1$mat,  "RNA", "1uM"),  per_drug_metrics(pred_prot_m_1,  d1$mat,  "Protein", "1uM"))
matched$n_matched_features <- length(common_genes)
write.csv(matched, file.path(outdir, "matched_per_drug_metrics.csv"), row.names = FALSE)
write.csv(filter(matched, Conc == "10uM"), file.path(outdir, "matched_per_drug_metrics_10uM.csv"), row.names = FALSE)  # legacy name
matched_obs <- bind_rows(
  per_observation(pred_rna_m_10, d10$mat, "RNA", "10uM"), per_observation(pred_prot_m_10, d10$mat, "Protein", "10uM"),
  per_observation(pred_rna_m_1,  d1$mat,  "RNA", "1uM"),  per_observation(pred_prot_m_1,  d1$mat,  "Protein", "1uM"))
write.csv(matched_obs, file.path(outdir, "matched_per_observation.csv"), row.names = FALSE)
message("Matched-feature per-observation r:")
print(as.data.frame(matched_obs %>% group_by(Modality, Conc) %>% summarise(r = cor(Predicted, Observed), n = n(), .groups = "drop")))

if (!is.null(pred_tail_10))
  write.csv(bind_rows(
    per_drug_metrics(pred_prot_10, d10$mat, "Protein_KNNk2", "10uM"), per_drug_metrics(pred_tail_10, d10$mat, "Protein_tail", "10uM"),
    per_drug_metrics(pred_prot_1,  d1$mat,  "Protein_KNNk2", "1uM"),  per_drug_metrics(pred_tail_1,  d1$mat,  "Protein_tail", "1uM")),
    file.path(outdir, "imputation_sensitivity_per_drug.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("DONE. Tables written to: ", outdir)