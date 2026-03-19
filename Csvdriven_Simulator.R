set.seed(123)

# -----------------------------
# User-defined settings
# -----------------------------
workdir <- *********************"  ## Change working dir
param_csv <- file.path(workdir, "*********.csv")   # Change csv file

n_other_samples    <- 140
n_total_proteins   <- 64
sample_col         <- "Sample"

# -----------------------------
# Checks for working directory / input file
# -----------------------------
if (!dir.exists(workdir)) {
  stop("Working directory does not exist: ", workdir)
}

if (!file.exists(param_csv)) {
  stop("Parameter CSV not found: ", param_csv)
}

# -----------------------------
# Read parameter CSV
# -----------------------------
param_df <- read.csv(param_csv, stringsAsFactors = FALSE, check.names = FALSE)

# Trim spaces from column names
colnames(param_df) <- trimws(colnames(param_df))

# Check required columns
required_cols <- c(
  "Label",
  "n_target_samples",
  "n_shifted_proteins",
  "others_sd",
  "others_mean",
  "shift_multiplier"
)

missing_cols <- setdiff(required_cols, colnames(param_df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in CSV: ", paste(missing_cols, collapse = ", "))
}

# -----------------------------
# Loop through each row of CSV
# -----------------------------
for (r in seq_len(nrow(param_df))) {
  
  # -----------------------------
  # Read one parameter set
  # -----------------------------
  label_value         <- as.character(param_df$Label[r])
  n_target_samples    <- as.integer(param_df$n_target_samples[r])
  n_shifted_proteins  <- as.integer(param_df$n_shifted_proteins[r])
  others_sd           <- as.numeric(param_df$others_sd[r])
  others_mean         <- as.numeric(param_df$others_mean[r])
  shift_multiplier    <- as.numeric(param_df$shift_multiplier[r])
  
  # Derived values
  target_mean <- (shift_multiplier * others_mean) + others_mean
  target_sd   <- others_sd
  
  # -----------------------------
  # Output files
  # -----------------------------
  file_raw_out    <- file.path(workdir, paste0(label_value, ".csv"))
  file_minmax_out <- file.path(workdir, paste0(label_value, "M.csv"))
  file_stats_out  <- file.path(workdir, paste0(label_value, "Stat.csv"))
  file_map_out    <- file.path(workdir, paste0(label_value, "Map.csv"))
  
  # -----------------------------
  # Checks
  # -----------------------------
  if (is.na(n_target_samples) || n_target_samples < 1) {
    stop("Invalid n_target_samples at row ", r)
  }
  if (n_other_samples < 1) {
    stop("n_other_samples must be >= 1")
  }
  if (n_total_proteins < 1) {
    stop("n_total_proteins must be >= 1")
  }
  if (is.na(n_shifted_proteins) || n_shifted_proteins < 0) {
    stop("Invalid n_shifted_proteins at row ", r)
  }
  if (n_shifted_proteins > n_total_proteins) {
    stop("n_shifted_proteins cannot be greater than n_total_proteins at row ", r)
  }
  if (is.na(others_mean) || is.na(others_sd) || is.na(shift_multiplier)) {
    stop("Invalid numeric parameters at row ", r)
  }
  
  # -----------------------------
  # Generate sample names
  # -----------------------------
  n_total_samples <- n_target_samples + n_other_samples
  sample_names <- paste0("Sample_", seq_len(n_total_samples))
  
  target_labels <- c(
    rep("Y", n_target_samples),
    rep("N", n_other_samples)
  )
  
  sample_info <- data.frame(
    Sample = sample_names,
    Target = target_labels,
    stringsAsFactors = FALSE
  )
  
  is_target <- sample_info$Target == "Y"
  
  # -----------------------------
  # Generate protein names
  # -----------------------------
  protein_cols <- paste0("Protein_", seq_len(n_total_proteins))
  
  if (n_shifted_proteins > 0) {
    target_proteins <- protein_cols[seq_len(n_shifted_proteins)]
  } else {
    target_proteins <- character(0)
  }
  
  other_proteins <- setdiff(protein_cols, target_proteins)
  
  # -----------------------------
  # Empty output data frame
  # -----------------------------
  tree_sim <- data.frame(
    Sample = sample_names,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  n_total  <- n_total_samples
  n_target <- sum(is_target)
  n_other  <- sum(!is_target)
  
  # =========================================================
  # 1) Simulate shifted proteins
  # =========================================================
  if (n_shifted_proteins > 0) {
    for (p in target_proteins) {
      
      vals <- numeric(n_total)
      
      # non-target samples
      vals[!is_target] <- rnorm(n_other, mean = others_mean, sd = others_sd)
      
      # target samples
      vals[is_target] <- rnorm(n_target, mean = target_mean, sd = target_sd)
      
      tree_sim[[p]] <- vals
    }
  }
  
  # =========================================================
  # 2) Simulate non-shifted proteins
  # =========================================================
  if (length(other_proteins) > 0) {
    for (p in other_proteins) {
      tree_sim[[p]] <- rnorm(n_total, mean = others_mean, sd = others_sd)
    }
  }
  
  # =========================================================
  # 3) Row-wise MinMax normalization
  # =========================================================
  tree_sim_minmax <- tree_sim
  
  numeric_mat <- as.matrix(tree_sim[, protein_cols, drop = FALSE])
  
  row_min <- apply(numeric_mat, 1, min, na.rm = TRUE)
  row_max <- apply(numeric_mat, 1, max, na.rm = TRUE)
  row_range <- row_max - row_min
  
  scaled_mat <- matrix(NA_real_, nrow = nrow(numeric_mat), ncol = ncol(numeric_mat))
  colnames(scaled_mat) <- colnames(numeric_mat)
  
  for (i in seq_len(nrow(numeric_mat))) {
    if (is.na(row_range[i]) || row_range[i] == 0) {
      scaled_mat[i, ] <- 0
    } else {
      scaled_mat[i, ] <- (numeric_mat[i, ] - row_min[i]) / row_range[i]
    }
  }
  
  tree_sim_minmax[, protein_cols] <- scaled_mat
  
  # =========================================================
  # 4) Save stats table
  # =========================================================
  stats_df <- data.frame(
    Label = label_value,
    Protein = protein_cols,
    Is_Shifted_Protein = protein_cols %in% target_proteins,
    n_target_samples = n_target_samples,
    n_other_samples = n_other_samples,
    n_shifted_proteins = n_shifted_proteins,
    n_total_proteins = n_total_proteins,
    Others_Mean = others_mean,
    Others_SD = others_sd,
    Shift_Multiplier = shift_multiplier,
    Target_Mean_For_ShiftedProteins = ifelse(
      protein_cols %in% target_proteins,
      target_mean,
      NA
    ),
    Target_SD_For_ShiftedProteins = ifelse(
      protein_cols %in% target_proteins,
      target_sd,
      NA
    ),
    stringsAsFactors = FALSE
  )
  
  # =========================================================
  # 5) Write outputs
  # =========================================================
  write.csv(tree_sim, file_raw_out, row.names = FALSE)
  write.csv(tree_sim_minmax, file_minmax_out, row.names = FALSE)
  write.csv(stats_df, file_stats_out, row.names = FALSE)
  write.csv(sample_info, file_map_out, row.names = FALSE)
  
  # -----------------------------
  # Console summary per row
  # -----------------------------
  cat("\nDone for Label:", label_value, "\n")
  cat("Generated:\n")
  cat(" - ", file_raw_out, "\n", sep = "")
  cat(" - ", file_minmax_out, "\n", sep = "")
  cat(" - ", file_stats_out, "\n", sep = "")
  cat(" - ", file_map_out, "\n", sep = "")
  
  cat("Counts:\n")
  cat("Target samples: ", n_target, "\n", sep = "")
  cat("Non-target samples: ", n_other, "\n", sep = "")
  cat("Total samples: ", n_total, "\n", sep = "")
  cat("Shifted proteins: ", n_shifted_proteins, "\n", sep = "")
  cat("Total proteins: ", n_total_proteins, "\n", sep = "")
  
  cat("Simulation settings:\n")
  cat("Others distribution: N(", others_mean, ", ", others_sd, ")\n", sep = "")
  cat("Target distribution for shifted proteins in target samples: N(",
      round(target_mean, 4), ", ", target_sd, ")\n", sep = "")
}

cat("\nAll simulations completed successfully.\n")
