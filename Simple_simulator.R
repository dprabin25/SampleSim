# =========================================================
# SELF-CONTAINED SIMPLE TREE-ONLY SIMULATION
#
# No input files needed
#
# You can change:
#   1. Number of target samples
#   2. Number of other samples
#   3. Number of shifted proteins
#   4. Number of total proteins
#
# Naming rules:
#   Samples  -> Sample_1, Sample_2, ...
#   Proteins -> Protein_1, Protein_2, ...
#
# Simulation rules:
#   - All non-shifted proteins: rnorm(mean = others_mean, sd = others_sd)
#   - Shifted proteins in non-target samples: same background
#   - Shifted proteins in target samples: shifted distribution
#
# Added:
#   - Row-wise MinMax normalization like your Python approach
# =========================================================

set.seed(123)

# -----------------------------
# User-defined settings
# -----------------------------
workdir <- "**************************************"

n_target_samples    <- 5
n_other_samples     <- 140
n_shifted_proteins  <- 3
n_total_proteins    <- 64

sample_col <- "Sample"

# -----------------------------
# Simulation parameters
# -----------------------------
others_mean <- 2.96
others_sd   <- 0.94
shift_multiplier <- 2.28

# Your original formula
target_mean <- (shift_multiplier * others_mean) + others_mean
target_sd   <- others_sd

# -----------------------------
# Output files
# -----------------------------
file_raw_out    <- file.path(workdir, "Pro1log10_Simulated317.csv")
file_minmax_out <- file.path(workdir, "Pro1log10MinMax_Simulated317.csv")
file_stats_out  <- file.path(workdir, "Simulation_Stats_Tree_SimpleShift317.csv")
file_map_out    <- file.path(workdir, "Sample_Target_Map317.csv")

# -----------------------------
# Checks
# -----------------------------
if (n_target_samples < 1) stop("n_target_samples must be >= 1")
if (n_other_samples < 1) stop("n_other_samples must be >= 1")
if (n_total_proteins < 1) stop("n_total_proteins must be >= 1")
if (n_shifted_proteins < 0) stop("n_shifted_proteins must be >= 0")
if (n_shifted_proteins > n_total_proteins) {
  stop("n_shifted_proteins cannot be greater than n_total_proteins")
}

# Create folder if needed
if (!dir.exists(workdir)) dir.create(workdir, recursive = TRUE)

# -----------------------------
# Generate sample names
# -----------------------------
n_total_samples <- n_target_samples + n_other_samples

sample_names <- paste0("Sample_", seq_len(n_total_samples))

# First target samples = target
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

# First n_shifted_proteins are shifted
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
    
    # non-target samples use background distribution
    vals[!is_target] <- rnorm(n_other, mean = others_mean, sd = others_sd)
    
    # target samples use shifted distribution
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
# 3) MAKE ROW-WISE MINMAX VERSION
# Same logic as your Python script:
# For each row:
#   (x - row_min) / (row_max - row_min)
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
  Protein = protein_cols,
  Is_Shifted_Protein = protein_cols %in% target_proteins,
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
# Console summary
# -----------------------------
cat("Done.\n")
cat("Generated:\n")
cat(" - ", file_raw_out, "\n", sep = "")
cat(" - ", file_minmax_out, "\n", sep = "")
cat(" - ", file_stats_out, "\n", sep = "")
cat(" - ", file_map_out, "\n", sep = "")

cat("\nCounts:\n")
cat("Target samples: ", n_target, "\n", sep = "")
cat("Non-target samples: ", n_other, "\n", sep = "")
cat("Total samples: ", n_total, "\n", sep = "")
cat("Shifted proteins: ", n_shifted_proteins, "\n", sep = "")
cat("Total proteins: ", n_total_proteins, "\n", sep = "")

cat("\nShifted proteins:\n")
if (length(target_proteins) > 0) {
  cat(paste(target_proteins, collapse = ", "), "\n")
} else {
  cat("None\n")
}

cat("\nSimulation settings:\n")
cat("Others distribution: N(", others_mean, ", ", others_sd, ")\n", sep = "")
cat("Target distribution for shifted proteins in target samples: N(",
    round(target_mean, 4), ", ", target_sd, ")\n", sep = "")
