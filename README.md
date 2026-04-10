# SampleSim

## Overview

**SampleSim** is an R-based simulation framework for generating synthetic protein datasets under a controlled **target-versus-non-target** design. The repository provides **two user options** built on the same simulation logic:

1. **Self-contained simple simulator** — all parameters are defined directly inside the script. No input file is required.
2. **CSV-driven simulator** — parameters are read from a CSV file, allowing multiple independent simulation scenarios to be executed in a single run.

In both versions, a fixed number of proteins is simulated across two sample groups. Most proteins follow a common background distribution, whereas a user-defined subset is intentionally shifted in the target group. This design creates transparent synthetic datasets in which the true signal structure is known in advance.

Such datasets are useful for benchmarking analytical methods, testing normalization strategies, validating feature-selection procedures, evaluating machine-learning workflows, and demonstrating computational pipelines under controlled conditions.

---

## Two user options

### Option 1: Self-contained simple simulator

This version is intended for users who want a quick simulation without preparing any external input file. All settings are edited directly inside the R script.

It is most useful when you want to:

- run a single simulation immediately,
- test the framework quickly,
- adjust one scenario manually, or
- generate a teaching or demonstration dataset with minimal setup.

#### Editing the parameters

Open the script and locate the **User-defined settings** block near the top:

```r
set.seed(123)

# -----------------------------
# User-defined settings
# -----------------------------
workdir <- "path/to/your/output/folder"   # <-- change this

n_target_samples    <- 5
n_other_samples     <- 140
n_shifted_proteins  <- 3
n_total_proteins    <- 64

others_mean      <- 2.96
others_sd        <- 0.94
shift_multiplier <- 2.28
```

| Parameter | Description |
|---|---|
| `workdir` | Path to the folder where all output files will be saved. The folder is created automatically if it does not exist. |
| `n_target_samples`  | Number of samples that belong to the **target group** (labeled `Y`). These are the samples in which shifted proteins will show elevated values. |
| `n_other_samples` | Number of samples that belong to the **non-target group** (labeled `N`). These samples draw all protein values from the background distribution. |
| `n_shifted_proteins` | How many proteins should be **intentionally elevated** in the target group. These are always the first `n` proteins (e.g., `Protein_1`, `Protein_2`, ...). Setting this to `0` produces a null dataset with no true signal. |
| `n_total_proteins` | Total number of proteins to simulate. Must be greater than or equal to `n_shifted_proteins`. The remaining proteins after the shifted set all follow the background distribution in every sample. |
| `others_mean` | Mean of the **background (non-target) distribution**, applied to all non-shifted proteins in all samples and to shifted proteins in non-target samples. Values are on a log₁₀ scale. |
| `others_sd` | Standard deviation of the background distribution. Controls how much natural variation exists across all samples and proteins. |
| `shift_multiplier` | Controls the magnitude of the signal in shifted proteins within target samples. The target mean is computed as `(shift_multiplier × others_mean) + others_mean`. A value of `0` produces no shift; larger values produce a stronger, more detectable signal. |

The output file names are defined just below and can also be changed if desired.

#### Running the script

From the R console or RStudio:

```r
source("path/to/simulate_simple.R")
```

From the command line:

```bash
Rscript simulate_simple.R
```

---

### Option 2: CSV-driven simulator

This version is intended for users who wish to run multiple simulations automatically from a parameter table. Each row of the CSV defines one independent simulation scenario, and all scenarios are executed in a single run.

It is most useful when you want to:

- generate multiple datasets in one run,
- vary sample size or signal strength systematically across scenarios,
- conduct benchmarking studies, or
- manage simulation settings in a reproducible and scalable way.

#### Preparing the parameter CSV

Create a CSV file with one row per scenario. The following columns are required:

| Column | Type | Description |
|---|---|---|
| `Label` | string | A short unique identifier for this scenario. Used as the filename prefix for all four output files (e.g., `Scenario_A.csv`, `Scenario_AM.csv`). |
| `n_target_samples` | integer | Number of target samples (labeled `Y`). See description in Option 1 above. |
| `n_shifted_proteins` | integer | Number of proteins elevated in target samples. Must not exceed `n_total_proteins`. See Option 1 above. |
| `others_mean` | numeric | Mean of the background distribution (log₁₀ scale). See Option 1 above. |
| `others_sd` | numeric | Standard deviation of the background distribution. See Option 1 above. |
| `shift_multiplier` | numeric | Magnitude of signal in shifted proteins. See Option 1 above. |

The number of non-target samples (`n_other_samples`) and total proteins (`n_total_proteins`) are fixed inside the script itself and apply to all rows.

An example CSV of input csv :

```
Label,n_target_samples,n_shifted_proteins,others_mean,others_sd,shift_multiplier

<img width="478" height="350" alt="image" src="https://github.com/user-attachments/assets/4e6ed4a9-29bc-4724-8186-6caf75c1e22c" />

```

An example file is also provided in this repository.

#### Editing the script

Open the CSV-driven script and set the two paths near the top:

```r
workdir   <- "path/to/your/output/folder"
param_csv <- file.path(workdir, "your_parameters.csv")
```

If needed, also update `n_other_samples` and `n_total_proteins`, which are fixed for all scenarios in this version.

#### Running the script

From the R console or RStudio:

```r
source("path/to/simulate_csv.R")
```

From the command line:

```bash
Rscript simulate_csv.R
```

The script will process every row in the CSV sequentially and print a summary to the console after each scenario completes.

---

## Output files

Both versions produce four output files per simulation run. In the simple simulator the file names are defined directly in the script. In the CSV-driven simulator they are derived automatically from the `Label` column.

| File suffix | Contents |
|---|---|
| *(base)* `.csv` | Raw simulated protein values on the original (log₁₀) scale, one row per sample, one column per protein. |
| `M.csv` | The same data after **row-wise min-max normalization**: each value is rescaled to [0, 1] relative to the minimum and maximum across all proteins within that sample. |
| `Stat.csv` | Simulation metadata: one row per protein indicating whether it is a shifted protein, the background distribution parameters, the shift multiplier, and the target distribution mean and SD used for shifted proteins. |
| `Map.csv` | Sample-to-group mapping: one row per sample with its name and group label (`Y` for target, `N` for non-target). |

### Raw output (base `.csv`)

Samples are named `Sample_1`, `Sample_2`, ... and proteins are named `Protein_1`, `Protein_2`, ...

The first `n_target_samples` rows correspond to the target group. The first `n_shifted_proteins` columns correspond to the proteins that were intentionally elevated. All values are on a log₁₀ scale.

### Min-max normalized output (`M.csv`)

Row-wise normalization is applied independently for each sample:

```
normalized_value = (x − row_min) / (row_max − row_min)
```

This rescales every sample's protein profile to the [0, 1] interval based on that sample's own range. The structure of the data (column names, row order) is otherwise identical to the raw file.

### Statistics file (`Stat.csv`)

Contains one row per protein with the following fields:

- `Is_Shifted_Protein` — `TRUE` if this protein was shifted in target samples, `FALSE` otherwise.
- `Others_Mean` / `Others_SD` — the background distribution parameters used for all non-shifted proteins and for shifted proteins in non-target samples.
- `Shift_Multiplier` — the multiplier used to derive the target distribution mean.
- `Target_Mean_For_ShiftedProteins` / `Target_SD_For_ShiftedProteins` — the distribution parameters actually used to draw values for shifted proteins in target samples. These are `NA` for non-shifted proteins.

### Sample map (`Map.csv`)

A two-column table with `Sample` (e.g., `Sample_1`) and `Target` (`Y` or `N`). Useful for downstream analysis steps that require group labels as a separate file.

---

## Shared simulation logic

Both options implement the same core statistical design.

**Sample groups.** Samples are divided into target (`Y`) and non-target (`N`). Target samples are listed first, followed by non-target samples.

**Protein panel.** A fixed total number of proteins is simulated. The first `n_shifted_proteins` are designated as shifted proteins; all remaining proteins are non-shifted.

**Background distribution.** All non-shifted proteins in all samples, and shifted proteins in non-target samples, are drawn from:

```
N(others_mean, others_sd)
```

**Target distribution.** Shifted proteins in target samples are drawn from:

```
N(target_mean, others_sd)
where target_mean = (shift_multiplier × others_mean) + others_mean
```

The standard deviation is the same in both distributions; only the mean is shifted.

**Reproducibility.** Both scripts call `set.seed(123)` at the top. Results are fully reproducible as long as the seed and parameters are unchanged.

---

## Reference

[1] Prabin Dawadi, Ryan M Tobin, Jorge Frias-Lopez, Alpdogan Kantarci, Sayaka Miura, Flavia Teles. Uncovering Periodontal Ecosystem Complexity with Sample Trees. (2026) Under Review.

---

## License
Copyright 2025, Authors and University of Mississippi.
BSD 3-Clause "New" or "Revised" License. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
- Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
