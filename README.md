# SampleSim

## Overview

**SampleSim** is an R-based simulation framework for generating synthetic protein datasets under a controlled **target-versus-non-target** design. The repository provides **two user options** built on the same simulation logic:

1. **Self-contained simple simulator**  
   All parameters are defined directly inside the script. No input file is required.

2. **CSV-driven simulator**  
   Parameters are read from a CSV file, allowing multiple independent simulation scenarios to be executed in a single run.

In both versions, a fixed number of proteins is simulated across two sample groups. Most proteins follow a common background distribution, whereas a user-defined subset is intentionally shifted in the target group. This design creates transparent synthetic datasets in which the true signal structure is known in advance.

Such datasets are useful for benchmarking analytical methods, testing normalization strategies, validating feature-selection procedures, evaluating machine-learning workflows, and demonstrating computational pipelines under controlled conditions.

---

## Usage

In real biological datasets, the true signal is rarely known with certainty. By contrast, simulation allows the investigator to define exactly:

- how many samples belong to the target group,
- how many proteins should carry signal,
- how strong that signal should be, and
- how much background variation should be present.

SampleSim was written for precisely this purpose: to generate simple, reproducible synthetic datasets in which **shifted** and **non-shifted** proteins can be known beforehand.

---

## Two user options

### Option 1: Self-contained simple simulator

This version is intended for users who want a quick simulation without preparing any external input file. All settings are edited directly inside the R script.

It is most useful when one wishes to:

- run a single simulation immediately,
- test the framework quickly,
- adjust one scenario manually, or
- generate a teaching or demonstration dataset with minimal setup.

Typical parameters edited directly in the script include:

- number of target samples,
- number of non-target samples,
- number of shifted proteins,
- total number of proteins,
- background mean,
- background standard deviation,
- shift multiplier,
- working directory, and
- output file names.

### Option 2: CSV-driven simulator

This version is intended for users who wish to run multiple simulations automatically from a parameter table. Each row of the CSV defines one independent simulation scenario.

It is most useful when one wishes to:

- generate multiple datasets in one run,
- vary sample size or signal strength across scenarios,
- conduct systematic benchmarking studies, or
- manage simulation settings in a reproducible and scalable way.

In this version, the value in the `Label` column is used as the prefix for all output files associated with that simulation run.

Note: Example csv file that could be used is in the SampleSim repository. 
---

## Shared simulation logic

Although the two versions differ in how parameters are supplied, both implement the same core statistical design.

### 1. Construct the sample groups

Samples are divided into two classes:

- **Target samples (`Y`)**
- **Non-target samples (`N`)**

Sample names are generated automatically as:

- `Sample_1`
- `Sample_2`
- `Sample_3`
- ...

### 2. Construct the protein panel

A fixed total number of proteins is simulated, named automatically as:

- `Protein_1`
- `Protein_2`
- `Protein_3`
- ...

The first `n_shifted_proteins` are designated as **shifted proteins**.  
All remaining proteins are treated as **non-shifted proteins**.

### 3. Simulate the protein values

All non-shifted proteins are sampled from the same background distribution in all samples

## Reference
[1] Prabin Dawadi, Ryan M Tobin, Jorge Frias-Lopez, Alpdogan Kantarci, Flavia Teles, Sayaka Miura. Uncovering Periodontal Ecosystem Complexity with Sample Trees. (2026) Under Review

Copyright 2025, Authors and University of Mississippi
BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
