# DEcodeTree Analysis with GEO Data

**DEcodeTree** is a machine learning model designed to predict differential expression statistics from gene regulatory features.

### 📜 Link to the Paper
Read our full research paper for more details:  
[**The YTHDF Proteins Shape the Brain Gene Signatures of Alzheimer's Disease**](https://www.biorxiv.org/content/10.1101/2024.10.23.619425v1)

## Complete Workflow: From GEO Data to DEcodeTree Analysis

This example demonstrates how to download experimental data from GEO and analyze it with DEcodeTree. We'll use GSE122069, which contains RNA-seq data from TDP-43 knockdown experiments in neuronal cell lines. TDP-43 (TAR DNA-binding protein 43) is a critical RNA-binding protein involved in RNA processing and transport, and its dysfunction is implicated in neurodegenerative diseases like ALS and frontotemporal dementia.

> **⚠️ Computational Requirements**: This analysis uses 16 CPU cores and can take anywhere from 1 hour to several hours to complete, depending on your CPU performance and availability. Plan accordingly and ensure your system has sufficient computational resources.

### Step 1: Install Required Packages
```r
# Install required packages
install.packages(c("devtools", "tidyverse", "data.table", "limma"))
devtools::install_github("stasaki/DEcodeTree")
```

### Step 2: Download and Load Data from GEO
```r
library(tidyverse)
library(data.table)
library(limma)
library(DEcodeTree)

# Download RNA-seq data from GEO
ftp_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122069/suppl/GSE122069%5FRSEM.gene.TPM.txt.gz"
destfile <- "GSE122069_RSEM.gene.TPM.txt.gz"
download.file(ftp_url, destfile, mode = "wb")
```

### Step 3: Process Expression Data
```r
# Load and process the expression matrix
expr <- fread("GSE122069_RSEM.gene.TPM.txt.gz") %>%
  mutate(V1 = gsub("\\|.+", "", V1)) %>%  # Remove gene symbols, keep only ENSEMBL IDs
  as.data.frame() %>%
  column_to_rownames("V1") %>%
  as.matrix()

# Filter for ShySy cell line samples
indx <- grepl("ShySy", colnames(expr))
expr_sub <- expr[, indx]

# Filter lowly expressed genes (keep genes with mean TPM > 1)
cat("Number of genes with mean TPM > 1:", sum(rowMeans(expr_sub) > 1), "\n")
indx <- rowMeans(expr_sub) > 1
expr_sub <- expr_sub[indx, ]

# Log2 transform
expr_sub <- log2(expr_sub + 0.5)
```

### Step 4: Perform Differential Expression Analysis
```r
# Set up experimental design (3 WT vs 3 KO samples)
design <- model.matrix(~ group, 
                      data = data.frame(group = factor(c(rep("WT", 3), rep("KO", 3)), 
                                                      levels = c("WT", "KO"))))

# Fit linear model using limma
fit <- lmFit(expr_sub, design)
fit <- eBayes(fit)

# Extract differential expression statistics
deg_stat_data <- topTable(fit, coef = "groupKO", number = Inf, adjust.method = "BH") %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  dplyr::select(gene_id, logFC, tstat = t)

cat("Generated DE statistics for", nrow(deg_stat_data), "genes\n")
head(deg_stat_data)
```

### Step 5: Run DEcodeTree Model

> **⏱️ Runtime Warning**: The following analysis uses 16 CPU cores and may take 1-4 hours to complete depending on your system's performance. The model performs extensive hyperparameter tuning (64 iterations) and cross-validation across multiple feature types.

```r
# Run DEcodeTree with comprehensive gene regulatory features
runDEcodeTreeModel(
  outcome = "tstat",                           # Predict t-statistics
  deg_stat = deg_stat_data,                   # Our DE results
  feature_type = c("base", "RBP", "miRNA", "TF"), # All regulatory features
  out_dir = "./results/",                     # Output directory
  n_splits = 5,                               # Cross-validation splits
  metric = c("COR", "RMSE"),                  # Performance metrics
  num_samples = 64,                           # Hyperparameter tuning iterations
  num_cpus = 16                               # Parallel processing cores
)
```

---

## Parameter Reference

| Parameter       | Description                                                                                     | Default          |
|------------------|-------------------------------------------------------------------------------------------------|------------------|
| `outcome`       | A character string specifying the outcome variable (e.g., "logFC", "tstat").                   | `"tstat"`       |
| `deg_stat`      | A data frame containing differential expression statistics.                                    | Required         |
| `feature_type`  | Features to include (`"base"`, `"RBP"`, `"miRNA"`, `"TF"`, `"cell_type"`).                     | `c("base")`     |
| `out_dir`       | Output directory for results.                                                                  | `"./results/"`  |
| `n_splits`      | Number of data splits for cross-validation.                                                    | `5`             |
| `metric`        | Performance metrics (`"COR"`, `"RMSE"`).                                                       | `c("COR", "RMSE")` |
| `measures`      | List of measures for performance evaluation compatible with `mlr`.                             | Default measures |
| `cv_folds`      | Number of folds for cross-validation.                                                          | `5`             |
| `num_samples`   | Number of hyperparameter tuning iterations.                                                    | `64`            |
| `num_cpus`      | Number of CPUs for parallel processing.                                                        | `4`             |
| `param_values`  | Default hyperparameter values for XGBoost models.                                              | Predefined list |
| `param_space`   | Parameter set defining the hyperparameter optimization space.                                  | Predefined set  |
| `save_model`    | Whether to save trained models.                                                                | `TRUE`          |

---

## Outputs

1. **Evaluation Results**: Summary of model performance metrics (e.g., RMSE, correlation).
2. **SHAP Values**: Feature importance values for interpretability.
3. **Plots**: Visualizations of model performance and SHAP scores.
4. **Trained Models**: Serialized XGBoost models for future use.

---

## License
The code for this project is licensed under the [BSD 3-Clause License](LICENSE). Gene regulatory features were prepared using the following databases. Please ensure compliance with the respective licenses for each data source:
- **Transcription factor binding targets**: [GTRD](https://doi.org/10.1093/nar/gky1128)  
- **RNA-binding protein binding targets**: [POSTAR2](https://doi.org/10.1093/nar/gky830)  
- **miRNA binding targets**: [TargetScan](https://doi.org/10.7554/eLife.05005)