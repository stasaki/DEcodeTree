**DEcodeTree** is a machine learning model designed to predict differential expression statistics from gene regulatory features.

### 📜 Link to the Paper
Read our full research paper for more details:  
[**The YTHDF Proteins Shape the Brain Gene Signatures of Alzheimer’s Disease**](https://www.biorxiv.org/content/10.1101/2024.10.23.619425v1)

### Installing
To install the DEcodeTree package, use the following commands:
``` r
install.packages("devtools")
devtools::install_github("stasaki/DEcodeTree")
```

### Running the model
```r
library(DEcodeTree)

# Load differential expression statistics for RNA-seq from DLPFC regions against cognitive decline in ROSMAP cohorts
deg_stat_data <- readRDS(system.file("data", "DLPFC_CogDec_DEG_stats.rds", package = "DEcodeTree"))

# Run the DEcodeTree model using cell type information as predictors
runDEcodeTreeModel(
  outcome = "tstat",
  deg_stat = deg_stat_data,
  feature_type = c("cell_type"),
  out_dir = "./results/",
  n_splits = 5,
  metric = c("COR", "RMSE"),
  num_samples = 10,
  num_cpus = 4
)
```


---

#### Parameters

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

#### Returns

Outputs saved to the specified directory include:

- **Evaluation Results**: Cross-validation metrics (e.g., COR, RMSE).
- **SHAP Scores**: Feature importance scores for model interpretability.
- **Plots**: Visualizations of evaluation metrics and SHAP scores.

---

### Outputs

1. **Evaluation Results**: Summary of model performance metrics (e.g., RMSE, correlation).
2. **SHAP Values**: Feature importance values for interpretability.
3. **Plots**: Visualizations of model performance and SHAP scores.
4. **Trained Models**: Serialized XGBoost models for future use.

---

### License
The code for this project is licensed under the [BSD 3-Clause License](LICENSE). Gene regulatory features were prepared using the following databases. Please ensure compliance with the respective licenses for each data source:
- **Transcription factor binding targets**: [GTRD](https://doi.org/10.1093/nar/gky1128)  
- **RNA-binding protein binding targets**: [POSTAR2](https://doi.org/10.1093/nar/gky830)  
- **miRNA binding targets**: [TargetScan](https://doi.org/10.7554/eLife.05005)  




