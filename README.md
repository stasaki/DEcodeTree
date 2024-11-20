**DEcodeTree** is a machine learning model designed to predict differential expression statistics from gene regulatory features.

## 📜 Link to the Paper
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

## License
The code for this project is licensed under the [BSD 3-Clause License](LICENSE). Gene regulatory features were prepared using the following databases. Please ensure compliance with the respective licenses for each data source:
- **Transcription factor binding targets**: [GTRD](https://doi.org/10.1093/nar/gky1128)  
- **RNA-binding protein binding targets**: [POSTAR2](https://doi.org/10.1093/nar/gky830)  
- **miRNA binding targets**: [TargetScan](https://doi.org/10.7554/eLife.05005)  
