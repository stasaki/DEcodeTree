 
#' Run DEcode-tree Model
#'
#' This function performs model training and evaluation for DEcode-tree models.
#'
#' @param outcome A character string specifying the outcome variable (default: "tstat").
#' @param deg_stat A data frame containing differential expression statistics. Must include `gene_id` and the specified `outcome`.
#' @param feature_type A character vector specifying feature types to include (e.g., "base", "RBP", "miRNA", "TF", "cell_type").
#' @param out_dir A character string specifying the output directory (default: "./result/").
#' @param n_splits Number of data splits for cross-validation (default: 5).
#' @param metric A character vector specifying performance metrics (default: c("COR", "RMSE")).
#' @param measures A list of measures for performance evaluation, compatible with mlr (default: \code{list(mlr::rmse, mlr::setAggregation(mlr::rmse, mlr::train.mean))}).
#' @param cv_folds Number of folds for cross-validation (default: 5).
#' @param num_samples Number of hyperparameter tuning iterations (default: 64).
#' @param num_cpus Number of CPUs for parallel processing (default: 4).
#' @param param_values A list of default hyperparameter values for XGBoost models.
#' @param param_space A parameter set defining the tuning space for hyperparameter optimization.
#' @param save_model Logical, whether to save the trained models (default: TRUE).
#'
#' @details
#' This function prepares data, performs cross-validation, trains models using XGBoost, evaluates them, and computes SHAP values.
#'
#' @return The function saves various outputs, including evaluation results, SHAP scores, and plots, to the specified output directory.
#'
#' @examples
#' \dontrun{
#' deg_stat_data <- readRDS(system.file("data", "DLPFC_CogDec_DEG_stats.rds", package = "DEcodeTree"))
#' runDEcodeTreeModel(
#'   outcome = "tstat",
#'   deg_stat = deg_stat_data,
#'   feature_type = c("cell_type"),
#'   out_dir = "./results/",
#'   n_splits = 5,
#'   metric = c("COR", "RMSE"),
#'   num_samples = 10,
#'   num_cpus = 4
#' )
#' }
#'
#' @export

runDEcodeTreeModel <- function(outcome = "tstat",
                               deg_stat,
                               feature_type = c("base","RBP","miRNA","TF","cell_type"),
                               out_dir = "./",
                               n_splits = 5,
                               metric = c("COR", "RMSE"),
                               measures = list(mlr::rmse, mlr::setAggregation(mlr::rmse, mlr::train.mean)),
                               cv_folds = 5,
                               num_samples = 1,
                               num_cpus = 1,
                               param_values = list(
                                 objective = "reg:linear",
                                 eval_metric = "rmse",
                                 nrounds = 100L,
                                 eta = 0.1,
                                 booster = "gbtree",
                                 nthread = 4,
                                 alpha = 0.5
                               ),
                               param_space = ParamHelpers::makeParamSet(
                                 ParamHelpers::makeIntegerParam("max_depth", lower = 3L, upper = 10L),
                                 ParamHelpers::makeNumericParam("min_child_weight", lower = 1L, upper = 10L),
                                 ParamHelpers::makeNumericParam("subsample", lower = 0.5, upper = 1),
                                 ParamHelpers::makeNumericParam("colsample_bytree", lower = 0.5, upper = 1),
                                 ParamHelpers::makeNumericParam("lambda", lower = -3, upper = 3, trafo = function(x) 10^x),
                                 ParamHelpers::makeNumericParam("gamma", lower = -3, upper = 0, trafo = function(x) 10^x)
                               ),
                               save_model = TRUE) {
  
  message("Saving configuration parameters...")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  # Save configuration parameters
  save(list = c("outcome", "out_dir", "n_splits", "metric",
                "cv_folds", "num_samples", "num_cpus", "param_values", "param_space", "save_model"),
       file = paste0(out_dir, "/config.RData"))
  
  # Define file paths
  model_score_file <- paste0(out_dir, "/summary/model.score.rds")
  shap_score_file <- paste0(out_dir, "/summary/shap.score.rds")
  if (file.exists(model_score_file) & file.exists(shap_score_file)) {
    message("Output files already exist. Skipping computation.")
    return(NULL)
  }
  
  message("Preparing data...")
  data = prep_data(deg_stat,outcome, feature_type)
    
  # Import required library
  library(caret)
  
  message("Checking if data splits already exist...")
  # Check if split files already exist
  split_exists <- TRUE
  for (i in 1:n_splits) {
    model_dir <- paste0(out_dir, "submodels/model", i, "/")
    split_exists <- split_exists & file.exists(paste0(model_dir, "inTrain.rds"))
  }
  
  # Create splits if not already created
  if (!split_exists) {
    message("Creating data splits...")
    set.seed(1234)
    fold_indices <- createMultiFolds(data[[outcome]], k = n_splits, times = 1)
    for (i in 1:length(fold_indices)) {
      model_dir <- paste0(out_dir, "submodels/model", i, "/")
      dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)
      inTrain <- fold_indices[[i]] %>% t() %>% t()
      saveRDS(inTrain, file = paste0(model_dir, "inTrain.rds"))
    }
  }
  
  message("Starting model building and evaluation...")
  # Apply model building and evaluation to each split
  lapply(1:n_splits, function(i) {
    message(paste("Building and evaluating model for split", i, "..."))
    buildAndEvaluateModel(
      outcome = outcome,
      out_dir = out_dir,
      data = data,
      metric = metric,
      measures = measures,
      cv_folds = cv_folds,
      num_samples = num_samples,
      num_cpus = num_cpus,
      param_values = param_values,
      param_space = param_space,
      save_model = save_model,
      n = i
    )
  })
  
  library(tidyverse)
  library(ggpubr)
  library(viridis)
  library(RColorBrewer)
  # Only run if the model score output does not already exist
  if (!file.exists(model_score_file)) {
    message("Combining performance results from all splits...")
    # Combine performance results from all splits
    model.score <- list.files(out_dir, pattern = "evalmodel.rds", recursive = TRUE) %>%
      lapply(., function(x) {
        evalmodel <- readRDS(paste0(out_dir, x))
        return(evalmodel$score %>% t() %>% as_tibble() %>% mutate(model = dirname(x) %>% basename()))
      }) %>%
      bind_rows()
    dir.create(paste0(out_dir, "/summary"), showWarnings = FALSE)
    saveRDS(model.score, file = model_score_file)
    
    print(model.score)
    
    test_prediction<- list.files(out_dir, pattern = "evalmodel.rds", recursive = TRUE,full.names = T) %>%
      lapply(., function(x) {
       readRDS(x)$predictions%>%
          mutate(model = dirname(x)%>%basename())%>%
          return()
      }) %>%
      bind_rows()
    test_prediction%>%
      ggplot(.,aes(predicted,actual))+
      geom_hex()+
      stat_cor(size=2)+
      theme_pubclean(base_size = 7,base_family = "Helvetica")+
      scale_fill_viridis( option = "D",
                          guide = guide_colorbar(
                            direction = "vertical",
                            title.position = "top"
                          ))+
      theme(legend.key.size = unit(2,units = "mm"),
            legend.key.height = unit(4,units = "mm"),
            legend.position  = c(0.8,0.2),
            legend.text=element_text(size=4),
            legend.title=element_text(size=5)) +
      ylab("Actual")+
      xlab("Prediction")+
      facet_wrap(~model)->p
  #  ggsave(p,file=paste0(out_dir,"/summary/prediction_vs_asctual.pdf"),width = 3,height = 3)
    ggsave(p,file=paste0(out_dir,"/summary/prediction_vs_asctual.png"),width = 4,height = 3)
  }
  
  # Only run if the SHAP score output does not already exist
  if (!file.exists(shap_score_file)) {
    message("Combining SHAP scores from all splits...")
    # Combine SHAP scores from all splits
    shap.score <- list.files(path = out_dir,
                             pattern = "shap.rds", recursive = TRUE, full.names = TRUE) %>%
      lapply(., function(x) {
        readRDS(x) %>% return()
      }) %>%
      bind_rows() %>%
      mutate(SHAP = replace(SHAP, is.na(SHAP), 0),
             Correlation = replace(Correlation, is.na(Correlation), 0)) %>%
      group_by(Variable) %>%
      summarise(SHAP = mean(SHAP),
                Correlation = mean(Correlation)) %>%
      arrange(desc(SHAP)) %>%
      mutate(Scaled_SHAP = scale(SHAP * sign(Correlation), center = FALSE, scale = TRUE) %>% as.numeric() %>% abs())
    
    saveRDS(shap.score, file = shap_score_file)
    
    
    shap.score%>%
      mutate(type="Base",
             type=replace(type,Variable%in%c("neurons","oligodendrocytes","astrocytes","OPC","microglia","endothelial"),"Cell type"),
             type=replace(type,grepl("RBP_",Variable),"RBP"),
             type=replace(type,grepl("TF_",Variable),"TF"),
             type=replace(type,grepl("miRNA_",Variable),"miRNA"),
             type=factor(type,levels=c("Base","Cell type","TF","miRNA","RBP")))%>%
      mutate(Variable=gsub("_",":",Variable))%>%
      ungroup()%>%
      arrange(Scaled_SHAP)%>%
      mutate(Variable=factor(Variable,levels=Variable))%>%
      arrange(desc(Scaled_SHAP))%>%.[1:min(20,nrow(.)),]%>%#.$type
      ggplot(.,aes(Variable,Scaled_SHAP,fill=type))+
      geom_bar(stat="identity")+coord_flip()+
      theme_pubclean(base_size = 7,base_family = "Helvetica")+
      scale_fill_manual(name="Type",values = brewer.pal(6, "Dark2")[-6],drop=FALSE)+
      theme(legend.key.size = unit(3,units = "mm"))+
      ylab("Normalized predictive contribution")+
      xlab("Input variable")+
      theme(legend.key.size = unit(2,units = "mm"),
            legend.key.height = unit(1,units = "mm"),
            legend.position  = c(0.7,0.2),
            legend.text=element_text(size=6))->p
    
    #ggsave(paste0(out_dir,"/summary/shap_bar.pdf"),p,width = 2.5,height = 0.5+2.5*p$data%>%nrow/20)
    ggsave(paste0(out_dir,"/summary/shap_bar.png"),p,width = 2.5,height = 0.5+2.5*p$data%>%nrow/20)
  }
  message("Model computation completed.")
}

prep_data = function(deg_stat,outcome, feature_type){
  
  # Check if deg_stat has the required columns
  required_columns <- c("gene_id", outcome)
  if (!all(required_columns %in% colnames(deg_stat))) {
    stop(paste("deg_stat must contain the following columns:", paste(required_columns, collapse = ", ")))
  }
  
  feature_data <- Reduce(function(x, y) merge(x, y, all=T,
                                              by=c("gene_id")), lapply(feature_type,function(x){readRDS(system.file("data", paste0(x,".rds"), package = "DEcodeTree"))}), accumulate=F)
  deg_stat[,c("gene_id",outcome)]%>%
    inner_join(.,feature_data,by=c("gene_id")) -> data
  # normalize features
  for (feat in intersect(colnames(data),c("intron_length","tr_length","exon_num"))){
    th=quantile(data[[feat]],prob=0.99)
    indx=data[[feat]]>th
    data[[feat]][indx]=th
    data[[feat]]=data[[feat]]/max(data[[feat]])
  }
  
  for (feat in colnames(data)%>%.[grepl("RBP_|TF_|miRNA_",.)]){
    if(sum(data[[feat]]>0)<30){
      data[[feat]]=NULL
      next()
    }
    data[[feat]]=data[[feat]]/max(data[[feat]])
  }
  
  for (feat in intersect(colnames(data),
                         c("astrocytes","endothelial","microglia",
                           "neurons","oligodendrocytes","OPC"))){
    data[[feat]]=scale(data[[feat]],center = F,scale = T)%>%as.vector()
  }
  
  
  
  data[[outcome]]=scale(data[[outcome]],center = F,scale = T)%>%as.vector()
  rm(feature_data)
  
  data$gene_id=NULL
  
  return(data)
}



buildAndEvaluateModel <- function(outcome, out_dir, data, metric, measures,
                                  cv_folds, num_samples, num_cpus, param_values, param_space, save_model = FALSE, n = 1) {
  library(caret)
  
  # Define model output directory
  model_dir <- paste0(out_dir, "submodels/model", n, "/")
  if (file.exists(paste0(model_dir, "evalmodel_train.rds"))) {
    return()
  }
  dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
  
  set.seed(n)
  
  # Load training indices
  inTrain <- readRDS(paste0(model_dir, "inTrain.rds"))
  
  # Prepare train and test datasets
  train_data <- data[inTrain, -1, drop = FALSE]
  test_data <- data[-inTrain, -1, drop = FALSE]
  train_data[[outcome]] <- data[[outcome]][inTrain]
  test_data[[outcome]] <- data[[outcome]][-inTrain]
  
  # Separate test set values and labels
  test_values <- test_data[, !colnames(test_data) %in% outcome, drop = FALSE]
  test_labels <- test_data[, outcome, drop = FALSE]
  
  # Build model
  model_fit <- xgbTreeModelBuild(
    outcome = outcome,
    data = train_data,
    param_values = param_values,
    cv_folds = cv_folds,
    num_samples = num_samples,
    num_cpus = num_cpus,
    measures = measures,
    param_space = param_space,
    output_path = model_dir
  )
  
  # Evaluate model on test data
  eval_test <- eval_model_mlr(
    model = model_fit,
    test_data = test_values,
    test_labels = test_labels[[1]],
    metric = metric,
    outcome
  )
  saveRDS(eval_test, file = paste0(model_dir, "evalmodel.rds"))
  
  # Evaluate model on train data
  eval_train <- eval_model_mlr(
    model = model_fit,
    test_data = train_data,
    test_labels = train_data[, outcome, drop = TRUE],
    metric = metric,
    outcome
  )
  saveRDS(eval_train, file = paste0(model_dir, "evalmodel_train.rds"))
  
  # Save the model if required
  if (save_model) {
    saveRDS(model_fit, file = paste0(model_dir, "model_fit.rds"))
  }
  
  # SHAP contributions
  data_matrix <- as.matrix(data[, -1, drop = FALSE])
  colnames(data_matrix) <- make.names(colnames(data_matrix))
  data_matrix <- data_matrix[, model_fit$learner.model$feature_names, drop = FALSE]
  shap_contrib <- predict(model_fit$learner.model, data_matrix, predcontrib = TRUE, approxcontrib = FALSE)
  
  # SHAP calculation function
  calc_shap <- function(shap_contrib) {
    rownames(shap_contrib) <- rownames(data)
    shap_contrib <- shap_contrib[, -which(colnames(shap_contrib) == "BIAS"), drop = FALSE]
    
    shap_data <- data.frame(
      Variable = colnames(shap_contrib),
      SHAP = colMeans(abs(shap_contrib)),
      Correlation = sapply(colnames(shap_contrib), function(var) {
        cor(data_matrix[, var], shap_contrib[, var], method = "spearman")
      })
    )
    return(shap_data)
  }
  
  # Apply SHAP calculation
  shap_values <- if (is.list(shap_contrib)) {
    lapply(shap_contrib, calc_shap)
  } else {
    calc_shap(shap_contrib)
  }
  
  saveRDS(shap_values, file = paste0(model_dir, "shap.rds"))
}




xgbTreeModelBuild <- function(outcome, data, param_values, cv_folds, num_samples,
                              num_cpus, measures, param_space, output_path) {
  library(mlr)
  library(xgboost)
  
  # Ensure column names are syntactically valid
  colnames(data) <- make.names(colnames(data))
  
  # Define the task based on the type of outcome
  if (is.factor(data[[outcome]]) && length(levels(data[[outcome]])) > 2) {
    task <- makeClassifTask(data = data, target = outcome)
    learner <- makeLearner("classif.xgboost", predict.type = "prob")
    param_values$objective <- "multi:softprob"
  } else if (is.factor(data[[outcome]]) && length(levels(data[[outcome]])) == 2) {
    task <- makeClassifTask(data = data, target = outcome)
    learner <- makeLearner("classif.xgboost", predict.type = "prob")
  } else {
    task <- makeRegrTask(data = data, target = outcome)
    learner <- makeLearner("regr.xgboost", predict.type = "response")
  }
  
  # Assign parameter values to learner
  learner$par.vals <- param_values
  
  # Set cross-validation resampling strategy
  resampling <- makeResampleDesc("CV", iters = cv_folds, predict = "both")
  
  # Set tuning control for random search
  tune_control <- makeTuneControlRandom(maxit = num_samples)
  
  # Start parallel processing
  library(parallel)
  library(parallelMap)
  parallelStartSocket(cpus = num_cpus)
  
  # Parameter tuning
  tuning_results <- mlr::tuneParams(
    learner = learner, task = task, resampling = resampling,
    measures = measures, par.set = param_space,
    control = tune_control, show.info = FALSE
  )
  parallelStop()
  
  # Save tuning results
  saveRDS(tuning_results, file = paste0(output_path, "tuning_results.rds"))
  
  # Generate and save hyperparameter effect data
  effect_data <- mlr::generateHyperParsEffectData(tuning_results, partial.dep = TRUE)
  saveRDS(effect_data, file = paste0(output_path, "effect_data.rds"))
  
  # Set learner with best parameters and train model
  learner <- mlr::setHyperPars(learner, par.vals = tuning_results$x)
  model <- mlr::train(learner, task)
  
  
  return(model)
}

eval_model_mlr <- function(model, test_data, test_labels, metric, outcome) {
  # Combine test data with outcome labels
  test_data[[outcome]] <- test_labels
  colnames(test_data) <- make.names(colnames(test_data))
  
  # Define task type based on outcome
  if (is.factor(test_data[[outcome]])) {
    model$learner$predict.type <- "prob"
    test_task <- makeClassifTask(data = test_data, target = outcome)
  } else {
    test_task <- makeRegrTask(data = test_data, target = outcome)
  }
  
  # Load rminer library for metric evaluation
  library(rminer)
  
  # Predict and evaluate the model
  predictions <- predict(model, test_task)
  prediction_results <- data.frame(
    predicted = predictions$data$response,
    actual = predictions$data$truth
  )
  rownames(prediction_results) <- rownames(predictions$data)
  
  # Calculate evaluation metrics
  score <- mmetric(predictions$data$response, predictions$data$truth, metric = metric)
  
  # Prepare output
  eval_results <- list(
    predictions = prediction_results,
    score = score
  )
  
  # Include probabilities for classification tasks
  if (is.factor(test_data[[outcome]])) {
    eval_results$probabilities <- predictions$data
    eval_results$probabilities$response <- NULL
    eval_results$probabilities$truth <- NULL
  }
  
  return(eval_results)
}

