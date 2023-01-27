#### Machine Learning PIPELINE
#### Weighted  Random forest models
###Load packages
library(optparse)
library(dplyr)
library(data.table)
library(caret)
library(ranger)
library(MLeval)
### Set seed
set.seed(10)
### Make options
option_list <- list(
  make_option(c("-a", "--train_dataset"), type="character", default="train_dataset.csv",help="Insert path to train dataset file (csv)"),
  make_option(c("-b", "--test_dataset"), type="character", default="test_dataset.csv",help="Insert path to test dataset file (csv)"),
  make_option(c("-m", "--metadata"), type="character", default="metadata.csv",help="Insert path to metadata file (csv)"),
  make_option(c("-d", "--drug"), type="character", default="rifampicin",help="Insert name of drug- column name must match to metadata"),
  make_option(c("-w", "--weight_file"), type="character", default="results/rifampicin_weights.csv",help="file containing weights (csv)"),
  make_option(c("-t", "--weight_type"), type="character", default="fitch",help="Weight you would like to use, containing two columns. First should be names Mutation.)"),
  make_option(c("-o", "--output"), type="character", default="results/rifampicin_weighted",help="prefix of file to output")
)

### Parse options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)

### Functions
## Get phenotype data
get_pheno <- function(geno){
  geno <- data.frame(geno)
  pheno <- geno$phenotype
  return(pheno)
}

## Function to remove samples with missing phenotype data
filter_missing <- function(metadata_file, geno, drug){
  meta <- read.csv(metadata_file, header=TRUE)
  filtered <- data.frame(meta$id, meta[,colnames(meta) %in% c(drug)])
  colnames(filtered) <- c("id", "phenotype")
  geno$id <-rownames(geno)
  geno <- left_join(geno, filtered)
  geno <- geno[!is.na(geno$phenotype),]
  sample_name <- geno$id
  geno <- geno[ , !(names(geno) %in% c("id"))]
  geno <- apply(geno, 2, function(x) as.numeric(as.character(x)))
  rownames(geno) <- sample_name
  return(geno)
}
                
## Function to normalize weights
norm <- function(x) {
  return ((x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE)))
}
                
## Function to remove columns with MAF of 0
filter_maf <- function(geno){
  maf <- colSums(geno[,-ncol(geno)])
  maf <- maf/(nrow(geno))
  maf <- data.frame(maf)
  maf <- maf %>% na.omit() %>% filter(maf>0) 
  geno <- data.frame(geno)
  geno <-geno[,names(geno) %in% c(rownames(maf), "phenotype")]
  return(geno)
}
                
## Function for converting parsimony score into probability for splitting
norm_fitch <- function(fitch_df){
  fitch <- fitch_df %>% filter(Mutation %in% colnames(train_geno[,-ncol(train_geno)]))
  fitch$norm <- norm(fitch[,2])
  return(fitch)
}
                
## Functions to perform prediction
# Training random forest
perf_pred <- function(geno_train, max_depth, select_weights, model_weights, tunegrid){
  train(
    phenotype ~ .
    ,data = geno_train
    ,method = "ranger"
    ,metric="ROC"
    ,trControl = trainControl(method="cv",  number = 5, allowParallel = TRUE, verbose = TRUE ,savePredictions = TRUE, classProbs = TRUE, summaryFunction=twoClassSummary)
    ,tuneGrid = tunegrid
    ,num.trees = 1000,  split.select.weights=select_weights, weights=model_weights, importance = "impurity", max.depth= max_depth, seed=10)
}
                
# Calculate class imbalance weights and fit random forest
train_model <- function(geno_train, fitch_weight, tunegrid){
  geno_train <- data.frame(geno_train)
  geno_train[,ncol(geno_train)] <- as.character(geno_train[,ncol(geno_train)])
  model_weights <- ifelse(geno_train[,ncol(geno_train)] == "0",(1/table(geno_train[,ncol(geno_train)])[1]) * 0.5,(1/table(geno_train[,ncol(geno_train)])[2]) * 0.5)
  geno_train$phenotype <- factor(geno_train$phenotype, levels = c("0","1"), labels = c("S", "R"))
  fit <- perf_pred(geno_train, max_depth=10, select_weights=fitch_weight, model_weights, tunegrid)
  return(fit)
}
                
# Evaluate Predicition
predict_eval <- function(geno_test, pheno_test, model){ 
  pheno_test<- factor(pheno_test, levels = c("0","1"), labels = c("S", "R"))
  geno_test <- geno_test[,-ncol(geno_test)]
  pred <- predict(model, test_geno)
  conf_matrix <- confusionMatrix(pred, pheno_test, positive="R")
  evaluate <- evalm(model, showplots=FALSE)
  evaluate <- evaluate$stdres
  return(list(conf_matrix, evaluate))
}
                
# Feature importance using Gini
var_imp <- function(model, importance_file_gini){
  imp <- varImp(model)
  return(imp)
}
                
# Fit second model using ranger for importance threshold and interactions only
train_model2 <- function(geno_train, weight){
  geno_train <- data.frame(geno_train)
  model_weights <- ifelse(geno_train[,ncol(geno_train)] == "0",(1/table(geno_train[,ncol(geno_train)])[1]) * 0.5,(1/table(geno_train[,ncol(geno_train)])[2]) * 0.5)
  geno_train$phenotype <- factor(geno_train$phenotype, levels = c("0","1"), labels = c("S", "R"))
  fit2 <- ranger(phenotype ~ ., data=geno_train,num.trees = 1000,case.weights=model_weights, split.select.weights=weight, importance = "impurity", max.depth= 10, seed=10, mtry = sqrt(ncol(geno_train)-1), splitrule = fit$final$tuneValue$splitrule)
  return(fit2)
}

# Calculate importance threshold for features
get_importance <- function(imp, fit,  train_geno, fw){
  subsets <- seq(1:nrow(imp))
  imp<- data.frame(imp)
  imp$feature <- rownames(imp)
  imp$fw <-fw
  imp <- imp[order(-imp$Overall),]
  results <- matrix(nrow=length(subsets), ncol=3)
  colnames(results) <- c("Features","ROC", "ROCSD")
  results[,1] <- subsets
  for (i in 1:length(subsets)){
    imp2 <- imp[1:subsets[i],]
    train_geno2 <- train_geno[,colnames(train_geno) %in% c(imp2$feature, "phenotype")]
    train_geno2 <- data.frame(train_geno2)
    train_geno2$phenotype <- as.factor(train_geno2$phenotype)
    if (subsets[i] ==1){
      fw2 <- imp$fw[1]
      mtry_optimum = 1
    }
    if (subsets[i] >1){
      fw2 <- imp$fw[1:i]
      n <- length(fw2[!fw2==0])
      mtry_optimum = c(sqrt(n))
    }
    fit_optimum <- expand.grid(
      mtry = mtry_optimum,
      min.node.size = fit$final$tuneValue$min.node.size, # default for classification
      splitrule = fit$final$tuneValue$splitrule)
    fit2 <- train_model(train_geno2, fw2, fit_optimum)
    results[i,2] <- fit2$results$ROC[1]
    results[i,3] <- fit2$results$ROCSD[1]
  }
  return(results)
}
                
# Funtion to get interactions  
get_interactions <- function(ranger_obj, n_trees){
  for (j in c(1:1000)){
    tree1 <- treeInfo(fit2, tree=j)
    tree1$leftChild <- tree1$leftChild + 1
    tree1$rightChild <- tree1$rightChild + 1
    child_parent <- matrix(ncol=2, nrow=nrow(tree1)*2)
    for (i in 1:nrow(tree1)){
      parent = tree1[i,5]
      child1_n = tree1[i,2]
      child2_n= tree1[i,3]
      child1 = tree1[child1_n,5]
      child2 = tree1[child2_n,5]
      row_1 = 2*i-1
      row_2 = 2*i
      child_parent[row_1,1] = parent
      child_parent[row_1,2] = child1
      child_parent[row_2,1] = parent
      child_parent[row_2,2] = child2
    }
    child_parent <- data.frame(child_parent)
    colnames(child_parent) <- c("parent","child")
    child_parent <- child_parent[!is.na(child_parent$child),]
    list_df[[j]] <- child_parent
  }
  child_parent_df <- do.call(rbind.data.frame, list_df)
  child_parent_count<- setDT(child_parent_df)[,list(Count=.N),names(child_parent_df)]
  return(child_parent_count)
}
                
# Stage 1- read in files and format
cat(" *Formatting training and testing data\n")
train_geno <- read.csv(opt$train_dataset, header=TRUE)
sample_name <- train_geno[,1]
rownames(train_geno) <- train_geno[,1]
train_geno<- train_geno[,-1]
test_geno <- read.csv(opt$test_dataset, header=TRUE)
rownames(test_geno) <- test_geno[,1]
test_geno<- test_geno[,-1]

# Remove samples with missing phenotype
train_geno <- filter_missing(opt$metadata, train_geno, opt$drug)
train_geno <- filter_maf(train_geno)
test_geno <- filter_missing(opt$metadata, test_geno, opt$drug)
train_pheno <- get_pheno(train_geno)
test_pheno <- get_pheno(test_geno)
if (opt$weight_type=="fitch"){
  fitch_df <- read.csv(opt$weight_file, header=TRUE)
  fitch_df <- norm_fitch(fitch_df)
  fw <- fitch_df$norm
}
                
# Ensure training and testing columns match
train_geno <- train_geno[,colnames(train_geno) %in% fitch_df$Mutation]
train_geno <- data.frame(train_geno)
train_geno$phenotype <- train_pheno
test_geno <- test_geno[,colnames(test_geno) %in% colnames(train_geno)]
test_geno <- test_geno[,!colnames(test_geno) %in% c("phenotype")]
                
# Tuning parameters
cat(" *Selecting tuning parameters\n")
n <- length(fw[!fw==0])
tunegrid <- expand.grid(
  mtry = c(sqrt(n)), #default is to use root
  min.node.size = c(1), # default for classification
  splitrule = c("extratrees","gini"))
n.cores <- 24 # choose number of cores
                
## Perform initial prediction
cat(" *Performing initial prediction\n")
fit <- train_model(train_geno, fw, tunegrid)
perf_list <- predict_eval(test_geno, test_pheno, fit)
                
## Calculate Importance
imp <- var_imp(fit)
imp <- imp$importance
write.csv(imp, paste0(opt$output,"_importance.csv"))
write.csv(perf_list[2], paste0(opt$output,"_AUC.csv"))
                
## Calculate importance threshold for features
results <-get_importance(imp, fit, train_geno, fw)
write.csv(results, paste0(opt$output,"_feature_threshold.csv"), row.names=FALSE)
                
## Get interactions
list_df <- list()
fit2 <- train_model2(train_geno,fw)
child_parent_count <- get_interactions(fit2, 1000)
write.csv(child_parent_count, paste0(opt$output,"_interactions.csv"), row.names=FALSE)







