## Data Splitting
## Load required packages
library(data.table)
library(dplyr)
library(optparse)
library(caret)
## options
option_list <- list(
  make_option(c("-a", "--dataset"), type="character", default="dataset.csv",help="Insert path to train dataset file (csv)"),
  make_option(c("-m", "--metadata"), type="character", default="metadata.csv",help="Insert path to metadata file (csv)"),
  make_option(c("-l", "--lineage"), type="character", default="global",help="lineage name, combined or global"),
  make_option(c("-o", "--output"), type="character", default="global",help="prefix of file to output")
)
### Parse options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)
## Functions
# Perform Split
perform_split <- function(meta, geno, output){
  train.rows <- createDataPartition(y= meta$lineage, p=0.8, list = FALSE)
  train.data <- meta[train.rows,] 
  train_data <- data.frame(train.data$id)
  colnames(train_data) <- "id"
  test.data <- meta[-train.rows,]
  test_data <- data.frame(test.data$id)
  colnames(test_data) <- "id"
  geno$id <- rownames(geno)
  train_geno <- left_join(train_data, geno)
  test_geno <- left_join(test_data, geno)
  write.csv(train_geno, paste0(output, "_train_dataset.csv"), row.names=FALSE)
  write.csv(test_geno, paste0(output, "_test_dataset.csv"), row.names=FALSE)
}
# Function to filter lineage
filter_meta_lineage <- function(meta, lineage, name_l, output){
  meta <- meta[grep(paste(lineage,collapse="|"), meta$lineage),]
  write.csv(meta, paste0(output, "_", name_l ,"_meta.csv"))
  return(meta)
}
## Make training and testing datasets
geno <- read.csv(opt$dataset, header=TRUE)
rownames(geno) <- geno$X
geno<- geno[,-1]
if(opt$lineage=="global"){
  metadata <- read.csv(opt$metadata, header=TRUE)
  perform_split(metadata, geno, opt$output)
}
if(opt$lineage=="lineage 2"){
  metadata <- read.csv(opt$metadata, header=TRUE)
  metadata_2 <- filter_meta_lineage(metadata, "lineage2.*", "lineage2", opt$output)
  perform_split(metadata_2, geno, opt$output)
}
if(opt$lineage=="lineage 4"){
  metadata <- read.csv(opt$metadata, header=TRUE)
  metadata_4 <- filter_meta_lineage(metadata, "lineage4.*", "lineage4",opt$output)
  perform_split(metadata_4, geno, opt$output)
}
if(opt$lineage=="combined"){
  metadata <- read.csv(opt$metadata, header=TRUE)
  metadata_comb <- filter_meta_lineage(metadata, c("lineage2.*", "lineage4.*"), "combined", opt$output)
  perform_split(metadata_comb, geno, opt$output)
}
