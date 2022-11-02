### Calculate Parsimony Score
## Load required packages
library(phytools)
library(phangorn)
library(data.table)
library(optparse)
library(dplyr)
library(ape)
### Set seed
set.seed(10)
### Make options
option_list <- list(
  make_option(c("-a", "--train_dataset"), type="character", default="train_dataset.csv",help="Insert path to train dataset file (csv)"),
  make_option(c("-m", "--metadata"), type="character", default="metadata.csv",help="Insert path to metadata file (csv)"),
  make_option(c("-d", "--drug"), type="character", default="rifampicin",help="Insert name of drug- column name must match to metadata"),
  make_option(c("-r", "--root_type"), type="character", default="mp",help="midpoint (mp) or outgroup (og)"),
  make_option(c("-l", "--root_label"), type="character", default="NULL",help="Name of outgroup label"),
  make_option(c("-t", "--treefile"), type="character", default="train_dataset.tree",help="tree in newick format"),
  make_option(c("-o", "--output"), type="character", default="results/rifampicin",help="prefix of file to output")
)
### Parse options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)
### Functions
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
## Function to remove columns with MAF of 0
filter_maf <- function(geno){
  maf <- colSums(geno[,-ncol(geno)])
  maf <- maf/(nrow(geno))
  maf <- data.frame(maf)
  maf <- maf %>% filter(maf>0)
  geno <- data.frame(geno)
  geno <- geno[,names(geno) %in% c(rownames(maf), "phenotype")]
  return(geno)
}
## pre-process the tree
preTree <- function(treefile, root_type=c("mp","og"), root_label) {
  tree <- read.tree(treefile)
  tree <- multi2di(tree)
  tree <- reorder.phylo(tree, order="pruningwise")
  if(is.null(root_type)){
    tree <- midpoint.root(tree)
  }
  if(root_type == "mp"){
    tree <- midpoint.root(tree)
  }
  if(root_type == "og"){
    tree <- root(tree, outgroup = root_label, resolve.root=TRUE)
  }
  return(tree)
}
## Filter tree to extract isolates with available phenotype data and calculate parsimony score
calc_fitch <- function(tree, geno){
  geno_1 <- geno[,-ncol(geno)]
  sample_name_1 <- rownames(geno_1)
  tree_filtered_1 <- keep.tip(tree, sample_name_1)
  phydat_1 <- as.phyDat(as.matrix(geno_1),type="USER", levels=c("0","1"))
  ind <- attr(phydat_1, "index") #get unique index number
  fitch.phangorn <- phangorn::fitch 
  fitch <- fitch.phangorn(tree_filtered_1, phydat_1 , site="site")
  fitch_score <- fitch[ind]
  fitch_score <- data.frame(colnames(geno_1), fitch_score)
  colnames(fitch_score) <- c("Mutation", "Fitch Score")
  return(fitch_score) 
}
### Read in and format files
cat(" *Formatting data\n")
train_geno <- read.csv(opt$train_dataset, header=TRUE)
sample_name <- train_geno[,1]
rownames(train_geno) <- sample_name
train_geno <- filter_missing(opt$meta, train_geno, opt$drug)
cat(" *Filtering data\n")
train_geno <- filter_maf(train_geno)
### Pre-process the tree
cat(" *Pre-processing tree\n")
phy <- preTree(treefile=opt$treefile, root_type=opt$root_type, root_label=opt$root_label)
### Calculate the parsimony score and Save file
cat(" *Calculating Fitch Score\n")
fitch_score <- calc_fitch(phy, train_geno)
write.csv(fitch_score, paste0(opt$output,"_weights.csv"), row.names=FALSE)
