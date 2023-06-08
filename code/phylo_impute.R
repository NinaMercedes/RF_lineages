##### Perform phylogenetic imputation according to nearest neighbour in tree
### Load required packages
library(phytools)
library(phangorn)
library(data.table)
library(reshape2)
library(dplyr)
library(caret)
library(ranger)
library(ape)

### Set seed
set.seed(10)

### Make options
option_list <- list(
  make_option(c("-a", "--train_dataset"), type="character", default="train_dataset.csv",help="Insert path to train dataset file (csv)"),
  make_option(c("-t", "--tree"), type="character", default="train_noppe.tree",help="Insert path to tree"),
  make_option(c("-s", "--sample_names"), type="character", default="train_names_noPPE.txt",help="Insert path to sample names(txt)"),
  make_option(c("-n", "--root"), type="character", default="root",help="Insert root sample name"),
  make_option(c("-r", "--root_type"), type="character", default="og",help="og or mp"),
  make_option(c("-o", "--output"), type="character", default="train_dataset_imputed.csv",help="Insert path to imputed output file (csv)")

)

### Import training data
train_geno <- read.csv(opt$train_dataset)

### Get SNP name
snp_name <- train_geno[,1]

### Read in sample names
sample_name <- read.table(opt$sample_names)

### Format
train_geno <- train_geno[,-c(1)]
is.na(train_geno)<- train_geno == "N" 
is.na(train_geno)<- train_geno == "NA"
train_geno <- t(train_geno)
rownames(train_geno)<- sample_name$V1
colnames(train_geno)<- snp_name

### Function to pre-process tree
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
  
imp_tree <- preTree(treefile=opt$tree, root_type=opt$root_type, root_label=opt$root)

### Function to perform imputation
geno_impute <- function(tree, geno){
  x <- cophenetic(tree)
  x <- melt(x)
  x <- x[!x$Var1==x$Var2,]
  y <- list()
  l <- data.frame()
  for (i in 1:ncol(geno)){
    for (j in rownames(geno)){
     if (is.na(geno[j,i])==TRUE){
       y[[i]] <- x %>% filter(Var1==j) %>% filter(value==min(value)) 
       l <- y[[i]][2]
       l <- as.character(l$Var2)
       l<- l[1]
       geno[j,i] <- geno[l, i]
       if (is.na(geno[l, i])==TRUE){
         geno[j,i] <- 0
         }
       }
     }
   }
  return(geno)
}

train_geno <- geno_impute(imp_tree, geno)

write.csv(train_geno, opt$output)
