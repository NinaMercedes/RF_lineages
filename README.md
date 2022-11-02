# Random Forest Lineages
Code for correcting for lineage-dependency in random forest prediction of drug-resistant phenotypes from *Mycobacterium tuberculosis* genome sequences.
Code is supplied in pipelines which can be run using the following examples. 
There are a few requirements:
* A csv file containing the binary genotype data. The first column must contain the sample id.
* A metadata file. The first column must contain the sample id with the name 'id', there should also be the column 'lineage' and name of drug. 
* A tree file in newick format. 

## Examples of how to run code on command line
### Data Splitting: training and testing data
```
Rscript Data_Splitting.R --dataset "dataset.csv" --metadata "metadata.csv" --lineage "global" --output "global"
```
### Unweighted Model
```
Rscript Unweighted_RF.R --train_dataset "global_train_dataset.csv" --test_dataset "global_train_dataset.csv" --metadata "metadata.csv" --drug "rifampicin" --output "global_unweighted"
```
### Calculate Weights
```
Rscript Parsimony_Score.R --train_dataset "global_train_dataset.csv"  --metadata "metadata.csv" --drug "rifampicin" --root_type "og" --root_label "ERR2510183" --treefile "noppe_train.tree" --output "global"
```
### Weighted Model
```
Rscript Unweighted_RF.R --train_dataset "global_train_dataset.csv" --test_dataset "global_train_dataset.csv" --metadata "metadata.csv" --drug "rifampicin" --weight_file "global_weights.csv" --weight_type "fitch" --output "global_weighted"
```
