# Random Forest Lineages
Code for correcting for lineage-dependency in random forest prediction of drug-resistant phenotypes from *Mycobacterium tuberculosis* genome sequences.
Code is supplied in pipelines which can be run using the following examples. 
There are a few requirements:
* A csv file containing the binary genotype data. The first column must contain the sample id.
* A metadata file. The first column must contain the sample id with the name 'id', there should also be the column 'lineage' and name of drug. 
* A tree file in newick format. 

## Examples of how to run code on command line

