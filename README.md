# DonorGenotypeFree
R Script for donor geneotype free analysis of donor fraction for SOT recipients

[![DOI](https://zenodo.org/badge/267396914.svg)](https://zenodo.org/badge/latestdoi/267396914)

This script was developed for the following article:

__*Cell free DNA in paediatric solid organ transplantation using a new detection method of separating donor-derived from recipient cell free DNA*__

Author List:

Evgenia Preka 1,2*, Drew Ellershaw 3*, Natalie Chandler 3, Helena Ahlfors 3, Helen Spencer 4,  Matthew J Fenton 4, Stephen D Marks 1,5

"*" These authors contributed equally to this work

1. Department of Paediatric Nephrology, Great Ormond Street Hospital for Children NHS Foundation Trust, Great Ormond Street, London, UK	
2. Paediatric Nephrology, University Hospital Southampton NHS Foundation Trust, Southampton, UK
3. London North Genomic Laboratory Hub, Great Ormond Street Hospital for Children NHS Foundation Trust, London, UK
4. Cardiothoracic Unit, Great Ormond Street Hospital for Children NHS Foundation Trust, Great Ormond Street, London, UK
5. University College London Great Ormond Street Institute of Child Health, NIHR Great Ormond Street Hospital Biomedical Research Centre, London, WC1N 1EH, UK

## Introduction
Solid organ transplant recipients require lifelong immunosupression to stop their bodies from rejecting the donor organ (allograft). However, allograft rejction can still occur and patients need regular monitoring to look for signs of rejection. Biopsies of the allograft are considered the gold standard to look for histological changes to the allograft but is costly, and carries a small risk. Non-invasive methods using DNA derived from the donor organ present the in recipient blood plasma are being developed, but some require a sample from the donor which is not always available. We have developed a SNP based assay which does not require a sample from the donor. This is termed a genotype free analysis. We use the R mixtools package normalmixEM function, previously described by [Gordon et al. (2016)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5031701/) to model for recipient and donor genotypes to determine the amount of donor-derived cell-free DNA (ddcfDNA) present in the recipient. Increased ddcfDNA can be detected before clinical signs of rejection and can allow for earlier intervention to improve patient outcomes.


### Quick Start using provided example data

1. Download this repository

2. Open MOAT-genotype-free.R in R studio

3. Install the mixtools package (you may be prompted to install it by R studio in a popup)

4. Run the script


This will create a "MOAT_Test" folder containing each run where the Results_Test.csv is, with a folder for each sample in the run. 
An overview of the results is output as a .pdf with 4 plots of processing steps.
Summary statistics are also saved as well as intermediate processing steps.


The Results_Test.csv can be saved into any subfolders. 
The filename MUST begin with "Results" to be found with the script.
This is example data only. ID30 and ID31 are mock examples representing approximately 10% and 5% donor fraction.
Columns can be added or deleted, and number of counts changed to experiment with the data. However the normalmixEM may fail if there is no variance in the data, or it does not include enough of each of the expected reciepient/donor genotype combinations.


Please refer to the published paper for study info
