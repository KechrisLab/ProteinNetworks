# Proteomic Networks and Related Genetic Variants Associated with Smoking and Chronic Obstructive Pulmonary Disease


**Note:** If you want to refer to the result in the proteomics network paper, please cite:

> Konigsberg, I. R., Vu, T., Liu, W., Litkowski, E. M., Pratte, K. A., Vargas, L. B., ... & Kechris, K. (2024). 
> Proteomic Networks and Related Genetic Variants Associated with Smoking and Chronic Obstructive Pulmonary Disease.
> medRxiv, 2024-02.

If you use SmCCNet in published research, please cite:

> Liu, W., Vu, T., Konigsberg, I., Pratte, K., Zhuang, Y., & Kechris, K. (2023). 
> SmCCNet 2.0: A Comprehensive Tool for Multi-omics Network Inference with Shiny Visualization. 
> bioRxiv, 2023-11.

> Shi, W. J., Zhuang, Y., Russell, P. H., Hobbs, B. D., Parker, M. M.,
> Castaldi, P. J., … & Kechris, K. (2019). Unsupervised discovery of
> phenotype-specific multi-omics networks. Bioinformatics, 35(21),
> 4336-4343.



This repository contains the R code for Koningsberg et al., Proteomic Networks and Related Genetic Variants Associated with Smoking and Chronic Obstructive Pulmonary Disease (submitted).

Two folders contain the source code and analysis code: 

- NetworkExtractCode: Analysis R script of extracting proteomics networks using SmCCNet for each combination of phenotype & cohort.
  - soma_adjusted_aa_fev1.R: SmCCNet pipeline analysis script with the phenotype of FEV1 within the African American cohort.
  - soma_adjusted_nhw_fev1.R: SmCCNet pipeline analysis script with the phenotype of FEV1 within the Non-Hispanic White cohort.
  - soma_adjusted_aa_emp.R: SmCCNet pipeline analysis script with the phenotype of percent emphysema within the African American cohort.
  - soma_adjusted_nhw_emp.R: SmCCNet pipeline analysis script with the phenotype of percent emphysema within the Non-Hispanic White cohort.
  - soma_adjusted_aa_smoking.R: SmCCNet pipeline analysis script with the phenotype of smoking status (1 yes, 0 no) within the African American cohort.
  - soma_adjusted_nhw_smoking.R: SmCCNet pipeline analysis script with the phenotype of smoking status (1 yes, 0 no) within the Non-Hispanic White cohort.
- SourceFunction: SmCCNet pipeline source functions. All the source R functions are required to be imported before running the analysis R script.
  - GSmCCNet_one.R: All necessary SmCCNet pipeline functions to run through the SmCCNet pipeline.
  - SCCAdiagTools.R: SmCCNet diagnosis function to assess the performance of SmCCNet algorithm.
 
Before running the SmCCNet analysis script, all source functions need to be imported by using 

```r
source('GSmCCNet_one.R')
source('SCCAdiagTools.R')
```

All data used in the analysis are controlled-access data that restricts the public access. If the user is interested in running single-omics SmCCNet, please refer to the SmCCNet package available in our [GitHub](https://github.com/KechrisLab/SmCCNet) repository. The package compiles all source codes and functions into an end-to-end pipeline with example data to run through. Below is the example implementation:


``` r
library(SmCCNet)
set.seed(123)
data("ExampleData")
Y_binary <- ifelse(Y > quantile(Y, 0.5), 1, 0)
# single-omics with binary phenotype
result <- fastAutoSmCCNet(X = list(X1), Y = as.factor(Y_binary), 
                          Kfold = 3, 
                          subSampNum = 100, DataType = c('Gene'),
                          saving_dir = getwd(), EvalMethod = 'auc', 
                          summarization = 'NetSHy', 
                          CutHeight = 1 - 0.1^10, ncomp_pls = 5)
# single-omics with quantitative phenotype
result <- fastAutoSmCCNet(X = list(X1), Y = Y, Kfold = 3, 
                          preprocess = FALSE,
                          subSampNum = 50, DataType = c('Gene'),
                          saving_dir = getwd(), summarization = 'NetSHy',
                          CutHeight = 1 - 0.1^10)
```
