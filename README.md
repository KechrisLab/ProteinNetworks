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

Two folders contains the source code and analysis code: 

- NetworkExtractCode: Analysis R script of extracting proteomics networks using SmCCNet for each combination of phenotype & cohort.
  - soma_adjusted_aa_fev1.R: SmCCNet pipeline analysis script with the phenotype of FEV1 within the African American cohort.
  - soma_adjusted_nhw_fev1.R: SmCCNet pipeline analysis script with the phenotype of FEV1 within the Non-Hispanic White cohort.
  - soma_adjusted_aa_emp.R: SmCCNet pipeline analysis script with the phenotype of percent emphysema within the African American cohort.
  - soma_adjusted_nhw_emp.R: SmCCNet pipeline analysis script with the phenotype of percent emphysema within the Non-Hispanic White cohort.
  - soma_adjusted_aa_smoking.R: SmCCNet pipeline analysis script with the phenotype of smoking status (1 yes, 0 no) within the African American cohort.
  - soma_adjusted_nhw_smoking.R: SmCCNet pipeline analysis script with the phenotype of smoking status (1 yes, 0 no) within the Non-Hispanic White cohort.
- SourceFunction: SmCCNet pipeline source functions. All the source R functions are required to be imported before running the analysis R script.
