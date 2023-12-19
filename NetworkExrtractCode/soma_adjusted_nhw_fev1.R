# import source code
source('SCCAdiagTools.R')
source("GSmCCNet_one.R")
# import required packages
library(dplyr)
library(purrr)
library(tidyr)
library(EnvStats)
library(pROC)
library(ggplot2)
library(reshape2)
library(dbplyr)
library(dplyr)
library(tidyverse)
library(WGCNA)
library(ggplot2)
library(igraph)
library(MASS)
library(lqmm)
library(pls)
library(igraph)
library(ggcorrplot)
library(dplyr)
library(purrr)
library(tidyr)
library(EnvStats)
library(RCy3)

######################################################Import and preprocess data
# load protein data (latest version)
load("FEV1FEV1_filtered_Weixuan.Rdata")
X1 <- read.table("COPDGeneSoma_SMP_5K_P2_16Jun20.txt", header = TRUE)
dat <- read.table("COPDGeneSoma_SMP_5K_P2_AA_NHW_Matched_5_25_21.txt", header = TRUE)
X1 <- merge(x = dat, y = X1, by = 'sid', all.y = FALSE, all.x = FALSE)
# calculate coefficient of variation
cv_value = apply(X1[,-(1:8)], 2, cv)
# load clinical data
clinical <- read.csv("COPDGene_P1P2P3_25SEP2020_SubjectFlattened.csv")
# filter protein data to have only Afican American
X1 <- X1 %>% 
  filter(aa == 0)
# set id name
colnames(clinical)[1] <- 'sid'
# preprocess clinical data
fev1_covariate_aa <- clinical %>%
  select(sid = sid,Platelets_P2,wbc_P2,BMI_P2, gender, smoking_status_P2,ccenter_P2, 
         Age_P2, FEV1_utah_P2) %>% # include all relevant clinical covariates, phenotypes
  merge(x = ., y = X1, by = "sid", all.y = TRUE) %>% # merge with protein data
  drop_na()%>% # drop na value
  filter(aa == 0) # include only African American
# subset covariate data
covariates <- fev1_covariate_aa[,2:8]
# subset protein data
X1_covariates <- fev1_covariate_aa[,-c(1:16)]
# center and sclae protein data
X1_covariates <- as.data.frame(scale(X1_covariates, center = TRUE, scale = TRUE))
# use coefficient of variation to filter out proteins 
X1_covariates <- X1_covariates[,cv_value >= 0.7]
# regress out sex, age, and clinical center
map_aa <- X1_covariates %>%   
  map(~lm(.x ~ 
            covariates$gender + 
            covariates$ccenter_P2 + covariates$Age_P2, 
          data = X1_covariates)) %>%
  map(resid)
# obtain adjusted data after regress-out
Adjusted_AA <- map_df(map_aa, ~as.data.frame(t(.)))
Adjusted_AA <- t(Adjusted_AA)
colnames(Adjusted_AA) <- colnames(X1_covariates)

#################################################################### Run SmCCNet
# define phenotype to FEV1 at P2
Y <- as.matrix(fev1_covariate_aa$FEV1_utah_P2)
# define protein data
X1 <- Adjusted_AA
# set phenotype name
colnames(Y)[1] <- "FEV1_NHW_Adjusted"
pheno <- "FEV1_NHW_Adjusted"
# set number of folds for k-fold cross-validation
K <- 5
# set directory where the data and result will be stored
CVDirs <- paste0(pheno,  "/", K, "foldCV", pheno, "_CCcoef_1_", 
                 "/")
# create directory
dir.create(CVDirs)
# create sub-directory to store result
resultDir <- paste0(pheno, "/")
CVDir <- CVDirs
# read in protein label data 
label <- read.csv("SOMAscan_Assay_v4_Annotations_version3.3.2.csv")
# transform protein sequence id to match format
newint <- intersect(gsub("_", "-",substring(colnames(X1), 2)), as.character(label$SeqId))
sum(gsub("_", "-",substring(colnames(X1), 2)) == newint)
label_data <- subset(label, SeqId %in% newint)
# fill in the protein target full name
label_data$TargetFullName <- ifelse(label_data$TargetFullName == "", 
                                    as.character(label_data$Target), 
                                    as.character(label_data$TargetFullName))
# transform sequence id to protein full name
label_data <- label_data %>%
  select(SeqId, TargetFullName)
colnames(X1) <- as.character(label_data$TargetFullName)

################################################## Source Code

# set random seed 
set.seed(12345)
# calculate correlation matrix between protein
bigCor2 <- cor(X1)
# store protein labels into a vector
AbarLabel <- colnames(bigCor2)
# store number of rows/columns
N <- nrow(X1); p1 <- ncol(X1)
# set subsampling proportion and sumsampling number
s1 <- 0.7; subSamp <- 100
# set a list of potential candidates
pen1 <- seq(.1, .5, by = .1)
# save data and parameter setup into local directory
save(X1, Y, s1, subSamp, pen1, 
     file = paste0(CVDir, "Data.Rdata"))
# split data into K fold and store them into local directory
foldIdx <- split(1:N, sample(1:N, K))
for(i in 1:K){
  iIdx <- foldIdx[[i]]
  x1.train <- scale(X1[-iIdx, ])
  yy.train <- scale(Y[-iIdx, ])
  x1.test <- scale(X1[iIdx, ])
  yy.test <- scale(Y[iIdx, ])
  
  if(is.na(min(min(x1.train), min(yy.train), min(x1.test), min(yy.test)))){
    stop("Invalid scaled data.")
  }
  
  subD <- paste0(CVDir, "CV_", i, "/")
  dir.create(subD)
  save(x1.train, yy.train, x1.test, yy.test,
       s1, pen1, p1, subSamp,
       file = paste0(subD, "Data.Rdata"))
}

# make parallel cluster for parallel computing
cl <- makeCluster(5, type = "PSOCK")
# export working directory to each parallel cluster
clusterExport(cl = cl, "CVDir")

# run 5-fold cross-validation
parSapply(cl, 1:5, function(CVidx){
  # the function need to be sourced again so that 
  source("GSmCCNet_one.R")
  # create subdirectory to load cross-validation result
  subD <- paste0(CVDir, "CV_", CVidx, "/")
  # load data
  load(paste0(subD, "Data.Rdata"))
  # create subdirectory to save cross-validation result
  dir.create(paste0(subD, "SmCCA/"))
  # create empty vector to save canonical correlation (train/test)
  RhoTrain <- RhoTest <- DeltaCor <- rep(0, length(pen1))
  # run SmCCA with each penalty candidate
  for(idx in 1:length(pen1)){
    # set penalty term for current iteration
    l1 <- pen1[idx]
    # run SmCCA
    Ws <- getRobustPseudoWeight_one(x1.train, yy.train, l1, s1,
                                    NoTrait = FALSE, FilterByTrait = FALSE, 
                                    SubsamplingNum = subSamp, CCcoef = CCcoef)
    # take average of all the subsamples
    meanW <- rowMeans(Ws)
    # store the canonical weight
    v <- meanW[1:p1]
    # calculate the training canonical correlation 
    rho.train <-  cor(x1.train %*% v, yy.train)
    # calculate the testing canonical correlation
    rho.test <- cor(x1.test %*% v, yy.test) 
    # store the training/testing canonical correlation and difference
    RhoTrain[idx] <- round(rho.train, digits = 5)
    RhoTest[idx] <- round(rho.test, digits = 5)
    DeltaCor[idx] <- abs(rho.train - rho.test)
    # save cross-validation result
    if(idx %% 10 == 0){
      save(pen1, RhoTrain, RhoTest, DeltaCor, idx,
           file = paste0(subD, "temp.Rdata"))
    }
    
  }
  # combine the penalty selection, training/testing canonical correlation and difference
  DeltaCor.all <- cbind(pen1, RhoTrain, RhoTest, DeltaCor)
  colnames(DeltaCor.all) <- c("l1", "Training CC", "Test CC", "CC Pred. Error")
  # save the result to local directory
  write.csv(DeltaCor.all,
            file = paste0(subD, "SmCCA/SCCA_", subSamp,"_allDeltaCor.csv"))
  # remove temporary file to reduce the memory issue
  system(paste0("rm ", subD, "temp.Rdata"))
  
  return(CVidx)
}
)

# close cluster
stopCluster(cl)
# aggregate cross-validation result and save it to the local directory
plotCVcontour(CVDir, "SmCCA", NumSubsamp = subSamp)
# set up local directory to store the final result 
plotD <- paste0(CVDir, "Figures/")
saveD <- paste0(CVDir, "Results/")
dataF <- paste0(CVDir, "Data.Rdata")
dir.create(plotD)
dir.create(saveD)
dir.create(dataF)
# set the method used
Method = "SmCCA"
# choose the optimal penalty and use it to run SmCCA on the whole dataset
for(Method in "SmCCA"){
  # read in aggregated cross-validation result
  T12 <- read.csv(paste0(CVDir, "Results/", Method, "CVmeanDeltaCors.csv"))[ , -1]
  # choose the optimal penalty term
  pen <- which.min(T12[ , 3]/T12[ ,2])
  # use SmCCA
  if(Method == "SmCCA"){
    FilterByTrait <- FALSE
  }else if(Method == "SsCCA"){
    FilterByTrait <- TRUE
  }
  # define optimal penalty term
  l1 <- pen;
  
  system.time({
    # run SmCCA on the complete dataset with optimal penalty term
    Ws <- getRobustPseudoWeight_one(X1, Y, l1, s1, NoTrait = FALSE,
                                    FilterByTrait = FilterByTrait,
                                    SubsamplingNum = subSamp)
    # obtain adjacency matrix
    Abar <- getAbars(Ws, P1 = p1, FeatureLabel = AbarLabel[1:p1])
    # save the result to local directory
    save(l1, X1, Y, s1, Ws, Abar,
         file = paste0(saveD, Method, K, "foldSamp", subSamp, "_", pen, ".Rdata"))
  })
  
  
}
############################## Obtain subnetwork from the global network
# load global network
load('D:/Research/SOMA/FEV1_NHW_Adjusted/5foldCVFEV1_NHW_Adjusted_CCcoef_1_/Results/SmCCA5foldSamp100_0.5.Rdata')
# clustering
mirGeneModule <- getMultiOmicsModules(Abar, p1, PlotTree = FALSE)
AbarLabel <- colnames(X1)
AbarLabel <- make.unique(AbarLabel)
row.names(Abar) <- colnames(Abar) <- AbarLabel
# global correlation matrix
bigCor2 <- cor(X1)
row.names(bigCor2) <- colnames(bigCor2) <- AbarLabel
colnames(X1) <- AbarLabel
p1 <- ncol(X1)
pheno <- 'FEV1_NHW_Adjusted'
cor_result <- list()
# store subnetwork with each specified size to the local directory
for (z in 1:length(mirGeneModule))
{
  if (length(mirGeneModule[[z]]) <= 10) next
  
  if (length(mirGeneModule[[z]]) > 200) 
    temp_size <- 200
  else 
    temp_size <- length(mirGeneModule[[z]])
  cor_result[[z]] <- data.frame(size = seq(from = 10, to = temp_size, by = 1),
                                pc1_cor = 0,
                                hybrid_cor = 0,
                                hybrid_zero_cor = 0)
  for (i in 1:nrow(cor_result[[z]]))
  {
    print(i)
    size <- cor_result[[z]]$size[i]  
    if(size == 0){
      netD <- paste0("D:/Research/SOMA/FEV1_NHW_Adjusted/5foldCVFEV1_NHW_Adjusted_CCcoef_1_/", "NetworkEdgeCutpt0/")
      titlePre <- paste0("Net ")
      savePre <- paste0(netD, "Net")
    }else{
      netD <- paste0("D:/Research/SOMA/FEV1_NHW_Adjusted/5foldCVFEV1_NHW_Adjusted_CCcoef_1_/", "NetworkEdgeCutpt", size, "/")
      titlePre <- paste0("Trimmed Net ")
      savePre <- paste0(netD, "TrimmedNet")
    }
    dir.create(netD)
    saveplot.z <- paste0(savePre, z, ".pdf")
    plottitle <- paste0(titlePre, z)
    #size <- cor_result[[z]]$size[i]  
    results <- plotMultiOmicsNetworks_PPR_Validation(Abar, bigCor2, mirGeneModule, 
                                                     data = X1, Pheno =Y, ModuleIdx = z,
                                                     type ='FEV1_NHW_Adjusted', 
                                                     id = pheno_data$sid,  P1 = p1,
                                                     folder = 'FEV1_NHW_Adjusted',
                                                     FeatureLabel = AbarLabel,
                                                     penalty = 0,
                                                     AddCorrSign = TRUE, 
                                                     SaveFile = saveplot.z,
                                                     ShowType1Label = TRUE, 
                                                     PlotTitle = plottitle, 
                                                     NetLayout = "circle",
                                                     ShowNodes = TRUE, 
                                                     VertexLabelCex = 1, 
                                                     VertexSize = 1,pheno = pheno, 
                                                     mod_size = size)
    cor_result[[z]][i,] <- results
  }
  pc_result <- cor_result[[z]]
  write.csv(pc_result,
            paste0('D:/Research/SOMA/FEV1_NHW_Adjusted/5foldCVFEV1_NHW_Adjusted_CCcoef_1_/pc_correlation_validation_', z, '.csv'),
            row.names = FALSE)
}


# network pruning
for (i in 10:200)
{
  temp_size <- i
  # load subnetwork with specified network size
  load(paste0('D:/Research/SOMA/FEV1_NHW_Adjusted/FEV1_NHW_Adjustednet', temp_size, '_1AA.Rdata'))
  if (i == 10)
    # create NetSHy summarization score data frame
    score_frame <- data.frame(pca_hybrid_zero[[1]][,1])
  else
    # store NetSHy summarization score data
    score_frame = cbind(score_frame, pca_hybrid_zero[[1]][,1])
}
# assign name to each column of NetSHy summarization score
colnames(score_frame) <- paste0('score_',size)
# calculate correlation matrix of the NetSHy summarization score among 
# different network size
cormat <- round(x = cor(score_frame), digits = 2)
# take absolute value
cormat <- abs(cormat)
# define vector of correlation between current network and baseline network size
corvec <- cormat[,1]
# define vector of correlation between current network and phenotype
corpheno <- cor(score_frame, Y)
# select the best network size
selection(corvec, corpheno)

# network visualization
# correlation filtering for subnetwork adjacency matrix
M[abs(correlation_sub) < 0.05] <- 0
# store direction of correlation
M_ind <- ifelse(correlation_sub > 0, 1, -1)
# apply direction of correlation to adjacency matrix
M_adj <- M * M_ind
# set diagonal value of adjacency matrix to 0
diag(M_adj) <- 0
# network visualization through cytoscape
graph <- igraph::graph_from_adjacency_matrix(M_adj, mode = 'undirected', weighted = TRUE,
                                             diag = TRUE, add.colnames = NULL, add.rownames = NA)


# define network node type
V(graph)$type <- sub_type
V(graph)$type
# create network through Cytoscape
createNetworkFromIgraph(graph,"fev1_nhw")