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
load("FEV1FEV1_filtered_Weixuan.Rdata")
# load protein data
X1 <- read.table("COPDGeneSoma_SMP_5K_P2_16Jun20.txt", header = TRUE)
dat <- read.table("COPDGeneSoma_SMP_5K_P2_AA_NHW_Matched_5_25_21.txt", header = TRUE)
X1 <- merge(x = dat, y = X1, by = 'sid', all.y = FALSE, all.x = FALSE)
# calculate coefficient of variation
cv_value = apply(X1[,-(1:8)], 2, cv)
#################### Read-in clinical data
clinical <- read.csv("COPDGene_P1P2P3_25SEP2020_SubjectFlattened.csv")
# filter protein data
X1 <- X1 %>% 
  filter(aa == 0)
colnames(clinical)[1] <- 'sid'
fev1_covariate_aa <- clinical %>%
  dplyr::select(sid = sid,Platelets_P2,wbc_P2,BMI_P2, gender, smoking_status_P2,ccenter_P2, Age_P2, FEV1_utah_P2) %>%
  merge(x = ., y = X1, by = "sid", all.y = TRUE) %>%
  drop_na()%>%
  filter(aa == 0)
# subset only the clinical covariates 
covariates <- fev1_covariate_aa[,2:8]
# subset only the protein data
X1_covariates <- fev1_covariate_aa[,-c(1:16)]
# center and scale protein data
X1_covariates <- as.data.frame(scale(X1_covariates, center = TRUE, scale = TRUE))
# filter out proteins with low CV
X1_covariates <- X1_covariates[,cv_value >= 0.7]
# regress out sex, age, and clinical center
map_aa <- X1_covariates %>%   
  map(~lm(.x ~ 
            covariates$gender + 
            covariates$ccenter_P2 + covariates$Age_P2, 
          data = X1_covariates)) %>%
  map(resid)
Adjusted_AA <- map_df(map_aa, ~as.data.frame(t(.)))
Adjusted_AA <- t(Adjusted_AA)
colnames(Adjusted_AA) <- colnames(X1_covariates)

##################################### Run SmCCNet
# set phentoype
Y <- as.matrix(fev1_covariate_aa$smoking_status-1)
X1 <- Adjusted_AA
colnames(Y)[1] <- "Smoking_AA_Adjusted_1d"
pheno <- "Smoking_AA_Adjusted_1d"
# set number of folds
K <- 5
# set result-saving directory
CVDirs <- paste0(pheno,  "/", K, "foldCV", pheno, "_CCcoef_1_", 
                 "/")
# create directory
dir.create(pheno)
dir.create(CVDirs)
resultDir <- paste0(pheno, "/")
CVDir <- CVDirs
# load protein label metadata
label <- read.csv("SOMAscan_Assay_v4_Annotations_version3.3.2.csv")
newint <- intersect(gsub("_", "-",substring(colnames(X1), 2)), as.character(label$SeqId))
sum(gsub("_", "-",substring(colnames(X1), 2)) == newint)
label_data <- subset(label, SeqId %in% newint)
# match protein label metadata and the features of protein data
label_data$TargetFullName <- ifelse(label_data$TargetFullName == "", as.character(label_data$Target), as.character(label_data$TargetFullName))
label_data <- label_data %>%
  dplyr::select(SeqId, TargetFullName)
# use target full name instead as the protein data column names
colnames(X1) <- as.character(label_data$TargetFullName)
################################################## Source Code
# Run 5fold CV locally
library(parallel)
# set random seed
set.seed(12345)
# calculate global correlation matrix
bigCor2 <- cor(X1)
# set feature labels
AbarLabel <- colnames(bigCor2)
# set number of rows/columns
N <- nrow(X1); p1 <- ncol(X1)
# set subsampling percentage and number
s1 <- 0.7; subSamp <- 100
# set penalty candidates
pen1 <- seq(.1, .9, by = .1)
K <- 5
# save data for reproducibility purpose
save(X1, Y, s1, subSamp, pen1, 
     file = paste0(CVDir, "Data.Rdata"))
# split data into multiple folds
foldIdx <- split(1:N, sample(1:N, K))
for(i in 1:K){
  iIdx <- foldIdx[[i]]
  x1.train <- scale(X1[-iIdx, ])
  yy.train <- Y[-iIdx, ]
  x1.test <- scale(X1[iIdx, ])
  yy.test <- Y[iIdx, ]
  
  if(is.na(min(min(x1.train), min(yy.train), min(x1.test), min(yy.test)))){
    stop("Invalid scaled data.")
  }
  
  subD <- paste0(CVDir, "CV_", i, "/")
  dir.create(subD)
  save(x1.train, yy.train, x1.test, yy.test,
       s1, pen1, p1, subSamp,
       file = paste0(subD, "Data.Rdata"))
}
# open parallel computing clusters
cl <- makeCluster(5, type = "PSOCK")
# export variables/functions to the clusters
clusterExport(cl = cl, "CVDir")
clusterExport(cl = cl, 'getRobustPseudoWeight_binary')
clusterExport(cl = cl, "caret")
# run SmCCNet binary phenotype with parallel computing
parSapply(cl, 1:5, function(CVidx){
  # define sub saving directory
  subD <- paste0("Smoking_AA_Adjusted_1d/5foldCVSmoking_AA_Adjusted_1d_CCcoef_1_/", "CV_", CVidx, "/")
  # load fold data
  load(paste0(subD, "Data.Rdata"))
  # create directory
  dir.create(paste0(subD, "SmCCA/"))
  # create empty vectors to store prediction accuracy
  AccTrain <- AccTest <- DeltaAcc <- rep(0, length(pen1))
  
  for(idx in 1:length(pen1)){
    # select a specified penalty parameter candidate
    l1 <- pen1[idx]
    # run splsda
    Ws <- spls::splsda(x = x1.train, y = yy.train, K = 3, eta = l1, kappa=0.5,
                       classifier=c('lda'), scale.x=FALSE)
    # create empty matrix to store the feature importance weights
    weight <- matrix(0,nrow = ncol(x1.train), ncol = 3)
    # store feature importance weight
    weight[Ws[["A"]],] <- Ws[["W"]]
    # define training data
    train_data <- data.frame(x = (x1.train %*% weight)[,1:3], y = as.factor(yy.train))
    # define testing data
    test_data <- data.frame(x = (x1.test %*% weight)[,1:3])
    # run lda classifier
    ldaFit <- caret::train(
      y ~.,
      data = train_data,
      method = "lda"
    )
    # make prediction on the test set
    test_pred <- predict(ldaFit, test_data)
    # obtain training accuracy
    acc.train <- ldaFit[["results"]][["Accuracy"]]
    # obtain testing accuracy
    acc.test <- mean(as.factor(test_pred) == as.factor(yy.test))  
    # store evaluation result
    AccTrain[idx] <- round(acc.train, digits = 5)
    AccTest[idx] <- round(acc.test, digits = 5)
    DeltaAcc[idx] <- abs(acc.train - acc.test)
    
    if(idx %% 10 == 0){
      save(pen1, AccTrain, AccTest, DeltaAcc, idx,
           file = paste0(subD, "temp.Rdata"))
    }
    
  }
  # combine all evaluation result into a dataframe
  DeltaAcc.all <- cbind(pen1, AccTrain, AccTest, DeltaAcc)
  colnames(DeltaAcc.all) <- c("l1", "Training ACC", "Test ACC", "CC Pred. Error")
  # store cross-validataion result into local directory
  write.csv(DeltaAcc.all,
            file = paste0(subD, "SmCCA/SCCA_", subSamp,"_allDeltaCor.csv"))
  
  system(paste0("rm ", subD, "temp.Rdata"))
  
  return(CVidx)
  
}
)



# Close cluster
stopCluster(cl)
# choose the best penalty term
plotCVcontour(CVDir, "SmCCA", NumSubsamp = subSamp)
# define result-saving directory
plotD <- paste0(CVDir, "Figures/")
saveD <- paste0(CVDir, "Results/")
dataF <- paste0(CVDir, "Data.Rdata")
# create result-saving directory
dir.create(plotD)
dir.create(saveD)
dir.create(dataF)
# set method to SmCCA
Method = "SmCCA"
# run SmCCNet with the best penalty term
for(Method in "SmCCA"){
  pen <- pen1[which.max(T12[ ,3])]
  pen <- 0.8
  
  if(Method == "SmCCA"){
    FilterByTrait <- FALSE
  }else if(Method == "SsCCA"){
    FilterByTrait <- TRUE
  }
  # set best penalty term
  l1 <- pen;
  system.time({
    # run SPLSDA and store the canonical weights
    Ws <- getRobustPseudoWeight_binary(X1 = X1, Trait = Y, Lambda1 = l1, 
                                       s1, SubsamplingNum = subSamp)
    # obtain adjacency matrix
    Abar <- getAbars(Ws, P1 = p1, FeatureLabel = AbarLabel[1:p1])
    # save adjacency matrix into the local directory
    save(l1, X1, Y, s1, Ws, Abar,
         file = paste0(saveD, Method, K, "foldSamp", subSamp, "_", pen, ".Rdata"))
  })
  
  
}
# load global network
load('D:/Research/SOMA/Smoking_NHW_Adjusted_1d/5foldCVSmoking_NHW_Adjusted_1d_CCcoef_1_/Results/SmCCA5foldSamp100_0.8.Rdata')
# clustering
mirGeneModule <- getMultiOmicsModules(Abar, p1, PlotTree = FALSE)
AbarLabel <- colnames(X1)
AbarLabel <- make.unique(AbarLabel)
row.names(Abar) <- colnames(Abar) <- AbarLabel
# global correlation matrix
bigCor2 <- cor(X1)
row.names(bigCor2) <- colnames(bigCor2) <- AbarLabel
colnames(X1) <- AbarLabel
# store subnetwork with each specified size to the local directory
cor_result <- list()
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
      netD <- paste0("D:/Research/SOMA/Smoking_NHW_Adjusted_1d/5foldCVSmoking_NHW_Adjusted_1d_CCcoef_1_/", "NetworkEdgeCutpt0/")
      titlePre <- paste0("Net ")
      savePre <- paste0(netD, "Net")
    }else{
      netD <- paste0("D:/Research/SOMA/Smoking_NHW_Adjusted_1d/5foldCVSmoking_NHW_Adjusted_1d_CCcoef_1_/", "NetworkEdgeCutpt", size, "/")
      titlePre <- paste0("Trimmed Net ")
      savePre <- paste0(netD, "TrimmedNet")
    }
    dir.create(netD)
    saveplot.z <- paste0(savePre, z, ".pdf")
    plottitle <- paste0(titlePre, z)
    #size <- cor_result[[z]]$size[i]  
    results <- plotMultiOmicsNetworks_PPR_Validation(Abar, bigCor2, mirGeneModule, 
                                                     data = X1, Pheno =Y, ModuleIdx = z,
                                                     type ='Smoking_NHW_Adjusted_1d', 
                                                     id = pheno_data$sid,  P1 = p1,
                                                     folder = 'Smoking_NHW_Adjusted_1d',
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
            paste0('D:/Research/SOMA/Smoking_AA_Adjusted_1d/5foldCVSmoking_NHW_Adjusted_1d_CCcoef_1_/pc_correlation_validation_', z, '.csv'),
            row.names = FALSE)
}
# network pruning
for (i in 10:34)
{
  temp_size <- i
  # load subnetwork with specified network size
  load(paste0('D:/Research/SOMA/Smoking_AA_Adjusted/Smoking_NHW_Adjustednet', temp_size, '_3AA.Rdata'))
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
createNetworkFromIgraph(graph,"smoking_nhw")