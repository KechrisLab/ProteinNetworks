######## Package that we need ########
# You only need to install once #
# install.packages("Matrix", repos="http://cran.r-project.org")
# install.packages("pbapply", repos="http://cran.r-project.org")
# install.packages("igraph", repos="http://cran.r-project.org")

requireNamespace("Matrix", quietly = TRUE)
requireNamespace("pbapply", quietly = TRUE)
requireNamespace("igraph", quietly = TRUE)



################################################################################
### Apply sparse multiple canonical correlation analysis to omics feature
### subsamples

#' Calculate the canonical correlation weights based on sparse multiple
#' canonical correlation analysis (SmCCA), sparse supervised canonical 
#' correlation analysis (SsCCA), or sparse canonical correlation analysis (SCCA).
#' 
#' Integrate two omics data type (and a quantitative phenotype), and calculate
#' the absolute canonical correlation weights for the omics features using SmCCA
#' SsCCA, or SCCA. SmCCA and SsCCA take into account a phenotype/trait. SmCCA 
#' maximizes the total (weighted or unweighted) pairwise canonical correlation 
#' weights between two omics data types and the trait. It requires the trait to 
#' be quantitative. SsCCA prioritizes omics features based on the trait, and 
#' assigns non-zero canonical weights to features that are more correlated to 
#' the trait. SCCA does not use any trait information for computing the 
#' canonical correlation weights. All of these three methods are included in 
#' this function, along with an omics feature subsampling scheme. 
#' 
#' To choose SmCCA, set \code{NoTrait = FALSE, FilterByTrait = FALSE}.  
#' To choose SsCCA, set \code{NoTrait = FALSE, FilterByTrait = TRUE}.
#' To choose SCCA, set \code{Trait = NULL, NoTrait = TRUE}.
#'
#' @param X1 An \eqn{n\times p_1} data matrix (e.g. mRNA) with \eqn{p_1} 
#'   features and \eqn{n} subjects.
#' @param X2 An \eqn{n\times p_2} data matrix (e.g. miRNA) with \eqn{p_2} 
#'   features and \eqn{n} subjects.
#' @param Trait An \eqn{n\times 1} trait data matrix for the same n subjects.
#' @param Lambda1 LASSO penalty parameter for \code{X1}. \code{Lambda1} needs
#'   to be between 0 and 1.
#' @param Lambda2 LASSO penalty parameter for \code{X2}. \code{Lambda2} needs 
#'   to be between 0 and 1.
#' @param s1 Proportion of mRNA features to be included, default at \code{s1 = 0.7}.
#'   \code{s1} needs to be between 0 and 1.
#' @param s2 Proportion of miRNA features to be included, default at \code{s1 = 0.7}.
#'   \code{s2} needs to be between 0 and 1.
#' @param NoTrait Logical, default is \code{FALSE}. Whether trait information is
#'   provided.
#' @param FilterByTrait Logical, default is \code{FALSE}. Whether only the top 
#'   (\eqn{80\%}) features with highest correlation to the trait will be assigned 
#'   nonzero weights. The choice of \eqn{80\%} is based on the PMA package.
#' @param SubsamplingNum Number of feature subsamples. Default is 1000. Larger
#'   number leads to more accurate results, but at a higher cost. 
#' @param CCcoef Optional coefficients for the SmCCA pairwise canonical 
#'   correlations. If \code{CCcoef = NULL} (default), then the objective 
#'   function is the total sum of all pairwise canonical correlations. It can 
#'   also be a coefficient vector that follows the column order of 
#'   \code{combn(K, 2)}.
#' @param trace Logical. Whether to display the CCA algorithm trace.
#' @return A canonical correlation weight matrix with \eqn{p_1+p_2} rows. Each 
#'   column is the canonical correlation weights based on subsampled \code{X1}
#'   and \code{X2} features. The number of columns is \code{SubsamplingNum}.
#'
#' @examples
#' 
#' ## For illustration, we only subsample 5 times.
#' set.seed(123)
#'
#' # Unweighted SmCCA
#' W1 <- getRobustPseudoWeights(X1, X2, Trait = Y, Lambda1 = 0.05,
#'   Lambda2 = 0.05, s1 = 0.7, s2 = 0.9, NoTrait = FALSE, FilterByTrait = FALSE,
#'   SubsamplingNum = 5, CCcoef = NULL, trace = FALSE)
#'   
#' # Weighted SmCCA
#' W2 <- getRobustPseudoWeights(X1, X2, Trait = Y, Lambda1 = 0.05,
#'   Lambda2 = 0.05, s1 = 0.7, s2 = 0.9, NoTrait = FALSE, FilterByTrait = FALSE,
#'   SubsamplingNum = 5, CCcoef = c(1, 5, 5), trace = FALSE)
#'   
#' # SsCCA
#' W3 <- getRobustPseudoWeights(X1, X2, Trait = Y, Lambda1 = .05, Lambda2 = 0.5,
#'   s1 = 0.7, s2 = 0.9, NoTrait = FALSE, FilterByTrait = TRUE,
#'   SubsamplingNum = 5, CCcoef = NULL, trace = FALSE)
#'
#' # SCCA
#' W4 <- getRobustPseudoWeights(X1, X2, Trait = NULL, Lambda1 = 0.05,
#'   Lambda2 = 0.05, s1 = 0.7, s2 = 0.9, NoTrait = TRUE, 
#'   SubsamplingNum = 5, CCcoef = NULL, trace = FALSE)
#'   
#'
#' @export
getRobustPseudoWeight_one <- function(X1, Trait, Lambda1, 
                                      s1 = 0.7, NoTrait = FALSE, 
                                      FilterByTrait = FALSE, SubsamplingNum = 1000,
                                      CCcoef = NULL, trace = FALSE){
  
  
  
  p1 <- ncol(X1); p <- p1 
  p1.sub <- ceiling(s1 * p1)
  X <- X1
  
  beta <- pbapply::pbsapply(seq_len(SubsamplingNum), function(x){
    # Subsample features
    samp1 <- sort(sample(seq_len(p1), p1.sub, replace = FALSE))
    
    
    x1.par <- scale(X1[ , samp1], center = TRUE, scale = TRUE)
    
    
    out <- getCCAout(x1.par, Trait, Lambda1,
                     NoTrait = NoTrait, FilterByTrait = FilterByTrait,
                     trace = trace)
    ### ???
    w <- rep(0, p)
    w[samp1] <- out$u
    coeff.avg <- w
    
    return(coeff.avg)
  })
  
  return(beta)
}




################################################################################
### Aggregate pseudo canonical weights and compute the network similarity matrix
### for all omics features.

#' Compute the similarity matrix based on one or more canonical correlation 
#' weight vectors.
#' 
#' Compute the similarity matrix based on the outer products of absolute 
#' canonical correlation weights.
#'
#'
#' @param Ws A canonical correlation weight vector or matrix. If \code{Ws} is a
#'   matrix, then each column corresponds to one weight vector.
#' @param P1 Total number of features for the first omics data type. 
#' @param FeatureLabel If \code{FeatureLabel = NULL} (default), the feature 
#'   names will be \eqn{\{TypeI_1, \cdots, TypeI_{p_1}, TypeII_1, \cdots, TypeII_{p-p_1}\}}, 
#'   where \eqn{p_1 = }\code{P1}, and \eqn{p} is the total number of omics features.
#' @return A \eqn{p\times p} symmetric non-negative matrix.
#'
#' @examples
#' w <- matrix(rnorm(6), nrow = 3)
#' Ws <- apply(w, 2, function(x)return(x/sqrt(sum(x^2))))
#' abar <- getAbar(Ws, P1 = 2, FeatureLabel = NULL)
#'
#' @export
getAbars <- function(Ws, P1 = NULL, P2 = NULL, FeatureLabel = NULL){
  
  
  if(is.null(dim(Ws))){
    Abar <- Matrix::Matrix(abs(Ws) %o% abs(Ws), sparse = TRUE)
  }else{
    b <- nrow(Ws)
    Abar <- matrix(0, nrow = b, ncol = b)
    for(ind in seq_len(ncol(Ws))){
      w <- abs(Ws[ , ind])
      A <- Matrix::Matrix(w %o% w, sparse = TRUE)
      Abar <- Abar + A
    }
  }
  
  diag(Abar) <- 0
  Abar <- Abar/max(Abar)
  
  if(is.null(colnames(Abar))){
    if(is.null(FeatureLabel)){
      if(is.null(P1)){
        stop("Need to provide FeatureLabel or the number of features 
                    for the first data type P1.")
      }else{
        p <- ncol(Abar)
        FeatureLabel <- c(paste0("TypeI_", seq_len(P1)), 
                          paste0("TypeII_", seq_len(P2)), 
                          paste0("TypeIII_", seq_len(p-P1-P2)))
      }
    }
    colnames(Abar) <- rownames(Abar) <- FeatureLabel
  }
  
  return(Abar)
}


################################################################################
### Extract multi-omics modules.

#' Extract multi-omics modules based on the similarity matrix.
#' 
#' Apply hierarchical tree cutting to the similarity matrix and extract
#' modules that contain both omics data types.
#'
#' @param Abar A similary matrix for all features (both omics data types).
#' @param P1 Total number of features for the first omics data type.
#' @param CutHeight Height threshold for the hierarchical tree cutting. Default 
#'   is \eqn{1-0.1^{10}}.
#' @param PlotTree Logical. Whether to create a hierarchical tree plot.
#' @return A list of multi-omics modules.
#'
#' @examples
#' set.seed(123)
#' w <- rnorm(5)
#' w <- w/sqrt(sum(w^2))
#' abar <- getAbar(w, P1 = 2, FeatureLabel = NULL)
#' modules <- getMultiOmicsModules(abar, P1 = 2, CutHeight = 0.5)
#'
#' @export
getMultiOmicsModules <- function(Abar, P1, CutHeight = 1-.1^10, PlotTree = TRUE){
  
  hc <- stats::hclust(stats::as.dist(1 - Abar))
  if(PlotTree){graphics::plot(hc)}
  cut.merge <- hc$merge[hc$height < CutHeight, ]
  lower.leaves <- sort(-cut.merge[cut.merge<0])
  
  grpID <- stats::cutree(hc, h = CutHeight)
  id <- grpID[lower.leaves]
  M <- lapply(seq_len(length(unique(id))), function(x){
    M.x <- lower.leaves[which(id == unique(id)[x])]
    return(M.x)
  })
  
  multiOmicsModule <- lapply(M, function(s){
    s.min <- min(s)
    s.max <- max(s)
    if(s.min <= P1 & s.max > P1)return(s)
  })
  
  if(length(multiOmicsModule) > 1){
    nullSet <- which(vapply(multiOmicsModule, is.null, logical(1)))
    if(length(nullSet) > 0){
      multiOmicsModule <- multiOmicsModule[-nullSet]
    }
  }
  
  return(multiOmicsModule)
}


################################################################################
### Visualize multi-omics subnetworks.

#' Plot multi-omics module networks.
#' 
#' Plot multi-omics modules based on similarity matrix derived from pseudo
#' canonical weights and pairwise feature correlations.
#'
#' @param Abar A \eqn{p\times p} similary matrix for both omics data types
#'   based on pseudo canonical correlation weights. \eqn{p} is the number of 
#'   total features for the two omics data types. All entries are non-negative.
#' @param CorrMatrix A \eqn{p\times p} correlation matrix that provides sign
#'   information for the network.
#' @param multiOmicsModule A list of multi-omics modules.
#' @param ModuleIdx Index for the module to be plotted. It can not exceed the
#'   length of \code{multiOmicsModule}.
#' @param P1 Total number of features for the first omics data type.
#' @param EdgeCut A numerical value between 0 and 1, indicating an edge
#'   threshold for the network. Any features (network nodes) without any edge
#'   strength that passes the threshold are excluded from the figure. If 
#'   \code{EdgeCut = 0} (default), then the full module network will be created. 
#' @param FeatureLabel A \eqn{1\times p} vector indicating feature names. If
#'   \code{FeatureLabel = NULL} (default), the feature names will be 
#'   \eqn{\{TypeI_1, \cdots, TypeI_{p_1}, TypeII_1, \cdots, TypeII_{p-p_1}\}}, 
#'   where \eqn{p_1 = }\code{P1}.
#' @param AddCorrSign Logical. Whether to add a positive or negative sign to
#'   each network edge based on pairwise feature correlations.
#' @param SaveFile A pdf file name for the figure output. 
#'   If \code{SaveFile = NULL} (default), the figure will not be saved.
#' @param ShowType1Label Logical. Whether to label the network nodes for the
#'   first omics data type.
#' @param ShowType2Label Logical. Whether to label the network nodes for the
#'   second omics data type.
#' @param PlotTitle A title for the figure. Default is without any title.
#' @param NetLayout Graphical layout for the network. Possible options are
#'   \code{circle} for circle layout, \code{sphere} for 3D sphere, \code{fr} for
#'   Fruchterman-Reinhold, and \code{lgl} for the LGL algorithm. Refer to igraph 
#'   manual for more details on the layout options.
#' @param ShowNodes Logical. Whether to show network nodes.
#' @param VertexLabelCex Scaling factor for the vertex labels.
#' @param VertexSize Size of the vertices.
#' @return A multi-omics network figure.
#'
#' @examples
#' set.seed(123)
#' w <- rnorm(5)
#' w <- w/sqrt(sum(w^2))
#' abar <- getAbar(w, P1 = 2, FeatureLabel = NULL)
#' modules <- getMultiOmicsModules(abar, P1 = 2, CutHeight = 0.5)
#' x <- cbind(X1[ ,seq_len(2)], X2[ , seq_len(3)])
#' corr <- cor(x)
#'
#' plotMultiOmicsNetwork(abar, corr, modules, ModuleIdx = 1, P1 = 2)
#'
#' @export
plotMultiOmicsNetworks <- function(Abar, CorrMatrix, multiOmicsModule,
                                   ModuleIdx, P1, P2, EdgeCut = 0, FeatureLabel = NULL,
                                   AddCorrSign = TRUE, SaveFile = NULL,
                                   ShowType1Label = TRUE, ShowType2Label = TRUE,ShowType3Label = TRUE,
                                   PlotTitle = "", NetLayout = "lgl",
                                   ShowNodes = TRUE,
                                   VertexLabelCex = 1, VertexSize = 1){
  
  p <- ncol(Abar)
  if(is.null(FeatureLabel)){
    FeatureLabel <- c(paste0("TypeI_", seq_len(P1)), paste0("TypeII_", seq_len(P2)),paste0("TypeII_", seq_len(p-P1-P2)))
  }
  colnames(Abar) <- rownames(Abar) <- FeatureLabel[seq_len(p)]
  colnames(CorrMatrix) <- rownames(CorrMatrix) <- FeatureLabel[seq_len(p)]
  
  grp <- multiOmicsModule[[ModuleIdx]]
  grp.memb <- colnames(Abar)[grp]
  M.node <- grp.memb
  
  # Trim the module by EdgeCut.
  M <- as.matrix(Abar[M.node, M.node])
  if(AddCorrSign){M <- M * sign(CorrMatrix[M.node, M.node])}
  M[which(abs(M) < EdgeCut)] <- 0
  newM.node <- M.node[which(apply(abs(M), 1, max) > 0)]
  
  if(length(newM.node) == 0){
    print("No edge passes threshold.")
  }else{
    M <- M[newM.node, newM.node]
    allidx <- matrix(seq_len(p), ncol = 1)
    rownames(allidx) <- rownames(Abar)
    
    NodeInfo <- data.frame(id = newM.node, idx = allidx[newM.node, ])
    net <- igraph::graph_from_adjacency_matrix(M, weighted = TRUE,
                                               diag = FALSE, mode = "undirected")
    
    # Define colors and shapes for vertex and label.
    k <- length(newM.node)
    type1 <- which(NodeInfo$idx <= P1)
    type2 <- which(NodeInfo$idx > P1 & NodeInfo$idx <= P1 + P2)
    type3 <- which(NodeInfo$idx > P1 + P2)
    vcol <- rep("purple", k); vcol[type2] <- "dark orange";vcol[type3] <- "green"
    vshape <- rep("square", k); vshape[type2] <- "circle"; vshape[type3] <- "circle"
    lcol <- vcol
    # If not show nodes.
    if(!ShowNodes){vshape <- "none"}
    # If no label, assign a place holder.
    if(!ShowType1Label){newM.node[type1] <- " "}
    if(!ShowType2Label){newM.node[type2] <- " "}
    if(!ShowType3Label){newM.node[type3] <- " "}
    
    # Define edge colors.
    ecol <- rep("gray80", igraph::ecount(net))
    ew <- abs(igraph::edge.attributes(net)$weight) * 5
    ecol[which(igraph::edge.attributes(net)$weight < 0)] <- "red"
    
    # Define network layout.
    if(NetLayout == "circle"){
      l <- igraph::layout_in_circle(net)
    }else if(NetLayout == "sphere"){
      l <- igraph::layout_on_sphere(net)
    }else if(NetLayout == "fr"){
      l <- igraph::layout_with_fr(net)
    }else if(NetLayout == "lgl"){
      l <- igraph::layout_with_lgl(net)
    }else{
      stop("Unrecognized NetLayout input. Acceptable options are 'circle',
                 'sphere', 'fr', and 'lgl'.")
    }
    
    
    if(!is.null(SaveFile)){
      grDevices::pdf(SaveFile)
      graphics::par(bg = "white")
      graphics::plot(net, vertex.color = vcol, vertex.shape = vshape,
                     vertex.label.cex = VertexLabelCex, layout = l,
                     vertex.size = VertexSize, vertex.label = newM.node,
                     edge.width = ew, vertex.label.color = lcol, edge.color = ecol,
                     vertex.label.font = 2)
      graphics::title(PlotTitle)
      grDevices::dev.off()
    }else{
      graphics::plot(net, vertex.color = vcol, vertex.shape = vshape,
                     vertex.label.cex = VertexLabelCex, layout = l,
                     vertex.size = VertexSize, vertex.label = newM.node,
                     edge.width = ew, vertex.label.color = lcol, edge.color = ecol,
                     vertex.label.font = 2, main = PlotTitle)
    }
    
  }
  return(M)
}









#######################################################
# Internal functions called by getRobustPseudoWeights #
#######################################################

# (INTERNAL)
getCCAout <- function(X1, Trait, Lambda1, CCcoef = NULL,
                      NoTrait = FALSE, FilterByTrait = FALSE, trace = FALSE){
  # Compute CCA weights.
  #
  # X1: An n by p1 mRNA expression matrix.
  # X2: An n by p2 miRNA expression matrix.
  # Trait: An n by k trait data for the same samples (k >=1).
  # Lambda1, Lambda2: LASSO pentalty parameters, need to be between 0 and 1.
  # CCcoef: A 3 by 1 vector indicating weights for each pairwise canonical
  #   correlation.
  # NoTrait: Logical. Whether trait information is provided.
  # FilterByTrait: Logical. Whether only the features with highest correlation
  #   to Trait will be assigned nonzero weights. The top 80% features are reserved.
  # trace: Logical. Whether to display CCA algorithm trace.
  
  if(abs(Lambda1 - 0.5) > 0.5){
    stop("Invalid penalty parameter. Lambda1 needs to be between zero and
           one.")}
  
  
  k <- ncol(Trait)    
  if(NoTrait | is.null(k)){
    out <- PMA::CCA(X1, X2, typex = "standard", typez = "standard", 
                    penaltyx = Lambda1, penaltyz = Lambda2, K = 1, 
                    trace = trace)
  }else{
    if(FilterByTrait){
      if(k > 1){
        stop("'FilterByTrait == TRUE' only allows one trait at a time.")
      }else{
        out <- PMA::CCA(X1, X2, outcome = "quantitative", y = Trait,
                        typex = "standard", typez = "standard", 
                        penaltyx = Lambda1, penaltyz = Lambda2, K = 1, 
                        trace = trace)
      }
    }else{
      xlist <- list(x1 = X1, y = scale(Trait))
      L1 <- max(1, sqrt(ncol(X1)) * Lambda1)
      out <- myMultiCCA(xlist, penalty = c(L1, sqrt(ncol(Trait))), trace = trace)
      out$u <- out$ws[[1]]
    }
  }
  
  return(out)
}

################################################# Additional source

requireNamespace("PMA", quietly = TRUE)
# library(PMA)



################################################################################
### Main function for performing SmCCA

# (INTERNAL)
myMultiCCA <- function(xlist, penalty=NULL, ws=NULL, niter=25,
                       type="standard", ncomponents=1, trace=TRUE,
                       standardize=TRUE, CCcoef = NULL){
  # Perform sparse multiple canonical correlation analysis (SmCCA) for a list
  #   of data inputs.
  #
  # xlist: A list of length K, where K is the number of data sets on which to
  #   perform SmCCA. Dataset k should be a matrix of dimension $n \times p_k$,
  #   where $p_k$ is the number of features in data set k. If some quantitative
  #   phenotype information is included in the list, it should be the last
  #   element of the list and of dimension $n \times 1$.
  # penalty: The penalty terms to be used. Can be a single value (if the same
  #   penalty term is to be applied to each data set) or a K-vector, indicating
  #   a different penalty term for each data set. There are 2 possible
  #   interpretations for the penalty term: If type = "standard" then this is
  #   an L1 bound on $w_k$, and it must be between 1 and $\sqrt(p_k)$ ($p_k$ is
  #   the number of features in matrix $X_k$). If type = "ordered" then this is
  #   the parameter for the fussed lasso penalty on $w_k$.
  # type: Either "standard" or "ordered" to indicate if the columns are
  #   unordered or ordered. If type = "standard", then a lasso penalty is aplied
  #   to v, to enforce sparsity. If type = "ordered" (generally used for CGH
  #   data), then a fused lasso penalty is applied, to enforce both sparsity and
  #   soothness. This argument can be a vector of lenght K (if different data
  #   sets are of different types) or it can be a single value "standard" or
  #   "ordered" (if all data sets are of the same type).
  # ncomponents: Number of factors to output. Default is 1.
  # niter: Number of iterations to be perfromed. Default is 25.
  # ws: A list of length K. The kth element contains the first ncomponents
  #   columns of the v matrix of the SVD of $X_k$. If NULL, then the SVD of
  #   $X_1, ..., X_K$ will be computed inside the MultiCCA function. However,
  #   if you plan to run this function multiple times, then save a copy of the
  #   argument so that it does not need to be re-computed.
  # trace: Logical. Whether to print out progress.
  # standardize: Whether to center and scale the columns of $X_1, ..., X_K$.
  #   Default is TRUE.
  # CCcoef: Optional coefficients for the pairwise canonical correlations (CC).
  #   If CCcoef = NULL (default), then the objective function is the total sum
  #   of all pairwise CC. It can also be a coefficient vector that follows the
  #   column order of combn(K, 2).
  
  if(ncol(xlist[[length(xlist)]]) > 1){
    out <- PMA::MultiCCA(xlist, penalty=penalty, ws=ws, niter=niter,
                         type=type, ncomponents=ncomponents, trace=ncomponents,
                         standardize=standardize)
    
  }else{
    # The Kth data set in xlist is a one dimensional phenotype
    call <- match.call()
    K <- length(xlist)
    pair_CC <- utils::combn(K, 2)
    num_CC <- ncol(utils::combn(K, 2))
    
    # Check canonical correlation coefficients.
    if(is.null(CCcoef)){CCcoef <- rep(1, num_CC)}
    if(length(CCcoef) != num_CC){
      stop(paste0("Invalid coefficients for pairwise canonical correlations.
                  Please provide ", num_CC, " numerical values for CCcoef.
                  Be sure to match combn(K,2) column order."))
    }
    
    # Check data type.
    if(length(type)==1){# Expand a single type to a vector of length(xlist).
      if(type != "standard"){
        stop("The phenotype data must be continuous and follow the type 'standard'. ")
      }
      type <- rep(type, K)}
    if(length(type)!=K){
      stop("Type must be a vector of length 1, or length(xlist)")}
    if(sum(type!="standard")>0){
      stop("Each element of type must be standard and not ordered.")}
    
    # Standardize or not.
    if(standardize){xlist <- lapply(xlist, scale)}
    
    # Initialize weights.
    if(!is.null(ws)){
      makenull <- FALSE
      for(i in seq_len(K-1)){
        if(ncol(ws[[i]])<ncomponents) makenull <- TRUE
      }
      if(makenull) ws <- NULL
    }
    if(is.null(ws)){
      ws <- list()
      for(i in seq_len(K-1)){
        ws[[i]] <- matrix(svd(xlist[[i]])$v[,seq_len(ncomponents)], ncol=ncomponents)
      }
      ws[[K]] <- 1
    }
    ws.init <- ws
    
    # Check penalties.
    if(is.null(penalty)){
      penalty <- rep(NA, K)
      penalty[type=="standard"] <- 4 # this is the default value of sumabs
      for(k in seq_len(K-1)){
        if(type[k]=="ordered"){
          stop("Current version requires all element types to be standard (not ordered).")
        }
      }
    }
    if(length(penalty)==1) penalty <- rep(penalty, K)
    if(sum(penalty<1 & type=="standard")){
      stop("Cannot constrain sum of absolute values of weights to be less than
           1.")
    }
    for(i in seq_len(K-1)){
      if(type[i]=="standard" && penalty[i]>sqrt(ncol(xlist[[i]]))){
        stop("L1 bound of weights should be no more than sqrt of the number of
             columns of the corresponding data set.", fill=TRUE)
      }
    }
    
    ws.final <- ws.init
    for(i in seq_len(K-1)){
      ws.final[[i]] <- matrix(0, nrow=ncol(xlist[[i]]), ncol=ncomponents)
    }
    cors <- NULL
    for(comp in seq_len(ncomponents)){
      ws <- list()
      for(i in seq_len(K-1)) ws[[i]] <- ws.init[[i]][,comp]
      ws[[K]] <- 1
      curiter <- 1
      crit.old <- -10
      crit <- -20
      storecrits <- NULL
      
      while(curiter<=niter && abs(crit.old-crit)/abs(crit.old)>.001 &&
            crit.old!=0){
        crit.old <- crit
        crit <- myGetCrit(xlist, ws, pair_CC, CCcoef)
        storecrits <- c(storecrits,crit)
        if(trace) cat(curiter, fill=FALSE)
        curiter <- curiter+1
        for(i in seq_len(K-1)){
          ws[[i]] <- myUpdateW(xlist, i, K, penalty[i], ws, type[i], ws.final,
                               pair_CC, CCcoef)
        }
      }
      
      for(i in seq_len(K-1)) ws.final[[i]][,comp] <- ws[[i]]
      cors <- c(cors, myGetCors(xlist, ws, pair_CC, CCcoef))
    }
    
    out <- list(ws=ws.final, ws.init=ws.init, K=K, call=call, type=type,
                penalty=penalty, cors=cors)
    class(out) <- "MultiCCA"
  }
  return(out)
}





############################################
# Internal functions called by myMultiCCA. #
############################################

# (INTERNAL)
myUpdateW <- function(xlist, i, K, sumabsthis, ws, type="standard", ws.final,
                      pair_CC, CCcoef){
  # Update canonical weights for the ith data set.
  #
  # xlist: Data list.
  # i: Data set index.
  # K: Total number of data sets.
  # sumabsthis: Penalty for the ith data set.
  # ws: First ncomponents columns of the v matrix of the SVD of $X_k$'s.
  # type: Type of data sets.
  # ws.final: Current weight matrix.
  # pair_CC: The output of combn(K, 2). A two-row table that include indices for
  #   pairwise canonical correlaltions between members of xlist.
  # CCcoef: Optional coefficients for the pairwise canonical correlations (CC).
  
  
  tots0 <- vapply(seq_len(length(CCcoef)), function(x){
    pairx <- pair_CC[ , x]
    Xi <- xlist[[i]]
    
    if(pairx[1] != i & pairx[[2]] != i){
      y <- rep(0, ncol(Xi))
    }else{
      if(pairx[1] == i){j <- pairx[2]}
      if(pairx[2] == i){j <- pairx[1]}
      Xj <- xlist[[j]]
      
      # diagmat is the diagonal correlation matrix calculated using previous
      #   canonical directions.
      # If phenotype is included, only the first canonical direction is used.
      #   diagmat is therefore a zero matrix.
      diagmat <- (t(ws.final[[i]])%*%t(Xi))%*%(Xj%*%ws.final[[j]])
      diagmat[row(diagmat)!=col(diagmat)] <- 0
      y <- t(Xi)%*%(Xj%*%ws[[j]]) -
        ws.final[[i]]%*%(diagmat%*%(t(ws.final[[j]])%*%ws[[j]]))
      y <- y * CCcoef[x]
    }
    
    return(y)
  }, numeric(ncol(xlist[[i]])))
  tots <- rowSums(tots0)
  
  if(type=="standard"){
    sumabsthis <- BinarySearch(tots, sumabsthis)
    w <- soft(tots, sumabsthis)/l2n(soft(tots, sumabsthis))
  } else {
    stop("Current version requires all element types to be standard (not ordered).")
  }
  return(w)
}


# (INTERNAL)
myGetCrit <- function(xlist, ws, pair_CC, CCcoef){
  # Compute the matrix form SmCCA objective function value for given weights.
  #
  # ws: First ncomponents columns of the v matrix of the SVD of $X_k$'s.
  
  
  
  crits <- apply(pair_CC, 2, function(x){
    i <- x[1]
    j <- x[2]
    y <- t(ws[[i]])%*%t(xlist[[i]])%*%xlist[[j]]%*%ws[[j]]
    return(y)
  })
  crit <- sum(crits * CCcoef)
  return(crit)
}


# (INTERNAL)
myGetCors <- function(xlist, ws, pair_CC, CCcoef){
  # Compute total weighted canonical correlations for given weights.
  #
  # xlist: Data list.
  # pair_CC: The output of combn(K, 2). A two-row table that include indices for
  #   pairwise canonical correlaltions between members of xlist.
  # CCcoef: Optional coefficients for the pairwise canonical correlations (CC).
  
  CCs <- apply(pair_CC, 2, function(x){
    i <- x[1]
    j <- x[2]
    y <- stats::cor(xlist[[i]]%*%ws[[i]], xlist[[j]]%*%ws[[j]])
    if(is.na(y)){y <- 0}
    return(y)
  })
  
  Cors <- sum(CCs * CCcoef)
  return(Cors)
}


# (INTERNAL)
BinarySearch <- function(argu,sumabs){
  #  Update sumabs so that the L1 norm of argu equals given penalty.
  
  if(l2n(argu)==0 || sum(abs(argu/l2n(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter < 150){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/l2n(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    if((lam2-lam1)<1e-6) return((lam1+lam2)/2)
    iter <- iter+1
  }
  warning("Didn't quite converge")
  return((lam1+lam2)/2)
}


# (INTERNAL)
l2n <- function(vec){
  # Computes the L2 norm. If the norm is zero, set it to 0.05.
  
  a <- sqrt(sum(vec^2))
  if(a==0) a <- .05
  return(a)
}

# (INTERNAL)
soft <- function(x,d){
  # Soft thresholding.
  
  return(sign(x)*pmax(0, abs(x)-d))
}


##################################################
# Additional internal functions for ordered data #
##################################################
# ChooseLambda1Lambda2()
# FLSA()

# source code for NetSHy summarization
hybrid_score = function(X, A, pc_id=1, is_alpha = TRUE){
  #X: data matrix (n, p)
  #A: corresponding adjacency matrix
  #pc_id: PC index
  g = graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  # laplacian
  L2 <- A %>% 
    graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE,diag = FALSE) %>% 
    graph.laplacian(normalized = F) %>% as.matrix
  # TOM
  TOM_A = TOMsimilarity(A, verbose = FALSE)
  if (is_alpha == TRUE)
  {  alpha = graph.density(g)}
  else
  {
    alpha = 0
  }
  X= scale(X, center = TRUE, scale = TRUE)
  # weighted approach
  temp  = (1-alpha)*(X %*% L2) + alpha*(X %*% TOM_A)
  temp = prcomp(temp)
  h_score = temp$x[,pc_id] 
  print(alpha)
  return(h_score)
  
}

# obtain subnetworks for each specified network size based on PageRank algorithm
plotMultiOmicsNetworks_PPR_Validation <- function(Abar, CorrMatrix, multiOmicsModule, data = X1,
                                                  Pheno, type, id, folder,
                                                  ModuleIdx, P1,  FeatureLabel = NULL,
                                                  penalty = 0,
                                                  AddCorrSign = TRUE, SaveFile = NULL,
                                                  ShowType1Label = TRUE,
                                                  PlotTitle = "", NetLayout = "lgl",
                                                  ShowNodes = TRUE,
                                                  VertexLabelCex = 1, VertexSize = 1, pheno, mod_size){
  
  p <- ncol(Abar)
  if(is.null(FeatureLabel)){
    FeatureLabel <- c(paste0("TypeI_", seq_len(P1)))
  }
  colnames(Abar) <- rownames(Abar) <- FeatureLabel[seq_len(p)]
  colnames(CorrMatrix) <- rownames(CorrMatrix) <- FeatureLabel[seq_len(p)]
  
  grp <- multiOmicsModule[[ModuleIdx]]
  grp.memb <- colnames(Abar)[grp]
  M.node <- grp.memb
  M <- as.matrix(Abar[M.node, M.node])
  
  # Trim module by PPR
  net_ppr <- igraph::graph_from_adjacency_matrix(M, weighted = TRUE,
                                                 diag = FALSE, mode = "undirected")
  # All parameter setups are based on the recommendation
  ranking <- igraph::page_rank(net_ppr, directed = FALSE, damping = 0.9, 
                               options = list(niter = 10^5,eps = 1e-06))
  # Obtain ranked protein names
  rank_names <- names(sort(ranking$vector, decreasing = TRUE))
  rank_value <- sort(ranking$vector, decreasing = TRUE)
  # Choose to include only the top proteins as we desired
  newM.node <- which(colnames(M) %in% rank_names[1:mod_size])
  
  if(length(newM.node) == 0){
    print("No edge passes threshold.")
  }else{
    M <- M[newM.node, newM.node]
    allidx <- matrix(seq_len(p), ncol = 1)
    rownames(allidx) <- rownames(Abar)
    NodeInfo <- data.frame(id = newM.node, idx = allidx[newM.node, ])
    net <- igraph::graph_from_adjacency_matrix(M, weighted = TRUE,
                                               diag = FALSE, mode = "undirected")
    
    # Define colors and shapes for vertex and label.
    k <- length(newM.node)
    type1 <- which(NodeInfo$idx <= P1)
    vcol <- rep("purple", k)
    vshape <- rep("square", k)
    lcol <- vcol
    # If not show nodes.
    if(!ShowNodes){vshape <- "none"}
    # If no label, assign a place holder.
    if(!ShowType1Label){newM.node[type1] <- " "}
    
    # Define edge colors.
    ecol <- rep("gray80", igraph::ecount(net))
    ew <- abs(igraph::edge.attributes(net)$weight) * 5
    ecol[which(igraph::edge.attributes(net)$weight < 0)] <- "red"
    
    # Define network layout.
    if(NetLayout == "circle"){
      l <- igraph::layout_in_circle(net)
    }else if(NetLayout == "sphere"){
      l <- igraph::layout_on_sphere(net)
    }else if(NetLayout == "fr"){
      l <- igraph::layout_with_fr(net)
    }else if(NetLayout == "lgl"){
      l <- igraph::layout_with_lgl(net)
    }else{
      stop("Unrecognized NetLayout input. Acceptable options are 'circle',
                 'sphere', 'fr', and 'lgl'.")
    }
    
    
    if(!is.null(SaveFile)){
      grDevices::pdf(SaveFile)
      graphics::par(bg = "white")
      graphics::plot(net, vertex.color = vcol, vertex.shape = vshape,
                     vertex.label.cex = VertexLabelCex, layout = l,
                     vertex.size = VertexSize, vertex.label = newM.node,
                     edge.width = ew, vertex.label.color = lcol, edge.color = ecol,
                     vertex.label.font = 2)
      graphics::title(PlotTitle)
      grDevices::dev.off()
    }else{
      graphics::plot(net, vertex.color = vcol, vertex.shape = vshape,
                     vertex.label.cex = VertexLabelCex, layout = l,
                     vertex.size = VertexSize, vertex.label = newM.node,
                     edge.width = ew, vertex.label.color = lcol, edge.color = ecol,
                     vertex.label.font = 2, main = PlotTitle)
    }
    
  }
  
  
  ###### section for principal component analysis
  X1_PC <- data[,which(colnames(data) %in% rank_names[1:mod_size])]
  # calculate individual protein correlation
  protein_correlation <- cor(X1_PC, Pheno)
  protein_correlation_data <- data.frame(name = colnames(X1_PC), correlation = protein_correlation)
  # PCA function
  pca_x1_summary <- prcomp(X1_PC)
  # Obtain summary
  summary_result <- summary(pca_x1_summary)
  # Extract the first three PC scores
  pca_x1_pc1 <- data.frame(pc1 = pca_x1_summary$x[,1], 
                           pc2 = pca_x1_summary$x[,2],
                           pc3 = pca_x1_summary$x[,3],
                           y = Pheno)
  correlation <- cor(pca_x1_pc1[,1:3], pca_x1_pc1[,4])
  
  pca_hybrid <- hybrid_score(data[,which(colnames(data) %in% rank_names[1:mod_size])], M)
  pca_hybrid_zero <- hybrid_score(data[,which(colnames(data) %in% rank_names[1:mod_size])], M, is_alpha = FALSE)
  correlation_hybrid <- cor(pca_hybrid, pca_x1_pc1[,4])
  correlation_hybrid_zero <- cor(pca_hybrid_zero, pca_x1_pc1[,4])
  # Save all the cc results into a same data file
  save(M = M, pca_score = pca_x1_pc1, rank_value = rank_value, pca_hybrid = pca_hybrid,
       pca_hybrid_zero = pca_hybrid_zero, correlation, protein_correlation_data,
       file = paste0(folder,'/',pheno,"net", mod_size,"_",ModuleIdx,"AA.Rdata"))
  
  return(c(nrow(M),correlation[1], correlation_hybrid, correlation_hybrid_zero))
}

#' Local Sub-network Pruning Algorithm
#' 
#' Local sub-network pruning to optimize the network size, this method contains
#' two steps: the first step is to find all the network size with correlation
#' between summarization score and phenotype that is within a certain range of
#' the maximum correlation across all the size; the second step is to first set
#' up a small baseline sub-network that contains only the top molecular
#' features based on the PageRank score, then setup a correlation cutoff such
#' that any network size with correlation between its summarization score and
#' the baseline summarization score below that will be excluded. The goal is
#' then to find the minimum/maximum of these two steps as the optimal network
#' depending on whether the user prefer a larger network or a smaller network.
#' 
#' 
#' @param corr_vec A correlation vector represents the correlation of
#' summarization score between a given size and the baseline size.
#' @param corr_pheno A correlation vector represents the summarization score
#' correlation with respect to phenotype for all the possible network size.
#' @param default_size The baseline network size.
#' @param cor_cut The correlation cut-off value for the second step.
#' @param network_preference With only the selection of 'small' or 'large',
#' whether to select a larger network or a smaller network.
#' @param tol the correlation cut-off value for the first step.

selection <- function (corr_vec, corr_pheno, default_size = 10, cor_cut = 0.8, 
                       network_preference = "small", tol = 0.9) {
  candidate_size_1 <- which(corr_pheno > (tol * max(corr_pheno)))
  candidate_size_2 <- which(corr_vec >= cor_cut)
  if (network_preference == "small") 
    candidate <- min(intersect(candidate_size_1, candidate_size_2))
  if (network_preference == "large") 
    candidate <- max(intersect(candidate_size_1, candidate_size_2))
  return(paste0("The optimal network size is: ", candidate + 
                  default_size - 1))
}

#' store subnetwork information with specified network size into the local directory
#' 
#' Based on global adjacency matrix, select the top omics features based on PageRank score,
#' and store the subnetwork into a specified directory, and plot the omics network using igraph
#' 
#' @param Abar global adjacency matrix
#' @param CorrMatrix corresponding correlation matrix of the global adjacency matrix
#' @param multiOmicsModule omics module after the hierarchical clustering
#' @param data omics data
#' @param Pheno phenotype data
#' @param type type of omcis data, e.g., protein
#' @param folder the folder directory to store the network
#' @param ModuleIdx network module index
#' @param P1 number of features in the global adjacency matrix4
#' @param FeatureLabel feature names of the global adjacency matrix
#' @param AddCorrSign add correlation direction sign to adjacency matrix
#' @param SaveFile save network file to local directory
#' @param pheno phenotype name
#' @param mod_size specified subnetwork size
    

plotMultiOmicsNetworks_PPR_Validation <- function(Abar, CorrMatrix, multiOmicsModule, data = X1,
                                                  Pheno, type, id, folder,
                                                  ModuleIdx, P1, penalty = 0, FeatureLabel = NULL,
                                                  AddCorrSign = TRUE, SaveFile = NULL,
                                                  ShowType1Label = TRUE,
                                                  PlotTitle = "", NetLayout = "lgl",
                                                  ShowNodes = TRUE,
                                                  VertexLabelCex = 1, VertexSize = 1, pheno, mod_size){
  
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
  
  hybrid_score = function(X, A, is_alpha = TRUE, pc_id=1){
    #X: data matrix (n, p)
    #A: corresponding adjacency matrix
    #pc_id: PC index
    g = graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
    # laplacian
    L2 <- A %>% 
      graph_from_adjacency_matrix(mode = "undirected", weighted = TRUE,diag = FALSE) %>% 
      graph.laplacian(normalized = F) %>% as.matrix
    # TOM
    TOM_A = TOMsimilarity(A, verbose = FALSE)
    
    alpha = graph.density(g)
    X= scale(X, center = TRUE, scale = TRUE)
    # weighted approach
    if (is_alpha == TRUE)
      temp  = (1-alpha)*(X %*% L2) + alpha*(X %*% TOM_A)
    else
      temp = (X %*% L2)
    temp = summary(prcomp(temp))
    h_score = temp$x[,1:pc_id] 
    importance = temp$importance[,1:pc_id]
    loading = temp$rotation[,1:pc_id]
    return(list(h_score, importance, loading))
    
  }
  
  
  p <- ncol(Abar)
  if(is.null(FeatureLabel)){
    FeatureLabel <- c(paste0("TypeI_", seq_len(P1)))
  }
  colnames(Abar) <- rownames(Abar) <- FeatureLabel[seq_len(p)]
  colnames(CorrMatrix) <- rownames(CorrMatrix) <- FeatureLabel[seq_len(p)]
  
  grp <- multiOmicsModule[[ModuleIdx]]
  grp.memb <- colnames(Abar)[grp]
  M.node <- grp.memb
  M <- as.matrix(Abar[M.node, M.node])
  
  # Trim module by PPR
  net_ppr <- igraph::graph_from_adjacency_matrix(M, weighted = TRUE,
                                                 diag = FALSE, mode = "undirected")
  # All parameter setups are based on the recommendation
  ranking <- igraph::page_rank(net_ppr, directed = FALSE, damping = 0.9, 
                               options = list(niter = 10^5,eps = 1e-06))
  # Obtain ranked protein names
  rank_names <- names(sort(ranking$vector, decreasing = TRUE))
  rank_value <- sort(ranking$vector, decreasing = TRUE)
  # Choose to include only the top proteins as we desired
  newM.node <- which(colnames(M) %in% rank_names[1:mod_size])
  
  if(length(newM.node) == 0){
    print("No edge passes threshold.")
  }else{
    M <- M[newM.node, newM.node]
    allidx <- matrix(seq_len(p), ncol = 1)
    rownames(allidx) <- rownames(Abar)
    NodeInfo <- data.frame(id = newM.node, idx = allidx[newM.node, ])
    net <- igraph::graph_from_adjacency_matrix(M, weighted = TRUE,
                                               diag = FALSE, mode = "undirected")
    
    # Define colors and shapes for vertex and label.
    k <- length(newM.node)
    type1 <- which(NodeInfo$idx <= P1)
    vcol <- rep("purple", k)
    vshape <- rep("square", k)
    lcol <- vcol
    # If not show nodes.
    if(!ShowNodes){vshape <- "none"}
    # If no label, assign a place holder.
    if(!ShowType1Label){newM.node[type1] <- " "}
    
    # Define edge colors.
    ecol <- rep("gray80", igraph::ecount(net))
    ew <- abs(igraph::edge.attributes(net)$weight) * 5
    ecol[which(igraph::edge.attributes(net)$weight < 0)] <- "red"
    
    # Define network layout.
    if(NetLayout == "circle"){
      l <- igraph::layout_in_circle(net)
    }else if(NetLayout == "sphere"){
      l <- igraph::layout_on_sphere(net)
    }else if(NetLayout == "fr"){
      l <- igraph::layout_with_fr(net)
    }else if(NetLayout == "lgl"){
      l <- igraph::layout_with_lgl(net)
    }else{
      stop("Unrecognized NetLayout input. Acceptable options are 'circle',
                 'sphere', 'fr', and 'lgl'.")
    }
    
    
    if(!is.null(SaveFile)){
      grDevices::pdf(SaveFile)
      graphics::par(bg = "white")
      graphics::plot(net, vertex.color = vcol, vertex.shape = vshape,
                     vertex.label.cex = VertexLabelCex, layout = l,
                     vertex.size = VertexSize, vertex.label = newM.node,
                     edge.width = ew, vertex.label.color = lcol, edge.color = ecol,
                     vertex.label.font = 2)
      graphics::title(PlotTitle)
      grDevices::dev.off()
    }else{
      graphics::plot(net, vertex.color = vcol, vertex.shape = vshape,
                     vertex.label.cex = VertexLabelCex, layout = l,
                     vertex.size = VertexSize, vertex.label = newM.node,
                     edge.width = ew, vertex.label.color = lcol, edge.color = ecol,
                     vertex.label.font = 2, main = PlotTitle)
    }
    
  }
  
  
  ###### section for principal component analysis
  X1_PC <- data[,which(colnames(data) %in% rank_names[1:mod_size])]
  # calculate individual protein correlation
  protein_correlation <- cor(X1_PC, Pheno)
  protein_correlation_data <- data.frame(name = colnames(X1_PC), correlation = protein_correlation)
  # PCA function
  pca_x1_summary <- prcomp(X1_PC)
  # Obtain summary
  summary_result <- summary(pca_x1_summary)
  
  # Obtain loadings
  pc_loading <- summary_result[["rotation"]]
  
  # Extract the first three PC scores
  pca_x1_pc1 <- data.frame(pc1 = pca_x1_summary$x[,1], 
                           pc2 = pca_x1_summary$x[,2],
                           pc3 = pca_x1_summary$x[,3],
                           y = Pheno)
  correlation <- cor(pca_x1_pc1[,1:3], pca_x1_pc1[,4])
  
  pca_hybrid <- hybrid_score(data[,which(colnames(data) %in% rank_names[1:mod_size])], M, pc_id = 3)
  pca_hybrid_zero <- hybrid_score(data[,which(colnames(data) %in% rank_names[1:mod_size])], M, pc_id = 3, is_alpha = FALSE)
  correlation_hybrid <- cor(pca_hybrid[[1]], pca_x1_pc1[,4])
  correlation_hybrid_zero <- cor(pca_hybrid_zero[[1]], pca_x1_pc1[,4])
  # Save all the cc results into a same data file
  save(M = M, pca_score = pca_x1_pc1, rank_value = rank_value, pca_hybrid = pca_hybrid, pca_hybrid_zero = pca_hybrid_zero, correlation, protein_correlation_data, pc_loading,
       file = paste0(folder,'/',pheno,"net", mod_size,"_",ModuleIdx,"AA.Rdata"))
  
  return(c(nrow(M),correlation[1], correlation_hybrid, correlation_hybrid_zero))
}

# source code for single-omics with quantitative phenotype
getRobustPseudoWeight_binary <- function(X1, Trait, Lambda1, 
                                         s1 = 0.7, SubsamplingNum = 1000
){
  
  
  # define the number of columns
  p1 <- ncol(X1); p <- p1 
  p1.sub <- ceiling(s1 * p1)
  X <- X1
  
  beta <- pbapply::pbsapply(seq_len(SubsamplingNum), function(x){
    # subsample features
    samp1 <- sort(sample(seq_len(p1), p1.sub, replace = FALSE))
    # center and scale the subsampled data 
    x1.par <- scale(X1[ , samp1], center = TRUE, scale = TRUE)
    # run splsda
    out <- spls::splsda(x = x1.par, y = Trait, K = 3, eta = Lambda1, kappa=0.5,
                        classifier=c('lda'), scale.x=FALSE)
    # store feature importance weight
    u <- rep(0, p1.sub)
    w <- rep(0, p)
    u[out[["spls.fit"]][["A"]]] <- out[["W"]]
    w[samp1] <- u
    coeff.avg <- w
    
    return(coeff.avg)
  })
  
  return(beta)
}
