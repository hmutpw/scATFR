setGeneric("inferGRNs", function(x, ...) standardGeneric("inferGRNs"))

#' Inferring Gene Regulatory Network (GRN) from Single-cell Expression Data
#'
#' @param x An object used for inferring GRNs. We support object with the
#'   following types: \code{\link{matrix}}, \code{\link{SingleCellExperiment}},
#'   \code{\link{Seurat}}, \code{\link{monocle3}}.
#' @param use_assay Assay used for GRN inferring, can be number or assay names.
#'   default: the first assay.
#' @param method Methods used for weighted matrix calculation. Current version only
#' support four types of methods. The "genie3" for \code{\link{GENIE3}} method from package \code{\link{GENIE3}}.
#' The "sclink" for package \code{\link{sclink}}. The "pcor" for \code{\link{ppcor}} package. The "spearman",
#' were common used method in \code{\link{cor()}} function.
#' @param use_regulator Regulators used for weighted matrix calculation, usually TFs. Default: NULL, stands for all
#' genes were used.
#' @param time_col A numeric \code{vector} or a one column \code{data.frame} contains time points for each sample.
#' If the input object x is a \code{\link{SingleCellExperiment}} or \code{\link{Seurat}} object with `colData`
#' information, the `time_col` could be the name of column with time information in `colData(x)`. The time points
#' must be numeric (eg. 0, 12, 24, 48) and t least 3 time points needed. Only useful for `sincerities` method.
#' @param ncores Number of cores used for parallel calculation. Using \code{\link{parallel::detectCores()}}
#' to detect the cores on your computer.Default: 1.
#' @param ... Other parameters passed to \code{\link{inferWeightExpMat}}
#'
#' @rdname inferGRNs
#' @export
setMethod("inferGRNs", "SingleCellExperiment", function(x,
                                                        regulon_list=NULL,
                                                        use_assay = 1L,
                                                        method = c("pcor", "genie3", "sclink", "sincerities"),
                                                        time_col = NULL,
                                                        ncores=1,
                                                        ...){
  exp_mat <- as.matrix(assay(x = x, i = use_assay))
  #---add time columns
  if(!is.null(time_col)){
    if(is.character(time_col) && length(time_col)==1){
      if(time_col %in% colnames(colData(x))){
        col_time <- data.frame(sample=row.names(colData(x)),time=colData(x)[[time_col]])
        row.names(col_time) <- row.names(colData(x))
      }else{
        stop("The 'time_col' is not found in colData(x)!")
      }
    }else if(is.numeric(time_col) && length(time_col)==nrow(colData(x))){
      if(!is.null(names(time_col))){
        colData(x)$time <- time_col[row.names(colData(x))]
      }else{
        colData(x)$time <- time_col
      }
      col_time <- data.frame(sample=row.names(colData(x)),time=colData(x)[["time"]])
      row.names(col_time) <- row.names(colData(x))
    }else if(is.data.frame(time_col) && nrow(time_col)==nrow(colData(x))){
      if(length(intersect(row.names(colData(x)),row.names(time_col)))==nrow(colData(x))){
        colData(x)$time <- time_col[row.names(colData(x)),1]
      }else{
        colData(x)$time <- time_col[,1]
      }
      col_time <- data.frame(sample=row.names(colData(x)),time=colData(x)[["time"]])
      row.names(col_time) <- row.names(colData(x))
    }else{
      stop("Current version of inferGRNs() only support 3 types of data.")
    }
  }else{
    col_time <- NULL
  }
  if(is.null(regulon_list)){
    use_regulator <- intersect(names(Regulon(x, i=1)),row.names(exp_mat))
    if(length(RegulonNames(x))==0){
      stop("No transcription factor regulons detected in your object. please add
           regulons `Regulon(x,i)<-value` into your object first!")
    }
    if(length(use_regulator)==0){
      stop("The names of each regulon list must be transcription factor symbols!")
    }
  }else{
    use_regulator <-Reduce(intersect,list(names(Regulon(x, i=1)), use_regulator, row.names(exp_mat)))
  }
  grn_mat <- inferWeightExpMat(exp_mat = exp_mat,
                               method = method,
                               use_regulator = use_regulator,
                               col_time = col_time,
                               ncores=ncores, ...)
  GRN(x, i = method) <- S4Vectors::DataFrame(grn_mat,check.names = FALSE)
  x
})

#' @rdname inferGRNs
#' @export
setMethod("inferGRNs", "Seurat",function(x,
                                         use_assay = 1L,
                                         method = c("sctenifoldnet","genie3", "sclink", "sincerities", "pcor", "spearman"),
                                         use_regulator = NULL,
                                         time_col = NULL, ...){
  exp_mat <- as.matrix(Seurat::GetAssayData(x, assay = use_assay, slot = "data"))
  #---add time columns
  if(!is.null(time_col)){
    if(is.character(time_col) && length(time_col)==1){
      if(time_col %in% colnames(x@meta.data)){
        col_time <- data.frame(sample=row.names(x@meta.data),time=x@meta.data[[time_col]])
        row.names(col_time) <- row.names(x@meta.data)
      }else{
        stop("The 'time_col' is not found in x@meta.data!")
      }
    }else if(is.numeric(time_col) && length(time_col)==nrow(x@meta.data)){
      if(!is.null(names(time_col))){
        x@meta.data$time <- time_col[row.names(x@meta.data)]
      }else{
        x@meta.data$time <- time_col
      }
      col_time <- data.frame(sample=row.names(x@meta.data),time=x@meta.data[["time"]])
      row.names(col_time) <- row.names(x@meta.data)
    }else if(is.data.frame(time_col) && nrow(time_col)==nrow(x@meta.data)){
      if(length(intersect(row.names(x@meta.data),row.names(time_col)))==nrow(x@meta.data)){
        x@meta.data$time <- time_col[row.names(x@meta.data),1]
      }else{
        x@meta.data$time <- time_col[,1]
      }
      col_time <- data.frame(sample=row.names(x@meta.data),time=x@meta.data[["time"]])
      row.names(col_time) <- row.names(x@meta.data)
    }else{
      stop("Current version of inferGRNs() only support 3 types of data.")
    }
  }else{
    col_time=NULL
  }

  grn_mat <- inferWeightExpMat(exp_mat = exp_mat,
                               method = method,
                               use_regulator = use_regulator,
                               col_time = col_time, ...)
  grn_assay <- Seurat::CreateAssayObject(data = grn_mat)
  Seurat::Key(grn_assay) <- "GRNs_"
  x[["GRNs"]] <- grn_assay
  x
})
#' @export
setMethod("inferGRNs", "matrix", function(x,
                                          method = c("sctenifoldnet","genie3", "sclink", "sincerities", "pcor", "spearman"),
                                          use_regulator = NULL,
                                          time_col = NULL, ...){
  exp_mat <- as.matrix(x)
  #---add time columns
  if(!is.null(time_col)){
    if(is.character(time_col)){
      stop("The 'time_col' can not be character when x is matrix!")
    }else if(is.numeric(time_col) && length(time_col)==ncol(x)){
      if(!is.null(names(time_col))){
        col_time <- data.frame(sample=colnames(x),time=time_col[colnames(x)])
      }else{
        col_time <- data.frame(sample=colnames(x),time=time_col)
      }
      row.names(col_time) <- colnames(x)
    }else if(is.data.frame(time_col) && nrow(time_col)==ncol(x)){
      if(length(intersect(colnames(x),row.names(time_col)))==ncol(x)){
        col_time <- data.frame(sample=colnames(x),time=time_col[colnames(x),1])
      }else{
        col_time <- data.frame(sample=colnames(x),time=time_col[,1])
      }
    }else{
      stop("Current version of inferGRNs() only support 3 types of data.")
    }
  }else{
    col_time=NULL
  }
  grn_mat <- inferWeightExpMat(exp_mat = exp_mat,
                               method = method,
                               use_regulator = use_regulator,
                               col_time = col_time, ...)
  grn_mat
})

#' @rdname inferGRNs
#' @export
setMethod("inferGRNs", "ExpressionSet",function(x,
                                                method = c("sctenifoldnet","genie3", "sclink", "sincerities", "pcor", "spearman"),
                                                use_regulator = NULL,
                                                time_col = NULL, ...){
  exp_mat <- as.matrix(Biobase::exprs(x))
  if(!is.null(time_col)){
    if(is.character(time_col)){
      stop("The 'time_col' can not be character when x is matrix!")
    }else if(is.numeric(time_col) && length(time_col)==ncol(x)){
      if(!is.null(names(time_col))){
        col_time <- data.frame(sample=colnames(x),time=time_col[colnames(x)])
      }else{
        col_time <- data.frame(sample=colnames(x),time=time_col)
      }
      row.names(col_time) <- colnames(x)
    }else if(is.data.frame(time_col) && nrow(time_col)==ncol(x)){
      if(length(intersect(colnames(x),row.names(time_col)))==ncol(x)){
        col_time <- data.frame(sample=colnames(x),time=time_col[colnames(x),1])
      }else{
        col_time <- data.frame(sample=colnames(x),time=time_col[,1])
      }
    }else{
      stop("Current version of inferGRNs() only support 3 types of data.")
    }
  }else{
    col_time=NULL
  }
  grn_mat <- inferWeightExpMat(exp_mat = exp_mat,
                               method = method,
                               use_regulator = use_regulator,
                               col_time = col_time, ...)
  Biobase::ExpressionSet(assayData = grn_mat)
})
#' Inferring weighted expression matrix
#'
#' @param exp_mat Matrix. Numeric matrix of gene expression data with rows were gene names and columns were samples.
#' @param method Methods used for weighted matrix calculation. Current version only
#' support four types of methods. The "genie3" for \code{\link{GENIE3}} method from package \code{\link{GENIE3}}.
#' The "sclink" for package \code{\link{sclink}}. The "pcor" for \code{\link{ppcor}} package. The "spearman",
#' were common used method in \code{\link{cor()}} function.
#' @param use_regulator Regulators used for weighted matrix calculation, usually TFs. Default: NULL, stands for all
#' genes were used.
#' @param col_time A two columns data.frame. The first column must be sample names which should be identical
#' with sample names in \code{exp_mat}. The second column must be time points for each sample which should
#' be numeric (eg. 0, 12, 24, 48). If you do not have exactly time information for your cell types, you can
#' just set time points from zero (eg. 0, 1, 2, 3, 4, ...). Note: at least 3 time points needed.
#' @param filter_genes Filtering low expressed genes before analysis. Default: TRUE.
#' @param ncores Number of cores used for parallel calculation. Using \code{\link{parallel::detectCores()}}
#' to detect the cores on your computer.Default: 1.
#' @param verbose Whether to show messages?
#' @param ... Other parameters passed to \pkg{GENIE3} \pkg{sincerities} and \pkg{sclink}.
#'
#' @return A matrix with weighted regulator-target correlation result.
#' @importFrom ppcor pcor
#' @importFrom GENIE3 GENIE3
#' @importFrom scLink sclink_cor sclink_norm
#' @importFrom parallel makeCluster clusterExport parSapply stopCluster
#'
#' @export
inferWeightExpMat <- function(exp_mat,
                              method = c("genie3", "pidc", "pcor","sincerities"),
                              use_regulator = NULL,
                              col_time = NULL,
                              filter_genes = TRUE,
                              ncores = 1,
                              verbose = interactive(),
                              ...){
  method <- match.arg(method)
  if(verbose) message("Inferring GRNs using method: ", method)
  #------check matrix format
  exp_mat <- as.matrix(exp_mat)
  if(filter_genes){
    exp_mat <- exp_mat[rowSums(exp_mat)>0,]
  }
  if(nrow(exp_mat)==0) stop("No gene names found in your input expression matrix!")
  if(!is.numeric(exp_mat)) stop("The value in your expression matrix is not numeric!")
  if(any(is.na(exp_mat))) stop("NA found in your expression matrix! Please check your input data!")
  #------check regulator
  if(is.null(use_regulator)){
    use_regulator <- row.names(exp_mat)
  }else{
    use_regulator <- intersect(use_regulator, row.names(exp_mat))
    if(length(use_regulator)==0) stop("No regulator detected while the 'use_regulator' is not NULL!")
  }
  #------check time_points
  if(is.null(col_time) && method=="sincerities"){
    stop("The 'col_time' can not be NULL when using 'sincerities' method!")
  }
  ##########
  #------start calculate GRNs
  if(method == "sctenifoldnet"){
    cor_mat <- run_scTenifoldNet(exp_mat = exp_mat,use_regulator = use_regulator, ...)
  }else if(method == "genie3"){
    #---GENIE3
    cor_mat <- suppressWarnings(GENIE3::GENIE3(exprMatrix=exp_mat, regulators=use_regulator, nCores=ncores,...))
  }else if(method == "sincerities"){
    #------sincerities
    sincerities_res <- run_sincerities(exp_mat = exp_mat, col_data = col_time, ...)
    cor_mat <- sincerities_res$adj_matrix
    cor_mat <- cor_mat[use_regulator,]
  }else if(method == "sclink"){
    #---scLink
    if(Sys.info()["sysname"]=="Windows" && ncores!=1){
      warning("When using 'sclink', ncores > 1 is not supported on Windows!")
      ncores = 1
    }
    cor_mat <- suppressWarnings(scLink::sclink_cor(expr=t(exp_mat),ncores=ncores,...))
    cor_mat <- cor_mat[use_regulator,]
  }else if(method == "pcor"){
    #---PPCOR
    pcor_res <-suppressWarnings(ppcor::pcor(t(exp_mat),method = "spearman"))
    cor_mat <- pcor_res$estimate
    dimnames(cor_mat) <- list(row.names(exp_mat),row.names(exp_mat))
    cor_mat <- cor_mat[use_regulator,]
  }else if(method == "spearman"){
    #---spearman
    use_regulator <- intersect(use_regulator,row.names(exp_mat))
    cor_mat <- .norm_methods(matrix=exp_mat, use_regulator=use_regulator, method=method, ncores = ncores)
  }else{
    stop("Method should be included in one of the methods.")
  }
  cor_mat
}

######
#------calculating cor with normal method
######
.norm_methods <- function(matrix,
                          use_regulator = NULL,
                          method = c("spearman", "pearson", "kendall"),
                          ncores = 1){
  #method <- match.arg(method)
  if(is.null(use_regulator)){
    use_regulator <- row.names(matrix)
  }
  use_regulator <- intersect(use_regulator, row.names(matrix))
  #------pa
  cl <- parallel::makeCluster(ncores)
  cor_mat <- pbapply::pbsapply(use_regulator,function(gname, matrix){
    exp_i <- matrix[gname,]
    apply(matrix,1,function(x){cor(x,exp_i,method = "spearman")})
  }, matrix=matrix,cl=cl)
  parallel::stopCluster(cl)
  t(cor_mat)
}

######
#------sincerities
######

#' Calculating GRNs with SINCERITIES
#'
#' @param exp_mat A gene by sample expression matrix.
#' @param col_data A two columns data.frame. The first column must be sample names which should be identical
#' with sample names in \code{exp_mat}. The second column must be time points for each sample which should
#' be numeric (eg. 0, 12, 24, 48). If you do not have exactly time information for your cell types, you can
#' just set time points from zero (eg. 0, 1, 2, 3, 4, ...). Note: at least 3 time points needed.
#' @param ... Other parameters passed to \code{SINCERITIES}, \code{SINCERITIES_PLUS} and
#' \code{final_ranked_predictions}.
#'
#' @return A list containing the following information:
#' \itemize{
#' \item adj_matrix: \emph{m} by \emph{m} matrix containing the weights of
#' regulatory edges. The larger \code{adj_matrix[i,j]} indicates higher confidence
#' that the corresponding edge exists (i.e., gene \emph{i} regulating gene \emph{j}).
#' \item final_table: a \code{data.frame} containing the regulation relationships
#' between source genes and targets and its regulatory directions.}
#'
#' @author \code{SINCERITIES} were first created by Nan Papili Gao and R version implemented
#' by Ziyi Hua. Institute for Chemical and Bioengineering.ETH Zurich. E-mail: nanp@ethz.ch.\cr
#' Puwen Tan Integrated the main functions from \code{SINCERITIES} into \code{run_sincerities()}, which
#' was included in \pkg{scATFR} package.
#'
#' @export
run_sincerities <- function(exp_mat, col_data, ...){
  #------check exp_mat
  exp_mat <- as.matrix(exp_mat)
  if(!identical(sort(colnames(exp_mat)),sort(as.character(col_data[,1])))){
    stop("The samples in exp_mat is not identical in col_data!")
  }
  #------check col_data
  if(is.null(col_data)) stop("The 'col_data' can't be NILL!")
  if(!is.data.frame(col_data)) stop("The 'col_data' must be a data.frame!")
  if(!is.numeric(col_data[,2])) stop("The second column of 'col_data' must be numeric!")
  exp_mat <- t(exp_mat)
  gene_num <- ncol(exp_mat)
  gene_names <- colnames(exp_mat)
  col_data <- col_data[order(col_data[,2]),]
  time_line <- col_data[,2]
  time <- sort(unique(time_line))
  exp_mat <- exp_mat[as.character(col_data[,1]),]
  single_cell_data <- by(exp_mat,time_line,identity)
  DATA <- list(time=time,num_time_points=length(time),
               totDATA=exp_mat,timeline=time_line,numGENES=gene_num,
               genes=gene_names,singleCELLdata=single_cell_data)
  #------get gene regulatory network
  if(length(time)>=5){
    sincerities_res <- SINCERITIES(DATA = DATA, ...)
  }else if(length(time)<5 && length(time)>2){
    sincerities_res <- SINCERITIES_PLUS(DATA = DATA, ...)
  }else{
    stop("SINCERITIES method must have at least 3 time points.")
  }
  #---normalization
  adj_matrix <- sincerities_res$adj_matrix/max(sincerities_res$adj_matrix)
  final_table <- final_ranked_predictions(adj_matrix,DATA$genes,SIGN=1, ...)
  list(adj_matrix = adj_matrix, final_table = final_table)
}




#' Run scTenifoldNet to infer Gene Regulatory network
#'
#' @param exp_mat A gene by sample expression matrix.
#' @param use_regulator Regulators used for weighted matrix calculation, usually TFs. Default: NULL, stands for all
#' genes were used.
#' @param nCells An integer value. The number of cells to subsample each time to generate a network. passed to
#' \code{\link{makeNetworks}}.
#' @param nNet An integer value. The number of networks based on principal components regression to generate.
#' passed to \code{\link{makeNetworks}}.
#' @param nComp An integer value. The number of principal components in PCA to generate the networks.
#' Should be greater than 2 and lower than the total number of genes. passed to \code{\link{makeNetworks}}.
#' @param ... Other parameters passed to \code{\link{makeNetworks}} \code{\link{tensorDecomposition}}
#'
#' @return A weighted matrix for gene regulation relationship.
#' @export
#'
run_scTenifoldNet <- function(exp_mat, use_regulator=NULL, nCells=500, nNet=100, nComp=10, ...){
  X <- as.matrix(exp_mat)
  #X <- scTenifoldNet::scQC(X=exp_mat)
  X <- scTenifoldNet::cpmNormalization(X=X)
  nCells <- min(nCells,ncol(X)/4)
  xList <- scTenifoldNet::makeNetworks(X = X,
                                       nCells = nCells,
                                       nNet = nNet,
                                       nComp = nComp, ...)
  tensorOut <- scTenifoldNet::tensorDecomposition(xList, nDecimal=10,...)
  tX <- as.matrix(tensorOut$X)
  tX <- (tX + t(tX))/2
  if(!is.null(use_regulator)){
    use_regulator <- intersect(use_regulator,row.names(tX))
    if(length(use_regulator)==0) stop("There is no overlap regulators!")
    tX <- tX[use_regulator,]
  }
  tX
}
