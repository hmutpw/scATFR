#' @export
setGeneric("runReduceDim", function(x, ...) standardGeneric("runReduceDim"))

#' Performing dimension reduction analysis
#'
#' Perform principal components analysis (PCA), t-stochastic neighbour embedding (t-SNE) and for the cells,
#' uniform manifold approximation and projection (UMAP) based on the data in a \code{matrix} \code{SingleCellExperiment}
#' or \code{Seurat} object.
#'
#' @param x A numeric matrix, \code\link{SingleCellExperiment}} object or \code\link{Seurat}} object contains
#' expression or regulon activity matrix.
#' @param alt_exp String or integer scalar specifying an alternative experiment containing the input data.
#' Can be obtained by \code{\link{altExpNames(x)}}.
#' @param exprsVal Integer scalar or string indicating which assay of x contains the expression values.
#' Can be obtained by \code{\link{assayNames(x)}}. For regulon activate matrix, use \code{\link{assayNames(altExp(x))}}
#' @param method String or integer scalar specifying the dimensionality reduction method to use.
#' @param nPCs Numeric scalar indicating the number of principal components to obtain. Same as \code{ncomponents} in
#' \pkg{scater}.
#' @param scale Logical scalar, should the expression values be standardized?
#' @param ... Other parameters passed to \code{\link{calculatePCA()}}, \code{\link{calculateTSNE()}} and
#' \code{\link{calculateUMAP()}}.
#'
#'
#' @return A \code{SingleCellExperiment} object or \code{matrix} with reduced dimension information.
#'
#'
#' @importFrom scater calculatePCA calculateTSNE calculateUMAP

#' @rdname runReduceDim
#' @export
setMethod("runReduceDim", "SingleCellExperiment",function(x,
                                                          alt_exp=NULL,
                                                          exprsVal=1L,
                                                          method = c("PCA","TSNE","UMAP"),
                                                          nPCs = 50,
                                                          scale=FALSE,
                                                          ...){
  method <- match.arg(method)
  if(is.null(alt_exp)){
    atfr <- x
  }else{
    atfr <- altExp(x,alt_exp)
  }
  if(method=="PCA"){
    reduced_dims <- scater::calculatePCA(atfr, exprs_values=exprsVal, ncomponents=nPCs, scale=scale, ...)
  }else if(method=="TSNE"){
    reduced_dims <- scater::calculateTSNE(atfr, exprs_values=exprsVal, ncomponents=nPCs, scale=scale, ...)
  }else if(method=="UMAP"){
    reduced_dims <- scater::calculateUMAP(atfr, exprs_values=exprsVal, ncomponents=nPCs, scale=scale, ...)
  }
  if(is.null(alt_exp)){
    reducedDim(x,method) <- reduced_dims
  }else{
    atfr <- altExp(x,alt_exp)
    reducedDim(altExp(x,alt_exp),method) <- reduced_dims
  }
  x
})

#' @rdname runReduceDim
#' @export
setMethod("runReduceDim", "Seurat",function(x,
                                            method = c("PCA","TSNE","UMAP"),
                                            nPCs = 50,
                                            scale=FALSE,
                                            ...){
  method <- match.arg(method)
  if(method=="PCA"){
    reduced_dims <- scater::calculatePCA(x, ncomponents=nPCs, scale=scale, ...)
  }else if(method=="TSNE"){
    reduced_dims <- scater::calculateTSNE(x, ncomponents=nPCs, scale=scale, ...)
  }else if(method=="UMAP"){
    reduced_dims <- scater::calculateUMAP(x, ncomponents=nPCs, scale=scale, ...)
  }
  reduced_dims
})


#' @rdname runReduceDim
#' @export
setMethod("runReduceDim", "matrix", function(x,
                                             method = c("PCA","TSNE","UMAP"),
                                             nPCs = 50,
                                             scale=FALSE,
                                             ...){
  method <- match.arg(method)
  x <- as.matrix(x)
  if(method=="PCA"){
    reduced_dims <- scater::calculatePCA(x, ncomponents=nPCs, scale=scale, ...)
  }else if(method=="TSNE"){
    reduced_dims <- scater::calculateTSNE(x, ncomponents=nPCs, scale=scale, ...)
  }else if(method=="UMAP"){
    reduced_dims <- scater::calculateUMAP(x, ncomponents=nPCs, scale=scale, ...)
  }
  reduced_dims
})









