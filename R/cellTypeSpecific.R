setGeneric("findCellTypeSignature", function(x, ...) standardGeneric("findCellTypeSignature"))
#' Identifying cell-type specific transcription factor regulons
#'
#' @param x An object with gene expression data used for regulon activity evaluation,
#' which can be \code{matrix}, \code{SingleCellExperiment}, \code{Seurat} and \code{ExpressionSet}.
#' @param cell_anno A one-column \code{data.frame}. Rows is cell name and columns is cell cluster or stages.
#' @param pos_clust Character, cluster used to calculate the cell-type specific score (CSS).
#' @param neg_clust Other clusters used as controls. Default: NULL (indicating all clusters expect for 'pos_clust').
#' @param nperm Number of permutations used for calculating p-values against CSS. Default:1000.
#' @param ncores Number of threads used for parallel calculation. Default: 1.
#'
#' @return A \code{data.frame} with differential CSS for each TFRs.
#' @import data.table
#' @export
#'
#'
#' @rdname findCellTypeSignature
#' @export
setMethod("findCellTypeSignature", "SingleCellExperiment", function(x,
                                                                    pos_clust,
                                                                    neg_clust = NULL,
                                                                    diff_by="cluster",
                                                                    use_altExp=TRUE,
                                                                    use_assay=1L,
                                                                    nperm = 1000,
                                                                    ncores = 1){
  if(use_altExp){
    mat <- as.matrix(assay(altExp(x=x,e=use_assay)))
    assay_name <- assayNames(altExp(x,e=use_assay))
  }else{
    mat <- as.matrix(assay(x=x,i=use_assay))
    assay_name <- assayNames(x,e=use_assay)
  }
  cell_anno <- as.data.frame(SummarizedExperiment::colData(x))
  if(!(diff_by %in% colnames(cell_anno))) stop("The 'diff_by' is not found in colData(x), ",
                                               "please check the names of colData() in your object!")
  cell_anno <- cell_anno[,diff_by,drop=FALSE]
  css_res <- calcCSS(mat = mat,
                     cell_anno = cell_anno,
                     pos_clust = pos_clust,
                     neg_clust = neg_clust,
                     nperm = nperm,
                     ncores = ncores)
  regulon_list <- Regulon(x)[css_res$regulon]
  css_res$targets <- sapply(regulon_list,function(x){paste0(x,collapse=";")})
  col_name_order <- c("regulon","targets",setdiff(colnames(css_res),c("regulon","targets")))
  css_res[,col_name_order]
})

#' @rdname findCellTypeSignature
#' @export
setMethod("findCellTypeSignature","matrix",function(x,
                                                    cell_anno,
                                                    pos_clust,
                                                    neg_clust = NULL,
                                                    nperm = 1000,
                                                    ncores = 1){
  mat <- as.matrix(x)
  css_res <- calcCSS(mat = mat,
                     cell_anno = cell_anno,
                     pos_clust = pos_clust,
                     neg_clust = neg_clust,
                     nperm = nperm,
                     ncores = ncores)
  css_res
})



#' @rdname findCellTypeSignature
#' @export
setMethod("findCellTypeSignature", "ExpressionSet", function(x,
                                                             cell_anno,
                                                             pos_clust,
                                                             neg_clust = NULL,
                                                             nperm = 1000,
                                                             ncores = 1){

  mat <- as.matrix(Biobase::exprs(x))
  css_res <- calcCSS(mat = mat, cell_anno = cell_anno,
                     pos_clust = pos_clust,
                     neg_clust = neg_clust,nperm = nperm,
                     ncores = ncores)
  css_res
})

#' @rdname findCellTypeSignature
#' @export
setMethod("findCellTypeSignature", "Seurat", function(x,
                                                 cell_anno,
                                                 pos_clust,
                                                 neg_clust = NULL,
                                                 nperm = 1000,
                                                 ncores = 1){
  if(!(names(x) %in% "atfr_mat"))
    stop("Activate Transcription Factor Regulon matrix not found in you object, please run regulonActivity() first!")

  mat <- as.matrix(x[["atfr_mat"]]@data)
  css_res <- calcCSS(mat = mat, cell_anno = cell_anno,
                     pos_clust = pos_clust,
                     neg_clust = neg_clust,nperm = nperm,
                     ncores = ncores)


  css_res
})

################################################
# find all cell type specific regulons
setGeneric("findAllCTSRegulons", function(x, ...) standardGeneric("findAllCTSRegulons"))

#' Calculating the significant enriched regulons for each cluster
#'
#' @param x An object with gene expression data used for regulon activity evaluation,
#' which can be \code{matrix}, \code{SingleCellExperiment}, \code{Seurat} and \code{ExpressionSet}.
#' @param cell_anno A two-column \code{data.frame} with cell name and cell cluster or cell (or stage) information.
#' @param fdr FDR cutoff used for filtering significant regulons.
#' @param css CSS cutoff used for filtering significant regulons.
#'
#' @return A \code{data.frame} with differential CSS for each TFRs.
#'
#' @import data.table
#' @export
#'
#' @rdname findAllCTSRegulons
#' @examples
setMethod("findAllCTSRegulons","SingleCellExperiment",function(x,
                                                               use_altExp=TRUE,
                                                               use_assay=1L,
                                                               diff_by="cluster",
                                                               fdr = 0.05,
                                                               css = 0.3,
                                                               nperm=1000,
                                                               ncores=1, ...){
  require(data.table)
  if(use_altExp){
    mat <- as.matrix(assay(altExp(x=x,e=use_assay)))
  }else{
    mat <- as.matrix(assay(x=x,i=use_assay))
  }
  regulon_list <- as.list(Regulon(x))
  col_data <- SummarizedExperiment::colData(x)
  if(!(diff_by %in% colnames(col_data))) stop("The 'diff_by' is not found in colData(x), ",
                                               "please check the names of colData() in your object!")
  cell_anno <- data.frame(sample=row.names(col_data),
                          cell_type=col_data[,diff_by,drop=TRUE],
                          row.names=row.names(col_data))
  pos_clusters <- unique(as.character(cell_anno$cell_type))
  #------par
  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl, c("calcCSS",".getCSS",".calcJSsp",".H"))
  parallel::clusterEvalQ(cl,library(data.table))
  all_res <- pbapply::pblapply(pos_clusters,function(clust, mat, cell_anno, nperm){
    calcCSS(mat=mat, cell_anno=cell_anno, pos_clust=clust,nperm=nperm)
  },mat=mat,cell_anno=cell_anno,nperm=nperm,cl=cl)
  parallel::stopCluster(cl)
  #------end par
  names(all_res) <- pos_clusters
  all_res <- data.table::rbindlist(all_res, use.names = TRUE, idcol = diff_by)
  all_res$targets <- regulon_list[all_res$regulon]
  all_res <- dplyr::tibble(all_res)
  diff_res <- all_res %>% filter(padj<fdr & CSS > css)
  diffRegulons(x) <- diff_res
  diff_res
})


#' @rdname findAllCTSRegulons
#' @examples
setMethod("findAllCTSRegulons","matrix",function(x,
                                                 cell_anno,
                                                 fdr = 0.05,
                                                 css = 0.5,
                                                 ...){
  mat <- as.matrix(x)
  pos_clusts <- unique(as.character(cell_anno[,2]))
  getSigCSS <- function(pos_clust, mat, cell_anno,
                        fdr=0.05,css=0.3,...){
    css_res <- calcCSS(mat = mat, cell_anno = cell_anno, pos_clust = pos_clust,...)
    css_res <- css_res[css_res$padj<fdr,]
    css_res[css_res$CSS>=css,]
  }
  all_res <- lapply(pos_clusts,getSigCSS,
                    mat=mat, cell_anno=cell_anno,
                    fdr=fdr,css=css,...)
  names(all_res) <- pos_clusts
  data.table::rbindlist(all_res,use.names = TRUE,idcol = "cluster")
})

################################################
#---calculate cell-type specific score
################################################
#' Identifying cell-type specific transcription factor regulons
#'
#' @param mat A \code{matrix} contains transcription factor regulons activate score.
#' @param cell_anno A two-column \code{data.frame} with cell name and cell cluster or cell (or stage) information.
#' @param pos_clust Character, cluster used to calculate the cell-type specific score (CSS).
#' @param neg_clust Other clusters used as controls. Default: NULL (indicating all clusters expect for 'pos_clust').
#' @param nperm Number of permutations used for calculating p-values against CSS. Default:1000.
#' @param ncores Number of threads used for parallel calculation. Default: 1.
#'
#' @return A \code{data.frame} with differential CSS for each TFRs.
#'
calcCSS <- function(mat, cell_anno, pos_clust, neg_clust = NULL, nperm=1000, ncores = 1){
  mat <- as.matrix(mat)
  pos_cells <- row.names(cell_anno[cell_anno$cell_type %in% pos_clust,])
  if(is.null(neg_clust)){
    neg_cells <- row.names(cell_anno)
  }else{
    neg_cells <- row.name(cell_anno[cell_anno$cell_type %in% c(pos_clust,neg_clust),])
  }
  mat <- mat[,neg_cells]
  if(nrow(mat)==0 || ncol(mat)==0) stop("No activate score found in matrix!")
  #---normalization
  mat <- t(apply(mat,1,function(x){(x-min(x))/(max(x)-min(x))}))
  #-----calculate average activate score
  mean_score <- t(apply(mat,1,function(x){
    pos_mean <- mean(x[pos_cells])
    neg_mean <- mean(x[setdiff(names(x),pos_cells)])
    c(pos_score = pos_mean, other_score = neg_mean)}))
  mean_score <- as.data.frame(mean_score)
  mean_score$ave_diff <- mean_score[,1]-mean_score[,2]
  #---get p values
  CSS <- as.data.frame(t(pbapply::pbapply(mat,1,.getCSS, pos_cells = pos_cells, nperm=nperm,ncores=ncores)))
  CSS$padj <- p.adjust(CSS$pval, method = "BH")
  CSS <- data.frame(regulon = row.names(mean_score),mean_score, CSS[row.names(mean_score),], check.names=FALSE)
  CSS[order(CSS$padj,CSS$CSS,decreasing = c(FALSE, TRUE)),]
}

#---calculating css score
.getCSS <- function(regulon, pos_cells, nperm=1000, ncores = 1){
  #normalization
  regulon_const <- rep(0,length(regulon))
  names(regulon_const) <- names(regulon)
  regulon_const[pos_cells] <- rep(1,length(pos_cells))
  # generate probility
  regulon <- regulon/sum(regulon)
  regulon_const <- regulon_const/sum(regulon_const)
  regulon_css <- .calcJSsp(regulon, regulon_const)
  #------permutation
  #cl <- parallel::makeCluster(ncores)
  #parallel::clusterExport(cl,c(".calcJSsp",".H","regulon","regulon_const"))
  nperm_regulon <- sapply(X=1:nperm,function(x, regulon, regulon_const){
    x <- sample(regulon,length(regulon))
    names(x) <- names(regulon_const)
    .calcJSsp(x, regulon_const)
  }, regulon, regulon_const)
  #parallel::stopCluster(cl)
  #------
  pval <- length(nperm_regulon[which(nperm_regulon>regulon_css)])/nperm
  c(CSS = regulon_css, pval = pval)
}

#---calculate cell-type specific score using JS divergence
.calcJSsp <- function(p1, p2, unit = c("log","log2","log10")){
  unit <- match.arg(unit)
  #---calculate JSD
  JSD <-.H(p=0.5*(p1+p2), unit = unit) - 0.5*(.H(p1, unit = unit)+.H(p2, unit = unit))
  1-sqrt(JSD)
}

#---calculate entropy from probability
.H <- function(p, unit = c("log","log2","log10")){
  unit <- match.arg(unit)
  if(unit == "log"){
    logp <- log(p)
  }else if(unit == "log2"){
    logp <- log2(p)
  }else{
    logp <- log10(p)
  }
  logp[is.infinite(logp)] <- 0
  -sum(p*(logp))
}
