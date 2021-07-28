# Getter/setter functions for filterRegulons.
setGeneric("filterRegulons", function(x, ...) standardGeneric("filterRegulons"))

#' Filtering Transcription factor regulons using Gene regulon networks
#'
#' @param x A \code{SingleCellExperiment} or \code{Seurat} object.
#' @param gene_list The regulon list used for filtering. default: NULL(use GRNs in your object).
#' @param grn_mat The gene regulatory network used for filtering.
#' @param use_regulon Regulons used for filtering.
#' @param use_grn GRNs used for regulon filtering.
#' @param score_type This parameter defines the GSEA score type. Possible options are ("std", "pos", "neg")
#' @param normalize Whether to normalize the ranking score?
#' @param withWeight Whether to keep weight in output regulons?
#' @param minSize The minimum number of genes in one TFR.
#' @param maxSize The maximum number of genes in one TFR.
#' @param filter_pval Set the cutoff of significant enriched TFRs. Default:1(no filtering).
#' @param ncores number of cores used.
#' @return An atfr object.
#' @export
#'
#' @examples
#'
#' @rdname filterRegulons
setMethod("filterRegulons", "SingleCellExperiment", function(x,
                                                             gene_list = NULL,
                                                             grn_mat = NULL,
                                                             use_regulon = 1L,
                                                             use_grn = 1L,
                                                             score_type = c("std", "pos", "neg"),
                                                             normalize=FALSE,
                                                             minSize = 5L,
                                                             maxSize = 1000L,
                                                             filter_pval=1L,
                                                             ncores=1){
  score_type <- match.arg(score_type)
  if(is.null(gene_list)){
    if(is.null(Regulon(x, use_regulon))) stop("The GRNs named ",use_grn," is not found in your object!")
    regulon <- Regulon(x, use_regulon)
  }else{
    regulon <- as.list(gene_list)
  }
  if(is.null(grn_mat)){
    if(is.null(GRN(x, use_grn))) stop("The GRNs named ",use_grn," is not found in your object!")
    grn_tab <- GRN(x,use_grn)
  }else{
    grn_tab <- as.matrix(grn_mat)
  }
  filter_regulon <- filterRegulonByGSEA(regulons = regulon,
                                        grn_tab = grn_tab,
                                        score_type = score_type,
                                        normalize = normalize,
                                        minSize = minSize,
                                        maxSize = maxSize,
                                        filter_pval = filter_pval,
                                        ncores = ncores)
  Regulon(x, i = 1L) <- filter_regulon
  x
})

#' Filtering Transcription factor regulons using GSEA
#'
#' @param regulons Gene list of regulons, the names of list must be transcription
#' factor, the values of each list should be gene symbols.
#' @param grn_tab A matrix contains weighted gene regulatory network.
#' @param score_type This parameter defines the GSEA score type. Possible options are ("std", "pos", "neg")
#' @param normalize Whether to normalize the ranking score?
#' @param withWeight Whether to keep weight in output regulons?
#' @param minSize The minimum number of genes in one TFR.
#' @param maxSize The maximum number of genes in one TFR.
#' @param filter_pval Set the cutoff of significant enriched TFRs. Default:1(no filtering).
#' @param ncores number of cores used.
#'
#' @return A gene list with regulon information.
#' @export
#'
#' @examples
filterRegulonByGSEA <- function(regulons,
                                grn_tab,
                                score_type = c("std", "pos", "neg"),
                                normalize=FALSE,
                                withWeight=TRUE,
                                minSize = 5L,
                                maxSize = 1000L,
                                filter_pval=1L,
                                ncores=1){
  score_type <- match.arg(score_type)
  TF_genes <- intersect(names(regulons),row.names(grn_tab))
  if(length(TF_genes)==0){
    stop("There is no overlap TFs between regulons and grn_tab, please check the gene names!")
  }else if(length(TF_genes)<50){
    warning("Only ",length(TF_genes)," TFs were both found in regulons and grn_tab.")
  }
  regulons <- as.list(regulons[TF_genes])
  regulons_type <- sapply(regulons, is.numeric)
  if(all(regulons_type)){
    withWeight <- TRUE
  }else{
    if(withWeight) warning("The input regulons is not weighted, set the 'withWeight=FALSE' ")
    withWeight <- FALSE
  }
  #---filter genes without expression.
  regulon_sbl <- lapply(regulons,function(x,grn_tab){
    if(is.numeric(x)){
      if(is.null(names(x))) stop("The weighted regulons must have target names for each weight!")
      targets <- intersect(names(x),colnames(grn_tab))
    }else if(is.character(x)){
      targets <- intersect(x,colnames(grn_tab))
    }
    targets
  },grn_tab=grn_tab)
  #---filter TF genes
  grn_tab <- as.matrix(grn_tab[TF_genes,,drop=FALSE])
  if(normalize) grn_tab <- t(apply(grn_tab,1,function(x){(x-min(x))/(max(x)-min(x))}))
  row.names(grn_tab) <- TF_genes
  score_type <- "pos"
  #---par
  cl <- parallel::makeCluster(ncores)
  parallel::clusterEvalQ(cl, library(fgsea))
  filter_res <- pbapply::pblapply(TF_genes, function(gname,
                                                     gene_list,
                                                     stat_mat,
                                                     score_type= c("std", "pos", "neg"),
                                                     minSize = 5L,
                                                     maxSize = 1000L){
    require(fgsea)
    score_type <- match.arg(score_type)
    suppressMessages(
      suppressWarnings(fgsea::fgseaMultilevel(pathways = gene_list[gname],
                                              stats = stat_mat[gname,],
                                              scoreType=score_type,
                                              minSize=minSize,
                                              maxSize=maxSize)))
  },gene_list=regulon_sbl, stat_mat=grn_tab, score_type=score_type, minSize=minSize, maxSize=maxSize, cl=cl)
  parallel::stopCluster(cl)
  #get result
  filter_tab <- data.table::rbindlist(filter_res, use.names = TRUE)
  #filter_tab <- filter_tab[filter_tab$NES>0,]
  regulon_targets <- filter_tab$leadingEdge
  names(regulon_targets) <- filter_tab$pathway
  if(withWeight){
    overlap_tfs <- intersect(names(regulons),names(regulon_targets))
    regulon_targets <- purrr::map2(.x = regulons[overlap_tfs],.y = regulon_targets[overlap_tfs],
                                   .f = function(x, y){x[intersect(names(x),y)]})
  }
  regulon_targets
}


