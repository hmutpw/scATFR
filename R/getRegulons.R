# Getter/setter functions for getRegulons.
setGeneric("getRegulons", function(x, ...) standardGeneric("getRegulons"))

#' Getting Transcription Factors Regulons using motif cisTarget information
#'
#'
#' @param x A matrix or data.frame containing gene regulatory network.
#' @param motif_ranking A \code{\link{rankingRcisTarget}} object or a \code{\link{matrix}} or
#' \code{\link{data.frame}} containing motif ranking matrix, in which rows
#' corresponding to motif ids and columns corresponding to genes.
#' @param motif_tf_anno A \code{\link{data.table}} with motif binding TF information.
#' @param method Methods used for TFR identification.
#' @param ... Others parameter passed to \code{\link{grnToLists}}, \code{\link{grnRcisTarget}}
#' \code{\link{doMotifRankGSEA}}
#'
#' @return A list of TF regulons.
#' @import data.table
#' @importFrom RcisTarget getRanking cisTarget
#'
#'
#' @rdname getRegulons
#' @export
setMethod("getRegulons", "matrix", function(x,
                                            motif_ranking,
                                            motif_tf_anno,
                                            method = c("cisTarget", "GSEA"),
                                            use_regulator = NULL,
                                            ncores = 1,
                                            ...){
  method <- match.arg(method)
  message("Transforming GRN to regulons...")
  grnList <- grnToLists(grn_mat = x, use_regulator = use_regulator, ...)
  message("Totally ",length(grnList)," regulotors!\nFiltering regulons using motif information...")
  if(method=="cisTarget"){
    grnFinal <- grnRcisTarget(gene_list = grnList,
                              motif_rank = motif_ranking,
                              motif_tf_anno = motif_tf_anno,...)
  }else if(method=="GSEA"){
    motif_rank_tab <- RcisTarget::getRanking(motif_ranking)
    motif_rank <- as.data.frame(motif_rank_tab)
    row.names(motif_rank) <- motif_rank_tab$features
    grnFinal <- doMotifRankGSEA(motif_rank = motif_rank,
                                gene_list = grnList,
                                motif_tf_anno = motif_tf_anno,
                                ncores = ncores, ...)
  }else{
    stop("Current version of 'getRegulons' only support 'cisTarget' and 'GSEA' methods!")
  }
  grnFinal
})

#' @rdname getRegulons
#' @export
setMethod("getRegulons", "data.frame", function(x,
                                                motif_ranking,
                                                motif_tf_anno,
                                                method = c("cisTarget", "GSEA"),
                                                use_regulator = NULL,
                                                ncores = 1,
                                                ...){
  method <- match.arg(method)
  message("Transforming GRN to regulons...")
  grnList <- grnToLists(grn_mat = x, use_regulator = use_regulator, ...)
  message("Totally ",length(grnList)," regulotors!\nFiltering regulons using motif information...")
  if(method=="cisTarget"){
    grnFinal <- grnRcisTarget(gene_list = grnList,
                              motif_rank = motif_ranking,
                              motif_tf_anno = motif_tf_anno, ...)
  }else if(method=="GSEA"){
    motif_rank_tab <- RcisTarget::getRanking(motif_ranking)
    motif_rank <- as.data.frame(motif_rank_tab)
    row.names(motif_rank) <- motif_rank_tab$features
    grnFinal <- doMotifRankGSEA(motif_rank = motif_rank,
                                gene_list = grnList,
                                motif_tf_anno = motif_tf_anno,
                                ncores = ncores, ...)
  }else{
    stop("Current version of 'getRegulons' only support 'cisTarget' and 'GSEA' methods!")
  }
  grnFinal
})

######
#------get top ranked genes for each regulator
#' Transform Gene Regulatory Network to regulons
#'
#'
#' @param grn_mat A matrix included weighted Gene Regulatory Network (GRN). Rows
#' were regulators (eg. TFs) and columns were target genes.
#' @param use_regulator The regulator used for regulons.
#' @param cutoff Percentage of top ranked genes used for constructing regulons. Default: 0.05.
#' @param use_target Keep wanted target genes in regulons.
#' @param rank_type  Get positive enriched or negative enriched targets? Default: pos.
#'
#' @return A list of TF regulons.
#' @export
grnToLists <- function(grn_mat,
                       use_regulator = NULL,
                       cutoff = 0.05,
                       use_target = NULL,
                       rank_type = c("pos","neg")){
  #---check arguments
  rank_type <- match.arg(rank_type)
  if(is.matrix(grn_mat) || is.data.frame(grn_mat)){
    grn_mat <- as.matrix(grn_mat)
  }else{
    stop("The 'grn_mat' should be matrix or data.frame!")
  }
  if(is.null(use_regulator)){
    use_regulator <- row.names(grn_mat)
  }else{
    use_regulator <- intersect(use_regulator,row.names(grn_mat))
  }
  if(length(use_regulator)==0) stop("No regulator found in your 'grn_mat'!")
  if(is.null(use_target)){
    use_target <- colnames(grn_mat)
  }else{
    use_target <- intersect(use_target,colnames(grn_mat))
  }
  if(length(use_target)==0) stop("No target found in your 'grn_mat'!")
  ######
  top_gene_num <- round(cutoff*ncol(grn_mat))
  module_list <- list()
  for(i in use_regulator){
    order_rank <- sort(grn_mat[i,setdiff(colnames(grn_mat),i)])
    if(rank_type=="pos"){
      module <- intersect(rev(tail(names(order_rank),top_gene_num)),use_target)
    }else if(rank_type=="neg"){
      module <- intersect(head(names(order_rank),top_gene_num),use_target)
    }else{
      stop("Wrong 'rank_type', must be 'pos' or 'neg'!")
    }
    module_list[[i]] <- module
  }
  module_list
}
##############
#RcisTarget methods
##############
#' Filtering Gene Regulatory Networks using motif information
#'
#'
#' @param gene_list A list of TF regulons used for filtering.
#' @param motif_rank motif ranking object import from RcisTarget.
#' @param motif_tf_anno A data.table including motif and TF annotation.
#' @param min_target_num Minimum number of targets in your final regulons.
#' @param max_target_num Maximum number of targets in your final regulons.
#' @param ncores Number of cores used for processing.
#' @param ... Others parameters passed to \code{\link{cisTarget}}
#'
#' @return A list of TF regulons
#' @importFrom  RcisTarget cisTarget
#'
#'
#' @export
grnRcisTarget <- function(gene_list,
                          motif_rank,
                          motif_tf_anno,
                          min_target_num = 5,
                          max_target_num = 500,
                          ncores = 1, ...){
  #---perform cis-target analysis
  cisTarget_res <- RcisTarget::cisTarget(geneSets = gene_list,
                                         motifRankings = motif_rank,
                                         motifAnnot = motif_tf_anno,
                                         highlightTFs = names(gene_list),
                                         nesThreshold = 1,
                                         nCores = ncores, ...)
  data.table::setnames(cisTarget_res,"geneSet","TF")
  #---get enriched TFs for each motif
  cisTarget_res$TF_highConf <- lapply(cisTarget_res$TF_highConf,function(x){
    genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
    genesSplit <- unique(unlist(strsplit(genes, "; ")))})
  cisTarget_res$TF_lowConf <- lapply(cisTarget_res$TF_lowConf,function(x){
    genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
    genesSplit <- unique(unlist(strsplit(genes, "; ")))})
  cisTarget_filter <-apply(cisTarget_res,1,function(x){
    if(x['TF_highConf'] %in% x['TF'] || x['TF_lowConf'] %in% x['TF']){
      return(TRUE)
    }else{
      return(FALSE)
    }})
  cisTarget_filtered <- cisTarget_res[cisTarget_filter,]
  #---get max NES motifs
  cisTarget_max_NES <- cisTarget_filtered[,.SD[which.max(abs(NES))],by=TF]
  cisTarget_regulon <- lapply(cisTarget_max_NES$enrichedGenes,function(x){
    unlist(strsplit(x,split = ';'))})
  names(cisTarget_regulon) <- cisTarget_max_NES$TF
  #---filter regulon by target gene number
  regulon_gene_num <- sapply(cisTarget_regulon, function(x){
    length(x[!is.na(x)])>=min_target_num && length(x[!is.na(x)])<=max_target_num})
  cisTarget_regulon[regulon_gene_num]
}

##############
#GSEA methods
##############

#' Performing GSEA with motif ranking data
#'
#'
#' @param motif_rank A \code{\link{data.frame}} containing motif ranking matrix,
#' in which rows corresponding to motifs and columns corresponding to genes. Note, The
#' first rows were motif ids
#' @param gene_list A list of pre-defined modules. Names of each mdule should be TFs.
#' @param motif_tf_anno A \code{\link{data.table}} with motif binding TF
#' information download from \pkg{\link{RcisTarget}}.
#' @param anno_type The annotation type of motifs with TFs, should be some of:
#' "directAnnotation", "inferred_Orthology", "inferred_MotifSimil". Default:
#' all above 3 types.
#' @param p_value set the cutoff of significant enriched motifs. Default: 0.05.
#' @param revValue Should reverse the values before perform GSEA?
#' @param verbose No messages should print while running.
#' @param ncores Number of cores used for parallel calculation. Default: 1.
#' @param ... Others parameters passed to \code{\link{fgsea}}.
#'
#' @return A \code{\link{data.table}} with motif GSEA result.
#' @import data.table
#' @importFrom RcisTarget getRanking
#' @importClassesFrom RcisTarget rankingRcisTarget
#' @importFrom fgsea fgseaMultilevel
#'
#' @export
doMotifRankGSEA <- function(motif_rank,
                            gene_list,
                            motif_tf_anno,
                            anno_type = c("directAnnotation","inferred_Orthology","inferred_MotifSimil"),
                            p_value = 0.05,
                            revValue = TRUE,
                            ncores = 1,
                            verbose = interactive(), ...){
  #---check 'motif_tf_anno' anno_type
  anno_type <- match.arg(anno_type, several.ok = TRUE)
  anno_type_table <- motif_tf_anno[, ..anno_type]
  motif_tf_anno[, type := apply(anno_type_table,1,any)]
  motif_tf_anno <- motif_tf_anno[type==TRUE ,]
  #------check arguments
  tf_genes <- intersect(names(gene_list),as.character(motif_tf_anno$TF))
  motif_ids <- intersect(as.character(motif_rank$features), motif_tf_anno[TF %in% tf_genes, ]$motif)
  motif_rank <- motif_rank[motif_rank$features %in% motif_ids,]
  motif_rank_mat <- as.matrix(motif_rank[,-1])
  row.names(motif_rank_mat) <- motif_rank$features

  message("Performing GSEA analysis for ",nrow(motif_rank),
          " motifs with ", length(gene_list), " TF modules")
  #------parallel
  cl <- parallel::makeCluster(ncores)
  parallel::clusterExport(cl,".do_fgsea")
  parallel::clusterEvalQ(cl, library(fgsea))
  parallel::clusterEvalQ(cl, library(data.table))
  fgsea_res <- parallel::parLapply(cl, motif_ids, .do_fgsea,
                                   gene_list = gene_list,
                                   motif_rank = motif_rank_mat,
                                   motif_tf_anno = motif_tf_anno,
                                   revValue = revValue, ...)
  parallel::stopCluster(cl)
  #------merge fgsea result
  fgsea_out <- data.table::rbindlist(fgsea_res)
  data.table::setnames(fgsea_out,c("pathway"),c("TF"))
  fgsea_out <- fgsea_out[fgsea_out$pval < p_value,]
  data.table::setkey(fgsea_out, TF, motif)
  fgsea_out_filter <- fgsea_out[,.SD[which.max(abs(NES)),],by=TF]
  tf_regulons <- fgsea_out_filter$leadingEdge
  names(tf_regulons) <- fgsea_out_filter$TF
  tf_regulons
}

#======do fgsea analysis on one motif rank for gene lists
.do_fgsea <- function(motif_id, gene_list, motif_rank, motif_tf_anno,
                      revValue = FALSE, scoreType = c("std", "pos", "neg"),
                      minSize = 5L, maxSize = 1000L,
                      verbose = interactive(), ...){
  scoreType <- match.arg(scoreType)
  #------check input motif rank
  motif_rank_vec <- motif_rank[motif_id,]
  if(revValue){
    motif_rank_vec <- rank(-motif_rank_vec)
    scoreType <- "pos"
  }
  #------check input gene list
  tf_genes <- unique(as.character(motif_tf_anno[motif_tf_anno$motif %in% motif_id,]$TF))
  tf_genes <- intersect(tf_genes, names(gene_list))
  gene_list <- gene_list[tf_genes]
  #------gsea analysis
  fgsea_res <- fgsea::fgseaMultilevel(pathways = gene_list,
                                      stats = motif_rank_vec,
                                      minSize = minSize,
                                      maxSize = maxSize,
                                      eps = 0.0,
                                      scoreType = scoreType, ...)
  fgsea_res$motif=rep(motif_id,nrow(fgsea_res))
  return(fgsea_res)
}


