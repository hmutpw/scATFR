setGeneric("inferRegulonFun", function(x, ...) standardGeneric("inferRegulonFun"))

#' Inferring the regulon function using pathway information
#'
#' This method estimate the function similarities of TFRs by using
#' Jaccard/Tanimoto similarity coefficients.
#'
#' @param x Object used for function inferring.
#' @param pathway_list A list of pathway gene sets. Names of the list were pathway names,
#' while the value of each pathway is genes.
#' @param method Method used compute p-values in \pkg{jaccard} package. The "asymptotic"
#' is the fastest, while the "bootstrap" is better but a little bit slow. The "exact"
#' solution is slow for a moderately large vector. Default: "asymptotic".
#' @param nperm Number of . Used only for "bootstrap" method.
#' @param ncores Number of cores used for parallel calculation.
#' @param sig_fdr Set the cutoff of similarity significance between \code{regulon_list} and \code{pathway_list}.
#'
#' @return a list of similarity test \code{data.frame} or a \code{matrix} of
#' normalized similarity scores.
#' @export
#'
#' @examples
#' @rdname inferRegulonFun
setMethod("inferRegulonFun", "SingleCellExperiment",function(x,
                                                             pathway_list,
                                                             regulon_list=NULL,
                                                             return_items = c("all","similarity"),
                                                             method = c("asymptotic","bootstrap","mca","exact"),
                                                             nperm=1000,
                                                             ncores=1,
                                                             sig_fdr=0.05){
  method <- match.arg(method)
  #---get regulon list
  diff_reg <- dplyr::tibble(diffRegulons(x)) %>% filter(padj<0.05 & CSS>0.3)
  if(is.null(diff_reg) && is.null(regulon_list)){
    stop("No regulon list found in your object 'x', and no 'regulon_list' input.\n",
         "Please run `findAllCTSRegulons()` or offer your own `regulon_list`!")
  }
  if(is.null(regulon_list) && !is.null(diff_reg)){
    message("Using regolon list in your object!")
    regulon_tbl <- diff_reg %>% select(regulon,targets) %>% distinct()
    regulon_list <- regulon_tbl$targets
    names(regulon_list) <- regulon_tbl$regulon
  }
  regulon_list <- as.list(regulon_list)
  #---pathway list
  pathway_list <- as.list(pathway_list)
  similarity_test <- calcuModuleSimilar(tfr_list= regulon_list,
                                        pathway_list = pathway_list,
                                        return_items = return_items,
                                        method = method,
                                        nperm=nperm,
                                        ncores=ncores)
  similarity_test <- similarity_test %>% filter(padj<sig_fdr)
  RegulonPathway(x) <- similarity_test
  similarity_test
})

#' Estimating the jaccard similarity between two gene lists
#'
#' This method estimate the similarities between two gene lists by using
#' Jaccard/Tanimoto similarity coefficients.
#'
#' @param tfr_list A list of transcription factor regulons. Names of the list were TF genes,
#' while the value of each list is its target genes.
#' @param pathway_list A list of pathway gene sets. Names of the list were pathway names,
#' while the value of each pathway is genes.
#' @param return_items The return value of similarity values. Default: "all".
#' @param method Method used compute p-values in \pkg{jaccard} package. The `asymptotic`
#' is the fastest, while the `bootstrap` is better but a little bit slow. The `exact`
#' solution is slow for a moderately large vector. Default: `asymptotic`.
#' @param nperm Number of . Used only for `bootstrap` method.
#' @param ncores Number of cores used for parallel calculation.
#'
#' @return a list of similarity test \code{data.frame} or a \code{matrix} of
#' normalized similarity scores.
#' @export
calcuModuleSimilar <- function(tfr_list,
                               pathway_list,
                               return_items = c("all","similarity"),
                               method = c("asymptotic","bootstrap","mca","exact"),
                               nperm=1000,
                               ncores=1){
  method <- match.arg(method)
  return_items <- match.arg(return_items)
  if(!is.list(tfr_list) || !is.list(pathway_list)){
    stop("The 'tfr_list' and the 'pathway_list' must be list!")
  }
  set.seed(10086)
  #------check args
  total_genes <- union(unique(unlist(tfr_list)),unique(unlist(pathway_list)))
  tfr_name <- names(tfr_list)
  pathway_name <- names(pathway_list)
  message("Perform similarity analysis using ",length(tfr_name)," regulons with ",
          length(pathway_name)," pathways. Using method: ",method)
  #---generate binary matrix for jaccard similarity analysis
  binary_mat <- matrix(0L,nrow=length(total_genes),ncol=length(c(tfr_name,pathway_name)),
                       dimnames = list(total_genes, c(tfr_name,pathway_name)))
  for(i in tfr_name){
    genes <- tfr_list[[i]]
    binary_mat[genes,i] <- 1L
  }
  for(j in pathway_name){
    genes <- pathway_list[[j]]
    binary_mat[genes,j] <- 1L
  }
  #---jaccard test
  #------par
  cl <- parallel::makeCluster(ncores)
  parallel::clusterEvalQ(cl,library(pbapply))
  parallel::clusterEvalQ(cl,library(jaccard))
  jaccard_list <- pbapply::pblapply(tfr_name,function(x,
                                                      pathway_name,
                                                      binary_mat,
                                                      method = c("asymptotic","bootstrap","mca","exact"),
                                                      nperm=1000){
    method <- match.arg(method)
    jaccard_res <- apply(binary_mat[,pathway_name,drop=FALSE],2,function(p){
      if(method=="bootstrap"){
        test_res <- suppressMessages(
          jaccard::jaccard.test(x=p,
                                y= binary_mat[names(p),x,drop=TRUE],
                                method = method, B=nperm))
      }else{
        test_res <- jaccard::jaccard.test(x=p,
                                          y= binary_mat[names(p),x,drop=TRUE],
                                          method = method)
      }
      c(similarity=test_res$statistics,pval=test_res$pvalue)
    })
    jaccard_res <- as.data.frame(t(jaccard_res))
    jaccard_res$pathway <- row.names(jaccard_res)
    jaccard_res$padj <- p.adjust(jaccard_res$pval,method = "BH")
    jaccard_res[,c("pathway","similarity","pval","padj")]
  }, pathway_name=pathway_name, binary_mat=binary_mat,method=method,nperm=nperm,cl=cl)
  names(jaccard_list) <- tfr_name
  parallel::stopCluster(cl)
  if(return_items=="similarity"){
    jaccard_list <- sapply(jaccard_list,function(x){x[,"similarity",drop=T]})
    names(jaccard_list) <- pathway_name
  }else{
    jaccard_list <- data.table::rbindlist(jaccard_list,use.names = TRUE, idcol = "regulon")
    jaccard_list[,targets:=tfr_list[jaccard_list$regulon]]
  }
  dplyr::tibble(jaccard_list)
}

