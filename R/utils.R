######################
# 1. transform data.frame to gene list
######################
#' Transform data.frame to gene list
#'
#' @param df A data.frame with TF and it's target information.
#' @param tf_col The column names of TFs.
#' @param target_col The column names of targets.
#'
#' @return A list of Transcription factor regulons. The names of list is TFs
#' @export
#'
#' @examples
#' data(dorothea_regulons)
#' mouse_regulons <- dfToList(df=dorothea_regulons$mm_regulons,
#' tf_col="tf", target_col="target")
#'
dfToList <- function(df, tf_col, target_col){
  regulons <- split(df[[target_col]],df[[tf_col]])
  lapply(regulons, unique)
}

######################
# 2. Get homologene from local database
######################
#' Get homologene from local database
#'
#' @param input_gene Input genes, which can be gene symbols or entrez gene ids.
#' @param input_tax Input gene taxonomy id or names. To browse the available organisms,
#' please use \code{\link{browserTaxonomy()}}
#' @param output_tax Output gene taxonomy id or names. To browse the available organisms,
#' please use \code{\link{browserTaxonomy()}}
#' @param output_type Output types of genes. Default: gene_symbol.
#' @param db Database used in our homologene transformation.
#'
#' @return A data.frame include homologene gene information
#' @export
#'
#' @examples
getHomologGene <- function(input_gene=NULL,
                           input_tax="human",
                           output_tax="mouse",
                           output_type=c("gene_symbol","gene_id","both"),
                           db=scATFR::homologene_data){
  require(dplyr)
  output_type <- match.arg(output_type)
  #---names
  short_name <- c("human","mouse","rat","zebrafish","chicken")
  names(short_name) <- c("Homo sapiens","Mus musculus","Rattus norvegicus","Danio rerio","Gallus gallus")
  tax_tbl <- browserTaxonomy()
  tax_tbl$short_name <- short_name[tax_tbl$Taxonomy_Name]
  tax_tbl[is.na(tax_tbl)] <- ""
  #------check input id
  input_tax_id <- .check_tax_id(input_tax=input_tax, tax_tbl = tax_tbl)
  output_tax_id <- .check_tax_id(input_tax=output_tax, tax_tbl = tax_tbl)
  #------check input gene
  if(is.null(input_gene)) input_gene <- db$Gene_Symbol
  input_gene <- unique(input_gene)
  #---input
  input_tax_tab <- db %>%
    dplyr::filter(Taxonomy_ID %in% input_tax_id & (Gene_Symbol %in% input_gene | Gene_ID %in% input_gene)) %>%
    dplyr::select(HID, Gene_ID,Gene_Symbol)
  names(input_tax_tab)[2:3] <- paste(c("Gene_ID","Gene_Symbol"),input_tax,sep="_")
  #---output
  output_tax_tab <- db %>%
    dplyr::filter(Taxonomy_ID %in% output_tax_id & HID %in% input_tax_tab$HID) %>%
    dplyr::select(HID, Gene_ID, Gene_Symbol)
  names(output_tax_tab)[2:3] <- paste(c("Gene_ID","Gene_Symbol"),output_tax,sep="_")

  output_tab <- dplyr::inner_join(input_tax_tab,output_tax_tab,by="HID",keep=FALSE,na_matches="never") %>%
    dplyr::select(2,3,4,5)
  #------arrage output
  if(output_type=="gene_id"){
    output_tab <- output_tab[,grep("Gene_ID",colnames(output_tab))] %>% dplyr::distinct()
  }else if(output_type=="gene_symbol"){
    output_tab <- output_tab[,grep("Gene_Symbol",colnames(output_tab))] %>% dplyr::distinct()
  }
  output_tab$sortBy <- factor(output_tab[,1], levels = input_gene)
  output_tab <- dplyr::arrange(output_tab, sortBy)
  output_tab$sortBy <- NULL
  output_tab
}

######
#---check output
.check_tax_id <- function(input_tax, tax_tbl=browserTaxonomy()){
  if(is.integer(input_tax)){
    if(!(input_tax %in% tax_tbl$Taxonomy_ID)) stop("No 'Taxonomy_ID' found for your input: ", input_tax)
  }else if(is.character(input_tax)){
    if(input_tax %in% tax_tbl$Taxonomy_Name){
      input_tax <- tax_tbl[tax_tbl$Taxonomy_Name==input_tax,"Taxonomy_ID"]
    }else if(input_tax %in% tax_tbl$short_name){
      input_tax <- tax_tbl[tax_tbl$short_name==input_tax,"Taxonomy_ID"]
    }else{
      if(as.integer(input_tax) %in% tax_tbl$Taxonomy_ID){
        input_tax <- as.integer(input_tax)
      }else{
        stop("No 'Taxonomy_ID' found for your input: ", input_tax)
      }
    }
  }else{
    stop("The input_tax is a ", class(input_tax), ", it must be intiger or character!")
  }
  input_tax
}

#---browser
browserTaxonomy <- function(db=scATFR::homologene_data){
  require(dplyr)
  as.data.frame(db %>% select(Taxonomy_ID, Taxonomy_Name) %>%
                  distinct())
}

######
#---get homolog gene list
#' Title
#'
#' @param input_gene Input gene list
#' @param input_tax Input gene taxonomy id or names. To browse the available organisms,
#' please use \code{\link{browserTaxonomy()}}
#' @param output_tax Output gene taxonomy id or names. To browse the available organisms,
#' please use \code{\link{browserTaxonomy()}}
#' @param isListNameGene Whether the gene names is gene symbols?
#' @param ... Others parameters passed to \code{\link{getHomologGene})
#'
#' @return A gene list with homolog transformed list.
#' @export
#'
#' @examples
mapHomologList <- function(input_gene,
                           input_tax="human",
                           output_tax="mouse",
                           isListNameGene=FALSE, ...){
  input_gene <- as.list(input_gene)
  if(isListNameGene){
    all_genes <- unique(c(names(input_gene),unlist(input_gene)))
  }else{
    all_genes <- unique(unlist(input_gene))
  }
  gene_tbl <- getHomologGene(input_gene = all_genes,
                             input_tax=input_tax,
                             output_tax=output_tax,...)
  colnames(gene_tbl) <- c("human","mouse")
  output_list <- lapply(input_gene,function(x, gene_tbl){
    out <- gene_tbl %>% dplyr::filter(human %in% x)
    out <- out[['mouse']]
    out[!is.na(out)]},gene_tbl)
  if(isListNameGene){
    gene_tbl_filter <- gene_tbl %>% dplyr::filter(human %in% names(output_list))
    gene_tbl_filter$list_data <- output_list[gene_tbl_filter[["human"]]]
    output_list <- gene_tbl_filter$list_data
    names(output_list) <- gene_tbl_filter[["mouse"]]
  }
  output_list
}

######################
# 3. Obtaining TFRs from DoRoThEA
######################

#' Get TFRs From DoRoThEA
#'
#' @param input The regulons \code{data.frame} from \pkg{dorothea} package
#' @param levels The levels of data used as regulons.
#' @param min_target_num Miniumn target gene bumber.
#' @param max_target_num Maxiumn target gene bumber.
#'
#' @export
tfrFromDoRoThEA <- function(input,
                            levels = c("A","B","C","D","E"),
                            min_target_num = 5,
                            max_target_num = 1000){
  levels <- match.arg(levels, several.ok = TRUE)
  input <- as.data.frame(input[input$confidence %in% levels,])
  input$tf <- as.factor(input$tf)
  regulons <- split(input$target, input$tf)
  regulon_gene_num <- sapply(regulons, function(x){
    length(x[!is.na(x)])>=min_target_num && length(x[!is.na(x)])<=max_target_num})
  regulons[regulon_gene_num]
}

######################
# 4. Obtaining TFRs from cisTarget
######################

#' extract TF regulons from Rcistarget database
#'
#' @param motifRnking The \code{\link{Rcistarget}} motif ranking file
#' @param motif_anno The motif annotation file
#' @param cutoff The cutoff of motif enrichment.
#' @param use_TFs User defined TFs.
#' @param motif_level The levels of motif used for filtering.
#' @param levels The level of motif used.
#'
#' @return
#'
#' @examples
#' @export
tfrFromRcisTarget <- function(motifRnking,
                              motif_anno,
                              cutoff=0.01,
                              use_TFs = NULL,
                              motif_level=c("high", "low", "custom"),
                              levels = c("directAnnotation","inferredBy_Orthology",
                                         "inferredBy_MotifSimilarity",
                                         "inferredBy_MotifSimilarity_n_Orthology")){
  motif_level <- match.arg(motif_level)
  if(motif_level=="high"){
    levels <- c("directAnnotation","inferredBy_Orthology")
  }else if(motif_level=="low"){
    levels <- c("directAnnotation","inferredBy_Orthology","inferredBy_MotifSimilarity",
                "inferredBy_MotifSimilarity_n_Orthology")
  }else{
    levels <- match.arg(levels,several.ok = TRUE)
  }
  ranking_tbl <- RcisTarget::getRanking(object = motifRnking)
  motif_anno_tbl <- dplyr::as_tibble(motif_anno) %>% dplyr::filter(annotationSource %in% levels)
  if(!is.null(use_TFs)){
    use_TFs <- intersect(use_TFs,motif_anno_tbl$TF)
    if(length(use_TFs)==0) stop("The TFs you input is not found!")
    motif_anno_tbl <- motif_anno_tbl %>% dplyr::filter(TF %in% use_TFs)
  }
  overlap_motif_id <- intersect(ranking_tbl$features, motif_anno_tbl$motif)
  if(length(overlap_motif_id)==0) stop("The motifRnking and motif_anno is not matching!")
  ranking_tbl <- ranking_tbl %>% dplyr::filter(features %in% overlap_motif_id)
  motif_anno_tbl <- motif_anno_tbl %>% dplyr::filter(motif %in% overlap_motif_id)

  #---get top ranked targets for each motif
  ranking_motif <- as.matrix(ranking_tbl[,-1])
  row.names(ranking_motif) <- ranking_tbl[['features']]
  message("Extracting motif list from motif ranking...")
  ranking_motif_list <- pbapply::pbapply(ranking_motif,1,function(x, cutoff){
    list(names(head(sort(x),round(length(x)*cutoff))))
  },cutoff=cutoff)
  ranking_motif_list <- lapply(ranking_motif_list,unlist)
  #---get enriched motif ids for each TF
  tf_motif_list <- split(motif_anno_tbl$motif,motif_anno_tbl$TF)
  message("Mergeing TF regulon list...")
  tf_regulon_list <- pbapply::pblapply(tf_motif_list,function(x,motif_list){
    unique(unlist(motif_list[intersect(x,names(motif_list))]))
  },motif_list=ranking_motif_list)
  tf_regulon_list
}


######################
# 5. re-level the cell-types
######################
#' @export
reLevel <- function(vect, level_order=unique(vect)){
  factor(vect, levels=level_order)
}


######################
# 6. filtering gene expression matrix
######################








