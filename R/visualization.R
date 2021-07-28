
##############
# 1. heatmap
##############
setGeneric("heatmapPlot", function(x, ...) standardGeneric("heatmapPlot"))

setMethod("heatmapPlot", "SingleCellExperiment",function(x,
                                                         use_feature = NULL,
                                                         regulon_assay=TRUE,
                                                         use_assay=1L,
                                                         use_sample=NULL,
                                                         anno_col_item = "cluster",
                                                         merge_cells_by = NULL,
                                                         ...){
  if(regulon_assay){
    input_mat <- as.matrix(assay(altExp(x), use_assay))
  }else{
    input_mat <- as.matrix(assay(x, use_assay))
  }
  if(is.null(use_feature)){
    use_feature <- diffRegulons(x)[['regulon']]
  }
  use_feature <- intersect(use_feature,row.names(input_mat))
  if(length(use_feature)==0) stop("No feature found from your input!")
  if(!is.null(use_sample)){
    use_sample <- intersect(use_sample, colnames(input_mat))
    if(length(use_sample)==0) stop("No samples found from your input!")
  }else{
    use_sample <- colnames(input_mat)
  }
  feature_mat <- input_mat[use_feature,use_sample]
  #---column annotation
  if(!is.null(anno_col_item)){
    anno_col_item <- intersect(anno_col_item, colnames(colData(x)))
    if(length(anno_col_item)==0) stop("No items forund from in colData!")
    anno_col <- as.data.frame(colData(x)[use_sample,anno_col_item,drop=FALSE])
  }else{
    anno_col <- NA
  }
  if(!is.null(merge_cells_by)){
    if(!(merge_cells_by %in% colnames(colData(x)))) stop("The ",merge_cells_by,
                                                         "colums is not found in colData(x)!")
    merge_cells <- as.data.frame(colData(x)[use_sample,merge_cells_by,drop=FALSE])
    input_mat <- t(apply(input_mat,1,function(x){tapply(x, merge_cells[names(x),1],mean)}))
    anno_col <- NA
  }
  pheatmap::pheatmap(mat = input_mat, annotation_col = anno_col, scale="row",...)
})

###########
# 2. expression plot
##########





