#'  Create an object for transcription factor regulon analysis
#'
#' Create an \code{\link{SingleCellExperiment}} object for transcription factor regulon (TFR) analysis
#'
#' @param exp_data Expression matrix contains expression data for individual samples.
#' Could be non-normalized data such as counts and UMIs or normalized data such as TPM or FPKM.
#' @param col_data A \code{data.frame} contains sample annotation information, such as cell types,
#' tissues and conditions. The row names of `col_data` must be the same with the column names of
#' `exp_data`. Default: null.
#' @param row_data A \code{data.frame} contains feature annotation information. Could be gene names
#' or gene types. The row names of `row_data` must be the same with the row names of `exp_data`.
#' Default: null.
#' @param regulons A list of transcription factor regulons (TFRs).
#'
#' @return A \code{\link{SingleCellExperiment}} object.
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors SimpleList
#' @export
#'
#' @examples
#' library(scATFR)
#' #load example data
#' expression_data <- readRDS(system.file("extdata", "mouse_HSC_formation_expression.rds", package = "scATFR"))
#' colData <- readRDS(system.file("extdata", "mouse_HSC_formation_colData.rds", package = "scATFR"))
#' #create atfr object.
#' atfr <- scATFR(exp_data = expression_data, col_data=colData)
#'
#'
scATFR <- function(exp_data, col_data=NULL, row_data=NULL, regulons=NULL){
  if(!(is.matrix(exp_data) || is(exp_data, "dgCMatrix"))){
    warnings("Your input exp_data is not matrix!")
  }
  if(is.null(col_data)){
    col_data=data.frame(sample=colnames(exp_data),
                        row.names = colnames(exp_data),
                        check.names = FALSE)
  }

  exp_data <- as.matrix(exp_data)
  if(!is.numeric(exp_data)) stop("There is non-numeric values in your input 'exp_data'!")
  if(any(is.na(exp_data))) stop("NA exists in your input 'exp_data'!")
  exp_data <- as(exp_data,"dgCMatrix")
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts=exp_data),
                                                    rowData = row_data,
                                                    colData = col_data)
  metadata(sce)$scATFR <- S4Vectors::SimpleList(regulons=S4Vectors::SimpleList(),
                                                GRNs=S4Vectors::SimpleList(),
                                                atfrMetadata=S4Vectors::SimpleList())
  if(!is.null(regulons)){
    if(!is(regulons,"list")) stop("The regulons should be list!")
    Regulon(x = sce, i = "regulon") <- regulons
  }
  metadata(sce)$atfr_version <- Biobase::package.version("scATFR")
  sce
}
