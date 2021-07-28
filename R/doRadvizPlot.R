#Setter/Getter
setGeneric("RadvizEnrichPlot",function(x, ...) standardGeneric("RadvizEnrichPlot"))

#' RadvizEnrichPlot
#'
#'
#' @param x A \code{\link{data.frame}} or objects used for \code{\link{Radviz}} plot.
#' @param anchors Anchors used to plot.
#' @param color_by Groups used to label the color.
#' @param return_rv Only return \code{\link{getRadviz}} object, not generate figures.
#' @param ... Other parameters passed to \code{\link{getRadviz}}.
#' @rdname RadvizEnrichPlot
#' @export
setMethod("RadvizEnrichPlot","data.frame",function(x,
                                                   anchors=colnames(x),
                                                   color_by=NULL,
                                                   return_rv=FALSE,
                                                   ...){
  rv <- getRadviz(input_frm = x,anchors = anchors, ...)
  rv_plot <- do_RadvizEnrichPlot(rv = rv, color_by = color_by, ...)
  if(return_rv) rv_plot <- rv
  rv_plot
})


####
#RadvizFeaturePlot
setGeneric("RadvizFeaturePlot",function(x, ...) standardGeneric("RadvizFeaturePlot"))

#' RadvizFeaturePlot
#'
#'
#' @param x A \code{\link{data.frame}} or objects used for \code{\link{Radviz}} plot.
#' @param anchors Anchors used to plot.
#' @param use_feature Features used to split the screens.
#' @param return_rv Only return \code{\link{getRadviz}} object, not generate figures.
#' @param ... Other parameters passed to \code{\link{getRadviz}}.
#' @rdname RadvizFeaturePlot
#' @export
setMethod("RadvizFeaturePlot","data.frame",function(x,
                                                    anchors = colnames(x),
                                                    use_feature,
                                                    return_rv = FALSE,
                                                    ...){
  rv <- getRadviz(input_frm = x,anchors = anchors, ...)
  rv_plot <- do_RadvizFeaturePlot(rv = rv, use_feature = use_feature, ...)
  if(return_rv) rv_plot <- rv
  rv_plot
})

#######################
# getRadviz object
#######################
#' getRadviz object
#'
#'
#' @param input_frm Matrix used to generate the Radviz object.
#' @param anchors Dimensions used to generate the map.
#' @param trim_ratio Whether to trim the top and bottom of each anchor which correspond to outliers.
#' @param optim_anchor Whether to optimize the anchor order?
#' @param recenter User defined start anchor.
#' @param label_col Label color, inherited from \code{\link{do.radviz}}.
#' @param label_size Label size, inherited from \code{\link{do.radviz}}.
#' @param ... Other parameters passed to \code{\link{do.radviz}}.
#' @return
#' @export
#' @importFrom Radviz do.L make.S cosine do.optimRadviz get.optim recenter do.radviz
getRadviz <- function(input_frm,
                      anchors=colnames(input_frm),
                      trim_ratio=0.005,
                      optim_anchor=TRUE,
                      recenter=NULL,
                      text_color="orangered4",
                      text_size=8,
                      ...){
  anchors <- intersect(anchors,colnames(input_frm))
  if(length(anchors)==0) stop("No anchors found from colnames of input_frm!")
  exp_mat <- as.matrix(input_frm[,anchors])
  if(!is.numeric(exp_mat)){
    stop("The 'anchors' columns in 'input_frm' must be numeric!")
  }
  #---normalization
  trans <- function(vec, cutoff=trim_ratio){
    Radviz::do.L(vec,fun = function(x) quantile(x,c(cutoff,1-cutoff)))
  }
  #---generate springs
  mat.S <- Radviz::make.S(anchors)
  if(optim_anchor){
    mat.sim <- Radviz::cosine(exp_mat)
    optim.mat <- Radviz::do.optimRadviz(mat.S,mat.sim,iter=10,n=100)
    mat.S <- Radviz::make.S(Radviz::get.optim(optim.mat))
  }
  if(!is.null(recenter) && length(recenter)==1 && (recenter %in% anchors)){
    mat.S <- Radviz::recenter(mat.S, newc = recenter)
  }
  Radviz::do.radviz(x = input_frm, mat.S, trans = trans,
                    label.color = text_color,
                    label.size = text_size, ...)
}

################################
#------feature plot
################################
#' do_RadvizEnrichPlot
#'
#' @param rv A \code{\link{getRadviz}} object generate from \code{\link{do.radviz}}.
#' @param color_by Groups used to label the color.
#' @param outline_circle Whether to add circle in figure?
#' @param plot_type Figure types to plot.
#' @param text_size Text size.
#' @param label_size Label size of point in \code{\link{geom_point()}}.
#' @param label_shape Label shape of point in \code{\link{geom_point()}}.
#' @param ... Other input parameters.
#'
#' @return a ggplot object
#'
#' @export
do_RadvizEnrichPlot <- function(rv,
                                color_by=NULL,
                                outline_circle=TRUE,
                                plot_type=c("point","density","hexagonal","buble"),
                                text_size=8,
                                label_size = 0.5,
                                label_shape=16,
                                ...){
  plot_type <- match.arg(plot_type)
  rv_data_frm <- as.data.frame(rv$proj$data)
  if(outline_circle) rv <- addRadvizCircle(rv)
  if(!is.null(color_by)){
    color_by <- intersect(color_by,colnames(rv_data_frm))
    if(length(color_by)!=1) stop("The 'color_by' shold only have one element!")
  }
  p <- rv$proj
  if(plot_type=="point"){
    p <- p + geom_point(aes_string(color=color_by),shape=label_shape, size=label_size)
  }else if(plot_type=="density"){
    p <- Radviz::smoothRadviz(rv)+geom_point(shape='.',alpha=0.5)
  }else if(plot_type=="hexagonal"){
    p <- Radviz::hexplot(rv)
  }else if(plot_type=="buble"){
    p <- Radviz::bubbleRadviz(rv,group=color_by)
  }
  p
}

################################
#------feature plot
################################
#' do_RadvizFeaturePlot
#'
#'
#' @param rv A \code{\link{getRadviz}} object generate from \code{\link{do.radviz}}
#' @param use_feature Features used to split the screens.
#' @param outline_circle Whether to add circle in figure?
#' @param plot_type Figure types to plot.
#' @param colors Colors for gradient.
#' @param text_size Text size.
#'
#' @return a ggplot object
#'
#' @export
do_RadvizFeaturePlot <- function(rv,
                                 use_feature,
                                 outline_circle=TRUE,
                                 plot_type=c("point","hexagonal"),
                                 colors=scale_color_gradient(low='grey80',high="dodgerblue4"),
                                 text_size=8){
  plot_type <- match.arg(plot_type)
  rv_data_frm <- as.data.frame(rv$proj$data)
  use.feature <- intersect(use_feature, colnames(rv_data_frm))
  if(length(use.feature)==0) stop("The 'use.features' in your input is not found!")
  feature_mat <- as.matrix(rv_data_frm[,use.feature])
  if(!is.numeric(feature_mat)) stop("All the values in your features should be numeric!")

  #---start plot
  if(outline_circle) rv <- addRadvizCircle(rv)
  p <- rv$proj
  plots <- lapply(use.feature,function(fea){

    if(plot_type=="point"){
      geoms <- geom_point(aes_string(color=fea),shape=16)
    }else if(plot_type=="hexagonal"){
      geoms <- stat_summary_hex(aes_string(z=fea), fun=median)
    }
    p+geoms+colors+labs(title=fea)+
      theme(text = element_text(size=text_size,hjust = 1,color="black"),
            plot.title=element_text(size=text_size, hjust=0.5,vjust=0.5))
  })
  patchwork::wrap_plots(plots, ncol = ceiling(sqrt(length(use_feature))))
}

#------add out.liner for circle
addRadvizCircle <- function(rv, color="black"){
  circle <- annotate("path", x=cos(seq(0,2*pi,length.out=100)),y=sin(seq(0,2*pi,length.out=100)),color=color)
  rv$proj <- rv$proj+circle
  rv
}




