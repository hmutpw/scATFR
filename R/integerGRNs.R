
setGeneric("integerGRNs", function(x, ...) standardGeneric("integerGRNs"))

#' Integer the gene regulatory network with pre-defined regulators and
#' targets regulatory network
#'
#' @param x A gene regulatory \code{matrix} or \code{SingleCellExperiment} object.
#' For a matrix, rows represent for regulators (eg, TFs) and columns represent
#' for targets.
#' @param baseNets A data.frame including the background regulatory network used
#' for integration. There are three columns in this data.frame: tf, target and weight.
#' Default: NULL. No background regulatory network.
#' @param use_assay The GRN assays used for integration.
#' @param cutoff The cutoff used for filtering the targets for each TF regulon.
#' @param weight The weight of each data.frame.Default: 1. Equal for each data.frame.
#' @param norm_weight Whether to normalize the weight for each data.frame? Default:FALSE.
#' @param keep_neg_regu Whether to keep negative weight for output? Default:FALSE
#' @param keepWeight Whether to keep merged weight in the final output. Default: TRUE.
#' @param beta0 Numeric. The sparsity parameter.
#' @param ncores Numeric. The number of cores used.
#' @param verbose Whether to print messages on screen?
#' @param ... Other parameters passed to \code{\link{integraNets}}.
#'
#' @return A list or an object with regulon list.
#' @export
#'
#' @examples
#' @rdname integerGRNs
setMethod("integerGRNs","matrix",function(x,
                                          baseNets=NULL,
                                          cutoff=NULL,
                                          weight = c(1,1),
                                          norm_weight = FALSE,
                                          keep_neg_regu = FALSE,
                                          keepWeight=TRUE,
                                          beta0 = 0,
                                          ncores=1,
                                          verbose=interactive(),
                                          ...){
  weightMat <-reshape2::melt(as.matrix(x), value.name = "weight")
  colnames(weightMat) <- c("tf","target","weight")
  #---no baseNets input, set cutoff directly
  if(is.null(baseNets)){
    if(verbose) message("No baseNets input, get regulons by setting cutoff...")
    if(is.null(cutoff)){
      cutoff <- round(ncol(x)*0.05)
      weightMat <- weightMat %>% dplyr::group_by(tf) %>% dplyr::top_n(.,n=cutoff, wt = weight)
    }else{
      weightMat <- weightMat[weightMat$weight>=cutoff,]
    }
    if(keepWeight){
      weightMat_list <- split(weightMat,weightMat$tf)
      regulons <- purrr::map(.x = weightMat_list,.f = function(x){
        targets <- x$weight
        names(targets) <- x$target
      })
    }else{
      regulons <- split(weightMat$target,weightMat$tf)
    }
  }else{
    #---merge baseNets input, you can also set cutoff directly
    if(verbose) message("Mergeing weighted tf-target matrix...")
    weightMat <-reshape2::melt(as.matrix(x), value.name = "weight")
    merged_net <- integraNets(grnNets = list(weightMat, baseNets), norm_weight=norm_weight,
                              keep_neg_regu=keep_neg_regu, beta0=beta0)
    if(is.null(cutoff)){
      if(verbose) message("Inferring weight cutoff for each regulon...")
      merged_net <- split(merged_net,merged_net$tf)
      cl <- parallel::makeCluster(ncores)
      #parallel::clusterExport(cl=cl, varlist = c("keepWeight",".getWeightThreshold"))
      regulons <- pbapply::pblapply(X = merged_net, FUN = .getWeightThreshold, keepWeight=keepWeight, cl=cl, ...)
      parallel::stopCluster(cl)
      #regulons <- purrr::map(.x = merged_net, .f = .getWeightThreshold, keepWeight=keepWeight,...)
    }else{
      merged_net <- merged_net[merged_net$weight>=cutoff,]
      if(keepWeight){
        merged_net_list <- split(merged_net, merged_net$tf)
        regulons <- purrr::map(.x = merged_net_list,.f = function(x){
          targets <- x$weight
          names(targets) <- x$target
        })
      }else{
        regulons <- split(merged_net$target,merged_net$tf)
      }
    }
  }
  regulons
})

#' @export
#' @rdname integerGRNs
setMethod("integerGRNs","SingleCellExperiment",function(x,
                                                        baseNets=NULL,
                                                        use_assay=NULL,
                                                        cutoff=NULL,
                                                        weight = c(1,1),
                                                        norm_weight = FALSE,
                                                        keep_neg_regu = TRUE,
                                                        beta0 = 0,
                                                        ncores=1,
                                                        verbose=interactive(),
                                                        ...){
  weightMat <- GRN(x = x, i = use_assay)
  weightMat <-reshape2::melt(as.matrix(x), value.name = "weight")
  colnames(weightMat) <- c("tf","target","weight")
  #---no baseNets input, set cutoff directly
  if(is.null(baseNets)){
    if(verbose) message("No baseNets input, get regulons by setting cutoff...")
    if(is.null(cutoff)){
      cutoff <- round(ncol(x)*0.05)
      weightMat <- weightMat %>% dplyr::group_by(tf) %>% dplyr::top_n(.,n=cutoff, wt = weight)
    }else{
      weightMat <- weightMat[weightMat$weight>=cutoff,]
    }
    if(keepWeight){
      weightMat_list <- split(weightMat,weightMat$tf)
      regulons <- purrr::map(.x = weightMat_list,.f = function(x){
        targets <- x$weight
        names(targets) <- x$target
      })
    }else{
      regulons <- split(weightMat$target,weightMat$tf)
    }
  }else{
    #---merge baseNets input, you can also set cutoff directly
    if(verbose) message("Mergeing weight tf-target matrix...")
    weightMat <-reshape2::melt(as.matrix(x), value.name = "weight")
    merged_net <- integraNets(grnNets = list(weightMat, baseNets), norm_weight=norm_weight,
                              keep_neg_regu=keep_neg_regu, beta0=beta0)
    if(is.null(cutoff)){
      if(verbose) message("Inferring weight cutoff for each regulon...")
      merged_net <- split(merged_net,merged_net$tf)
      cl <- parallel::makeCluster(ncores)
      parallel::clusterExport(cl=cl, varlist = c("keepWeight",".getWeightThreshold"))
      regulons <- pbapply::pblapply(X = merged_net, FUN = .getWeightThreshold, keepWeight=keepWeight, cl=cl, ...)
      parallel::stopCluster(cl)
      #regulons <- purrr::map(.x = merged_net, .f = .getWeightThreshold, keepWeight=keepWeight,...)
    }else{
      merged_net <- merged_net[merged_net$weight>=cutoff,]
      if(keepWeight){
        merged_net_list <- split(merged_net, merged_net$tf)
        regulons <- purrr::map(.x = merged_net_list,.f = function(x){
          targets <- x$weight
          names(targets) <- x$target
        })
      }else{
        regulons <- split(merged_net$target,merged_net$tf)
      }
    }
  }
  regulons
})

#' Integrating weighted regulator-target regulatory networks
#'
#' @param grnNets A list of weighted regulator-target regulatory networks from different
#' resources, for example: motif prediction, Chip-seq and Knockout experiment. For each
#' network in the list, there must three columns corresponding to regulator, target and
#' weight, respectively. If one of your network do not have weight, you can set 1 to all
#' the regulator-target pairs in `weight` column.
#' @param weight Positive number. Specific the weight of resources. default: 1
#' @param norm_weight Whether to normalize the weight to [0,1] of regulator-target pairs
#' from each resource. Default: FALSE
#' @param keep_neg_regu Whether to keep negative regulatory relationships between regulators
#' and targets. Default: TRUE.
#' @param join_type The merge types among different data frames. Default:full.
#' @param beta0 Numeric. The sparsity parameter.
#'
#' @return A data.frame with integrated weighted regulator-target regulatory networks
#' @importFrom dplyr full_join
#' @importFrom purrr reduce pmap_dbl
#'
#' @export
#'
#' @examples
integraNets <- function(grnNets, weight = rep(1,length(grnNets)), norm_weight = FALSE,
                        keep_neg_regu = TRUE, join_type = c("full","inner","left","right"),
                        beta0 = 0, verbose = interactive()){
  #------check arguments
  if(!is.list(grnNets)) stop("The parameter grnNets must be a list!")
  if(length(grnNets)<2) stop("There need more than two (included) networks for integration!")
  if(verbose) message("[1] Checking the input list...")
  grnNets <- .check_df(grnNets)
  net_names <- names(grnNets)
  if(is.null(net_names)){
    net_names <- paste("Net",1:length(grnNets),sep = "")
  }
  names(grnNets) <- net_names
  #---merge all tables
  if(verbose) message("[2] Mergeing the weighted data.frame ...")
  join_type <- match.arg(join_type)
  if(join_type=="full"){
    grn_tbl <- grnNets %>% purrr::reduce(dplyr::full_join,by=c("tf","target"))
  }else if(join_type=="inner"){
    grn_tbl <- grnNets %>% purrr::reduce(dplyr::inner_join,by=c("tf","target"))
  }else if(join_type=="left"){
    grn_tbl <- grnNets %>% purrr::reduce(dplyr::left_join,by=c("tf","target"))
  }else if(join_type=="right"){
    grn_tbl <- grnNets %>% purrr::reduce(dplyr::right_join,by=c("tf","target"))
  }
  grn_tbl[is.na(grn_tbl)] <- 0
  #---weight direction
  if(keep_neg_regu){
    regu_direc <- purrr::pmap_dbl(.l = grn_tbl[,-c(1:2)],.f = function(x){
      ifelse(any(x<0),-1,1)})
  }
  #--- check weights and normalization
  if(!is.numeric(weight)) stop("The weight must be numeric!")
  if(length(weight)!=length(grnNets)) stop("The weight must be the same length with grnNets!")
  names(weight) <- net_names
  if(verbose) message("[3] Normalizing the merged weight...")
  grnNets <- pbapply::pblapply(net_names, function(x, grnNets, weight, norm_weight){
    tbl_one <- grnNets[[x]]
    colnames(tbl_one) <- c("regulator","target", x)
    if(norm_weight){
      tbl_one[[x]] <- (tbl_one[[x]]-min(tbl_one[[x]]))/(max(tbl_one[[x]])-min(tbl_one[[x]]))
    }
    tbl_one[[x]] <- tbl_one[[x]]*weight[x]
    tbl_one
  }, grnNets=grnNets, weight=weight, norm_weight=norm_weight)
  names(grnNets) <- net_names
  if(verbose) message("[4] Calculating the merged weight...")
  merge_weight <- purrr::pmap_dbl(.l = grn_tbl[,-c(1:2)],
                                  .f = function(..., beta0=0){
                                    x <- as.numeric(unlist(list(...)))
                                    1/(1+exp(-1*sum(beta0,abs(x))))
                                    }, beta0=beta0)
  final_tbl <- dplyr::bind_cols(grn_tbl[,1:2],weight=merge_weight) %>% unique(.)
  final_tbl[order(final_tbl$tf,final_tbl$weight,final_tbl$target,decreasing=c(TRUE, TRUE, FALSE),method = "radix"),]
}

#---check the GRNs in your data
.check_df <- function(df_list){
  df_class <- sapply(df_list, is.data.frame)
  if(!all(df_class)) stop("Not all the data stored in you list is data.frame!")
  df_cols <- sapply(df_list, ncol)
  if(any(df_cols < 2)) stop("There must be at least 2 columns in your data.frame!")
  if(any(df_cols > 3)){
    warning("There are more than 3 columns in your data.frame, only using the first three columns!")
  }
  df_list <- lapply(df_list,function(x){
    if(ncol(x)>3){
      x <- x[,1:3]
      if(!is.numeric(x[[3]])) warning("The third column is not numeric")
      x[[3]] <- as.numeric(x[[3]])
    }else if(ncol(x)==2){
      x[[3]] <- rep(1,nrow(x))
    }
    colnames(x) <- c("tf","target","weight")
    x
  })
  df_list
}


#---calculate the cutoff of score using method from AUCell package with some modification
.getWeightThreshold <- function(df, keepWeight=TRUE, plotHist=FALSE,
                                smallestPopPercent=.1, densAdjust=2, thrP=0.01, nBreaks=100){
  #---progress bar
  if(length(unique(df$tf))!=1) stop("Only support one tf input!")
  auc <- df[["weight"]]
  names(auc) <- df[["target"]]
  gSetName <- unique(df[["tf"]])

  nCells <- length(auc)
  skipGlobal <- TRUE
  skipRed <- FALSE
  skipSmallDens <- FALSE
  commentMsg <- ""
  aucThrs <- c()

  notPopPercent <- 1 - smallestPopPercent
  if(sum(auc==0) > (nCells*notPopPercent))
  {
    skipGlobal <- FALSE
    commentMsg <- paste(commentMsg,
                        round((sum(auc==0)/nCells)*100),
                        "% (more than ", notPopPercent,"%) of AUC are zero. ", sep="")
  }

  meanAUC <- mean(auc)
  sdAUC <- sd(auc)
  maybeNormalDistr <- !suppressWarnings(
    ks.test(auc, rnorm(max(100,length(auc)),mean=meanAUC, sd = sdAUC),
            alternative = "less")$p.value < .01)
  if(maybeNormalDistr){
    commentMsg <- paste0(commentMsg,
                         "The AUC might follow a normal distribution (random gene-set?). ")
    skipGlobal <- FALSE

    # aucThrs["outlierOfGlobal"] <- meanAUC + 2*sdAUC
    aucThrs["outlierOfGlobal"] <- qnorm(1-(thrP/nCells), mean=meanAUC, sd=sdAUC)
  }

  #V6
  histogram <- hist(c(0, auc/max(auc)), breaks=100, plot=FALSE)$count
  if((sum(histogram[1:5]) / sum(histogram)) >= notPopPercent*.75) {
    skipGlobal <- FALSE
    skipRed <- TRUE
    skipSmallDens <- TRUE
  }
  if((sum(histogram[1:10]) / sum(histogram)) >= notPopPercent*.50) {
    skipSmallDens <- TRUE
    skipGlobal <- FALSE
    # skipRed <- TRUE ?
    aucThrs["tenPercentOfMax"] <- max(auc)*.10
  }
  # print(skipRed)

  densCurve <- density(auc, adjust=densAdjust, cut=0)
  maximumsDens <- NULL
  inflPoints <- diff(sign(diff(densCurve$y)))
  maximumsDens <- which(inflPoints==-2)
  globalMax <- maximumsDens[which.max(densCurve$y[maximumsDens])]
  minimumDens <- which(inflPoints==2)
  smallMin <- NULL
  if(!skipSmallDens)
    smallMin <- data.table::last(minimumDens[which(minimumDens < globalMax)]) #1prev to max
  minimumDens <- c(smallMin,
                   minimumDens[which(minimumDens > globalMax)]) # all after maximum

  # Density-based threshold (V4):
  # First minimum after the biggest maximum   (adjust=2)
  densTrh <- NULL
  if(length(minimumDens)>0) # && (!skipMinimumDens))
  {
    densTrh <- densCurve$x[min(minimumDens)]
    # Commented on V6
    # Only keep if it is a real inflextion point
    # (i.e. next max at least 5% of the global max)
    if(length(maximumsDens)>0)
    {
      nextMaxs <- maximumsDens[which(densCurve$x[maximumsDens] > densTrh)]
      if((max(densCurve$y[nextMaxs])/max(densCurve$y))<.05)
      {
        densTrh <- NULL
        # print(gSetName)
      }
    }
  }

  ## TO DO: Check special cases with many zeroes
  auc <- sort(auc)
  distrs <- list()
  distrs[["Global_k1"]] <- list(mu=c(meanAUC, NA), sigma=c(sdAUC, NA), x=auc)


  if("mixtools" %in% rownames(installed.packages()))
  {
    na <- capture.output(distrs[["k2"]] <-
                           tryCatch(mixtools::normalmixEM(auc, fast=FALSE, k=2, verb=FALSE),
                                    # With fast, if there are many zeroes, it fails quite often
                                    error = function(e) {
                                      return(NULL)
                                    }))

    na <- capture.output(distrs[["k3"]] <-
                           tryCatch(mixtools::normalmixEM(auc, fast=FALSE, k=3, verb=FALSE),
                                    error = function(e) {
                                      return(NULL)
                                    }))

    if(is.null(distrs[["k2"]]) && is.null(distrs[["k3"]]))
    {
      if(sum(auc==0)<(nCells*notPopPercent*.5))
        skipGlobal <- FALSE    # only if not too many zeroes??
      # qpois(1-(thrP/nCells), 1, log = FALSE)
      # plot(sort(rpois(auc, lambda=var(auc)), decreasing=TRUE))

      # commented V6 why was it here??
      # qPop <- quantile(auc, 1-smallestPopPercent)
      # if(sum(auc<qPop) >0)
      #   distrs[["k2"]] <- list(mu=c(mean(auc[auc<qPop]), NA),
      #                          sigma=c(sd(auc[auc<qPop]), NA),
      #                          lambda=c(1,NA), x=auc)
    }
    # if(!skipGlobal) print(gSetName) Warning?

    if(!is.null(distrs[["k2"]]))
    {
      compL <- which.min(distrs[["k2"]][["mu"]])
      compR <- which.max(distrs[["k2"]][["mu"]])
      ### Check distributions
      # Second distribution is "taller" than first one
      height1 <- .4/distrs[["k2"]][["sigma"]][compL]*
        distrs[["k2"]][["lambda"]][compL]
      height2 <- .4/distrs[["k2"]][["sigma"]][compR]*
        distrs[["k2"]][["lambda"]][compR]
      taller <- height1 < height2
      # Use global distr:
      # Mean of the global distr is included within the SD of the first
      # & Both means are included within the mean+SD of the Global distribution
      globalInclInFirst <-
        (distrs[["Global_k1"]]$mu[1] <
           (distrs[["k2"]][["mu"]][compL]+(1.5*distrs[["k2"]][["sigma"]][compL])))
      includedInGlobal <-
        ((distrs[["k2"]][["mu"]][compL] >
            (distrs[["Global_k1"]]$mu[1]-distrs[["Global_k1"]]$sigma[1])) &&
           (distrs[["k2"]][["mu"]][compR] <
              (distrs[["Global_k1"]]$mu[1]+distrs[["Global_k1"]]$sigma[1])))
      if(taller || (globalInclInFirst && includedInGlobal))
      {
        skipGlobal <- FALSE

        if(globalInclInFirst && includedInGlobal)
          commentMsg <- paste(commentMsg,
                              "The global distribution overlaps the partial distributions. ")
        if(taller && !includedInGlobal)
          commentMsg <- paste(commentMsg, "The right distribution is taller. ")
      }
    }
  }else{
    warning("Package 'mixtools' is not available to calculate the sub-distributions.")
  }

  glProb <- 1-(thrP/nCells + smallestPopPercent)   ## CORRECT?!?!
  aucThrs["Global_k1"] <- qnorm(glProb,# qnorm(1-(thrP/nCells),
                                mean=distrs[["Global_k1"]][["mu"]][1],
                                sd=distrs[["Global_k1"]][["sigma"]][1])
  if(!is.null(distrs[["k2"]]))
  {
    k2_L <- which.min(distrs[["k2"]][["mu"]]) # (sometimes the indexes are shifted)
    aucThrs["L_k2"] <- qnorm(1-(thrP/nCells),
                             mean=distrs[["k2"]][["mu"]][k2_L],
                             sd=distrs[["k2"]][["sigma"]][k2_L])
  }

  if(!is.null(distrs[["k3"]]))
  {
    k3_R <- which.max(distrs[["k3"]][["mu"]]) # R: right distribution
    k3_R_threshold <- qnorm(thrP,
                            mean=distrs[["k3"]][["mu"]][k3_R],
                            sd=distrs[["k3"]][["sigma"]][k3_R])
    if(k3_R_threshold > 0) aucThrs["R_k3"] <- k3_R_threshold
  }

  if(!is.null(densTrh))
  {
    aucThrs["minimumDens"] <- densTrh
  }

  aucThr <- aucThrs
  if(skipGlobal)
    aucThr <- aucThrs[which(!names(aucThrs) %in% "Global_k1")]
  # TO DO: Decide when to merge with GLOBAL

  if(skipRed)
    aucThr <- aucThrs[which(!names(aucThrs) %in% "L_k2")]
  # TO DO: Decide when to merge with GLOBAL

  aucThr <- aucThr[which.max(aucThr)] # to keep name
  if((length(aucThr)>0) && (names(aucThr) == "minimumDens"))
  {
    maximumsDens <- maximumsDens[which(densCurve$y[maximumsDens]>1)]
    if(length(maximumsDens) > 2)
    {
      tmp <- cbind(minimumDens[seq_len(length(maximumsDens)-1)],
                   maximumsDens[-1])
      FCs <- densCurve$y[tmp[,2]]/densCurve$y[tmp[,1]]
      if(any(FCs > 1.5))
        warning(gSetName,
                ":\tCheck the AUC histogram. ",
                "'minimumDens' was selected as the best threshold, ",
                "but there might be several distributions in the AUC.")
    }
  }

  if("minimumDens" %in% names(aucThrs))
    aucThr <- aucThrs["minimumDens"]
  if(length(aucThr)==0)
    aucThr <-  aucThrs[which.max(aucThrs)]

  if(plotHist)
  {
    histInfo <- AUCell_plotHist(aucRow,
                                aucThr=aucThr,
                                nBreaks=nBreaks)
    histMax <- max(histInfo[[gSetName]]$counts)

    # Plot density
    densCurve$y <- densCurve$y*(histMax/max(densCurve$y))
    thisLwd <- ifelse(
      (aucThrs["minimumDens"]==aucThr) &&
        (!is.null(aucThr) && !is.null(aucThrs["minimumDens"])),
      3,
      1)
    lines(densCurve, lty=1, lwd=thisLwd, col="blue")
    if(!is.null(minimumDens))
      points(densCurve$x[minimumDens], densCurve$y[minimumDens],
             pch=16, col="darkblue")

    ### Plot distributions
    scalFact <- 1
    # if(!skipGlobal)
    # {
    aucDistr <- dnorm(distrs[["Global_k1"]][["x"]],
                      mean=distrs[["Global_k1"]][["mu"]][1],
                      sd=distrs[["Global_k1"]][["sigma"]][1])
    scalFact <- (histMax/max(aucDistr))*.95

    thisLwd <- ifelse(aucThrs["Global_k1"]==aucThr, 3, 1)
    lines(distrs[["Global_k1"]][["x"]],
          scalFact * aucDistr,
          col="darkgrey", lwd=thisLwd, lty=2)

    if(!is.null(distrs[["k2"]]))
    {
      aucDistr <- dnorm(distrs[["k2"]][["x"]],
                        mean=distrs[["k2"]][["mu"]][k2_L],
                        sd=distrs[["k2"]][["sigma"]][k2_L])
      scalFact <- (histMax/max(aucDistr))*.95


      thisLwd <- ifelse(aucThrs["k2"]==aucThr, 3, 1)
      lines(distrs[["k2"]][["x"]],
            scalFact * aucDistr,
            col="red", lwd=thisLwd, lty=2)

      rect(distrs[["k2"]][["mu"]][k2_L]-distrs[["k2"]][["sigma"]][k2_L],
           histMax-(histMax*.02),
           distrs[["k2"]][["mu"]][k2_L]+distrs[["k2"]][["sigma"]][k2_L],
           histMax, col="#70000030", border="#00009000")
    }

    # print(aucThrs)
    if((!is.null(distrs[["k3"]])) && ("R_k3" %in% names(aucThrs)))
    {
      k3_L <- which.min(distrs[["k3"]][["mu"]]) # (index position not constant)

      aucDistr2 <- dnorm(distrs[["k3"]][["x"]],
                         mean=distrs[["k3"]][["mu"]][k3_R],
                         sd=distrs[["k3"]][["sigma"]][k3_R])
      scalFact2 <- scalFact *
        (distrs[["k3"]][["lambda"]][k3_R]/distrs[["k3"]][["lambda"]][k3_L])

      thisLwd <- ifelse(aucThrs["k3"]==aucThr, 3, 1)
      lines(distrs[["k3"]][["x"]],
            scalFact2*aucDistr2,
            col="magenta", lwd=thisLwd, lty=2)

      rect(distrs[["k3"]][["mu"]][k3_R]-distrs[["k3"]][["sigma"]][k3_R],
           histMax-(histMax*.02),
           distrs[["k3"]][["mu"]][k3_R]+distrs[["k3"]][["sigma"]][k3_R],
           histMax, col="#80808030", border="#80808030")
    }

    ## Add threshold lines
    aucThrs <- aucThrs[!is.na(aucThrs)]
    if(length(aucThrs)>0)
    {
      pars <- list()
      pars[["Global_k1"]] <- c(col1="#909090", col2="black", pos=.9)
      pars[["L_k2"]] <- c(col1="red", col2="darkred", pos=.8)
      # pars[["Max"]] <- c(col1="grey", col2="black", pos=.4)
      pars[["R_k3"]] <- c(col1="magenta", col2="magenta", pos=.6)
      pars[["minimumDens"]] <- c(col1="blue", col2="darkblue", pos=.4)
      pars[["tenPercentOfMax"]] <- c(col1="darkgreen", col2="darkgreen", pos=.9)
      pars[["outlierOfGlobal"]] <- c(col1="darkgreen", col2="darkgreen", pos=.9)

      for(thr in names(aucThrs))
      {
        thisLwd <- ifelse(aucThrs[thr]==aucThr, 5, 2)
        thisLty <- ifelse(aucThrs[thr]==aucThr, 1, 3)

        abline(v=aucThrs[thr], col=pars[[thr]][1], lwd=thisLwd, lty=thisLty)
        xPos <- aucThrs[thr]*1.01
        if(aucThrs[thr] > (max(auc)*.8))
          xPos <- 0
        if(aucThrs[thr]==aucThr)
          text(xPos, histMax*as.numeric(pars[[thr]][3]),
               pos=4, col=pars[[thr]][2], cex=.8,
               paste("AUC > ", signif(aucThrs[thr],2),
                     "\n(",sum(auc>aucThrs[thr])," cells)", sep=""))
      }
    }
  }
  df <- df[df$weight>=aucThr,]
  if(keepWeight){
    targets <- df[["weight"]]
    names(targets) <- df$target
  }else{
    targets <- df$target
  }
  targets
}


