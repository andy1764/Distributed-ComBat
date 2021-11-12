# Author: Andrew Chen, andrewac@sas.upenn.edu
# Date: September 7, 2020

#' Distributed ComBat step at each site
#'
#' @param dat A \emph{p x n} matrix (or object coercible by
#'   \link[base]{as.matrix} to a numeric matrix) of observations where \emph{p}
#'   is the number of features and \emph{n} is the number of subjects.
#' @param batch Factor indicating batch. Needs to have the same levels across
#'   all individual sites, but can have multiple batches per site (i.e.
#'   multiple levels in each site)
#' @param mod Optional design matrix of covariates to preserve, usually from 
#'    \link[stats]{model.matrix}. This matrix needs to have the same columns
#'    across sites. The rows must be in the same order as the data columns.
#' @param ref.batch Optional, reference batch used to determine target mean and
#'   variance. Must be specified for all sites.
#' @param central.out Output list from \code{distributedCombat_central}. Output
#'   of \code{distributedCombat_site} will depend on the values of 
#'   \code{central.out}. If \code{NULL}, then the output will be sufficient for
#'   estimation of \code{B.hat}. If \code{B.hat} is provided, then the output
#'   will be sufficient for estimation of \code{sigma} or for harmonization if
#'   \code{mean.only} is \code{TRUE}. If \code{sigma} is provided, then
#'   harmonization will be performed.
#' @param eb If \code{TRUE}, the empirical Bayes step is used to pool
#'   information across features, as per the original ComBat methodology. If
#'   \code{FALSE}, adjustments are made for each feature individually.
#'   Recommended left as \code{TRUE}.
#' @param parametric If \code{TRUE}, parametric priors are used for the
#'   empirical Bayes step, otherwise non-parametric priors are used. See
#'   neuroComBat package for more details. 
#' @param mean.only If \code{TRUE}, distributed ComBat does not harmonize the
#'   variance of features.
#' @param verbose If \code{TRUE}, print progress updates to the console.
#' @param file File name of .Rdata file to export
#' 
distributedCombat_site <- function(dat, 
                                   batch, 
                                   mod=NULL,
                                   ref.batch=NULL,
                                   central.out=NULL,
                                   eb=TRUE, 
                                   parametric=TRUE,
                                   mean.only=FALSE,
                                   verbose=TRUE, 
                                   file=NULL
){
  if (!is.character(file)) {
    warning("Must specify filename to output results as a file. Currently
            saving output to current workspace only.")
  }
  
  if (is.character(central.out)) {
    load(central.out)
    central.out <- central_out
  }
  
  dat <- as.matrix(dat)
  .checkConstantRows(dat)
  .checkNARows(dat)
  ## Check for missing values
  hasNAs <- any(is.na(dat))
  if (hasNAs & verbose){
    cat(paste0("[neuroCombat] Found ", sum(is.na(dat)), " missing data values. \n"))
  }
  if(mean.only){
    if (verbose) cat("[neuroCombat] Performing ComBat with mean only\n")
  }
  
  ##################### Getting design ############################
  dataDict <- getDataDictDC(batch, mod, verbose=verbose, mean.only=mean.only, ref.batch=ref.batch)
  design <- dataDict[["design"]]
  ####################################################################
  
  
  ############### Site matrices for standardization #################
  # W^T W used in LS estimation
  ls_site <- NULL
  ls_site[[1]] <- crossprod(design)
  ls_site[[2]] <- tcrossprod(t(design), dat)
  
  dataDictOut <- dataDict
  dataDictOut$design <- NULL
  
  # new dataDict with batches within current site
  inclBat <- dataDict$n.batches > 0
  dataDictSite <- dataDict
  dataDictSite$batch <- droplevels(dataDict$batch)
  dataDictSite$batches <- dataDict$batches[inclBat]
  dataDictSite$n.batch <- sum(inclBat)
  dataDictSite$n.batches <- dataDict$n.batches[inclBat]
  dataDictSite$batch.design <- as.matrix(dataDict$batch.design[,inclBat])
  
  # remove reference batch information if reference batch is not in site
  if (!is.null(ref.batch)) {
    if (dataDictSite$ref %in% dataDictSite$batch) {
      dataDictSite$ref <- which(levels(as.factor(dataDictSite$batch))==ref.batch)
    } else {
      dataDictSite$ref <- NULL
      dataDictSite$ref.batch <- NULL
    }
  }
  
  if (is.null(central.out)) {
    site_out <- list(
      ls.site = ls_site,
      dataDict = dataDict,
      sigma.site = NULL
    )
    if (is.character(file)) {
      save(site_out, file = file)
      return(invisible())
    } else {
      return(site_out)
    }
  } 
  
  # If beta.estimates given, get summary statistics for sigma estimation
  
  if (is.null(central.out$var.pooled)) {
    sigma_site <- getSigmaSummary(dat, dataDict, design, hasNAs, central.out)
    
    site_out <- list(
      ls.site = ls_site,
      dataDict = dataDict,
      sigma.site = sigma_site
    )
    if (is.character(file)) {
      save(site_out, file = file)
      return(invisible())
    } else {
      return(site_out)
    }
  }
  
  stdObjects <- getStandardizedDataDC(dat=dat, 
                                      dataDict=dataDict,
                                      design=design,
                                      hasNAs=hasNAs,
                                      central.out=central.out
  )
  s.data <- stdObjects[["s.data"]]
  ####################################################################
  
  
  
  ##################### Getting L/S estimates #######################
  if (verbose) cat("[distributedCombat] Fitting L/S model and finding priors\n")
  naiveEstimators <- getNaiveEstimators(s.data=s.data,
                                        dataDict=dataDictSite, 
                                        hasNAs=hasNAs,
                                        mean.only=mean.only
  )
  ####################################################################
  
  
  ######################### Getting final estimators ####################
  if (eb){
    if (parametric){
      if (verbose) cat("[distributedCombat] Finding parametric adjustments\n")}else{
        if (verbose) cat("[distributedCombat] Finding non-parametric adjustments\n")
      }
    estimators <- getEbEstimators(naiveEstimators=naiveEstimators, 
                                  s.data=s.data, 
                                  dataDict=dataDictSite,
                                  parametric=parametric,
                                  mean.only=mean.only
    )
  } else {
    estimators <- getNonEbEstimators(naiveEstimators=naiveEstimators, dataDict=dataDict)
  }
  ####################################################################
  
  
  
  ######################### Correct data #############################
  if (verbose) cat("[distributedCombat] Adjusting the Data\n")
  bayesdata <- getCorrectedData(dat=dat,
                                s.data=s.data,
                                dataDict=dataDictSite,
                                estimators=estimators,
                                naiveEstimators=naiveEstimators,
                                stdObjects=stdObjects,
                                eb=eb
  )
  ####################################################################
  
  
  # List of estimates:
  estimates <- list(gamma.hat=naiveEstimators[["gamma.hat"]], 
                    delta.hat=naiveEstimators[["delta.hat"]], 
                    gamma.star=estimators[["gamma.star"]],
                    delta.star=estimators[["delta.star"]], 
                    gamma.bar=estimators[["gamma.bar"]], 
                    t2=estimators[["t2"]], 
                    a.prior=estimators[["a.prior"]], 
                    b.prior=estimators[["b.prior"]], 
                    stand.mean=stdObjects[["stand.mean"]], 
                    mod.mean=stdObjects[["mod.mean"]], 
                    var.pooled=stdObjects[["var.pooled"]],
                    beta.hat=stdObjects[["beta.hat"]],
                    mod=mod, 
                    batch=batch, 
                    ref.batch=ref.batch, 
                    eb=eb, 
                    parametric=parametric, 
                    mean.only=mean.only
  )
  
  site_out <- list(dat.combat=bayesdata, estimates=estimates)
  if (is.character(file)) {
    save(site_out, file = file)
    return(invisible())
  } else {
    return(site_out)
  }
}

#' Distributed ComBat step at analysis core
#' 
#' @param site.outs List or vector of filenames containing site outputs.
#' @param file File name of .Rdata file to export
#' @param ref.batch Optional, reference batch used to determine target mean and
#'   variance
#' @param verbose Whether to print messages to console
distributedCombat_central <- function(site.outs,
                                      file = NULL,
                                      ref.batch = NULL,
                                      verbose = FALSE) {
  if (!is.character(file)) {
    warning("Must specify filename to output results as a file. Currently
            saving output to current workspace only.")
  }
  
  if (is.character(site.outs)) {
    fnames <- site.outs
    site.outs <- lapply(fnames, function(file) {
      load(file)
      site_out
    })
  }
  m <- length(site.outs) # number of sites
  
  # get n.batches and n.array from sites
  batch_levels <- levels(site.outs[[1]]$dataDict$batch)
  n.batches <- Reduce("+", lapply(site.outs, function(x) x$dataDict$n.batches))
  n.batch <- length(n.batches)
  n.array <- sum(n.batches)
  n.arrays <- lapply(site.outs, function(x) x$dataDict$n.array)
  
  # # get reference batch if specified
  if (!is.null(ref.batch)){
    if (!(ref.batch%in%levels(batch))) {
      stop("reference level ref.batch is not found in batch")
    }
    if (verbose){
      cat(paste0("[combat] Using batch=",ref.batch, " as a reference batch \n"))
    }
    ref <- which(batch_levels==ref.batch) # find the reference
  } else {
    ref <- NULL
  }

  # check if beta estimates have been given to sites
  step1s <- sapply(site.outs, function(x) is.null(x$sigma.site))
  if (length(unique(step1s)) > 1) {
    stop("Not all sites are at the same step, please confirm with each site.")
  }
  step1 <- all(step1s)
  
  #### Step 1: Get LS estimate across sites ####
  ls1 <- Reduce("+", lapply(site.outs, function(x) x$ls.site[[1]]))
  ls2 <- Reduce("+", lapply(site.outs, function(x) x$ls.site[[2]]))
  B.hat <- crossprod(solve(ls1), ls2)
  
  if (!is.null(ref.batch)) {
    grand.mean <- t(B.hat[ref, ])
  } else {
    grand.mean <- crossprod(n.batches/n.array, B.hat[1:n.batch,])
  }
  stand.mean <- crossprod(grand.mean, t(rep(1,n.array)))
  
  if (step1) {
    central_out <- list(
      B.hat = B.hat,
      stand.mean = stand.mean,
      var.pooled = NULL
    )
    if (is.character(file)) {
      save(central_out, file = file)
      return(invisible())
    } else {
      return(central_out)
    }
  }
  
  # #### Step 2: Get standardization parameters ####
  vars <- lapply(site.outs, function(x) x$sigma.site)
  
  # if ref.batch specified, use estimated variance from reference site
  if (!is.null(ref.batch)){
    var.pooled = vars[[ref]]
  } else {
    var.pooled = rep(0, length(vars[[1]]))
    for (i in 1:m) {
      var.pooled = var.pooled + n.arrays[[i]]*vars[[i]]
    }
    var.pooled = var.pooled/n.array
  }
  
  central_out <- list(
    B.hat = B.hat,
    stand.mean = stand.mean,
    var.pooled = var.pooled
  )
  if (is.character(file)) {
    save(central_out, file = file)
    return(invisible())
  } else {
    return(central_out)
  }
}

# modified to not check design matrix
getDataDictDC <- function(batch, mod, verbose, mean.only, ref.batch=NULL){
  batch <- as.factor(batch)
  n.batch <- nlevels(batch)
  batches <- lapply(levels(batch), function(x)which(batch==x))
  n.batches <- sapply(batches, length)
  n.array  <- sum(n.batches)
  batchmod <- model.matrix(~-1+batch)  
  if (verbose) cat("[combat] Found",nlevels(batch),'batches\n')
  if(any(n.batches==1) & mean.only==FALSE){
    stop("Found one site with only one sample. Consider using the mean.only=TRUE option")
  }
  if (!is.null(ref.batch)){
    if (!(ref.batch%in%levels(batch))) {
      stop("reference level ref.batch is not found in batch")
    }
    if (verbose){
      cat(paste0("[combat] Using batch=",ref.batch, " as a reference batch \n"))
    }
    ref <- which(levels(as.factor(batch))==ref.batch) # find the reference
    batchmod[,ref] <- 1
  } else {
    ref <- NULL
  }
  #combine batch variable and covariates
  design <- cbind(batchmod,mod)
  # check for intercept in covariates, and drop if present
  # check  <- apply(design, 2, function(x) all(x == 1))
  # if(!is.null(ref)){
  #   check[ref] <- FALSE
  # }
  # design <- as.matrix(design[,!check])
  # design <- .checkDesign(design, n.batch)
  n.covariates <- ncol(design)-ncol(batchmod)
  if (verbose) cat("[combat] Adjusting for ",n.covariates,' covariate(s) or covariate level(s)\n')
  out <- list()
  #Making sure to keep track of names:
  names(batches)   <- names(n.batches) <- levels(batch)
  colnames(design) <- gsub("batch", "", colnames(design))
  out[["batch"]] <- batch
  out[["batches"]] <- batches
  out[["n.batch"]] <- n.batch
  out[["n.batches"]] <- n.batches
  out[["n.array"]] <- n.array
  out[["n.covariates"]] <- n.covariates
  out[["design"]] <- design
  out[["batch.design"]] <- design[,1:n.batch]
  out[["ref"]] <- ref
  out[["ref.batch"]] <- ref.batch
  return(out)
}

getSigmaSummary <- function(dat, dataDict, design, hasNAs, central.out){
  batches=dataDict$batches
  n.batches=dataDict$n.batches
  n.array=dataDict$n.array
  n.batch=dataDict$n.batch
  ref.batch=dataDict$ref.batch
  ref=dataDict$ref
  B.hat <- central.out$B.hat
  stand.mean <- central.out$stand.mean[,1:n.array]
  
  if (!hasNAs){
    if (!is.null(ref.batch)){
      ref.dat <- dat[, batches[[ref]]]
      factors <- (n.batches[ref]/(n.batches[ref]-1))
      var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat), na.rm=TRUE)/factors
    } else {
      factors <- (n.array/(n.array-1))
      var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)/factors
    }
  } else {
    if (!is.null(ref.batch)){
      ref.dat <- dat[, batches[[ref]]]  
      ns <- rowSums(!is.na(ref.dat))
      factors <- (ns/(ns-1))
      var.pooled <- rowVars(ref.dat-t(design[batches[[ref]], ]%*%B.hat), na.rm=TRUE)/factors
    } else {
      ns <- rowSums(!is.na(dat))
      factors <- (ns/(ns-1))
      var.pooled <- rowVars(dat-t(design %*% B.hat), na.rm=TRUE)/factors
    }
  }
  
  return(var.pooled)
}

getStandardizedDataDC <- function(dat, dataDict, design, hasNAs, central.out){
  batches=dataDict$batches
  n.batches=dataDict$n.batches
  n.batch=dataDict$n.batch
  n.array=dataDict$n.array
  ref.batch=dataDict$ref.batch
  ref=dataDict$ref
  
  B.hat <- central.out$B.hat
  stand.mean <- central.out$stand.mean[,1:n.array]
  var.pooled <- central.out$var.pooled
  
  if(!is.null(design)){
    tmp <- design
    tmp[,c(1:n.batch)] <- 0
    mod.mean <- t(tmp%*%B.hat)
    #stand.mean <- stand.mean+t(tmp%*%B.hat)
  } else {
    mod.mean <- 0
  }
  s.data <- (dat-stand.mean-mod.mean)/(tcrossprod(sqrt(var.pooled), rep(1,n.array)))
  names(var.pooled) <- rownames(dat)
  rownames(stand.mean) <- rownames(mod.mean) <- rownames(dat)
  colnames(stand.mean) <- colnames(mod.mean) <- colnames(dat)
  return(list(s.data=s.data, 
              stand.mean=stand.mean,
              mod.mean=mod.mean, 
              var.pooled=var.pooled,
              beta.hat=B.hat
  )
  )
}
