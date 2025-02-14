#Contains code for:
#stepFLX4sim, stepKM4sim, stepDaisy4sim
#These functions recreate stepFlexmix and stepFlexclust (and in the last case,
#just adapt cluster::pam output to the same format as the other two), but with
#additional error catching, parallelization control, and reformatting the output
#to data.tables

#' Recreating stepFlexmix
#'@param data input data
#'@param k number of clusters
#'@param model object of class FLXM or list of FLXM objects, which drivers to be
#'       used in the mixture?
#'@param nrep how often to restart flexmix (solution with max. likelihood is returned)
#'@param verbose show progress information?
#'@param multicore use multiple cores?
#'@param set_seed when using multiple cores?
#'
stepFLX4sim <- function(data, k, model, nrep=4L, verbose=FALSE, multicore=FALSE, set_seed=TRUE, ...) {
  fn = \(k) {
    flexmix(data~1, k=k, model=model, ...)
  }
  
  .stepGeneric(fn, data, k=k, nrep=nrep, multicore=multicore, max_by=logLik,
               verbose=verbose, set_seed=set_seed)
}

#' Recreating stepFlexclust
#'@param x input data
#'@param k number of clusters
#'@param FUN which function to use for clustering? Default: flexclust::kcca
#'      (alternatively: flexclust::cclust)
#'@param nrep how often to restart flexclust (solution with min. within cluster
#'       distance is returned)
#'@param multicore use multiple cores?
#'@param ... other parameters. of specific interest:
#'       family (enter object of class 'kccaFamily' here to influence distance
#'       and centroid functions used)
stepKM4sim <- function(x, k, FUN=kcca, nrep=10, multicore=FALSE, ...) {
  fn = \(k) {
    FUN(x, k=k, ...)
  }
  
  max_by = \(mod) {
    -info(mod, "distsum")
  }
  
  .stepGeneric(fn, data, k=k, nrep=nrep, multicore=multicore, max_by=max_by)
}

#' Fitting cluster::daisy output into the same format as from the above functions
#'@param x input data
#'@param k number of clusters
#'@param diss either a dissimilarity matrix, or a distance function to calculate
#'       one. (Optional. If NULL, daisy(metric='Gower') will be used.)
#'@param ordinal logical. Is data ordinal? (x is expected to be coded numerically,
#'       thus, this information is not available in the data set.)
stepDaisy4sim <- function(x, k, diss=NULL, ordinal=TRUE) {
  requireNamespace("cluster")
  
  dissmat.x = .dissmat_pam(x, diss, ordinal)
  
  stepfun = \(k) {
    mod = cluster::pam(dissmat.x, diss=TRUE, k=k)
    rownames_med = rownames(x)[mod$id.med] %>%
      sub("\\..+$", "", .)
    mod$rownames.med = rownames_med
    mod
  }
  
  .stepGeneric(stepfun, x, k, nrep=1, multicore=FALSE, \(x) 42, clusterfun=.clusters_pam)
  
}