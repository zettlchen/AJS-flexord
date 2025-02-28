# Code for the original Gower's distance that was used in the simulation study
# (`flexord::distGower` makes use of new capabilities in `flexclust` V. 1.5.0, that were
# not yet implemented at the time of the simulation study)

#' Gower's distance between arbitrary data objects for purely ordinal data without
#' missing values
#' @param x A numeric matrix or data frame. Ordinal variables need to be coded as
#'        `1:length(levels(x[,i]))` in steps of one.
#' @param centers A numeric matrix with the same coding scheme as in `x`,
#'       `ncol(centers)==ncol(x)`, and `nrow(centers)<=nrow(x)`.
#' @param xrange response level range of the variables. Implemented options:
#'       'all': range of x, same for all variables (Default for xrange=NULL)
#'       'columnwise': range of each variable in the data
#'       c(lower, upper): range vector specified by user, upper and lower limit for all variables
#'       list(x1= c(lower, upper), x2=c(lower, upper), ...): list of range vectors specified by user, upper and lower limits are variable specific

distGower.ordinal <- function(x, centers, xrange=NULL) {
  
  if (ncol(x) != ncol(centers))
    stop(sQuote('x'), ' and ', sQuote('centers'), ' must have the same number of columns')
  z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
  
  if(is.null(xrange)) xrange <- 'all'
  if(sum(xrange=='all')) {
    rng <- rep(range(x), ncol(x)) |>
      matrix(nrow=2)
  } else if(sum(xrange=='columnwise')) {
    rng <- apply(x, 2, range)
  } else if(is.vector(xrange, mode='numeric')) {
    rng <- rep(xrange, ncol(x)) |> 
      matrix(nrow=2)
  } else {
    if(length(xrange) != ncol(x))
      stop('Either supply 1 range vector, or list of ranges for all variables')
    rng <- unlist(xrange) |> 
      matrix(nrow=2)
  }
  
  scl <- apply(rng, 2, diff)
  scl <- ifelse(scl==0, 1, scl)
  
  xr <- scale(x, center=rng[1,], scale=scl)
  centr <- scale(centers, center=rng[1,], scale=scl)
  
  for(k in 1:nrow(centers)) {
    z[, k] <- colMeans(abs(t(xr) - centr[k, ]))
  }
  z  
}