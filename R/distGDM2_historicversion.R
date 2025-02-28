# Code for the original GDM2 distance for arbitrary data objects that was used
# in the simulation study (`flexord::distGDM2` makes use of new capabilities in
# `flexclust` V. 1.5.0, that were not yet implemented at the time of the simulation study)


#' @param x A numeric matrix or data frame. Ordinal variables need to be coded as
#'        `1:length(levels(x[,i]))` in steps of one.
#' @param centers A numeric matrix with the same coding scheme as in `x`,
#'       `ncol(centers)==ncol(x)`, and `nrow(centers)<=nrow(x)`.
#' @param fmly object of class `kccaFamilyGDM2` created by `kccaFamilyGDM2()`.
#'        Default: NULL (if applied within kcca, it will extract the family object
#'        from there.)
distGDM2 <- function(x, centers, fmly=NULL, ...) {
  
  if (ncol(x) != ncol(centers))
    stop(sQuote('x'), ' and ', sQuote('centers'), ' must have the same number of columns')
  z <- matrix(0, nrow=nrow(x), ncol=nrow(centers))
  
  if(is.null(fmly)) {
    whereisfamily <- sapply(sys.frames(), ls)
    whereisfamily <- grep('family', whereisfamily)[1]
    fmly <- get('family', envir=sys.frame(whereisfamily),
                inherits=F) #necessary cause else it finds stats::family
  }
  
  N <- nrow(x)
  fx <- fmly@freqDist$proj.x
  hats <- fmly@freqDist$freq.dists
  fc <- .projectIntofx(centers, hats)
  
  for(k in 1:nrow(centers)) {
    for(i in 1:N) {
      deltaeq <- x[i,]==centers[k,]
      numneq <- (!deltaeq)*(1/N + 2*abs(fx$Ftildex[i,]-fc$Ftildex[k,]))
      numeq <- deltaeq*fx$epdf[i,]
      denom <- sqrt(sum(1-fx$epdf[i,])*sum(1-fc$epdf[k,]))
      
      z[i,k] <- (1-(sum(1-numneq-numeq)/denom))/2
    }
  }
  z 
  
}


#' kccaFamilyGDM2: class and filler function that add an additional slot to `flexclust::kcca`
#' and that conduct the necessary preproc steps for `distGDM2`. Superseded by
#' `flexord::kccaExtendedFamily` and `flexclust` V. 1.5.0
setClass(
  "kccaFamilyGDM2",
  contains = "kccaFamily",
  slots = c(freqDist = "list")  # Additional slot for storing frequency distribution data
)

#' @param x A numeric matrix or data frame. Ordinal variables need to be coded as
#'        `1:length(levels(x[,i]))` in steps of one.
#' @param centers A numeric matrix with the same coding scheme as in `x`,
#'       `ncol(centers)==ncol(x)`, and `nrow(centers)<=nrow(x)`.
#' @param xrange response level range of the variables. Implemented options:
#'       'all': range of x, same for all variables (Default for xrange=NULL)
#'       'columnwise': range of each variable in the data
#'       c(lower, upper): range vector specified by user, upper and lower limit for all variables
#'       list(x1= c(lower, upper), x2=c(lower, upper), ...): list of range vectors specified by user, upper and lower limits are variable specific
kccaFamilyGDM2 <- function(x, cent=NULL, xrange='all') {
  
  freq.projected.x <- .xAsFreqdist(x, xrange=xrange)
  
  baseFamily <- kccaFamily(dist=distGDM2,
                           cent=cent)
  
  extendedFamily <- new("kccaFamilyGDM2",
                        baseFamily,
                        freqDist = freq.projected.x)
  
  return(extendedFamily)
}

.projectIntofx <- function(x, hats) {
  fxs <- sapply(c('epdf', 'ecdf'), \(type) {
    z <- matrix(0, nrow=nrow(x), ncol=ncol(x),
                dimnames=dimnames(x))
    for(j in 1:ncol(x)) {
      for(i in 1:nrow(x)) {
        z[i,j] <- hats[[j]][[type]](x[i,j])
      }
    }
    z
  }, simplify = F)
  fxs$Ftildex <- fxs$ecdf - (fxs$epdf/2)
  return(fxs)
}

.xAsFreqdist <- function(x, xrange){
  
  if(sum(xrange=='all')) {
    rng <- rep(range(x, na.rm=T), ncol(x)) |>
      matrix(nrow=2)
  } else if(sum(xrange=='columnwise')) {
    rng <- apply(x, 2, range, na.rm=T)
  } else if(is.vector(xrange, mode='numeric')) {
    rng <- rep(xrange, ncol(x)) |> 
      matrix(nrow=2)
  } else {
    if(length(xrange) != ncol(x))
      stop('Either supply 1 range vector, or list of ranges for all variables')
    rng <- unlist(xrange) |> 
      matrix(nrow=2)
  }
  
  hats <- lapply(1:ncol(x), \(y) {
    level <- factor(x[,y], levels=seq(rng[1,y], rng[2,y]))
    pdf <- table(level)/nrow(x)# |> 
    pdf <- as.data.frame.table(pdf)
    pdf$level <- as.numeric(pdf$level)
    epdf <- function(i) {
      ind <- which(pdf[,'level']<=i)
      if(length(ind)==0) {
        return(0)
      } else {
        return(pdf$Freq[ind[length(ind)]])
      }
    }
    list(epdf=epdf, ecdf=ecdf(x[,y]))
  })
  names(hats) <- colnames(x) #pdf tables remain accesible via get('pdf', environment(hats[[1]]$epdf))
  
  fxs <- .projectIntofx(x, hats)
  
  return(list(freq.dists=hats,
              proj.x=fxs))
}