#For the documentation of this - what is it, or better, what's the difference to FLXMCbinomial?

#' @export
FLXMCbinomial2 = function(formula=.~., size=NULL, alpha2=0, eps=1e-5)
{
  #size <- NULL
  
  z <- new("FLXMC", weighted=TRUE, formula=formula,
           name="binomial")
  
  z@preproc.y <- function(y) {
    if (any(y < 0, na.rm=TRUE))
      stop("negative values are not allowed for the binomial family")
    
    if(is.null(size)) {
      size <<- apply(y, 2, max, na.rm=TRUE)
    }
    
    y
  }
  
  defineComponent <- function(probs, df) {
    predict <- function(x, ...) {
      stop("not implemented")
    }
    
    logLik <- function(x, y) {
      ty <- t(y)
      bc <- log(choose(size, ty))
      
      l <- bc + ty*log(probs) + (size-ty)*log(1-probs)
      cs <- colSums(l, na.rm=TRUE)
      if(any(!is.finite(cs))) stop("non-finite values while calculating log-likelihood")
      cs
    }
    
    new("FLXcomponent",
        parameters=list(coef=probs),
        logLik=logLik, predict=predict,
        df=df)
  }
  
  z@fit <- function(x, y, w) {
    has_na <- any(is.na(y))
    which_na <- which(is.na(y))
    
    ymarg <- colMeans(y, na.rm=TRUE)/size
    
    b_alpha <- ymarg*alpha2
    b_beta <- (1-ymarg)*alpha2
    
    if(has_na) {
      w <- matrix(w, nrow=length(w), ncol=ncol(y))
      w[which_na] <- NA
      p <- (b_alpha + colSums(w*y, na.rm=TRUE)) / (b_alpha+b_beta+size*colSums(w, na.rm=TRUE))
    } else {
      p <- (b_alpha + colSums(w*y, na.rm=TRUE)) / (b_alpha+b_beta+size*sum(w, na.rm=TRUE))
    }
    
    p <- ifelse(p >= 1-eps, 1-eps, ifelse(p <= eps, eps, p))
    defineComponent(p, ncol(y))
  }
  
  z
}