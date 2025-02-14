#' Simulate data from parameters of a fitted `flexmix` model
#' @param mod flexmix object
#' @param size size Parameter of a binomial distribution
#' @param N sample size (defaults to sample size from original data)
#
# NOTE TO SELVES: removed the p_noise part because we didn't end up using it (alpha
#is instead introduced when fitting the binomial model on the input data)
#
data_sim_from_model = function(mod, size, N=NULL) {
  prop = prior(mod)
  param = parameters(mod)
  
  k = ncol(param)
  nvar = nrow(param)
  
  if(is.null(N)) {
    N = nrow(posterior(mod))
  }
  
  ni = round(N*prop)
  
  # true clusters
  z = rep(seq(k), ni)
  
  # data
  x = lapply(seq(k), \(i) {
    lapply(seq(nvar), \(var) {
      p = param[var,i]
      rbinom(ni[i], size, p)
    }) %>% do.call(cbind, .)
  }) %>% do.call(rbind, .)
  
  
  shuffle = sample(sum(ni), size=sum(ni), replace=FALSE)
  dat = x[shuffle, ]
  
  res = list(group = z[shuffle], data = dat)
  res
}