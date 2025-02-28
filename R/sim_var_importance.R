#' Compute 'variable importance' of the variables in a given dataset based on decreasing log-Likelihood
#' 
#' Dataset is a simulated binary dataset based on backpain data.
#' For each variable individually a multinomial regression model is computed
#' where the dependent variable is the true clustering.
#' Variables are then sorted by logLik decreasing (1st = highest importance)
#' @param inputdata binary data set that will be 'seed data set' for simulating
#'        data sets. In the paper, it was run on `data('lowbackpain')`.
#' @importFrom nnet multinom
# 
sim_var_importance = function(inputdata) {
  
  mod1 = stepFLX4sim(data=inputdata, model=FLXMCbinomial2(alpha2=0), k=3,
                  nrep=10,
                  verbose=FALSE, multicore=TRUE)
  
  sim_data = data_sim_from_model(mod1$mod[[1]], 1, 500)
  
  cls = as.factor(sim_data$group)
  mods = lapply(seq(ncol(sim_data$data)), \(i) {
    sink("/dev/null") # wtf rly?
    on.exit(sink(NULL))
    multinom(cls ~ sim_data$data[,i])
  })
  
  
  logliks = sapply(mods, logLik)
  ord = order(logliks, decreasing=TRUE)
  
  test = sort(logliks, decreasing=TRUE) == logliks[order(logliks, decreasing=TRUE)]
  stopifnot(all(test))
  
  ord
}