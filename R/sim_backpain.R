#copied from SimulationStudy.Rnw
#fix: back to print to file option (adding mkdir step to simulation script)
#general TODO: add scripts with old distGDM2 and distGower
#general TODO: add script for FLXMCbinomial2 (it isn't actually the same, at least alpha2 is
#row count, not proportion)

#' 'The big simulation wrapper' where we apply the simulation and the clustering
#' algorithms
#' @param inputdata data set that will be used as the basis for simulating data
#'        sets
#' @param nIter number of data sets to create for each combination of data
#'        characteristics
#' @param alpha2 goal simulated data characteristic:
#'        numeric vector of regularization parameter to be used in data set
#'        simulations. Possible values: 0-`nrow(sample_size)`, i.e. the number of
#'        observations in the simulated data set that shall be replaced by the mean value
#' @param n_vars goal simulated data characteristic:
#'        numeric vector of the numbers of variables in the simulated data sets.
#'        (The maximum value here is `ncol(inputdata)`).
#' @param sample_size goal simulated data characteristic:
#'        numeric vector of the numbers of observations in the simulated data
#'        sets. (The maximum value here is `nrow(inputdata)`).
#' @param size goal simulated data characteristic:
#'        numeric vector of the response level lengths of the ordinal variables
#'        in the simulated data sets (equal across all variables of one data set).
#
#NOTE TO SELVES: removed the `algos` argument
sim_backpain <- function(inputdata,
                         nIter, alpha2,
                         n_vars, sample_size,
                         size) {
  with_seed(1234, {
    
    var_imp = sim_var_importance(inputdata=inputdata$data)
    
    # set default parameters
    simfn = \(...) {
      sim_backpain_apply2grid(sizes=size,
                              sample_size=sample_size,
                              n_vars=n_vars,
                              var_imp=var_imp,
                              Niter=nIter, ...)
    }
    
    # running the base model for simulation (within seed,
    # so the resulting data sets will be the same)
    mod1 = stepFLX4sim(data=inputdata$data,
                       model=FLXMCbinomial(alpha2=alpha2), k=3, #is this the same as ajs.ordinal.clustering::FLXMCbinomial2?
                       nrep=10,
                       verbose=FALSE, multicore=FALSE)
  })
  
  res <- list()
  
  # clustering steps (+application of simulation)
  #binomial model
  res[[1]] <- sim_res_binom = simfn(mod1, \(data, k, size) {
    stepFLX4sim(data=data, model=FLXMCbinomial(size=size), k=k,
                nrep=10,
                verbose=FALSE, multicore=FALSE)
  }) %>% .[, algo := "binomial"]
  
  #regularized multivariate normal model
  res[[2]] <- sim_res_regnorm = simfn(mod1, \(data, k, size) {
    stepFLX4sim(data=data, model=FLXMCregnorm(G=k), k=k,
                nrep=10,
                verbose=FALSE, multicore=FALSE)
  }) %>% .[, algo := "regnorm"]
  
  # kmeans model
  res[[3]] <- sim_res_kmeans = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE)
  }) %>% .[, algo := "kmeans"]
  
  # betabinomial model
  res[[4]] <- sim_res_betabinom = simfn(mod1, \(data, k, size) {
    stepFLX4sim(data=data, model=FLXMCbetabinom2(size=size), k=k,
                nrep=10,
                verbose=FALSE, multicore=FALSE)
  }) %>% .[, algo := "betabinom"]
  
  # multinomial model
  res[[5]] <- sim_res_multinom = simfn(mod1, \(data, k, size) {
    stepFLX4sim(data=data+1, model=FLXMCregmultinom(size=size+1, alpha2=1), k=k,
                nrep=10,
                verbose=FALSE, multicore=FALSE)
  }) %>% .[, algo := "multinom"]
  
  # kmedians modeL
  res[[6]] <- sim_res_kmedians = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE,
               family=kccaFamily('kmedians'))
  }) %>% .[, algo := "kmedians"]
  
  # kmodes model
  res[[7]] <- sim_res_kmodes = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE,
               family=kccaFamily(dist=distSimMatch,
                                 cent=centMode))
  }) %>% .[, algo := "kmodes"]
  
  # PAM:Gower model
  res[[8]] <- sim_res_PAMGower = simfn(mod1, \(data, k, size) {
    stepDaisy4sim(data, k=k, ordinal = T)
  }) %>% .[, algo := "PAM:Gower"]
  
  # PAM:GDM2 model
  res[[9]] <- sim_res_PAMGDM2 = simfn(mod1, \(data, k, size) {
    stepDaisy4sim(data, k=k, diss = clusterSim::GDM2)
  }) %>% .[, algo := "PAM:GDM2"]
  
  # kgower model
  res[[10]] <- sim_res_kgower = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE,
               family=kccaFamily(dist=\(y, centers) distGower.ordinal(y,
                                                                      centers,
                                                                      range=range(data)),
                                 cent=\(y) centMin(y, distGower.ordinal,
                                                   range=range(data))))
  }) %>% .[, algo := "kgower"]
  
  # kGower model with medians as centroids
  res[[11]] <- sim_res_kgowerMed = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE,
               family=kccaFamily(dist=\(y, centers) distGower.ordinal(y,
                                                                      centers,
                                                                      range=range(data)),
                                 cent=centMedian))
  }) %>% .[, algo := "kgowerMed"]
  
  # kGDM2 model
  res[[12]] <- sim_res_kGDM2 = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE,
               family=kccaFamilyGDM2(data, xrange=range(data))) #picks centOptim as default
  }) %>% .[, algo := "kGDM2centOpt"]
  
  # kGDM2 model with modes as centroids  
  res[[13]] <- sim_res_kGDM2_1 = simfn(mod1, \(data, k, size) {
    stepKM(data, k=k, nrep=10, multicore=FALSE,
           family=kccaFamilyGDM2(data, cent=centMode, xrange=range(data)))
  }) %>% .[, algo := "kGDM2centMode"]
  
  bind_rows(res)    
}