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
#'        
sim_backpain <- function(inputdata,
                         nIter, alpha2,
                         n_vars, sample_size,
                         size) {
  withr::with_seed(1234, {
    
    var_imp = sim_var_importance(inputdata=inputdata)
    n_debug = NA
    
    # set default parameters
    simfn = \(...) {
      sim_backpain_apply2grid(sizes=size,
                              sample_size=sample_size,
                              n_vars=n_vars,
                              var_imp=var_imp,
                              Niter=nIter,
                              n_debug=NA, ...)
    }
    
    # running the base model for simulation (within seed,
    # so the resulting data sets will be the same)
    mod1 = stepFLX4sim(data=inputdata,
                       model=FLXMCbinomial2(alpha2=alpha2), k=3,
                       nrep=10,
                       verbose=FALSE, multicore=FALSE)
  })
  
  # clustering steps (+application of simulation)
  # binomial model
  sim_res_binom = simfn(mod1, \(data, k, size) {
    stepFLX4sim(data=data, model=FLXMCbinomial(size=size), k=k,
                nrep=10,
                verbose=FALSE, multicore=FALSE)
  }) %>% .[, algo := "binomial"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_binomial_alpha%d.rds", alpha2) |> 
    saveRDS(sim_res_binom, file=_)
  
  # multivariate normal model
  sim_res_regnorm = simfn(mod1, \(data, k, size) {
    stepFLX4sim(data=data, model=FLXMCregnorm(G=k), k=k,
                nrep=10,
                verbose=FALSE, multicore=FALSE)
  }) %>% .[, algo := "regnorm"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_regnorm_alpha%d.rds", alpha2)  |> 
    saveRDS(sim_res_regnorm, file=_)
  
  # kmeans model
  sim_res_kmeans = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE)
  }) %>% .[, algo := "kmeans"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_kmeans_alpha%d.rds", alpha2) |> 
    saveRDS(sim_res_kmeans, file=_)
  
  # betabinomial model
  sim_res_betabinom = simfn(mod1, \(data, k, size) {
    stepFLX4sim(data=data, model=FLXMCbetabinom2(size=size), k=k,
                nrep=10,
                verbose=FALSE, multicore=FALSE)
  }) %>% .[, algo := "betabinom"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_betabinom_alpha%d.rds", alpha2) |> 
    saveRDS(sim_res_betabinom, file=_)
  
  # multinomial model
  sim_res_multinom = simfn(mod1, \(data, k, size) {
    stepFLX4sim(data=data+1, model=FLXMCregmultinom(size=size+1, alpha2=1), k=k,
                nrep=10,
                verbose=FALSE, multicore=FALSE)
  }) %>% .[, algo := "multinom"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_multinom_alpha%d.rds", alpha2) |> 
    saveRDS(sim_res_multinom, file=_)
  
  # kmedians modeL
  sim_res_kmedians = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE,
               family=kccaFamily('kmedians'))
  }) %>% .[, algo := "kmedians"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_kmedians_alpha%d.rds", alpha2) |> 
    saveRDS(sim_res_kmedians, file=_)
  
  # kmodes model
  sim_res_kmodes = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE,
               family=kccaFamily(dist=distSimMatch,
                                 cent=centMode))
  }) %>% .[, algo := "kmodes"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_kmodes_alpha%d.rds", alpha2) |> 
    saveRDS(sim_res_kmodes, file=_)
  
  # PAM:Gower model
  sim_res_PAMGower = simfn(mod1, \(data, k, size) {
    stepDaisy4sim(data, k=k, ordinal = T)
  }) %>% .[, algo := "PAM:Gower"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_PAMGower_alpha%d.rds", alpha2) |> 
    saveRDS(sim_res_PAMGower, file=_)
  
  # PAM:GDM2 model
  sim_res_PAMGDM2 = simfn(mod1, \(data, k, size) {
    stepDaisy4sim(data, k=k, diss = clusterSim::GDM2)
  }) %>% .[, algo := "PAM:GDM2"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_PAMGDM2_alpha%d.rds", alpha2) |> 
    saveRDS(sim_res_PAMGDM2, file=_)
  
  # kgower model
  sim_res_kgower = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE,
               family=kccaFamily(dist=\(y, centers) distGower.ordinal(y,
                                                                      centers,
                                                                      range=range(data)),
                                 cent=\(y) centMin(y, distGower.ordinal,
                                                   xrange=range(data))))
  }) %>% .[, algo := "kgower"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_kgower_alpha%d.rds", alpha2) |> 
    saveRDS(sim_res_kgower, file=_)
  
  # kGower model with medians as centroids
  sim_res_kgowerMed = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE,
               family=kccaFamily(dist=\(y, centers) distGower.ordinal(y,
                                                                      centers,
                                                                      range=range(data)),
                                 cent=centMedian))
  }) %>% .[, algo := "kgowerMed"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_kgowerMed_alpha%d.rds", alpha2) |> 
    saveRDS(sim_res_kgowerMed, file=_)
  
  # kGDM2 model
  sim_res_kGDM2 = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE,
               family=kccaFamilyGDM2(data, xrange=range(data)))
  }) %>% .[, algo := "kGDM2centOpt"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_kGDM2centOpt_alpha%d.rds", alpha2) |> 
    saveRDS(sim_res_kGDM2, file=_)
  
  # kGDM2 model with modes as centroids  
  sim_res_kGDM2_1 = simfn(mod1, \(data, k, size) {
    stepKM4sim(data, k=k, nrep=10, multicore=FALSE,
           family=kccaFamilyGDM2(data, cent=centMode, xrange=range(data)))
  }) %>% .[, algo := "kGDM2centMode"]
  
  sprintf("../data/sim_backpain_sample_size_ncat_kGDM2centMode_alpha%d.rds", alpha2) |> 
    saveRDS(sim_res_kGDM2_1, file=_)
  
}