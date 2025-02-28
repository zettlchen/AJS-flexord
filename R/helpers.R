# helper functions for the simulation study

#' Apply simulation and clustering to specified parameters.
#' 
#' Conducts the simulation step for the specified parameters (`alpha`, `N`, `r`,
#' and `m`), shapes the output into a grid, and applies the desired stepwise
#' clustering function to the grid.
#' 
#'@param mod flexmix model fitted onto the input data, input model from which
#'       simulated data is generated.
#'@param stepfn desired stepwise clustering function to be applied to the
#'       simulated data set
#'@param sizes the number of response levels in the simulated data sets (equal for
#'       all variables), defined as `r` in the paper
#'@param sample_size the number of observations in the simulated data sets,
#'       defined as `N` in the paper
#'@param n_vars the number of variables in the simulated data sets, defined as `m`
#'       in the paper
#'@param Niter number of simulation+clustering runs. 100 in the original study
#'@param var_imp the function with which to calculate 'variable importance' for
#'       selecting the `m` variables. Default is `sim_var_importance()`
#'@param seed seed on which to create the 100 data sets each (only on that right??
#'       Else all 27,000 results would be equal for each algo right??)
#'@param multicore use several cores? Default=TRUE
#
#NOTE TO SELVES: in the orig code there was an n_debug option, but I think it is unnecessary
#here, as we didn't use it in the final sim function (was always n_debug=NA, and
#the condition was if(!is.na(n_debug)))
sim_backpain_apply2grid = function(mod, stepfn,
                                   sizes, sample_size,
                                   n_vars, Niter, n_debug,
                                   var_imp = sim_var_importance(),
                                   seed = 0xBA0BAB, multicore=TRUE) {
  
  full_sample_size <- c(50, 100, 200, 300, 400, 500)
  full_n_vars <- 2:11
  #NOTE TO SELVES: leaving in the original params here (i.e. the 'original original' ones, before
  #we selected 3 Ns and 3 ms), cause in the code we filtered that afterward to
  #maintain the same seeds. And I guess for a 'true reproduction' that is necessary
  #here as well if we wanted other people to get the same results on the same seed
  #data sets
  
  with_seed(seed, {
    sim1 = \(size, N, m, iter, seed, index) {
      set.seed(seed)
      progr(index)
      sim_data = data_sim_from_model(mod, size, N)
      data = sim_data$data
      group = sim_data$group
      
      which_vars = head(var_imp, m)
      data = data[, which_vars]
      checksum = digest::sha1(sim_data$data)
      
      m1 = try(stepfn(data, k=3, size))
      
      if(!is(m1, "try-error")) {
        if(any(is.na(m1$clusters[[1]]))) {
          ari = NA
        } else {
          ari = randIndex(m1$clusters[[1]], group)
        }
        data.table(size, N, n_vars=m, iter, ari, checksum=checksum)
      } else {
        rr = \(debug=FALSE) {
          if(debug) debugonce(sim1)
          sim1(size, N, iter, seed, index)
        }
        data.table(size, N, n_vars=m, iter, error=m1, retry=rr, seed, checksum=checksum)
      }
    }
    if(is.data.table(mod)) {
      mod = mod$mod[[1]]
    }
    
    force(var_imp)
    
    grid = expand.grid(size=sizes, N=full_sample_size, m=full_n_vars, iter = seq(Niter)) %>%
      as.data.table %>%
      .[, seed := .get_seeds(.N)] %>%
      .[sample(nrow(.), size=nrow(.), replace=FALSE),] %>%
      .[, index := seq(.N)]
    
    grid = grid[N %in% sample_size & m %in% n_vars]
    
    if(!is.na(n_debug)) {
      grid = head(grid, n_debug)
    }
    
    grid[, index := seq(.N)]
    
    applyfn = if(multicore) mclapply else lapply
    
    progr = .progress_eta(nrow(grid))
    
    res = applyfn(seq(nrow(grid)), \(i) {
      do.call(sim1, as.list(grid[i,]))
    }) %>% rbindlist(fill=TRUE)
    print(Sys.time())
    
    res
  })
}

#helper for progress display on console
.progress_eta = function(ntotal) {
  start = Sys.time()
  
  function(i) {
    rest = ntotal-i
    
    now = Sys.time()
    eta = start + (now - start)/i * ntotal
    
    cat(sprintf("%d\t%d\t%s\n", i, ntotal,
                strftime(eta, "%Y-%m-%d %H:%M:%S")))
  }
}

#Conduct the 'stepping-step' for the different algorithm types
#
#recreates 'stepping'-step provided in stepFlexmix/stepFlexclust, but with
#error catching, parallelization and reformatting the output to a data.table.
.stepGeneric = function(fn, data, k=2:10, nrep=10, multicore=TRUE,
                       clusterfun = clusters,
                       max_by,
                       verbose=FALSE,
                       set_seed=multicore) {
  
  progr = .progress_eta(length(k)*nrep)
  
  step1 = \(ki, iter, seed) {
    if(verbose) { cat(sprintf("k=%d\t", ki)); progr(iter) }
    if(set_seed) set.seed(seed)
    
    tt = system.time({
      res = try(fn(ki))
    })
    
    if(!is(res, "try-error")) {
      newcls = clusterfun(res)
      data.table(k=ki, time=tt["elapsed"], nrep=nrep[1],
                 clusters=list(newcls),
                 mod = list(res))
    } else {
      data.table(k=ki, time=tt["elapsed"], nrep=nrep[1],
                 clusters=list(NA),
                 mod = list(NA), error=res)
      
    }
  }
  
  map = if(multicore) mcmapply else mapply
  grid = expand.grid(k=k, rep=seq.int(nrep)) %>%
    as.data.table %>%
    .[order(k, rep)] %>%
    .[, iter := seq(.N)] %>%
    .[, seeds := .get_seeds(.N)]
  
  steps_each = map(\(k, iter, seed) {
    step1(k, iter, seed)
  }, grid$k, grid$iter, grid$seeds, SIMPLIFY=FALSE) %>%
    rbindlist(fill=TRUE)
  
  nrep1 = nrep
  steps_each[, time := sum(time), by=.(k)]
  steps_each[, nrep := nrep1]
  steps_each[, n_success := sum(!is.na(mod)), by=.(k)]
  steps_each[, crit := NA_real_]
  steps_each[!is.na(mod), crit := sapply(mod, max_by)]
  #steps_each[, max_crit := crit == max(crit), by=.(k)]
  
  steps = steps_each[order(k, -crit)] %>%
    unique(by="k") %>%
    .[, -c("crit")]
  
  steps
}

#Helper to increase flexibility in daisy
#
#If user introduces a distance matrix, it is returned.
#If user introduces a distance function to be used, it is applied on x.
#If user does not specify a `diss` argument, `cluster::daisy(metric='Gower')` is
#used as default. Furthermore, if `ordinal=TRUE`, x is transformed to a data.frame
#with all ordinal variables first. (This is done because most of the clustering
#algorithms applied here expect `x` to be a numeric matrix, and information on 
#the coding of variables is thus lost. Here, the ordinal coding is reinstated for
#use in algorithms that rely on aptly coded variables, such as `cluster::daisy`)
.dissmat_pam = function(x, diss, ordinal) {
  
  if(data.class(diss) %in% c('dist', 'matrix')) {
    dissmat = diss
  } else {
    
    if(is.function(diss)) {
      dissfun <- \(data, ...) diss(data, ...)
    } else {
      dissfun <- \(data, ...){
        if(ordinal) {
          cluster::daisy(data, metric='gower', ...)
        } else {cluster::daisy(data, ...)}
      }
    }
    
    if(ordinal) {
      x = do.call(data.frame, lapply(as.data.frame(x), as.ordered))
    }
    dissmat = dissfun(x)
  }    
  
  dissmat
}

#Returns the clusters in a pam object, OR clusters a new data set by assigning
#observations to the nearest medoid of the pam object (equivalent to
#`flexclust::clusters` for kcca objects)
.clusters_pam = function(mod, newdata=NULL, dissmat) {
  if(is.null(newdata)) {
    return(mod$clustering)
  }
  
  as.matrix(dissmat)[,mod$rownames.med] |> 
    apply(1, which.min)
}

#generate different seeds for different starting values
.get_seeds = function(n) {
  round(2^31 * runif(n, -1, 1))
}