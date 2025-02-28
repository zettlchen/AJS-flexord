#Simulation script, last checked: 28.02.25

#open points:

#fix: back to print to file option (adding mkdir step to simulation script) [X]
#general TODO: add scripts with old distGDM2 and distGower [X]
#general TODO: add script for FLXMCbinomial2 [~]
#general TODO: add script for FLXMCbetabinom2 [~]
#are the other FLXs the same as in flexord? [ ]

library(flexclust)
library(flexmix)
list.files('../R', full.names = TRUE) |> 
  sapply(source)
library(flexord)
library(parallel)
library(data.table)
library(here)
library(nnet, include.only = 'multinom')
library(magrittr, include.only = '%>%')
library(clusterSim, include.only = 'GDM2')
library(withr, include.only = 'with_seed')

(num_cores <- detectCores(logical = TRUE)) 
(num_cores_use <- max(1, floor(num_cores * 0.75)))

data('lowbackpain', package='flexord')

alpha2 <- c(0, 75, 150) |> 
  setNames(c('easy', 'medium', 'hard'))
m <- c(3, 6, 11)
N <- c(50, 200, 500)
r <- 2:11
nIter <- 100

sapply('../data', \(dir) {
  if(!dir.exists(dir)) {
    dir.create(dir, recursive=T, showWarnings=F)
  }
})

mclapply(alpha2,
         \(a2) sim_backpain(inputdata = lowbackpain$data,
                            nIter = nIter,
                            alpha2 = a2,
                            n_vars = m,
                            sample_size = N,
                            size=r-1),
         mc.cores = num_cores_use)

