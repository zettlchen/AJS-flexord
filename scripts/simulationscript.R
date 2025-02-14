#copied from ajs-ordinal-clustering/inst/scripts/sim4nonas_v2.R
#cleanup required - change saving space, and function names to new structure
#furthermore: text stuff from .Rmd vignettes needs to move to the README

list.files('../R') |> sapply(source) #gonna test whether that works in a different .Rproj

library(flexclust)
library(flexord)
#do we not need flexmix??
library(magrittr) #do we need anything else from there?
library(parallel)
library(data.table)

backpain <- data('lowbackpain', package='flexord')$data

alpha2 <- c(0, 75, 150) |> 
  setNames(c('easy', 'medium', 'hard'))
m <- c(3, 6, 11)
N <- c(50, 200, 500)
r <- 2:11 #do I need to do the -1 or not?

#removed algos argument

sapply('../data', \(dir) { #again, test in new .Rproj
  if(!dir.exists(dir)) {
    dir.create(dir, recursive=T, showWarnings=F)
  }
})

mclapply(alpha2,
         \(a2) sim_backpain(inputdata = lowbackpain$data,
                            nIter = 100,
                            alpha2 = a2,
                            n_vars = m,
                            sample_size = N,
                            size=r),
         mc.cores = 4)

