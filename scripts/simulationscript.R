#copied from ajs-ordinal-clustering/inst/scripts/sim4nonas_v2.R
#cleanup required - change saving space, and function names to new structure
#furthermore: text stuff from .Rmd vignettes needs to move to the README

r <- unclass(lsf.str(envir = asNamespace("ajs.ordinal.clustering"), all = T))
for(name in r) eval(parse(text=paste0(name, '<-ajs.ordinal.clustering:::', name)))

library(flexclust)
library(tidyverse)
library(parallel)
library(data.table)

#add data loading part here

alpha2 <- c(0, 75, 150) |> 
  setNames(c('easy', 'medium', 'hard'))

algs <- c('regnormal', 'multinomial', 'kmedians',
          'kmodes', 'PAM:Gower', 'PAM:GDM2') #leaving the kgs still out as WIP

#add makedir thing here

mclapply(alpha2,
         \(a2) sim_backpain(a2, algos = algs),
         mc.cores = 4)

