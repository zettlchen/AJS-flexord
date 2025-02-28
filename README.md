# AJS-flexord

## About 
In this repository, we provide reproduction code for the simulation study published
in *Ordinal Clustering with the flex-scheme* [@ernst_ordinal_2025]. An introduction
to the `flexord` package can be found [here](Intro2Flexord.html).

In the study, we have compared the performance of 12 different algorithms against
the *true clustermembership* via de Adjusted Rand Index [@hubert_arabie_1985]
on simulated ordinal data sets with equal response level lengths and no missing values.

We varied the following aspects of the data sets:

* the sample size $N \in \{50,200,500\}$;
* the number of response levels of each ordinal variable
  $r \in \{2,\ldots,11\}$;
* the number of variables $m \in \{3, 6, 11\}$; and
* the regularization parameter $\alpha \in \{0, 75, 150\}$, where
  more regularization results in more diffuse clusters.
  
For each configuration combination, we created `nIter=100` data sets, on
which we ran the algorithms.

To simulate our `nIter` data sets, we

1) took the binary data on presence/absence of 11 symptoms of low backpain  in 464 patients and 
 their respective diagnoses provided by @fop_variable_2017 as input data,
2) on which we fit finite mixtures of multivariate independent Bernoulli distributions with
three^[The value of `k=3` is taken from the original study by @fop_variable_2017.]
components (regularized by `alpha`, in order to obtain moderately well separated clusters).
3) From these fitted models, we then generated data by drawing `N` times using binomial distributions
for components, where the number of trials is set to `r`.
4) Finally, we selected the `m` variables of the simulated data set with the highest 'variable importance'.
We calculated this 'variable importance' by fitting a mixture of three multinomial regressions
explaining the expert diagnoses (and thus *true cluster memberships*) separately by each of the 11
original variables, and ordering them by decreasing Log-Likelihood.

## Contents of this repository

- `scripts`: Scripts for the simulation and visualizations published in the paper
- `R`: R code for additional functions
- `session_info.txt` and `system_info.txt`: Information of the R session, package versions
    and system architecture on which the simulation was run.

The main simulation functions (in `R/`) are:
`data_sim_from_model()`
: takes a flexmix model that has been fitted onto the input data while using the
regulation parameter `alpha`, and generates data with desired `N` and `r`,  from
it by drawing from a binomial distribution using the model priors and parameters.
Furthermore, it returns the model partitions as *true clusters* of the simulated
data sets.

`sim_backpain_apply2grid()`
: applies `data_sim_from_model()` and the tested clustering algorithms to the
parameter grid (`N`, `r`, `m`, and `alpha`). Plus some error catching for cases
where combinations of simulated data set and clustering algorithm cannot converge.

`sim_var_importance()`
: Computes the 'variable importance' of each variable in the input data set in
order to select the `m` variables of highest importance (=most information on
cluster structure). This is done by fitting independent multinomial models to
each variable explaining the *true clustering*, and then sorting variables by
decreasing log-Likelihood.

`step___4sim()`
: replicate `flexclust::stepFlexclust`/`flexmix::stepFlexmix` with additional
error catching, parallelization customization, and reformatting to data.tables.

## Required packages and system architecture

The packages required for running this simulation (and plotting its results) are:

-methods
-cluster
-clusterSim
-withr
-data.table
-flexclust
-flexmix
-flexord
-ggplot2
-here
-magrittr
-nnet
-parallel

Detailed information on the versions used can be found in `session_info.txt`.

However, we do want to point out, that the simulation was written depending on `flexclust` version 1.4.2.
While it also runs with `flexclust` 1.5.0, the new capabilities implemented would allow for more elegant
options in some of the partitioning cluster algorithms. While we have implemented these new versions in
our package `flexord`, this changes are not implemented in this simulation. For this reason, we also provide
to scripts on distance measures that show the 'historical' code as it was used in the simulation, marked
by `_historicversion.R`.