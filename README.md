# Overview
This is the code repository for "A General Machine Learning Framework for Survival Analysis" published at ECML 2020

The analyses in the publication were based on a prototype implementation of
Piece-wise-exponential models (PEMs) using XGBoost in R. This prototyp implementation
is available as an R package `pem.xgb` and must be installed to run the
benchmarks. The package is available from a [separate repository](https://github.com/adibender/pem.xgb). If you are only interested in how to make XGBoost
estimate PEMs head over to that repository. Note that the package is a prototype and should not be used in production. In the near future we are planing to integrate the general framework, including GBT (PEM) into **`mlr3`** via **`mlr3proba`** and
**`mlr3pipelines`**.

If you experience any problems or need advice on how to run the models, don't
hesitate to open an issue or contact the first author.

# Benchmark experiments
To run the code for the benchmark experiments you first need to install
two packages:


```r
# install pem.xgb
devtools::install_github("adibender/pem.xgb") # PEM via XGBoost
devtools::install_github("adibender/pammtools", ref = "ecml")# Data trafo
```

## Folder structure
- Code for benchmarks based on real data sets is contained within folder
**`benchmark`**
  - `prep_data.R`: preprocesses the different data sets used for benchmarking
  - `problems.R`: contains functions that generate an instances during each
  iteration of the benchmark
  - `algorithms.R`: contains wrapper functions for the different algorithms that
  - are compared to each other
  - `benchmark.R`: contains the code that will performs the benchmark experiments
  and is based on package `batchtools`.
  - `bench-time.R`: benchmark w.r.t. different strategies w.r.t. cut-point selection (scaling experiment)
  - `benchmark-metabric.R`: benchmarks for `metabric` data set for comparison with DeepHit
  - `benchmark-pbc-tvf`: benchmark of the PBC data sets (with time-dependent covariates)
  - `benchmark-mgus2.CR`: benchmark of the MGUS data with competing risks for
  comparison with DeepHit
  - `extract-benchmark-instances.R` script that extracts results from benchmark
  as well as instances to facilitate comparison with DeepHit

- Code for benchmarks based on synthetic data sets is contained within
folder **`simulation`**
  - `sim_funs.R` contains functions that simulate data for right-censored and
  competing risks data
  - `tv-effects.R`: contains code to run benchmark w.r.t. to performance in presence of time-varying effects/non-proportional hazards (right-censored and competing
  risks)

- Code for the aggregation and evaluation of benchmark studies are contained in folder **`paper`**
  - `results.R`: reads in the results from different experiments, combines results, creates raw tables (manually processed for publication)
