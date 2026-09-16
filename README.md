# nsddm

Maximum likelihood and Bayesian estimation of the parameters of a hierarchical 4-parameter drift diffusion model. Functions to estimate the parameters of single-person models are also implemented. The package also contains functions for the 7-parameter model. However, these have not been tested.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("nestlersteffen/nsddm")
```

To update, simply rerun the installation command.

## Example 1: single-person data

In the first example, we fit a single-person 4-parameter DDM with ML. 

``` r

# get some data

data(sim1)

# a control-argument list, see ddm_args() for further information

# here we use the dao-parametrization, the standard ddm density 
# (see Tuerlinckx, 2004)

ctrl <- list( type_alpha="dao", type_ddm="std" )

# fit the model 

fit <- ddm(rt=df$rt, xs=df$resp, estimator="ML", control=ctrl, verbose=TRUE )

# we have not implemented a summary function, but you can use:

fit$parm_table

```

In the next example, we fit the same model with Bayes using two chains

``` r

# some control arguments

ctrl <- list( type_alpha="dao", type_ddm="std", biter=5000, burnin=2500, nchain=2 )

# fit the model 

fit <- ddm(rt=df$rt, xs=df$resp, estimator="Bayes", control=ctrl, verbose=TRUE )

# and now the parm_table:

fit$parm_table

```

## Example 2: hierarchical data

In this example, we fit a hierarchical 4-parameter DDM with ML. 

...

## Contributing

Issues and pull requests are not actively monitored. For questions, 
suggestions, or bug reports, please contact me directly via mail.

## Status

Work in progress.