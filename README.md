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

fit <- ddm(rt=sim1$rt, xs=sim1$resp, estimator="ML", control=ctrl, verbose=TRUE )

# we have not implemented a summary function, but you can use:

fit$parm_table

```

In the next example, we fit the same model with Bayes using two chains

``` r

# some control arguments

ctrl <- list( type_alpha="dao", type_ddm="std", biter=5000, burnin=2500, nchain=2 )

# fit the model 

fit <- ddm(rt=sim1$rt, xs=sim1$resp, estimator="Bayes", control=ctrl, verbose=TRUE )

# and now the parm_table:

fit$parm_table

```

## Example 2: hierarchical data

We first fit a hierarchical 4-parameter DDM with ML. 

``` r

data(sim2)

# again some control arguments, see ddm_args() for further information

# we use five points in AGH, we want to see the inner optimization and 10 cores are used during fitting 

ctrl <- list( type_alpha="dao", type_ddm="std", nPoints=5, verbose_inn=TRUE, n_threads=10L )

# fit the model, note that the function is called hddm, and that you have to provide an id-column 

fit <- hddm(rt=sim2$rt, xs=sim2$resp, id=sim2$id, estimator="ML", control=ctrl, verbose=TRUE )

# and now the parm_table:

fit$parm_table

```

And here is how the same data can be fit with Bayes using the PmWG - algorithm described in Dao et al. (2025)

``` r

# the control arguments, see ddm_args() for further information

# the chain length is 3,000 with 1,500 used as burnin, we draw 20 particles 

ctrl <- list( type_alpha="dao", type_ddm="std", biter=500, burnin=250, nchain=1, 
    type_proposal="pmwg", type_sigma_prior="huang_wand", R=20 )

# fit the model 

fit <- hddm(rt=sim2$rt, xs=sim2$resp, id=sim2$id, estimator="Bayes", control=ctrl, verbose=TRUE )

# the result as parm_table:

fit$parm_table

```

## Contributing

Issues and pull requests are not actively monitored. For questions, 
suggestions, or bug reports, please contact me directly via mail.

## Status

Work in progress.