
# An example SEIR compartmental model implemented in stan

> [There are many like it, but this one is
> mine](https://en.wikipedia.org/wiki/Rifleman%27s_Creed)

An example SEIR model implemented in stan for the CMMID computational
and inference theme meeting.

## What

> Stan is a state-of-the-art platform for statistical modeling and
> high-performance statistical computation. Thousands of users rely on
> Stan for statistical modeling, data analysis, and prediction in the
> social, biological, and physical sciences, engineering, and business.

> Users specify log density functions in Stan’s probabilistic
> programming language and get:

>   - full Bayesian statistical inference with MCMC sampling (NUTS, HMC)
>   - approximate Bayesian inference with variational inference (ADVI)
>   - penalized maximum likelihood estimation with optimization (L-BFGS)

> Stan’s math library provides differentiable probability functions &
> linear algebra (C++ autodiff). Additional R packages provide
> expression-based linear modeling, posterior visualization, and
> leave-one-out cross-validation.

See [here](https://mc-stan.org) for more.

## Why

  - Community, documention and ongoing development.
  - A good enough tool for a broad range of problems. Useful if wanting
    to do anything other than compartmental models.
  - Clean DSL which can be extended with functions written in C++.
  - Largely automated and optimised MCMC reduces the cognitive load when
    modelling.

## Why not

  - MCMC often not a great choice for complex compartmental model
    systems. For these models tools that support PMCMC or SMC^2 are
    likely more optimal (Libbi/Birch, ODIN/Dust, etc).

  - 
## Resouces (a small selection):

  - [Bayesian workflow for disease transmission modelling in
    stan](https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html)
  - [Fitting Bayesian models using Stan and
    R](https://www.weirdfishes.blog/blog/fitting-bayesian-models-with-stan-and-r/)
  - [Fitting a basic SIR model in
    stan](https://www.generable.com/blog/2020/04/fitting-a-basic-sir-model-in-stan/)
  - [Contemporary statistical inference for infectious disease models
    using Stan](https://arxiv.org/abs/1903.00423v3)
  - [covidseir](https://github.com/seananderson/covidseir)

## Example

Code below is heavily based on [Bayesian workflow for disease
transmission modelling in
stan](https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html).

  - Load some packages

<!-- end list -->

``` r
library(rstan)
```

  - Load the model

<!-- end list -->

``` r
model <- stan_model("model.stan")
```
