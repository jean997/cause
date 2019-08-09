CAUSE: Causal Analysis Using Summary Effect Estimates
======

[![Travis build status](https://travis-ci.com/jean997/cause.svg?branch=master)](https://travis-ci.com/jean997/cause)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/jean997/cause?branch=master&svg=true)](https://ci.appveyor.com/project/jean997/cause)
[![CircleCI build status](https://circleci.com/gh/jean997/cause.svg?style=svg)](https://circleci.com/gh/jean997/cause)

This R package implements the CAUSE method described in Morrison et al 2018 (pre-print coming soon).

Get started with an example analysis: https://jean997.github.io/cause/ldl_cad.html


To install:
```{r}
devtools::install_github("jean997/cause")
```
 You will need both `mixsqp` and `ashr`. Please use `mixsqp-0.1-97` which is currently the version available through CRAN. 
 

To install development versions you can use 
 
```{r, eval=FALSE}
devtools::install_github("stephenslab/mixsqp")
devtools::install_github("stephens999/ashr")
```

However, right now it is most advisable to use the versions available through CRAN.
