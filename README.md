CAUSE: Causal Analysis Using Summary Effect Estimates
======

This R package implements the CAUSE method described in Morrison et al 2019 (BioRxiv https://www.biorxiv.org/content/10.1101/682237v4).
Get started with an example analysis: https://jean997.github.io/cause/ldl_cad.html


To install:
```{r}
devtools::install_github("jean997/cause")
```

_________________________________________
Install notes beginning 10-17-2019:

CAUSE is only compatible with `mixsqp-0.1-97` and `ashr-2.2-32`. Please use the following commands to install older versions.
```{r}
devtools::install_version("mixsqp", version = "0.1-97", repos = "http://cran.us.r-project.org")
devtools::install_version("ashr", version = "2.2-32", repos = "http://cran.us.r-project.org")
```
_________________________________________
