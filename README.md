CAUSE: Causal Analysis Using Summary Effect Estimates
======

This R package implements the CAUSE method described in Morrison et al 2019 (BioRxiv https://www.biorxiv.org/content/10.1101/682237v3).

_________________________________________
Install notes beginning 10-17-2019:

Two package dependencies `mixsqp` and `ashr` have recently been updated on CRAN. CAUSE has only been tested with the previous versions and there are reasons to believe that I might have to make some modifications for compatibility with the new ones. In the mean time, please use the following commands to install older versions. This will hopefully be sorted out in a few weeks.
```{r}
devtools::install_version("mixsqp", version = "0.1-97", repos = "http://cran.us.r-project.org")
devtools::install_version("ashr", version = "2.2-37", repos = "http://cran.us.r-project.org")
```
_________________________________________

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
