CAUSE: Causal Analysis Using Summary Effect Estimates
======

This R package implements the CAUSE method described in Morrison et al 2019 (BioRxiv https://www.biorxiv.org/content/10.1101/682237v4).
Get started with an example analysis: https://jean997.github.io/cause/ldl_cad.html


### Important announcement for v1.1.0 and earlier

We have discovered that there is an error introduced in cause elpd computations using newer versions of the loo package with cause version 1.1.0 due to a change in the order of the `loo_compare` output. 

If you are new to CAUSE please use v1.2.0 or the most recent github version. We will issue a new release soon. 

If you have used cause version 1.1.0 or used cause with loo version 2.3.0 or 2.3.1 please do the following steps. 

1. Insall v1.2.0 or the most recent version from github. 
2. Recompute the elpd model table using the `recompute_elpd_table` function. 
```
elpd_table <- recompute_elpd_table(res)
res$elpd <- elpd_table
```
Here `res` is an object produced by running the cause function. You may find that the new table and the old table are the same or the sign of the z-score comparing models will flip, changing the p-value. The error only occurred in situations when the sharing model is better than the causal model. 


### Release Notes:

### Latest

+ The `pval_cols` option in `gwas_merge` is available as of version `1.2.0.0320`.

#### v1.2.0:

+ Corrected bug that came from re-ordering of loo_compare output in >= 2.3.0.
+ Added function to recompute z-score for previously run analysis.
+ Minor updates to parameter checking in gwas_merge.

#### v1.1.0:

+ Compatible with current versions of mixsqp (0.3-43) and ashr (2.2-47). This version is no longer compatible with versions of mixsqp prior to 0.3.
+ Compatible with loo 2.3.0 (no more warning about `compare` being depricated).
+ elpd object no longer contains comparisons with null model.
+ Improved documentation and arguments for cause function. Added the pval_thresh argument.



### Installation Instructions:

The original version of the `cause` R package is only compatible with earlier versions of `mixsqp` and `ashr`. The latest version is compatible with newer versions of those packages. If you want to exactly replicate the results in the paper you should use version 1.0.0. It is possible that the newer version is slightly faster. This gives two installation options:

#### 1. Version 1.0.0 with older `mixsqp` and `ashr`:

```{r}
devtools::install_version("mixsqp", version = "0.1-97", repos = "http://cran.us.r-project.org")
devtools::install_version("ashr", version = "2.2-32", repos = "http://cran.us.r-project.org")
devtools::install_github("jean997/cause@v1.0.0")
```
Don't allow R to update `mixsqp` or `ashr`.

#### 2. Latest version or version 1.2.0

For the release version
```{r}
devtools::install_github("jean997/cause@v1.2.0")
```

For the development version
```{r}
devtools::install_github("jean997/cause")
```
Note that for CAUSE versions after 1.0.0, you must use the newer version of `mixsqp` (0.3.XX)
