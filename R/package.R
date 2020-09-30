#' CAUSE
#'
#' Implementation of CAUSE described in Morrison et al. (2020)
#'
#' @docType package
#' @author Jean Morrison <jvmorr@umich.edu>
#' @import purrr tidyr mixsqp dplyr ashr ggplot2
#' @importFrom scales muted
#' @importFrom matrixStats logSumExp
#' @importFrom intervals Intervals interval_overlap
#' @importFrom gridExtra arrangeGrob tableGrob
#' @importFrom loo loo loo_compare
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom readr read_tsv write_tsv
#' @useDynLib cause
#' @name cause
NULL
