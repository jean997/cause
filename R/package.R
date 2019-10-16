#' CAUSE
#'
#' Implementation of CAUSE described in Morrison et al. (2019)
#'
#' @docType package
#' @author Jean Morrison <jeanm@uchicago.edu>
#' @import Rcpp ashr dplyr ggplot2
#' @import purrr tidyr mixsqp
#' @importFrom scales muted
#' @importFrom matrixStats logSumExp
#' @importFrom numDeriv hessian
#' @importFrom intervals Intervals interval_overlap
#' @importFrom gridExtra arrangeGrob tableGrob
#' @importFrom loo loo compare
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom readr read_tsv write_tsv
#' @useDynLib cause
#' @name cause
NULL
