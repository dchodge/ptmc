# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

run_ptmc <- function(model, dataList, settings, update_ind, PTMCpar, i) {
    .Call('_ptmc_run_ptmc', PACKAGE = 'ptmc', model, dataList, settings, update_ind, PTMCpar, i)
}

run_ptmc_discrete <- function(model, dataList, settings, update_ind, PTMCpar, i) {
    .Call('_ptmc_run_ptmc_discrete', PACKAGE = 'ptmc', model, dataList, settings, update_ind, PTMCpar, i)
}

