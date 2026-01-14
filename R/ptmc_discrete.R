
#' Parallel Tempering Monte Carlo for Discrete Model
#'
#' This function runs the Parallel Tempering Monte Carlo (PTMC) simulation for a discrete model. It 
#' checks and validates the settings, and if the parameter list is empty, it initializes the parameters 
#' with default values. The function then calls the `get_discrete_output` function to perform the actual 
#' simulation and return the results.
#'
#' @param model A list representing the model, which contains necessary information for running the 
#'              PTMC. This could include the model structure, parameter names, and other model-specific 
#'              settings.
#'
#' @param data A list containing the data required for running the model. This might include observed 
#'             data, priors, or any other inputs that the model requires to perform the simulation.
#'
#' @param settings A list of settings for the PTMC simulation. The list must include the following:
#'   - \code{numberChainRuns}: The number of chains to run in parallel.
#'   - Other settings relevant for running the PTMC simulation, which will be validated using the 
#'     \code{check_settings_discete} function.
#'
#' @param par (Optional) A list of parameters for the chains. If provided, each element should correspond 
#'            to the parameters for one chain. If not provided (or empty), default parameters will be used.
#'
#' @return A list containing the results of the PTMC simulation. The structure of the returned list will
#'         depend on the result of the `get_discrete_output` function, which includes:
#'   - \code{mcmc}: An MCMC object with posterior samples of the model parameters.
#'   - \code{discrete}: A list of discrete states generated during the simulation for each chain.
#'   - \code{lpost}: A data frame of log-posterior values, with columns representing different chains 
#'                  and rows representing samples.
#'   - \code{temp}: A data frame of temperatures for each chain at each sample.
#'   - \code{acc}: A data frame of acceptance rates for each chain at each sample.
#'   - \code{outPTpar}: A list containing the parameter values for each chain.
#'
#' @importFrom dplyr gather
#' @importFrom coda mcmc
#' @importFrom parallel mclapply
#' 
#' @export
ptmc_discrete_func <- function(model, data, settings, par = NULL) {
    settings <- check_settings_discete(settings, model)

  if (length(par) == 0) {
    par <- rep(list(list(type = "None")), settings[["numberChainRuns"]])
    output <- get_discrete_output(model, data, settings, FALSE, par)
  } else {
    output <- get_discrete_output(model, data, settings, TRUE, par)
  }
  output
}


#' Get Discrete Output from PTMC Model
#'
#' This function runs multiple Markov Chains in parallel or sequentially to generate output from a 
#' Parallel Tempering Monte Carlo (PTMC) model. The output includes the posterior parameter samples,
#' discrete state information, log posterior values, temperatures, and acceptance rates.
#'
#' @param model A list representing the model containing information such as parameter names and any 
#'              model-specific settings required for running the PTMC.
#' 
#' @param data_list A list containing the data necessary for running the model. This could include 
#'                  observed data, priors, and any other variables needed for the simulation.
#' 
#' @param settings A list containing the settings for the model execution. The list must contain:
#'                 - \code{numberChainRuns}: The number of chains to run in parallel.
#'                 - \code{runParallel}: A boolean indicating whether to run the chains in parallel (TRUE)
#'                   or sequentially (FALSE).
#'                 - \code{numberFittedPar}: The number of parameters to fit in the model.
#' 
#' @param update_ind An index or flag used to control which part of the model or data to update
#'                   during the PTMC simulation.
#' 
#' @param par A list of parameters or starting values for the chains. Each element corresponds to 
#'            one chain.
#'
#' @return A list containing the results of the PTMC simulation. The list includes:
#'   \item{mcmc}{An MCMC object containing the posterior samples of the model parameters.}
#'   \item{discrete}{A list of discrete states generated during the simulation for each chain.}
#'   \item{lpost}{A data frame of log-posterior values, with columns representing different chains 
#'                and rows representing samples.}
#'   \item{temp}{A data frame of temperatures for each chain at each sample.}
#'   \item{acc}{A data frame of acceptance rates for each chain at each sample.}
#'   \item{outPTpar}{A list containing the parameter values for each chain.}
#'
#' @importFrom dplyr gather
#' @importFrom coda mcmc
#' @importFrom parallel mclapply
#' @export
get_discrete_output <- function(model, data_list, settings, update_ind, par) {

  outPTpost <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTdiscrete <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTlp <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTtemp <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTacc <- vector(mode = "list", length = settings[["numberChainRuns"]])
  outPTpar <- vector(mode = "list", length = settings[["numberChainRuns"]])
  out_raw <- list()

  # Run the chains in parallel
  if (settings[["runParallel"]]) {
    out_raw <- mclapply(1:settings[["numberChainRuns"]], 
      function(i) {
        run_ptmc_discrete(model, data_list, settings, update_ind, par[[i]], i)
      },
      mc.cores = settings[["numberCores"]]
    )
  } else {
    for (i in 1:settings[["numberChainRuns"]]) {
      out_raw[[i]] <- run_ptmc_discrete(model, data_list, settings, update_ind, par[[i]], i)
    }
  }

  for(i in 1:settings[["numberChainRuns"]]) {
    out_post <- out_raw[[i]][["output"]][, 1:settings$numberFittedPar]
    outPTpar[[i]] <- out_raw[[i]][["PTMCpar"]]
    if (settings$numberFittedPar > 1){
        colnames(out_post) <- model[["namesOfParameters"]]
    }
    outPTpost[[i]] <- mcmc(out_post)
    outPTdiscrete[[i]] <- out_raw[[i]][["discrete"]]
    outPTlp[[i]] <- out_raw[[i]][["output"]][, settings$numberFittedPar + 1]
    outPTtemp[[i]] <- out_raw[[i]][["output"]][, settings$numberFittedPar + 2]
    outPTacc[[i]] <- out_raw[[i]][["output"]][, settings$numberFittedPar + 3]
  }

  outlpv <- data.frame(matrix(unlist(outPTlp), nrow = length(outPTlp[[1]])))
  colnames(outlpv) <- c(1:settings[["numberChainRuns"]])
  outlpv <- outlpv %>% gather(colnames(outlpv), key="chain_no",value="lpost")
  outlpv$sample_no <-rep(1:length(outPTlp[[1]]), settings[["numberChainRuns"]])

 # outdiscretev <- data.frame(matrix(unlist(outPTdiscrete), nrow = length(outPTdiscrete[[1]])))
 # colnames(outdiscretev) <- c(1:settings[["numberChainRuns"]])
  #outdiscretev <- outdiscretev %>% gather(colnames(outdiscretev), key="chain_no",value="lpost")
  #outdiscretev$sample_no <- rep(1:length(outdiscretev[[1]]), settings[["numberChainRuns"]])

  outltempv <- data.frame(matrix(unlist(outPTtemp), nrow=length(outPTtemp[[1]])))
  colnames(outltempv) <- c(1:settings[["numberChainRuns"]])
  outltempv <- outltempv %>% gather(colnames(outltempv), key="chain_no", value="temperature")
  outltempv$sample_no <- rep(1:length(outPTtemp[[1]]), settings[["numberChainRuns"]])
  
  outlaccv <- data.frame(matrix(unlist(outPTacc), nrow=length(outPTacc[[1]])))
  colnames(outlaccv) <- c(1:settings[["numberChainRuns"]])
  outlaccv <- outlaccv %>% gather(colnames(outlaccv), key="chain_no", value="acceptance rate")
  outlaccv$sample_no <- rep(1:length(outPTacc[[1]]), settings[["numberChainRuns"]])

  output <- list(
    mcmc = as.mcmc.list(outPTpost),
    discrete = outPTdiscrete,
    lpost = outlpv,
    temp = outltempv,
    acc = outlaccv,
    outPTpar = outPTpar
  )
  output
}

#' Check and Set Default Values for PTMC Settings
#'
#' This function checks the provided settings list for missing values and assigns default values where necessary. 
#' It ensures that all required settings for the Parallel Tempering Monte Carlo (PTMC) simulation are specified 
#' and provides appropriate defaults if any are missing. The function also validates specific parameters by checking 
#' the model for any necessary attributes.
#'
#' @param settings A list of settings for the PTMC simulation. The function checks and fills missing settings with
#'                 default values. The following settings can be provided:
#'   - `numberChainRuns`: The number of chains to run in parallel.
#'   - `numberCores`: The number of CPU cores to use for parallel processing (default is equal to `numberChainRuns`).
#'   - `numberTempChains`: The number of temperature chains for the Parallel Tempering (default is 10).
#'   - `iterations`: The number of iterations for each chain (default is 20,000).
#'   - `burninPosterior`: The burn-in period for posterior samples (default is 10,000).
#'   - `thin`: The thinning interval for samples (default is 100).
#'   - `consoleUpdates`: The frequency of console updates (default is 100).
#'   - `numberFittedPar`: The number of parameters to fit in the model. If not specified, it is set to the length 
#'                         of `model$namesOfParameters`.
#'   - `onAdaptiveCov`: Whether to use adaptive covariance updates (default is TRUE).
#'   - `updatesAdaptiveCov`: The frequency of updates for the adaptive covariance (default is 100).
#'   - `burninAdaptiveCov`: The burn-in period for adaptive covariance updates (default is 2,000).
#'   - `onAdaptiveTemp`: Whether to use adaptive temperature updates (default is TRUE).
#'   - `updatesAdaptiveTemp`: The frequency of updates for adaptive temperature (default is 10).
#'   - `onDebug`: Whether to enable debugging output (default is FALSE).
#'   - `lowerParBounds`: The lower bounds for the parameters. If not specified, it is set to `model$lowerParSupport_fitted`.
#'   - `upperParBounds`: The upper bounds for the parameters. If not specified, it is set to `model$upperParSupport_fitted`.
#'   - `covarInitVal`: The initial covariance value (default is 1e-10).
#'   - `covarInitValAdapt`: The initial adaptive covariance value (default is 1e-10).
#'   - `covarMaxVal`: The maximum allowed covariance value (default is 1).
#'   - `runParallel`: Whether to run the simulation in parallel (default is TRUE).
#'   - `lengthDiscreteVec`: The length of the discrete vector, which should be specified as `model$discrete_length`.
#'   - `updateDiscreteFreq`: The frequency at which to update the discrete vector (default is 0).
#'
#' @param model A list representing the model used in the PTMC simulation. The function checks if specific model-related
#'              attributes are available to fill in missing settings (e.g., `lowerParSupport_fitted`, `upperParSupport_fitted`, 
#'              and `discrete_length`).
#'
#' @return A list of validated and updated settings. If any setting was missing or invalid, the corresponding default 
#'         value is filled in, and a message is printed to the console.
#'
#' @examples
#' # Example usage of check_settings_discete function
#' model <- list(
#'   namesOfParameters = c("param1", "param2"),
#'   lowerParSupport_fitted = c(-5, -5),
#'   upperParSupport_fitted = c(5, 5),
#'   discrete_length = 10
#' )
#' settings <- list(
#'   numberChainRuns = 4,
#'   numberCores = 4,
#'   numberTempChains = 10,
#'   iterations = 20000,
#'   burninPosterior = 10000,
#'   thin = 100
#' )
#' validated_settings <- check_settings_discete(settings, model)
#' 
#' @importFrom utils cat
#' @export
check_settings_discete <- function(settings, model) {
  if (is.null(settings[["numberChainRuns"]])) {
    settings[["numberChainRuns"]] <- 4
    cat("`numberChainRuns` not specified in settings. Default value 4. \n")
  }

  if (is.null(settings[["numberCores"]])) {
    settings[["numberCores"]] <- settings[["numberChainRuns"]]
    cat("`numberCores` not specified in settings. Default value equal to `numberChainRuns`. \n")
  }

  if (is.null(settings[["numberTempChains"]])) {
    settings[["numberTempChains"]] <- 10
    cat("`numberTempChains` not specified in settings. Default value 10. \n")
  }  
  if (is.null(settings[["iterations"]])) {
    settings[["iterations"]] <- 20000
    cat("`iterations` not specified in settings. Default value 20,000. \n")
  }  
  if (is.null(settings[["burninPosterior"]])) {
    settings[["burninPosterior"]] <- 10000
    cat("`numberChainRuns` not specified in settings. Default value 10,000. \n")
  }  
  if (is.null(settings[["thin"]])) {
    settings[["thin"]] <- 100
    cat("`thin` not specified in settings. Default value 100. \n")
  }
  if (is.null(settings[["consoleUpdates"]])) {
    settings[["consoleUpdates"]] <- 100
    cat("`consoleUpdates` not specified in settings. Default value 100. \n")
  }
  if (is.null(settings[["numberFittedPar"]])) {
    if (is.null(model$namesOfParameters)) {
      stop("`numberFittedPar` not specified in settings. MUST be specified. \n")
    } 
    settings[["numberFittedPar"]] <- length(model$namesOfParameters)
    cat("`numberFittedPar` not specified in settings. Default value equal to the number of parameters in the model ", length(model$namesOfParameters), ". \n")
  }
  if (is.null(settings[["onAdaptiveCov"]])) {
        settings[["onAdaptiveCov"]] <- TRUE
    cat("`onAdaptiveCov` not specified in settings. Default value TRUE. \n")
  }
  if (is.null(settings[["updatesAdaptiveCov"]])) {
        settings[["updatesAdaptiveCov"]] <- 100
    cat("`updatesAdaptiveCov` not specified in settings. Default value 100. \n")
  }
  if (is.null(settings[["burninAdaptiveCov"]])) {
        settings[["burninAdaptiveCov"]] <- 2000
    cat("`burninAdaptiveCov` not specified in settings. Default value 2000. \n")
  }
  if (is.null(settings[["onAdaptiveTemp"]])) {
        settings[["onAdaptiveTemp"]] <- TRUE
    cat("`onAdaptiveTemp` not specified in settings.  Default value TRUE. \n")
  }
  if (is.null(settings[["updatesAdaptiveTemp"]])) {
        settings[["updatesAdaptiveTemp"]] <- 10
    cat("`updatesAdaptiveTemp` not specified in settings.  Default value 10. \n")
  }
  if (is.null(settings[["onDebug"]])) {
        settings[["onDebug"]] <- FALSE
  }
  if (is.null(settings[["lowerParBounds"]])) {
    if (is.null(model$lowerParSupport_fitted)) {
      stop("`lowerParBounds` not specified in settings. MUST be specified. \n")
    } 
    settings[["lowerParBounds"]] <- model$lowerParSupport_fitted
    cat("`lowerParBounds` not specified in settings. Defaults to lowerParSupport_fitted. \n")
  }
  if (is.null(settings[["upperParBounds"]])) {
    if (is.null(model$lowerParSupport_fitted)) {
      stop("`upperParBounds` not specified in settings. MUST be specified. \n")
    } 
    settings[["upperParBounds"]] <- model$upperParSupport_fitted
    cat("`upperParBounds` not specified in settings. Defaults to upperParSupport_fitted \n")
  }
  if (is.null(settings[["covarInitVal"]])) {
        settings[["covarInitVal"]] <- 1e-10
    cat("`covarInitVal` not specified in settings.  Default value 1e-10. \n")
  }
  if (is.null(settings[["covarInitValAdapt"]])) {
        settings[["covarInitValAdapt"]] <- 1e-10
    cat("`covarInitValAdapt` not specified in settings.  Default value 1e-10. \n")
  }
  if (is.null(settings[["covarMaxVal"]])) {
        settings[["covarMaxVal"]] <- 1
    cat("`covarMaxVal` not specified in settings. Default value 1. \n")
  }
  if (is.null(settings[["runParallel"]])) {

    settings[["runParallel"]] <- TRUE 
    cat("`runParallel` not specified in settings. Default value TRUE. \n")
  }

  if (is.null(settings[["lengthDiscreteVec"]])) {
    if (is.null(model$discrete_length)) {
      stop("`lengthDiscreteVec` not specified in settings. MUST be specified. \n")
    } 
    settings[["lengthDiscreteVec"]] <- model$discrete_length
    cat("`lengthDiscreteVec` not specified in settings. Defaults to ", model$discrete_length, ". \n")
  }
  if (is.null(settings[["updateDiscreteFreq"]])) {
        settings[["updateDiscreteFreq"]] <- 0
    cat("`updateDiscreteFreq` not specified in settings. Default value 0. \n")
  }

  settings
}