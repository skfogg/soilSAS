##
## Functions to compute fSAS function as a certain shape
##

#' beta distribution SAS shape
#'
#' @param P_s the cumulative storage age distribution (S_T/S)
#' @param parameters list of shape parameters for the beta distribution: 'a' and 'b'. Lower values of parameter 'a'
#' indicate a higher affinity for younger ages, whereas lower values of parameter 'b' indicate a higher affinity for
#' older ages.
#' @param ... unused arguments
#'
#' @returns returns the cumulative SAS function, Omega.
#' @export
#'

shape_beta <- function(P_s, parameters, ...){

  if (length(parameters) != 2) {
    stop("Wrong number of input parameters. The beta function need 2 parameters: 'a' and 'b'.")
  }
  if (any(parameters < 0)) {
    stop("Beta function parameters must be positive.")
  }
  if(!names(parameters) %in% c("a", "b")){
    stop("Beta function parameters must be named 'a' and 'b'.")
  }

  # distribution of beta function (cdf)
  Omega <- stats::pbeta(P_s, parameters$a, parameters$b)

  return(Omega)
}

#' inverse power law SAS shape
#'
#' @description Calculates an inverse power law aSAS function: Omega = (1-P_s)^a.
#'
#' @param P_s the cumulative storage age distribution (S_T/S)
#' @param parameters A list with either parameter 'a' for a fixed exponent or parameters
#' 'a_min' and 'a_max' corresponding to wetness index of 1 ('a_min') and wetness index of 0 ('a_max')
#' @param wetness The wetness index of the system. A value between 0 and 1.
#'
#' @returns Returns the cumulative SAS function, Omega.
#' @export
#'

shape_inv_pl <- function(P_s, parameters, wetness = NULL) {

  if(!names(parameters) == "a" | !names(parameters) %in% c("a_min", "a_max")){
    stop("Inverse power law parameters must be named 'a' for a constant exponent or 'a_min' and 'a_max', which produces a variable exponent based on the system's wetness index.")
  }
  if(wetness < 0 | wetness > 1){
    stop("The wetness index for the variable exponent inverse power law shape function must be between 0 and 1.")
  }

  ## Variable exponent option:
  if (length(parameters) == 2) {
    a <- parameters$a_min + (1 - wetness) * (parameters$a_max - parameters$a_min)
  }

  Omega <- 1 - (1 - P_s)^parameters$a

  return(Omega)
}

#' normal (truncated) SAS shape
#'
#' @description Calculates a truncated normal SAS function.
#'
#' @param P_s the cumulative storage age distribution (S_T/S)
#' @param parameters A list of parameters 'm' for the mean and 'sd' for the standard deviation of the normal distribution.
#' @param ... unused arguments
#'
#' @returns Returns the cumulative SAS function, Omega.
#' @export
#'
shape_norm_trunc <- function(P_s, parameters, ...) {

  # Define the error function
  # see https://stackoverflow.com/questions/29067916/error-function-erfz
  erf <- function(x) 2 * stats::pnorm(x * sqrt(2)) - 1

  if(!names(parameters) %in% c("m", "sd")){
    stop("The parameters of the normal SAS shape must be 'm' (mean) and 'sd' (standard deviation). ")
  }
  if (parameters$sd < 0) {
    stop("The parameter sd of the normal SAS shap must be positive.")
  }

  a <- 0
  b <- 1
  phi_a <- erf((a - parameters$m) / (parameters$sd * sqrt(2)))
  phi_b <- erf((b - parameters$m) / (parameters$sd * sqrt(2)))
  Z <- phi_b - phi_a

  Omega <- (erf((P_s - parameters$m) / (parameters$sd * sqrt(2))) - phi_a) / Z

  return(Omega)
}

#' power law SAS shape
#'
#' @description Function to compute a power law fSAS function: Om = P_s^k.
#'
#' @param P_s the cumulative storage age distribution (S_T/S)
#' @param parameters A list including the power law exponent, 'k'.
#' @param ... unused arguments
#'
#' @returns Returns the cumulative SAS function, Omega.
#' @export
#'

shape_pl <- function(P_s, parameters, ...) {

  if (!names(parameters) == "k") {
    stop("The power law parameter list must include the exponent 'k'.")
  }
  if (parameters$k < 0) {
    stop("The exponent of the power-law function, 'k', must be positive.")
  }

  Omega <- P_s^parameters$k

  return(Omega)
}

#' time-variant power law SAS shape
#'
#' @param P_s the cumulative storage age distribution (S_T/S)
#' @param parameters list of parameters: 'k_min' and 'k_max'. 'k_min' corresponds
#' to a wetness index of 1 and 'k_max' corresponds to a wetness index of 0.
#' @param wetness the wetness index of the system. A value between 0 and 1.
#'
#' @returns Returns the cumulative SAS function, Omega.
#' @export
#'

shape_pl_tv <- function(P_s, parameters, wetness) {
  # function to compute a pl (power) fSAS function: Om = P_s^k(wi)
  # Om = fSAS_pl(P_s,par=[kmin,kmax],wi)  %power function with variable exponent
  # between kmin and kmax, depending on the system state wi

  k <- parameters$k_min + (1 - wetness) * (parameters$k_max - parameters$k_min)

  # do some error check
  if (!names(parameters)  %in% c("k_min", "k_max")) {
    stop("The parameters of the time-variant power-law function must be: 'k_min' and 'k_max'.")
  }
  if (any(parameters < 0)) {
    stop("The parameters of the time-variant power-law function must be positive.")
  }

  Omega <- P_s^k

  return(Omega)
}

#' step SAS shape
#'
#' @description A function to compute a step fSAS function (a "left" fSAS function).
#'
#'
#' @param P_s the cumulative storage age distribution (S_T/S)
#' @param parameters list with parameter 'u', which indicates the top interval of
#' where Omega is linear. 'u' must be between 0 and 1.
#' @param ... unused arguments
#'
#' @returns Returns the cumulative SAS function, Omega.
#' @export
#'
#' @details When P_s < u, Omega = P_u/u. When P_s >= u, Omega = 1.
#'

shape_step <- function(P_s, parameters, ...) {

  if (names(parameters) != "u") { # error in parameter input
    stop("The parameter for the step SAS funciton must be named 'u'.")
  }
  if (parameters$u <= 0 | parameters$u > 1) {
    stop("The parameter, 'u', of the step function must be 0 < u <=1.")
  }

  Omega <- P_s / parameters$u
  Omega[P_s >= parameters$u] <- 1

  return(Omega)
}
