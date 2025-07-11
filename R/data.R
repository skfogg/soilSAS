#' Example 1 parameter data
#'
#' Example 1 shows a step function input of rain and conservative tracer. Q is the only outflow, with a
#' SAS Omega shape set to a step function with a value of 0.5.
#'
#' @format ## `example_1_params`
#' A list of 12 items:
#' \describe{
#'  \item{inputdata}{a data frame with 24,000 rows and 6 columns:
#'      'dates' Observation date-times;
#'      'J' Rain (mm/h);
#'    Q' Outflow (mm/h);
#'    'ET' Outflow due to Evapotranspiration (mm/h);
#'    'WI' Wetness Index (-);
#'    'C_J' Input concentration of tracer (mg/L)
#'  }
#'  \item{S0}{Initial Storage (mm)}
#'  \item{C_S0}{Initial Tracer Concentration in Storage (mg/L)}
#'  \item{dt}{Time Step}
#'  \item{SAS_Q_shape}{Shape of the SAS function (Omega) of Q}
#'  \item{SAS_ET_shape}{Shape of the SAS function (Omega) of ET}
#'  \item{shape_params_Q}{List of shape parameters of the SAS function of Q}
#'  \item{shape_params_ET}{List of shape parameters of the SAS function of ET}
#'  \item{f_thresh}{Fraction of rank storage after which the storage is sampled uniformly (1 indicates option not used)}
#'  \item{age_distributions}{Vector containng the date times of age distributions to output}
#'  \item{conservative}{Logical indicating whether or not the tracer is conservative}
#'  \item{decay_rate}{The decay rate of the non-conservative tracer}
#' }
#'
#' @source Authored by Katie Fogg-Tesar as a simple example
"example_1_params"
