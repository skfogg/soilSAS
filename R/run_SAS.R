#' run a SAS model
#'
#' @param model_parameters a list of the model parameters. See Details.
#'
#' @returns returns a list with model_parameters and model_output
#'
#' @seealso [example_1_params], [plot_SAS_input()], and [plot_SAS_output()]
#'
#' @examples
#' data("example_1_params")
#' ex1model <- run_SAS(example_1_params)
#' #plot_SAS_input(ex1model$model_parameters)
#' #plot_SAS_output(ex1model$model_output)
#'
#' @export
#'
#'
run_SAS <- function(model_parameters = example_1_params){

  requireNamespace("wesanderson")
  requireNamespace("lubridate")

  ##  INPUT WARNINGS
  if(is.list(model_parameters) == FALSE){
    stop("Model parameters are not a list. Model parameters must be in a list with an 'inputdata' dataframe, and the following named parameters: 'S0', 'C_S0', 'dt', 'SAS_Q_shape', 'SAS_ET_shape', 'shape_params_Q', 'shape_params_ET', f_thresh', 'age_distributions', 'conservative', 'decay_rate'.")
  }
  if(any(colnames(model_parameters$inputdata) != c("dates", "J", "Q", "ET", "WI", "C_J"))){
    stop("Input data columns missing or incorrect. Column names must be: 'dates', 'J', 'Q', 'ET', 'WI', 'C_J'.")
  }
  if(any(names(model_parameters) != c("inputdata", "S0", "C_S0", "dt", "SAS_Q_shape",
                                      "SAS_ET_shape", "shape_params_Q", "shape_params_ET",
                                      "f_thresh", "age_distributions", "conservative", "decay_rate"))){
    stop("Model parameters missing or named incorrectly. Model parameters must be in a list with an 'inputdata' dataframe, and the following named parameters: 'S0', 'C_S0', 'dt', 'SAS_Q_shape', 'SAS_ET_shape', 'shape_params_Q', 'shape_params_ET', f_thresh', 'age_distributions', 'conservative', 'decay_rate'.")
  }
  if(is.data.frame(model_parameters$inputdata) == FALSE | is.list(model_parameters$inputdata) == FALSE){
    stop("Model input data must be a data frame or named list.")
  }
  if(!any(is.numeric(model_parameters$S0),
         is.numeric(model_parameters$C_S0),
         is.numeric(model_parameters$dt),
         is.numeric(model_parameters$shape_params_Q),
         is.numeric(model_parameters$shape_params_ET),
         is.numeric(model_parameters$f_thresh),
         is.numeric(model_parameters$decay_rate),
         is.character(model_parameters$SAS_Q_shape),
         is.character(model_parameters$SAS_ET_shape),
         # is.POSIXt(model_parameters$age_distributions),
         is.logical(model_parameters$conservative))){
    stop("Model parameter class error.
         Numeric parameters: 'S0', 'C_S0', 'dt', 'shape_params_Q', 'shape_params_ET', f_thresh', decay_rate'.
         Character parameters: 'SAS_Q_shape', 'SAS_ET_shape'.
         POSIXt parameters: 'age_distributions'.
         Logical parameters: 'conservative'.")
  }
  if(!model_parameters$SAS_Q_shape %in% c("shape_beta", "shape_inv_pl", "shape_norm_truc", "shape_pl", "shape_pl_tv", "shape_step") |
     !model_parameters$SAS_ET_shape %in% c("shape_beta", "shape_inv_pl", "shape_norm_truc", "shape_pl", "shape_pl_tv", "shape_step")){
    stop("SAS shape not available. Choose from: 'shape_beta', 'shape_inv_pl', 'shape_norm_truc', 'shape_pl', 'shape_pl_tv', 'shape_step'.")
  }


  inputdata <- model_parameters$inputdata

  # assign CALIBRATION parameters
  parQ  <- model_parameters$shape_params_Q
  parET <- model_parameters$shape_params_ET
  S0 <- model_parameters$S0

  # set a few constants
  NN     <- length(model_parameters$inputdata$J) # length of the timeseries
  ndistr <- length(model_parameters$age_distributions) # number of age distributions that will be saved

  # preallocate variables
  S_T      <- numeric(NN) # rank storage (function of age T)
  C_ST     <- numeric(NN) # rank storage concentration (function of age T)
  C_Q      <- numeric(NN) # stream concentration (function of time t)
  age_matr <- matrix(0, nrow = NN, ncol = ndistr) # matrix to store the desired age distributions
  med      <- matrix(0, nrow = NN, ncol = 3) # selected percentile(s) of discharge age
  Fyw      <- matrix(0, nrow = NN, ncol = 3) # young water fraction(s)

  # initial conditions
  C_Q[1]   <- model_parameters$C_S0 # initial streamflow concentration equal to the initial storage
  length_s <- 1                                                        # rank storage vector length
  S_T[1]   <- S0                                                       # initial rank storage [mm]
  C_ST[1]  <- model_parameters$C_S0                                                # mean concentration of the initial rank storage
  Omega_Q  <- do.call(model_parameters$SAS_Q_shape, list(P_s = S_T[1:length_s]/S_T[1],
                                                         parameters = parQ, wetness = model_parameters$inputdata$WI[1]))   # [-]
  Omega_ET <- do.call(model_parameters$SAS_ET_shape, list(P_s = S_T[1:length_s]/S_T[1], parameters = parET, wetness = model_parameters$inputdata$WI[1]))

  #--------------------------------------------------------------------------
  # MODEL LOOPS
  #--------------------------------------------------------------------------

  for (j in 1:(NN - 1)) {
    # j = 2
    #------------------------------------------------------------------
    # SOLVE THE AGE BALANCE and evaluate the rank storage concentration
    #------------------------------------------------------------------

    # 0) define the domain for the SAS function evaluation (basically a shifted S_T with new water addition)
    age1 <- max(0, model_parameters$dt * (inputdata$J[j] - inputdata$Q[j]*Omega_Q[1] - inputdata$ET[j]*Omega_ET[1]))
    dom  <- (c(0, S_T[1:length_s]) + age1) / (S_T[length_s] + age1) # rescaled domain

    # 1) evaluate the SAS functions Omega over the domain 'dom'
    Omega_Q  <- do.call(model_parameters$SAS_Q_shape, list(dom, parQ, inputdata$WI[j])) # [-]
    Omega_ET <- do.call(model_parameters$SAS_ET_shape, list(dom, parET, inputdata$WI[j])) # [-]

    # 2) solve the water age balance (S_T = rank storage)
    S_T[1:(length_s + 1)] <- pmax(0, c(0, S_T[1:length_s]) +
                                    model_parameters$dt * inputdata$J[j] -
                                    model_parameters$dt * (inputdata$Q[j]*Omega_Q + inputdata$ET[j]*Omega_ET))
    for (i in 2:(length_s + 1)) {
      S_T[i] <- max(S_T[i], S_T[i - 1]) # ensure non-decreasing S_T
    }

    # 3) update solute concentration for each parcel
    if(model_parameters$conservative){
      # conservative transport:
      C_ST[2:(length_s + 1)] <- C_ST[1:length_s]
    }else{
      # non-conservative transport:
      C_ST[2:(length_s + 1)] <- C_ST[1:length_s] * model_parameters$decay_rate
    }



    C_ST[1] <- inputdata$C_J[j] # concentration of new input

    # 4) check if the vectors need to grow or not
    if (j == 1 || S_T[length_s] < model_parameters$f_thresh * S_T[length_s + 1]) {
      length_s <- length_s + 1                                                             # still need to grow
    } else {
      # oldest elements are merged into an old pool
      numer <- C_ST[length_s + 1] * (S_T[length_s + 1] - S_T[length_s]) +
        C_ST[length_s]     * (S_T[length_s]     - S_T[length_s - 1])
      denom <- S_T[length_s + 1] - S_T[length_s - 1]
      C_ST[length_s] <- max(0, numer / denom)                                              # update mean concentration
      S_T[length_s]  <- S_T[length_s + 1]                                                  # merge storage
      Omega_Q[length_s]  <- Omega_Q[length_s + 1]                                          # merge SAS values
      Omega_ET[length_s] <- Omega_ET[length_s + 1]
      Omega_Q  <- Omega_Q[1:length_s]
      Omega_ET <- Omega_ET[1:length_s]
    }

    #----------------------------------------
    # COMPUTE output: stream concentration
    #----------------------------------------

    pQ <- diff(c(0, Omega_Q)) # this is pQ(T)*dT and it is equivalent to omegaQ(S_T)*dS_T
    C_Q[j + 1] <- sum(C_ST[1:length_s] * pQ) # modeled streamflow concentration

    #-------------------------
    # COMPUTE output: other
    #-------------------------

    # for the selected dates, store discharge age distributions at selected dates
    if (inputdata$dates[j] %in% model_parameters$age_distributions) {
      idx <- which(inputdata$dates[j] == model_parameters$age_distributions)
      age_matr[1:length_s, idx] <- pQ
    }

    # ex2 - compute percentile ages
    pp <- c(0.1, 0.25, 0.5)                                                                # percentiles [-]
    for (i in 1:length(pp)) {
      med[j + 1, i] <- which(Omega_Q >= pp[i])[1]
    }

    # ex3 - compute young water fractions
    ywt <- c(7, 60, 90)                                                                    # thresholds [days]
    for (i in 1:length(ywt)) {
      idx <- min(length(pQ), round(ywt[i] * 24/model_parameters$dt))
      Fyw[j + 1, i] <- Omega_Q[idx]
    }
  }

  #--------------------------------------------------------------------------
  # RETURN
  #--------------------------------------------------------------------------
  output <- list(
    C_Q = C_Q,
    S_T = S_T,
    C_ST = C_ST,
    age_matr = age_matr,
    med = med,
    Fyw = Fyw,
    Omega_Q = Omega_Q,
    Omega_ET = Omega_ET
  )

  return(list(model_parameters = model_parameters,
              output = output))
}
