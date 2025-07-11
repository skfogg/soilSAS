##
## Code to prepare `example_1_params` dataset
##

library(lubridate)

## fake J [mm/h]
fJ <- c(rep(0.2, 7*24*10),
        rep(0, 22320))

## fake Q [mm/h]
fQ <- c(rep(0.08, 7*24*10),
        seq(0.08, 0.025, length.out = 22320))

## fake ET [mm/h]
fET <- rep(0, times = length(fQ))

## fake wi [mm/h]
fWI <- NA

## fake Cin for atmospheric rain tracer [mg/L]
fCin <-  c(rep(20, 7*24*10),
           rep(0, 22320))

## Dates
dates <- seq(mdy_hms("02-01-2020 00:00:00"), length.out = length(fJ), by = 3600)

## Create list of all model input data
example_1_params <- list(inputdata = data.frame(dates = dates,
                                                J = fJ,
                                                Q = fQ,
                                                ET = fET,
                                                WI = fWI,
                                                C_J = fCin
                                                ),
                         S0 = 1000, # initial storage [mm]
                         C_S0 = 0, # concentration of the initial storage [mg/L]
                         dt = as.numeric(dates[2] - dates[1]),
                         SAS_Q_shape ="shape_step", #power-function SAS
                         SAS_ET_shape = "shape_step", #step SAS
                         shape_params_Q = list(u = 0.5), # [-]
                         shape_params_ET = list(u = 0.5), # [-]
                         f_thresh = 1,
                         age_distributions = mdy_hms(c("03-01-2021 12:00:00",
                                                       "04-01-2021 12:00:00",
                                                       "05-01-2021 12:00:00",
                                                       "06-01-2021 12:00:00",
                                                       "07-01-2021 12:00:00",
                                                       "08-01-2021 12:00:00")), # choose dates for which TTDs are saved # fraction of rank storage [-] after which the storage is sampled uniformly (leave f_thresh = 1 if not interested in this option)
                         conservative = TRUE,
                         decay_rate = 1
)

usethis::use_data(example_1_params, overwrite = TRUE)

