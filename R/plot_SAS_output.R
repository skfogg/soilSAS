#' standardized SAS output plots
#'
#' @description Plots standardized output for trans-SAS soil models.
#'
#' @param model_output output list produced by running run_SAS()
#' @param save a logical of whether or not you want the plots to be saved to disk
#' @param dir_path path to the directory where you wan the plots to be saved. Used only if save = TRUE.
#'
#' @returns Returns four plots: Age quantiles, CDFs, Concentrations, PDFs, and Young Water Fractions.
#' @export
#' @import wesanderson
#' @importFrom grDevices dev.off hcl.colors png
#' @importFrom graphics abline axis legend lines mtext par segments
plot_SAS_output <- function(model_output, save = FALSE, dir_path = NULL){


  model_output <- ex1model

  requireNamespace("wesanderson")

  inputdata <- model_output$model_parameters$inputdata
  outputdata <- model_output$output
  model_params <- model_output$model_parameters
  dt <- model_params$dt

  input_pal <- wesanderson::wes_palette("AsteroidCity2")

  #### CONCENTRATION PLOT ####
  if(save){
    png(paste0(dir_path, "/out_concentrations.png"),
        width = 600*3,
        height = 450*3,
        res = 72*3)
  }
  par(mfrow = c(1,1))
  plot(inputdata$dates,
       inputdata$C_J,
       type = "l",
       col = input_pal[1],
       lwd = 2,
       ylab = "C [mg/L]",
       xlab = "Dates",
       ylim = c(0,max(inputdata$C_J)+(0.1*max(inputdata$C_J))))
  lines(inputdata$dates,
        outputdata$C_Q,
        col = "orange", lwd = 2)
  legend("topright",
         c("Input Concentration", "Output Concentration"),
         col = c(input_pal[1], "orange"),
         lwd = 2,
         cex = 0.6)
  if(save){
    dev.off()
  }

  ### AGE PDF ###
  # get last age of age matrix
  last_age <- numeric(ncol(outputdata$age_matr))
  for(i in 1:ncol(outputdata$age_matr)){
    last_age[i] <- max(which(outputdata$age_matr[,i] > 0))
  }

  #### PDF ####
  if(save){
    png(paste0(dir_path, "/pdfs.png"),
        width = 600*3,
        height = 450*3,
        res = 72*3)
  }
  plot(dt/24*1:length(outputdata$age_matr[1:last_age[1]-1, 1]),
       outputdata$age_matr[1:last_age[1]-1,1],
       xlim = c(0,max(last_age/24)),
       type = "n",
       xlab = "Age [d]",
       ylab = "Frequency [1/d]",
       main = "Selected Streamflow Age PDF")
  for(i in 1:length(last_age)){
    segments(x0 = dt/24*1:length(outputdata$age_matr[1:last_age[i]-1,i]),
             y0 = rep(0, times = length(outputdata$age_matr[1:last_age[i]-1,i])),
             x1 = dt/24*1:length(outputdata$age_matr[1:last_age[i]-1,i]),
             y1 = outputdata$age_matr[1:last_age[i]-1,i],
             col = hcl.colors(length(model_params$age_distributions), "Purp")[i],
             lwd = 1
    )
  }
  legend("topright",
         as.character(paste0(lubridate::month(model_params$age_distributions,
                                   label = T,
                                   abbr= T),
                             "-",
                             lubridate::day(model_params$age_distributions))),
         lty = 1,
         col = hcl.colors(length(model_params$age_distributions), "Purp"),
         lwd = 2,
         cex = 0.6)
  if(save){
    dev.off()
  }

  #### AGE CDFs ####
  if(save){
    png(paste0(dir_path, "/cdfs.png"),
        width = 600*3,
        height = 450*3,
        res = 72*3)
  }
  plot(dt/24*1:length(outputdata$age_matr[1:last_age[1]-1,1]),
       cumsum(outputdata$age_matr[1:last_age[1]-1,1]),
       type = "l",
       xlab = "Age [d]",
       ylab = "cumulative frequency [-]",
       main = "Selected Streamflow Age CDF",
       col = hcl.colors(length(model_params$age_distributions), "Purp")[1],
       ylim = c(0,1),
       xlim = c(0,max(last_age/24)),
       lwd = 2)
  for(i in 1:length(model_params$age_distributions)){
    lines(dt/24*1:length(outputdata$age_matr[1:last_age[i]-1,i]),
          cumsum(outputdata$age_matr[1:last_age[i]-1,i]),
          col = hcl.colors(length(model_params$age_distributions), "Purp")[i],
          lwd = 2)
  }
  legend("topleft",
         as.character(paste0(lubridate::month(model_params$age_distributions,
                                   label = T,
                                   abbr= T),
                             "-",
                             lubridate::day(model_params$age_distributions))),
         lty = 1,
         col = hcl.colors(length(model_params$age_distributions), "Purp"),
         lwd = 2,
         cex = 0.6)
  if(save){
    dev.off()
  }

  #### AGE QUANTILES ####
  agestat_pal <- wesanderson::wes_palette("FantasticFox1")

  if(save){
    png(paste0(dir_path, "/age_quantiles.png"),
        width = 600*3,
        height = 450*3,
        res = 72*3)
  }
  par(mfrow = c(2,1),
      mar = c(3,4,3,4))
  plot(inputdata$dates,
       outputdata$med[, 1]/(24/dt),
       type = "l",
       col = agestat_pal[1],
       ylim = c(0,max(outputdata$med[,1]/24)),
       ylab = "Age [d]",
       xlab = "",
       main = "Age Quantile Time Series",
       lwd = 2)
  for(i in 2:3){
    lines(inputdata$dates,
          outputdata$med[, i]/(24/dt),
          col = agestat_pal[i],
          lwd = 2)
  }
  legend("topleft",
         c("q 0.10", "q 0.25", "q 0.50"),
         lty = 1,
         col = agestat_pal,
         lwd = 2,
         cex = 0.6)

  if(max(inputdata$Q)*10 > max(inputdata$J)){
    plot(J ~ dates,
         data = inputdata,
         type = "l",
         lwd = 2,
         col = input_pal[1],
         pch = 16,
         ylab = "J, Q, & ET [mm/h]",
         xlab = "Date",
         main = "Model Input: Flows")
    lines(I(Q) ~ dates,
          data = inputdata,
          col = input_pal[3],
          lwd = 2)
    lines(I(ET) ~ dates,
          data = inputdata,
          col = input_pal[5],
          lwd = 2)
  }else{
    plot(J ~ dates,
         data = inputdata,
         type = "l",
         lwd = 2,
         col = input_pal[1],
         pch = 16,
         ylab = "J [mm/h]",
         xlab = "Date",
         main = "Model Input: Flows")
    axis(side = 4,
         at = seq(0, 2, by = 0.5),
         labels = seq(0, 2, by = 0.5)/10)
    mtext("Q & ET [mm/h]",
          side = 4,
          line = 2.5)
    lines(I(Q*10) ~ dates,
          data = inputdata,
          col = input_pal[3],
          lwd = 2)
    lines(I(ET*10) ~ dates,
          data = inputdata,
          col = input_pal[5],
          lwd = 2)
  }
  if(anyNA(inputdata$WI)){
    legend("topright",
           c("J", "Q", "ET"),
           col = input_pal[c(1,3,5)],
           lty = 1,
           cex = 0.6)
  }else{
    lines(WI ~ dates,
          data = inputdata,
          type = "l",
          col = input_pal[6],
          lwd = 2)
    legend("topright",
           c("J", "Q", "ET", "WI"),
           col = input_pal[c(1,3,5,7)],
           lty = 1,
           cex = 0.6)
  }
  if(save){
    dev.off()
  }



  # show the young water fraction (Fyw)
  youngwat_pal <- wesanderson::wes_palette("Darjeeling1")

  if(save){
    png(paste0(dir_path, "/young_water_frac.png"),
        WIdth = 600*3,
        height = 450*3,
        res = 72*3)
  }
  par(mfrow = c(2,1),
      mar = c(3,4,3,4))
  plot(inputdata$dates,
       outputdata$Fyw[, 1],
       type = "l",
       col = youngwat_pal[1],
       ylim = c(0,1),
       ylab = "Fraction of Young Water [-]",
       xlab = "",
       main = "Young Water Fraction Time Series", lwd = 2)
  for(i in 2:3){
    lines(inputdata$dates,
          outputdata$Fyw[, i],
          col = youngwat_pal[i],
          lwd = 2)
  }
  legend("topright",
         c("7 d", "60 d", "90 d"),
         lty = 1,
         col = youngwat_pal[1:3],
         title = "''Young Water'' Threshold",
         lwd = 2,
         cex = 0.6)
  if(max(inputdata$Q)*10 > max(inputdata$J)){
    plot(J ~ dates,
         data = inputdata,
         type = "l",
         lwd = 2,
         col = input_pal[1],
         pch = 16,
         ylab = "J, Q, & ET [mm/h]",
         xlab = "Date",
         main = "Model Input: Flows")
    lines(I(Q) ~ dates,
          data = inputdata,
          col = input_pal[3],
          lwd = 2)
    lines(I(ET) ~ dates,
          data = inputdata,
          col = input_pal[5],
          lwd = 2)
  }else{
    plot(J ~ dates,
         data = inputdata,
         type = "l",
         lwd = 2,
         col = input_pal[1],
         pch = 16,
         ylab = "J [mm/h]",
         xlab = "Date",
         main = "Model Input: Flows")
    axis(side = 4,
         at = seq(0, 2, by = 0.5),
         labels = seq(0, 2, by = 0.5)/10)
    mtext("Q & ET [mm/h]",
          side = 4,
          line = 2.5)
    lines(I(Q*10) ~ dates,
          data = inputdata,
          col = input_pal[3],
          lwd = 2)
    lines(I(ET*10) ~ dates,
          data = inputdata,
          col = input_pal[5],
          lwd = 2)
  }
  if(anyNA(inputdata$WI)){
    legend("topright",
           c("J", "Q", "ET"),
           col = input_pal[c(1,3,5)],
           lty = 1,
           cex = 0.6)
  }else{
    lines(WI ~ dates,
          data = inputdata,
          type = "l",
          col = input_pal[6],
          lwd = 2)
    legend("topright",
           c("J", "Q", "ET", ""),
           col = input_pal[c(1,3,5,7)],
           lty = 1,
           cex = 0.6)
  }
  if(save){
    dev.off()
  }
}
