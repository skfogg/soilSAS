#' Plots soil SAS inputs
#'
#' @description
#' Produces plots of the inputs to the SAS model
#'
#' @param model_parameters a list of SAS model input parameters
#' @param save a logical of whether or not you want the plot to be saved to disk as a .png
#' @param file_path used only if save = TRUE. The file path (including plot name) to be saved
#'


plot_SAS_input <- function(model_parameters, save = FALSE, file_path = NULL){

  require(wesanderson)

  pal <- wes_palette("AsteroidCity2")

  if(save){
    png(paste0(file_path, ".png"),
        height = 1000*3,
        width = 600*3,
        res = 72*3)
  }
  par(mar = c(4,5,3,4),
      mfrow = c(3,1),
      cex.axis = 1.5,
      cex.lab = 1.5,
      cex.main = 1.8)

  nf <- layout(matrix(c(1,2,
                        3,3,
                        4,4,
                        5,5), 4, 2, byrow = TRUE), respect = TRUE)
  layout.show(nf)

  ps <- seq(0,1,0.01)
  big_omega_q <- do.call(model_params$SAS_Q_shape,
                         list(ps, model_params$k_Q))
  big_omega_et <- do.call(model_params$SAS_ET_shape,
                          list(ps, model_params$k_ET))

  plot(ps,
       big_omega_q,
       type = "l",
       ylab = "CDF [-]",
       xlab = "normalized rank storage [-]",
       col = pal[3])
  lines(ps,
        big_omega_et,
        col = pal[5])
  legend("bottomright",
         legend = c(bquote(Omega[Q] == .(model_params$k_Q)),
                    bquote(Omega[ET] == .(model_params$k_ET))),
         lty = 1,
         col = pal[c(3,5)])

  lil_omega_q <- diff(big_omega_q)/diff(ps)
  lil_omega_et <- diff(big_omega_et)/diff(ps)

  plot(ps[2:length(ps)],
       lil_omega_q,
       type = "l",
       ylim = c(0,max(lil_omega_et, lil_omega_q)),
       ylab = "PDF [-]",
       xlab = "normalized rank storage [-]",
       col = pal[3])
  lines(ps[2:length(ps)],
        lil_omega_et,
        col = pal[5])
  legend("topright",
         legend = c(bquote(omega[Q] == .(model_params$k_Q)),
                    bquote(omega[ET] == .(model_params$k_ET))),
         lty = 1,
         col = pal[c(3,5)])

  if(max(model_params$inputdata$Q)*10 > max(model_params$inputdata$J)){
    plot(J ~ dates,
         data = model_params$inputdata,
         type = "o",
         lwd = 2,
         col = pal[1],
         pch = 16,
         ylab = "J, Q, & ET [mm/h]",
         xlab = "Date",
         main = "Model Input: Flows")
    lines(I(Q) ~ dates,
          data = model_params$inputdata,
          col = pal[3],
          lwd = 2,
          type = "o")
    lines(I(ET) ~ dates,
          data = model_params$inputdata,
          col = pal[5],
          lwd = 2,
          type = "o")
  }else{
    plot(J ~ dates,
         data = model_params$inputdata,
         type = "o",
         lwd = 2,
         col = pal[1],
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
          data = model_params$inputdata,
          col = pal[3],
          lwd = 2,
          type = "o")
    lines(I(ET*10) ~ dates,
          data = model_params$inputdata,
          col = pal[5],
          lwd = 2,
          type = "o")
  }
  if(anyNA(model_params$inputdata$wi)){
    legend("topright",
           c("J", "Q", "ET"),
           col = pal[c(1,3,5)],
           pch = 16,
           lty = 1,
           cex = 1.3)
  }else{
    lines(wi ~ dates,
          data = dates,
          type = "o",
          col = pal[6],
          pch = 16)
    legend("topright",
           c("J", "Q", "ET", "wi"),
           col = pal[c(1,3,5,7)],
           pch = 16,
           lty = 1,
           cex = 1.3)
  }
  abline(v = model_params$age_distributions,
         lty = 4,
         col = hcl.colors(length(model_params$age_distributions), "Purp"),
         lwd = 3)
  legend("right",
         as.character(paste0(month(model_params$age_distributions,
                                   label = T,
                                   abbr= T),
                             "-",
                             day(model_params$age_distributions))),
         lty = 4,
         lwd = 3,
         col = hcl.colors(length(model_params$age_distributions), "Purp"),
         title = "Water Ages Output:",
         cex = 1)


  plot(C_J ~ dates,
       data = model_params$inputdata,
       type = "o",
       col = pal[2],
       pch = 16,
       ylab = "C [mg/L]",
       xlab = "Date",
       main = "Model Input: Rain Concentration")
  abline(v = model_params$age_distributions,
         lty = 4,
         col = hcl.colors(length(model_params$age_distributions), "Purp"),
         lwd = 3)
  legend("right",
         as.character(paste0(month(model_params$age_distributions,
                                   label = T,
                                   abbr= T),
                             "-",
                             day(model_params$age_distributions))),
         lty = 4,
         lwd = 3,
         col = hcl.colors(length(model_params$age_distributions), "Purp"),
         title = "Water Ages Output:",
         cex = 1)

  ### STORAGE CHANGES
  delta_S <- model_params$inputdata$J - model_params$inputdata$Q - model_params$inputdata$ET

  storage <- numeric(length(model_params$inputdata$J))
  storage[1] <- model_params$S0 + delta_S[1]

  for(i in 2:length(model_params$inputdata$J)){
    storage[i] <- storage[i-1] + delta_S[i]
  }
  plot(model_params$inputdata$dates,
       storage,
       ylab = "S [mm]",
       xlab = "Date",
       col = pal[4],
       main = "Storage Change")
  abline(v = model_params$age_distributions,
         lty = 4,
         col = hcl.colors(length(model_params$age_distributions), "Purp"),
         lwd = 3)
  legend("right",
         as.character(paste0(month(model_params$age_distributions,
                                   label = T,
                                   abbr= T),
                             "-",
                             day(model_params$age_distributions))),
         lty = 4,
         lwd = 3,
         col = hcl.colors(length(model_params$age_distributions), "Purp"),
         title = "Water Ages Output:",
         cex = 1)
  if(save){
    dev.off()
  }

}
