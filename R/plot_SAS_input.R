#' standardized SAS input plot
#'
#' @description
#' Produces plots of the inputs to the SAS model.
#'
#' @param model_parameters a list of SAS model input parameters
#' @param save a logical of whether or not you want the plot to be saved to disk as a .png
#' @param file_path used only if save = TRUE. The file path (including plot name) to be saved
#'
#' @details
#' The input plots included are: (1) the CDFs of the transit time distributions for Q and ET,
#' (2) the PDFs of the transit time distribution for Q and ET, (3) Model flow inputs (including J, Q, and ET),
#' (4) concentration input (C), and (5) the change in storage (S) of the system.
#'
#' @seealso plot_sas_output
#'
#' @export
#' @import wesanderson
#' @importFrom grDevices dev.off hcl.colors png
#' @importFrom graphics abline axis legend lines mtext par segments
#'
#'


plot_SAS_input <- function(model_parameters, save = FALSE, file_path = NULL){

  requireNamespace("wesanderson")

  pal <- wesanderson::wes_palette("AsteroidCity2")

  # dev.off()
  if(save){
    grDevices::png(paste0(file_path, ".png"),
        height = 1000*3,
        width = 600*3,
        res = 72*3)
  }
  # }else{
  #   if(Sys.info()["sysname"] == "Windows"){
  #     windows(height = 1000,
  #             width = 700)
  #   }else if(Sys.info()["sysname"] == "Linux"){
  #     x11(height = 1000,
  #         width = 700)
  #   }else{
  #     quartz(height = 1000,
  #            width = 700)
  #   }
  # }
  graphics::par(mar = c(4,5,3,4),
      mfrow = c(3,1),
      cex.axis = 1.5,
      cex.lab = 1.5,
      cex.main = 1.8)

  nf <- graphics::layout(matrix(c(1,2,
                        3,3,
                        4,4,
                        5,5), 4, 2, byrow = TRUE), respect = TRUE)
  # layout.show(nf)

  ps <- seq(0,1,0.01)
  big_omega_q <- do.call(model_parameters$SAS_Q_shape,
                         list(ps, model_parameters$shape_params_Q))
  big_omega_et <- do.call(model_parameters$SAS_ET_shape,
                          list(ps, model_parameters$shape_params_ET))
  # CDF
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
         legend = c(bquote(Omega[Q] == .(model_parameters$shape_params_Q)),
                    bquote(Omega[ET] == .(model_parameters$shape_params_ET))),
         lty = 1,
         col = pal[c(3,5)])

  lil_omega_q <- diff(big_omega_q)/diff(ps)
  lil_omega_et <- diff(big_omega_et)/diff(ps)

  # PDF
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
         legend = c(bquote(omega[Q] == .(model_parameters$shape_params_Q)),
                    bquote(omega[ET] == .(model_parameters$shape_params_ET))),
         lty = 1,
         col = pal[c(3,5)])

  if(max(model_parameters$inputdata$Q)*10 > max(model_parameters$inputdata$J)){
    plot(J ~ dates,
         data = model_parameters$inputdata,
         type = "o",
         lwd = 2,
         col = pal[1],
         pch = 16,
         ylab = "J, Q, & ET [mm/h]",
         xlab = "Date",
         main = "Model Input: Flows")
    lines(I(Q) ~ dates,
          data = model_parameters$inputdata,
          col = pal[3],
          lwd = 2,
          type = "o")
    lines(I(ET) ~ dates,
          data = model_parameters$inputdata,
          col = pal[5],
          lwd = 2,
          type = "o")
  }else{
    plot(J ~ dates,
         data = model_parameters$inputdata,
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
          data = model_parameters$inputdata,
          col = pal[3],
          lwd = 2,
          type = "o")
    lines(I(ET*10) ~ dates,
          data = model_parameters$inputdata,
          col = pal[5],
          lwd = 2,
          type = "o")
  }
  if(anyNA(model_parameters$inputdata$WI)){
    legend("topright",
           c("J", "Q", "ET"),
           col = pal[c(1,3,5)],
           pch = 16,
           lty = 1,
           cex = 1.3)
  }else{
    lines(WI ~ dates,
          data = model_parameters$inputdata,
          type = "o",
          col = pal[6],
          pch = 16)
    legend("topright",
           c("J", "Q", "ET", "WI"),
           col = pal[c(1,3,5,7)],
           pch = 16,
           lty = 1,
           cex = 1.3)
  }
  abline(v = model_parameters$age_distributions,
         lty = 4,
         col = hcl.colors(length(model_parameters$age_distributions), "Purp"),
         lwd = 3)
  legend("right",
         as.character(paste0(lubridate::month(model_parameters$age_distributions,
                                              label = T,
                                              abbr= T),
                             "-",
                             lubridate::day(model_parameters$age_distributions))),
         lty = 4,
         lwd = 3,
         col = hcl.colors(length(model_parameters$age_distributions), "Purp"),
         title = "Water Ages Output:",
         cex = 1)

  # Input flows
  plot(C_J ~ dates,
       data = model_parameters$inputdata,
       type = "o",
       col = pal[2],
       pch = 16,
       ylab = "C [mg/L]",
       xlab = "Date",
       main = "Model Input: Rain Concentration")
  abline(v = model_parameters$age_distributions,
         lty = 4,
         col = hcl.colors(length(model_parameters$age_distributions), "Purp"),
         lwd = 3)
  legend("right",
         as.character(paste0(lubridate::month(model_parameters$age_distributions,
                                   label = T,
                                   abbr= T),
                             "-",
                             lubridate::day(model_parameters$age_distributions))),
         lty = 4,
         lwd = 3,
         col = hcl.colors(length(model_parameters$age_distributions), "Purp"),
         title = "Water Ages Output:",
         cex = 1)

  ### STORAGE CHANGES
  delta_S <- model_parameters$inputdata$J - model_parameters$inputdata$Q - model_parameters$inputdata$ET

  storage <- numeric(length(model_parameters$inputdata$J))
  storage[1] <- model_parameters$S0 + delta_S[1]

  for(i in 2:length(model_parameters$inputdata$J)){
    storage[i] <- storage[i-1] + delta_S[i]
  }
  plot(model_parameters$inputdata$dates,
       storage,
       ylab = "S [mm]",
       xlab = "Date",
       col = pal[4],
       main = "Storage Change")
  abline(v = model_parameters$age_distributions,
         lty = 4,
         col = hcl.colors(length(model_parameters$age_distributions), "Purp"),
         lwd = 3)
  legend("right",
         as.character(paste0(lubridate::month(model_parameters$age_distributions,
                                   label = T,
                                   abbr= T),
                             "-",
                             lubridate::day(model_parameters$age_distributions))),
         lty = 4,
         lwd = 3,
         col = hcl.colors(length(model_parameters$age_distributions), "Purp"),
         title = "Water Ages Output:",
         cex = 1)

    dev.off()


}
