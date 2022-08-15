ROIplot_SW <- function(spq_summary, ROI = 2.15, ROI.width = 0.05,
                       alpha = 0.5, color.code = NULL,
                       incl.plots = c(1,2,3,4)) {
  require(egg)

  tmp1 <- copy(spq_summary$spectra[ppm >= ROI - ROI.width*1.1 & ppm <= ROI + ROI.width*1.1])
  tmp1 <- merge(tmp1, spq_summary$color)

  p1 <- ggplot(tmp1, aes(ppm, value, group = sample, col = color)) +
    scale_x_reverse() +
    coord_cartesian(xlim = c(ROI + ROI.width, ROI - ROI.width),
                    ylim = c(min(tmp1$value), max(tmp1$value))) +
    geom_line(alpha = alpha) +
    facet_grid(.~color) +
    theme_bw() +
    ylab("Signal intensity") +
    ggtitle("Raw spectra") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none")
  if(!is.null(color.code))
    p1 <- p1 + scale_color_manual(values = color.code)


  tmp2 <- copy(spq_summary$peaks[ppm >= ROI - ROI.width*1.1 & ppm <= ROI + ROI.width*1.1])
  tmp2 <- merge(tmp2, spq_summary$color)

  p2 <- ggplot(tmp2, aes(ppm, value, group = sample, col = color)) +
    scale_x_reverse() +
    coord_cartesian(xlim = c(ROI + ROI.width, ROI - ROI.width),
                    ylim = c(min(tmp2$value), max(tmp2$value))) +
    geom_point(alpha = alpha, shape = 4) +
    theme_bw() +
    facet_grid(.~color) +
    ylab("AUC") +
    ggtitle("Peak detection and AUC quantification") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none")
  if(!is.null(color.code))
    p2 <- p2 + scale_color_manual(values = color.code)

  tmp3 <- copy(spq_summary$grouped[ppm >= ROI - ROI.width*1.1 & ppm <= ROI + ROI.width*1.1])
  tmp3 <- merge(tmp3, spq_summary$color)

  p3 <- ggplot(tmp3, aes(ppm, peakValue, group = sample, col = color)) +
    scale_x_reverse() +
    coord_cartesian(xlim = c(ROI + ROI.width, ROI - ROI.width),
                    ylim = c(min(tmp3$peakValue), max(tmp3$peakValue))) +
    geom_point(alpha = alpha, shape = 4) +
    theme_bw() +
    facet_grid(.~color) +
    ylab("AUC") +
    ggtitle("Peak grouping") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none")
  if(!is.null(color.code))
    p3 <- p3 + scale_color_manual(values = color.code)


  tmp4 <- copy(spq_summary$filled[ppm >= ROI - ROI.width*1.1 & ppm <= ROI + ROI.width*1.1])
  tmp4 <- merge(tmp4, spq_summary$color)

  p4 <- ggplot(tmp4, aes(ppm, peakValue, group = sample, col = color)) +
    scale_x_reverse() +
    coord_cartesian(xlim = c(ROI + ROI.width, ROI - ROI.width),
                    ylim = c(min(tmp4$peakValue), max(tmp4$peakValue))) +
    geom_point(alpha = alpha, shape = 4) +
    theme_bw() +
    facet_grid(.~color) +
    ylab("AUC") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "none")
  if(!is.null(color.code))
    p4 <- p4 + scale_color_manual(values = color.code)


  plist <- list(p1,p2,p3,p4)
  plist <- plist[incl.plots]
  
  p_comb <- egg::ggarrange(plots = plist, ncol = 1, draw = F)
  return(p_comb)
}
