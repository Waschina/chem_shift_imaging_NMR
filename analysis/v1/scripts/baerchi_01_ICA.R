
comp_comb <- combn(names(baer_H$ICA$importance),2)

#––––––––––––––––––#
# Bierschinken (H) #
#––––––––––––––––––#
p_tmp <- list()
for(i in 1:ncol(comp_comb)) {
  # IC_OI <- c("IC1","IC2")
  IC_OI <- comp_comb[,i]
  
  dt.ica.loading <- copy(baer_H$ICA$loadings)
  tmp_filled <- baer_H$filled[,.(ppm = median(ppm),
                                 peakVal = median(peakValue)),
                              by = peakIndex]
  tmp_filled[, peakIndex := as.character(peakIndex)]
  dt.ica.loading <- merge(dt.ica.loading, tmp_filled,
                          by.x = "feat", by.y = "peakIndex")
  dt.ica.loading[, scale := baer_H$PCA$pca.object$scale[feat]]
  n <- 10
  #relfeats <-dt.ica.loading[order(sqrt((get(IC_OI[1])*scale)^2 + (get(IC_OI[2])*scale)^2), decreasing = T)][1:n, feat]
  relfeats <-dt.ica.loading[order(sqrt((get(IC_OI[1]))^2 + (get(IC_OI[2]))^2), decreasing = T)][1:n, feat]
  #tmp_scores <- merge(baer_H$ICA$scores, baerchi$info[, .(sample = as.character(exp.H), color)])
  tmp_scores <- baer_H$ICA$scores
  
  p_scores <- ggplot(tmp_scores, aes(get(IC_OI[1]), get(IC_OI[2]), fill = color)) +
    geom_point(shape = 21) +
    scale_fill_manual(values = baerchi_colors) +
    labs(x = paste0(IC_OI[1]," (",
                    round(100*baer_H$ICA$importance[IC_OI[1]]),"%)"),
         y = paste0(IC_OI[2]," (",
                    round(100*baer_H$ICA$importance[IC_OI[2]]),"%)"),
         fill = "Color of sample") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          legend.position = "right")
  
  scaling_fac <- 1
  tmpmat <- as.matrix(dt.ica.loading[feat %in% relfeats,.(get(IC_OI[1]),get(IC_OI[2]))])
  for(isf in 2:5000) {
    tmptmpmat <- tmpmat * isf
    if(min(tmptmpmat[,1]) < min(tmp_scores[[IC_OI[1]]]) |
       max(tmptmpmat[,1]) > max(tmp_scores[[IC_OI[1]]]) |
       min(tmptmpmat[,2]) < min(tmp_scores[[IC_OI[2]]]) |
       max(tmptmpmat[,2]) > max(tmp_scores[[IC_OI[2]]]))
      break
    scaling_fac <- isf
  }
  
  p_loading <- ggplot(dt.ica.loading[feat %in% relfeats],
                      aes(xend = get(IC_OI[1]) * scaling_fac,
                          yend = get(IC_OI[2]) * scaling_fac)) +
    geom_segment(x = 0, y = 0, alpha = 0.4)  +
    geom_text_repel(aes(x=get(IC_OI[1]) * scaling_fac,
                        y=get(IC_OI[2]) * scaling_fac,
                        label = round(ppm, digits = 2)),
                    max.overlaps = 100, show.legend = FALSE,
                    min.segment.length = 100) +
    scale_x_continuous(limits = c(min(tmp_scores[[IC_OI[1]]]),max(tmp_scores[[IC_OI[1]]]))) +
    scale_y_continuous(limits = c(min(tmp_scores[[IC_OI[2]]]),max(tmp_scores[[IC_OI[2]]]))) +
    labs(x = paste0(IC_OI[1]," (",
                    round(100*baer_H$ICA$importance[IC_OI[1]]),"%)"),
         y = paste0(IC_OI[2]," (",
                    round(100*baer_H$ICA$importance[IC_OI[2]]),"%)"),
         col = "Phase") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          legend.position = "right")
  
  p_tmp[[i]] <- egg::ggarrange(p_scores, p_loading, nrow = 1, draw = FALSE)
}

p_out <- ggpubr::ggarrange(plotlist = p_tmp, ncol = 1)
ggsave("analysis/v1/plots/baerchen_H_ICA.pdf", plot = p_out,
       width = 7, height = 8)

pdf("analysis/v1/plots/baerchen_H_PCA_elbow.pdf", width = 6.5, height = 4)
plot(baer_H$PCA$importance["Proportion of Variance",1:20], type = "b",
     xlab = "Principal component index", ylab = "Proportion of Variance")
dev.off()

#––––––––––––––––––#
# Bierschinken (L) #
#––––––––––––––––––#
p_tmp <- list()
for(i in 1:ncol(comp_comb)) {
  # IC_OI <- c("IC1","IC2")
  IC_OI <- comp_comb[,i]
  
  dt.ica.loading <- copy(baer_L$ICA$loadings)
  tmp_filled <- baer_L$filled[,.(ppm = median(ppm),
                                 peakVal = median(peakValue)),
                              by = peakIndex]
  tmp_filled[, peakIndex := as.character(peakIndex)]
  dt.ica.loading <- merge(dt.ica.loading, tmp_filled,
                          by.x = "feat", by.y = "peakIndex")
  dt.ica.loading[, scale := baer_L$PCA$pca.object$scale[feat]]
  n <- 10
  #relfeats <-dt.ica.loading[order(sqrt((get(IC_OI[1])*scale)^2 + (get(IC_OI[2])*scale)^2), decreasing = T)][1:n, feat]
  relfeats <-dt.ica.loading[order(sqrt((get(IC_OI[1]))^2 + (get(IC_OI[2]))^2), decreasing = T)][1:n, feat]
  #tmp_scores <- merge(baer_L$ICA$scores, baerchi$info[, .(sample = as.character(exp.L), color)])
  tmp_scores <- baer_L$ICA$scores
  
  p_scores <- ggplot(tmp_scores, aes(get(IC_OI[1]), get(IC_OI[2]), fill = color)) +
    geom_point(shape = 21) +
    scale_fill_manual(values = baerchi_colors) +
    labs(x = paste0(IC_OI[1]," (",
                    round(100*baer_L$ICA$importance[IC_OI[1]]),"%)"),
         y = paste0(IC_OI[2]," (",
                    round(100*baer_L$ICA$importance[IC_OI[2]]),"%)"),
         fill = "Color of sample") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          legend.position = "right")
  
  scaling_fac <- 1
  tmpmat <- as.matrix(dt.ica.loading[feat %in% relfeats,.(get(IC_OI[1]),get(IC_OI[2]))])
  for(isf in 2:5000) {
    tmptmpmat <- tmpmat * isf
    if(min(tmptmpmat[,1]) < min(tmp_scores[[IC_OI[1]]]) |
       max(tmptmpmat[,1]) > max(tmp_scores[[IC_OI[1]]]) |
       min(tmptmpmat[,2]) < min(tmp_scores[[IC_OI[2]]]) |
       max(tmptmpmat[,2]) > max(tmp_scores[[IC_OI[2]]]))
      break
    scaling_fac <- isf
  }
  
  p_loading <- ggplot(dt.ica.loading[feat %in% relfeats],
                      aes(xend = get(IC_OI[1]) * scaling_fac,
                          yend = get(IC_OI[2]) * scaling_fac)) +
    geom_segment(x = 0, y = 0, alpha = 0.4)  +
    geom_text_repel(aes(x=get(IC_OI[1]) * scaling_fac,
                        y=get(IC_OI[2]) * scaling_fac,
                        label = round(ppm, digits = 2)),
                    max.overlaps = 100, show.legend = FALSE,
                    min.segment.length = 100) +
    scale_x_continuous(limits = c(min(tmp_scores[[IC_OI[1]]]),max(tmp_scores[[IC_OI[1]]]))) +
    scale_y_continuous(limits = c(min(tmp_scores[[IC_OI[2]]]),max(tmp_scores[[IC_OI[2]]]))) +
    labs(x = paste0(IC_OI[1]," (",
                    round(100*baer_L$ICA$importance[IC_OI[1]]),"%)"),
         y = paste0(IC_OI[2]," (",
                    round(100*baer_L$ICA$importance[IC_OI[2]]),"%)"),
         col = "Phase") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          legend.position = "right")
  
  p_tmp[[i]] <- egg::ggarrange(p_scores, p_loading, nrow = 1, draw = FALSE)
}

p_out <- ggpubr::ggarrange(plotlist = p_tmp, ncol = 1)
ggsave("analysis/v1/plots/baerchen_L_ICA.pdf", plot = p_out,
       width = 7, height = 8)

pdf("analysis/v1/plots/baerchen_L_PCA_elbow.pdf", width = 6.5, height = 4)
plot(baer_L$PCA$importance["Proportion of Variance",1:20], type = "b",
     xlab = "Principal component index", ylab = "Proportion of Variance")
dev.off()

