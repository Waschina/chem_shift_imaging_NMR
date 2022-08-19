library(data.table)
library(speaq)
library(stringr)
library(ggplot2)
library(ggrepel)

baerchi_colors <- c(red = "#b30000", lightred = "#fc8d59", white = "#fef0d9")

suppressMessages(library(bit64))
source("analysis/v1/scripts/ROIplot_SW.R")
source("analysis/v1/scripts/spectra_analysis.R")

# Let's call it baerchi
baerchi <- list()
baerchi$info <- fread("data/clean/baerchi_SI.csv")
baerchi$info[color == "", color := "white"]
baerchi$info[color == "rot", color := "red"]
baerchi$info[color == "halbrot" | color == "halbrot (eher weiß)",
             color := "lightred"]
baerchi$info[, newID := paste0("baer_",1:.N)]
tmp <- fread("data/raw/baerchi_points_all.gz")
tmp[[1]] <- str_match(tmp[[1]],"Slice-foodomics-Baerchi-600cryo1_\\s*([0-9]{4}?)\\s*_1$")[,2]
tmp_names <- tmp[[1]]
baerchi$spectra <- as.matrix(tmp[,-1])
rownames(baerchi$spectra) <- tmp_names
baerchi$ppm <- as.numeric(colnames(tmp)[-1])
rm(tmp)
rm(tmp_names)

# Normalisation
excl_from_norm <- with(baerchi,which(ppm < 0.4 | (ppm > 4.4 & ppm < 6.01)))
specta_sum <- apply(baerchi$spectra[,-excl_from_norm], 1, function(x) sum(ifelse(x < 0, 0, x)))
baerchi$spectra.norm <- baerchi$spectra / specta_sum
baerchi$spectra.norm <- baerchi$spectra.norm / max(baerchi$spectra.norm) * max(baerchi$spectra)


# Spectra analysis
baer_H <- spectra_analysis(baerchi$spectra.norm[as.character(baerchi$info$exp.H),],
                           baerchi$ppm,
                           baerchi$info$color)

baer_L <- spectra_analysis(baerchi$spectra.norm[as.character(baerchi$info$exp.L),],
                           baerchi$ppm,
                           baerchi$info$color)

baer_G <- spectra_analysis(baerchi$spectra.norm[as.character(baerchi$info$exp.G),],
                           baerchi$ppm,
                           baerchi$info$color)

# plots for response letter
source("analysis/v1/scripts/baerchi_ROI_plots.R")

# PCA
baer_H <- feat.PCA(baer_H, max.na = 0.05)
baer_L <- feat.PCA(baer_L, max.na = 0.05)
baer_G <- feat.PCA(baer_G, max.na = 0.05)


# visual inspection for PC selection (drop axes that explain only a signal outlier)
pdf("analysis/v1/plots/baerchen_PCA_12345.pdf")
pairs(baer_L$PCA$scores[,2:6],
      main = "Bärchen (L)", col = as.factor(baerchi$info$color))
pairs(baer_H$PCA$scores[,2:6],
      main = "Bärchen (H)", col = as.factor(baerchi$info$color))
pairs(baer_G$PCA$scores[,2:6],
      main = "Bärchen (G)", col = as.factor(baerchi$info$color))
dev.off()


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# PCA Plot only for (L) #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
PC_OI <- c("PC1","PC3")

dt.pca.loading <- copy(baer_L$PCA$loadings)
tmp_filled <- baer_L$filled[,.(ppm = median(ppm),
                               peakVal = median(peakValue)),
                            by = peakIndex]
tmp_filled[, peakIndex := as.character(peakIndex)]
dt.pca.loading <- merge(dt.pca.loading, tmp_filled,
                        by.x = "feat", by.y = "peakIndex")
dt.pca.loading[, scale := baer_L$PCA$pca.object$scale[feat]]
n <- 10
#relfeats <-dt.pca.loading[order(sqrt((get(PC_OI[1])*scale)^2 + (get(PC_OI[2])*scale)^2), decreasing = T)][1:n, feat]
relfeats <-dt.pca.loading[order(sqrt((get(PC_OI[1]))^2 + (get(PC_OI[2]))^2), decreasing = T)][1:n, feat]
#tmp_scores <- merge(baer_L$PCA$scores, baerchi$info[, .(sample = as.character(exp.L), color)])
tmp_scores <- baer_L$PCA$scores

p_scores <- ggplot(tmp_scores, aes(get(PC_OI[1]), get(PC_OI[2]), fill = color)) +
  geom_point(shape = 21) +
  scale_fill_manual(values = baerchi_colors) +
  labs(x = paste0(PC_OI[1]," (",
                  round(100*baer_L$PCA$importance["Proportion of Variance",PC_OI[1]]),"%)"),
       y = paste0(PC_OI[2]," (",
                  round(100*baer_L$PCA$importance["Proportion of Variance",PC_OI[2]]),"%)"),
       fill = "Color of sample") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "right")

scaling_fac <- 1
tmpmat <- as.matrix(dt.pca.loading[feat %in% relfeats,.(get(PC_OI[1]),get(PC_OI[2]))])
for(isf in 2:5350) {
  tmptmpmat <- tmpmat * isf
  if(min(tmptmpmat[,1]) < min(tmp_scores[[PC_OI[1]]]) |
     max(tmptmpmat[,1]) > max(tmp_scores[[PC_OI[1]]]) |
     min(tmptmpmat[,2]) < min(tmp_scores[[PC_OI[2]]]) |
     max(tmptmpmat[,2]) > max(tmp_scores[[PC_OI[2]]]))
    break
  scaling_fac <- isf
}

p_loading <- ggplot(dt.pca.loading[feat %in% relfeats],
                    aes(xend = get(PC_OI[1]) * scaling_fac,
                        yend = get(PC_OI[2]) * scaling_fac)) +
  geom_segment(x = 0, y = 0, alpha = 0.4)  +
  geom_text_repel(aes(x=get(PC_OI[1]) * scaling_fac,
                      y=get(PC_OI[2]) * scaling_fac,
                      label = round(ppm, digits = 2)),
                  max.overlaps = 100, show.legend = FALSE,
                  min.segment.length = 100) +
  scale_x_continuous(limits = c(min(tmp_scores[[PC_OI[1]]]),max(tmp_scores[[PC_OI[1]]]))) +
  scale_y_continuous(limits = c(min(tmp_scores[[PC_OI[2]]]),max(tmp_scores[[PC_OI[2]]]))) +
  labs(x = paste0(PC_OI[1]," (",
                  round(100*baer_L$PCA$importance["Proportion of Variance",PC_OI[1]]),"%)"),
       y = paste0(PC_OI[2]," (",
                  round(100*baer_L$PCA$importance["Proportion of Variance",PC_OI[2]]),"%)"),
       col = "Phase") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "right")

P_L_12 <- egg::ggarrange(p_scores, p_loading, nrow = 1)
ggsave("analysis/v1/plots/baerchen_L_PCA.pdf", plot = P_L_12,
       width = 7.5, height = 3.2)

# Bear plot
dt_scaled <- data.table(as.table(baer_L$feat.tab))
dt_scaled <- cbind(dt_scaled[,.(feat = V2, signal = N)], baer_L$color[,.(exp.L = as.integer(sample))])
#setnames(dt_scaled, c("newID","feat","signal"))
dt_scaled <- merge(dt_scaled, baerchi$info, by = "exp.L")
dt_scaled[, pos.y.new := pos.y]
dt_scaled[rack == 2, pos.y.new := pos.y.new + 12]
dt_scaled[, pos.x.new := pos.x]
#dt_scaled[, signal_z := scale(signal), by = "feat"]
dt_scaled[, signal_z := (signal-min(signal))/(max(signal)-min(signal)), by = "feat"]

dt_scaled <- merge(dt_scaled, baer_L$filled[!duplicated(peakIndex),
                                            .(ppm = paste(round(ppm, digits = 2),"ppm"),
                                              feat = as.character(peakIndex))],
                   by = "feat")

p_relfeat_dist <- ggplot(dt_scaled[feat %in% relfeats], aes(x = pos.x.new, y = pos.y.new,
                                                            fill = signal_z, col = signal_z)) +
  geom_tile() +
  facet_wrap("ppm", scales = "free", nrow = 2) +
  labs(x = "Index (x)", y = "Index (y)", fill = "range-scaled\nsignal",
       col = "range-scaled\nsignal") +
  scale_y_continuous(expand = c(0,0), trans = "reverse") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_viridis_c(option = "magma") + scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
ggsave("analysis/v1/plots/baerchen_L_relfeats.pdf", plot = p_relfeat_dist,
       width = 7.5, height = 4.25)

dt_pca_baer <- copy(tmp_scores[,1:10][, exp.L := as.integer(sample)])
dt_pca_baer <- merge(dt_pca_baer, baerchi$info, by = "exp.L")
dt_pca_baer[, pos.y.new := pos.y]
dt_pca_baer[rack == 2, pos.y.new := pos.y.new + 12]
dt_pca_baer[, pos.x.new := pos.x]

p_dist_pc1 <- ggplot(dt_pca_baer, aes(x = pos.x.new, y = pos.y.new,
                                      fill = get(PC_OI[1]), col = get(PC_OI[1]))) +
  geom_tile() +
  labs(x = "Index (x)", y = "Index (y)", title = PC_OI[1],
       fill = PC_OI[1], col = PC_OI[1]) +
  scale_y_continuous(expand = c(0,0), trans = "reverse") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_viridis_c(option = "magma") +
  scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))

p_dist_pc2 <- ggplot(dt_pca_baer, aes(x = pos.x.new, y = pos.y.new,
                                      fill = get(PC_OI[2]), col = get(PC_OI[2]))) +
  geom_tile() +
  labs(x = "Index (x)", y = "Index (y)", title = PC_OI[2],
       fill = PC_OI[2], col = PC_OI[2]) +
  scale_y_continuous(expand = c(0,0), trans = "reverse") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_viridis_c(option = "magma") +
  scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))

p_dist_pc12 <- egg::ggarrange(p_dist_pc1, p_dist_pc2, nrow = 1)
ggsave("analysis/v1/plots/baerchen_L_PC1_PC2.pdf", plot = p_dist_pc12,
       width = 5.25, height = 3.25)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# PCA Plot only for (H) #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
PC_OI <- c("PC1","PC2")

dt.pca.loading <- copy(baer_H$PCA$loadings)
tmp_filled <- baer_H$filled[,.(ppm = median(ppm),
                               peakVal = median(peakValue)),
                            by = peakIndex]
tmp_filled[, peakIndex := as.character(peakIndex)]
dt.pca.loading <- merge(dt.pca.loading, tmp_filled,
                        by.x = "feat", by.y = "peakIndex")
dt.pca.loading[, scale := baer_H$PCA$pca.object$scale[feat]]
n <- 10
#relfeats <-dt.pca.loading[order(sqrt((get(PC_OI[1])*scale)^2 + (get(PC_OI[2])*scale)^2), decreasing = T)][1:n, feat]
relfeats <-dt.pca.loading[order(sqrt((get(PC_OI[1]))^2 + (get(PC_OI[2]))^2), decreasing = T)][1:n, feat]
#tmp_scores <- merge(baer_H$PCA$scores, baerchi$info[, .(sample = as.character(exp.H), color)])
tmp_scores <- baer_H$PCA$scores

p_scores <- ggplot(tmp_scores, aes(get(PC_OI[1]), get(PC_OI[2]), fill = color)) +
  geom_point(shape = 21) +
  scale_fill_manual(values = baerchi_colors) +
  labs(x = paste0(PC_OI[1]," (",
                  round(100*baer_H$PCA$importance["Proportion of Variance",PC_OI[1]]),"%)"),
       y = paste0(PC_OI[2]," (",
                  round(100*baer_H$PCA$importance["Proportion of Variance",PC_OI[2]]),"%)"),
       fill = "Color of sample") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "right")

scaling_fac <- 1
tmpmat <- as.matrix(dt.pca.loading[feat %in% relfeats,.(get(PC_OI[1]),get(PC_OI[2]))])
for(isf in 2:50) {
  tmptmpmat <- tmpmat * isf
  if(min(tmptmpmat[,1]) < min(tmp_scores[[PC_OI[1]]]) |
     max(tmptmpmat[,1]) > max(tmp_scores[[PC_OI[1]]]) |
     min(tmptmpmat[,2]) < min(tmp_scores[[PC_OI[2]]]) |
     max(tmptmpmat[,2]) > max(tmp_scores[[PC_OI[2]]]))
    break
  scaling_fac <- isf
}

p_loading <- ggplot(dt.pca.loading[feat %in% relfeats],
                    aes(xend = get(PC_OI[1]) * scaling_fac,
                        yend = get(PC_OI[2]) * scaling_fac)) +
  geom_segment(x = 0, y = 0, alpha = 0.4)  +
  geom_text_repel(aes(x=get(PC_OI[1]) * scaling_fac,
                      y=get(PC_OI[2]) * scaling_fac,
                      label = round(ppm, digits = 2)),
                  max.overlaps = 100, show.legend = FALSE,
                  min.segment.length = 100) +
  scale_x_continuous(limits = c(min(tmp_scores[[PC_OI[1]]]),max(tmp_scores[[PC_OI[1]]]))) +
  scale_y_continuous(limits = c(min(tmp_scores[[PC_OI[2]]]),max(tmp_scores[[PC_OI[2]]]))) +
  labs(x = paste0(PC_OI[1]," (",
                  round(100*baer_H$PCA$importance["Proportion of Variance",PC_OI[1]]),"%)"),
       y = paste0(PC_OI[2]," (",
                  round(100*baer_H$PCA$importance["Proportion of Variance",PC_OI[2]]),"%)"),
       col = "Phase") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "right")

P_L_12 <- egg::ggarrange(p_scores, p_loading, nrow = 1)
ggsave("analysis/v1/plots/baerchen_H_PCA.pdf", plot = P_L_12,
       width = 7.5, height = 3.2)

# Bear plot
dt_scaled <- data.table(as.table(baer_H$feat.tab))
dt_scaled <- cbind(dt_scaled[,.(feat = V2, signal = N)], baer_H$color[,.(exp.H = as.integer(sample))])
#setnames(dt_scaled, c("newID","feat","signal"))
dt_scaled <- merge(dt_scaled, baerchi$info, by = "exp.H")
dt_scaled[, pos.y.new := pos.y]
dt_scaled[rack == 2, pos.y.new := pos.y.new + 12]
dt_scaled[, pos.x.new := pos.x]
#dt_scaled[, signal_z := scale(signal), by = "feat"]
dt_scaled[, signal_z := (signal-min(signal))/(max(signal)-min(signal)), by = "feat"]
dt_scaled <- merge(dt_scaled, baer_H$filled[!duplicated(peakIndex),
                                            .(ppm = paste(round(ppm, digits = 2),"ppm"),
                                              feat = as.character(peakIndex))],
                   by = "feat")

p_relfeat_dist <- ggplot(dt_scaled[feat %in% relfeats], aes(x = pos.x.new, y = pos.y.new,
                                                            fill = signal_z, col = signal_z)) +
  geom_tile() +
  facet_wrap("ppm", scales = "free", nrow = 2) +
  labs(x = "Index (x)", y = "Index (y)", fill = "range-scaled\nsignal",
       col = "range-scaled\nsignal") +
  scale_y_continuous(expand = c(0,0), trans = "reverse") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_viridis_c(option = "magma") + scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
ggsave("analysis/v1/plots/baerchen_H_relfeats.pdf", plot = p_relfeat_dist,
       width = 7.5, height = 4.25)

dt_pca_baer <- copy(tmp_scores[,1:10][, exp.H := as.integer(sample)])
dt_pca_baer <- merge(dt_pca_baer, baerchi$info, by = "exp.H")
dt_pca_baer[, pos.y.new := pos.y]
dt_pca_baer[rack == 2, pos.y.new := pos.y.new + 12]
dt_pca_baer[, pos.x.new := pos.x]

p_dist_pc1 <- ggplot(dt_pca_baer, aes(x = pos.x.new, y = pos.y.new,
                                      fill = get(PC_OI[1]), col = get(PC_OI[1]))) +
  geom_tile() +
  labs(x = "Index (x)", y = "Index (y)", title = PC_OI[1],
       fill = PC_OI[1], col = PC_OI[1]) +
  scale_y_continuous(expand = c(0,0), trans = "reverse") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_viridis_c(option = "magma") +
  scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))

p_dist_pc2 <- ggplot(dt_pca_baer, aes(x = pos.x.new, y = pos.y.new,
                                      fill = get(PC_OI[2]), col = get(PC_OI[2]))) +
  geom_tile() +
  labs(x = "Index (x)", y = "Index (y)", title = PC_OI[2],
       fill = PC_OI[2], col = PC_OI[2]) +
  scale_y_continuous(expand = c(0,0), trans = "reverse") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_viridis_c(option = "magma") +
  scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))

p_dist_pc12 <- egg::ggarrange(p_dist_pc1, p_dist_pc2, nrow = 1)
ggsave("analysis/v1/plots/baerchen_H_PC1_PC2.pdf", plot = p_dist_pc12,
       width = 5.25, height = 3.25)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# PCA Plot only for (G) #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
PC_OI <- c("PC1","PC2")

dt.pca.loading <- copy(baer_G$PCA$loadings)
tmp_filled <- baer_G$filled[,.(ppm = median(ppm),
                               peakVal = median(peakValue)),
                            by = peakIndex]
tmp_filled[, peakIndex := as.character(peakIndex)]
dt.pca.loading <- merge(dt.pca.loading, tmp_filled,
                        by.x = "feat", by.y = "peakIndex")
dt.pca.loading[, scale := baer_G$PCA$pca.object$scale[feat]]
n <- 10
#relfeats <-dt.pca.loading[order(sqrt((get(PC_OI[1])*scale)^2 + (get(PC_OI[2])*scale)^2), decreasing = T)][1:n, feat]
relfeats <-dt.pca.loading[order(sqrt((get(PC_OI[1]))^2 + (get(PC_OI[2]))^2), decreasing = T)][1:n, feat]
#tmp_scores <- merge(baer_G$PCA$scores, baerchi$info[, .(sample = as.character(exp.G), color)])
tmp_scores <- baer_G$PCA$scores

p_scores <- ggplot(tmp_scores, aes(get(PC_OI[1]), get(PC_OI[2]), fill = color)) +
  geom_point(shape = 21) +
  scale_fill_manual(values = baerchi_colors) +
  labs(x = paste0(PC_OI[1]," (",
                  round(100*baer_G$PCA$importance["Proportion of Variance",PC_OI[1]]),"%)"),
       y = paste0(PC_OI[2]," (",
                  round(100*baer_G$PCA$importance["Proportion of Variance",PC_OI[2]]),"%)"),
       fill = "Color of sample") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "right")

scaling_fac <- 1
tmpmat <- as.matrix(dt.pca.loading[feat %in% relfeats,.(get(PC_OI[1]),get(PC_OI[2]))])
for(isf in 2:50) {
  tmptmpmat <- tmpmat * isf
  if(min(tmptmpmat[,1]) < min(tmp_scores[[PC_OI[1]]]) |
     max(tmptmpmat[,1]) > max(tmp_scores[[PC_OI[1]]]) |
     min(tmptmpmat[,2]) < min(tmp_scores[[PC_OI[2]]]) |
     max(tmptmpmat[,2]) > max(tmp_scores[[PC_OI[2]]]))
    break
  scaling_fac <- isf
}

p_loading <- ggplot(dt.pca.loading[feat %in% relfeats],
                    aes(xend = get(PC_OI[1]) * scaling_fac,
                        yend = get(PC_OI[2]) * scaling_fac)) +
  geom_segment(x = 0, y = 0, alpha = 0.4)  +
  geom_text_repel(aes(x=get(PC_OI[1]) * scaling_fac,
                      y=get(PC_OI[2]) * scaling_fac,
                      label = round(ppm, digits = 2)),
                  max.overlaps = 100, show.legend = FALSE,
                  min.segment.length = 100) +
  scale_x_continuous(limits = c(min(tmp_scores[[PC_OI[1]]]),max(tmp_scores[[PC_OI[1]]]))) +
  scale_y_continuous(limits = c(min(tmp_scores[[PC_OI[2]]]),max(tmp_scores[[PC_OI[2]]]))) +
  labs(x = paste0(PC_OI[1]," (",
                  round(100*baer_G$PCA$importance["Proportion of Variance",PC_OI[1]]),"%)"),
       y = paste0(PC_OI[2]," (",
                  round(100*baer_G$PCA$importance["Proportion of Variance",PC_OI[2]]),"%)"),
       col = "Phase") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "right")

P_L_12 <- egg::ggarrange(p_scores, p_loading, nrow = 1)
ggsave("analysis/v1/plots/baerchen_G_PCA.pdf", plot = P_L_12,
       width = 7.5, height = 3.2)

# Bear plot
dt_scaled <- data.table(as.table(baer_G$feat.tab))
dt_scaled <- cbind(dt_scaled[,.(feat = V2, signal = N)], baer_G$color[,.(exp.G = as.integer(sample))])
#setnames(dt_scaled, c("newID","feat","signal"))
dt_scaled <- merge(dt_scaled, baerchi$info, by = "exp.G")
dt_scaled[, pos.y.new := pos.y]
dt_scaled[rack == 2, pos.y.new := pos.y.new + 12]
dt_scaled[, pos.x.new := pos.x]
#dt_scaled[, signal_z := scale(signal), by = "feat"]
dt_scaled[, signal_z := (signal-min(signal))/(max(signal)-min(signal)), by = "feat"]
dt_scaled <- merge(dt_scaled, baer_G$filled[!duplicated(peakIndex),
                                            .(ppm = paste(round(ppm, digits = 2),"ppm"),
                                              feat = as.character(peakIndex))],
                   by = "feat")

p_relfeat_dist <- ggplot(dt_scaled[feat %in% relfeats], aes(x = pos.x.new, y = pos.y.new,
                                                            fill = signal_z, col = signal_z)) +
  geom_tile() +
  facet_wrap("ppm", scales = "free", nrow = 2) +
  labs(x = "Index (x)", y = "Index (y)", fill = "range-scaled\nsignal",
       col = "range-scaled\nsignal") +
  scale_y_continuous(expand = c(0,0), trans = "reverse") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_viridis_c(option = "magma") + scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
ggsave("analysis/v1/plots/baerchen_G_relfeats.pdf", plot = p_relfeat_dist,
       width = 7.5, height = 4.25)

dt_pca_baer <- copy(tmp_scores[,1:10][, exp.G := as.integer(sample)])
dt_pca_baer <- merge(dt_pca_baer, baerchi$info, by = "exp.G")
dt_pca_baer[, pos.y.new := pos.y]
dt_pca_baer[rack == 2, pos.y.new := pos.y.new + 12]
dt_pca_baer[, pos.x.new := pos.x]

p_dist_pc1 <- ggplot(dt_pca_baer, aes(x = pos.x.new, y = pos.y.new,
                                      fill = get(PC_OI[1]), col = get(PC_OI[1]))) +
  geom_tile() +
  labs(x = "Index (x)", y = "Index (y)", title = PC_OI[1],
       fill = PC_OI[1], col = PC_OI[1]) +
  scale_y_continuous(expand = c(0,0), trans = "reverse") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_viridis_c(option = "magma") +
  scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))

p_dist_pc2 <- ggplot(dt_pca_baer, aes(x = pos.x.new, y = pos.y.new,
                                      fill = get(PC_OI[2]), col = get(PC_OI[2]))) +
  geom_tile() +
  labs(x = "Index (x)", y = "Index (y)", title = PC_OI[2],
       fill = PC_OI[2], col = PC_OI[2]) +
  scale_y_continuous(expand = c(0,0), trans = "reverse") +
  scale_x_discrete(expand = c(0,0)) +
  scale_fill_viridis_c(option = "magma") +
  scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))

p_dist_pc12 <- egg::ggarrange(p_dist_pc1, p_dist_pc2, nrow = 1)
ggsave("analysis/v1/plots/baerchen_G_PC1_PC2.pdf", plot = p_dist_pc12,
       width = 5.25, height = 3.25)


#–––––––––––––––––––––––––––#
# ICA Plots for supplements #
#–––––––––––––––––––––––––––#
source("analysis/v1/scripts/baerchi_01_ICA.R")






# ROI.ppm <- 3.03
# roiWidth.ppm <- 0.05
# 
# ROIplot_SW(baer_H, ROI = ROI.ppm, ROI.width = roiWidth.ppm,
#            color.code = MetBrewer::met.brewer("Isfahan2",3))
# 
# any(is.na(baer_H$feat.tab[1:10,1:10]))
# 
# erni <- dcast(baer_L$grouped[peakSNR >= 2], sample ~ peakIndex, value.var = "peakValue")
# tmpnames <- rownames(erni)
# erni <- as.matrix(erni[,-1])
# rownames(erni) <- tmpnames
# erni <- erni[,apply(erni,
#                     2,
#                     function(x) sum(is.na(x))) <= nrow(erni) * 0.05]
# # imputation with kNN (only features, which occur in at least 97.5% of samples)
# #peaks.features[which(peaks.features == 0, arr.ind = TRUE)] <- NA
# erni_ci <- apply(erni, 2, quantile, na.rm = TRUE, prob = c(0.025, 0.975))
# up_rm <- which(t(erni) > erni_ci[2,] + (erni_ci[2,] - erni_ci[1,]), arr.ind = TRUE)
# dw_rm <- which(t(erni) < erni_ci[1,] - (erni_ci[2,] - erni_ci[1,]), arr.ind = TRUE)
# 
# erni[up_rm[,c(2,1)]] <- NA
# erni[dw_rm[,c(2,1)]] <- NA
# 
# baer_L$grouped[peakIndex %in% rownames(up_rm)][order(-peakValue)][!duplicated(peakIndex)]
# 
# 
# erni2 <- missForest::missForest(erni)$ximp

