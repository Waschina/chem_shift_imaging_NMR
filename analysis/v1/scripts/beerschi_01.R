library(data.table)
library(speaq)
library(stringr)
library(ggplot2)
library(ggrepel)

beerschi_colors <- c(w = "#f7f7f7", wg = "#cccccc", g = "#969696",
                     gr = "#df65b0", r = "#980043", wr = "#f1eef6")
beerschi_colors2 <- c(w = "#f1eef6", r = "#980043", g = "#df65b0")

suppressMessages(library(bit64))
source("analysis/v1/scripts/ROIplot_SW.R")
source("analysis/v1/scripts/spectra_analysis.R")

# Let's call it beerschi
beerschi <- list()
beerschi$info <- fread("data/clean/bierschinken_SI.csv")
beerschi$info[color == "", color := "w"]
beerschi$info[, color2 := color]
beerschi$info[color2 == "gr", color2 := "r"]
beerschi$info[color2 == "wg", color2 := "w"]
beerschi$info[rack == 1 & pos.x == "F" & pos.y == 4, color2 := "g"]
beerschi$info[rack == 1 & pos.x == "D" & pos.y == 4, color2 := "g"]
beerschi$info[color2 == "wr", color2 := "r"]
beerschi$info[, newID := paste0("beer_",1:.N)]
tmp <- fread("data/raw/bierschinken_pointwise_all_v2")
tmp[[1]] <- str_match(tmp[[1]],"Slice-foodomics-Bierschinken-600cryo1_\\s*([0-9]{4}?)\\s*_1$")[,2]
tmp_names <- tmp[[1]]
beerschi$spectra <- as.matrix(tmp[,-1])
rownames(beerschi$spectra) <- tmp_names
beerschi$ppm <- as.numeric(colnames(tmp)[-1])
rm(tmp)
rm(tmp_names)

# Normalisation
excl_from_norm <- with(beerschi,which(ppm < 0.4 | (ppm > 4.4 & ppm < 6.01)))
specta_sum <- apply(beerschi$spectra[,-excl_from_norm], 1, function(x) sum(ifelse(x < 0, 0, x)))
beerschi$spectra.norm <- beerschi$spectra / specta_sum
beerschi$spectra.norm <- beerschi$spectra.norm / max(beerschi$spectra.norm) * max(beerschi$spectra)


beer_H <- spectra_analysis(beerschi$spectra.norm[as.character(beerschi$info$exp.H),],
                           beerschi$ppm,
                           beerschi$info$color)

beer_L <- spectra_analysis(beerschi$spectra.norm[as.character(beerschi$info$exp.L),],
                           beerschi$ppm,
                           beerschi$info$color)

beer_G <- spectra_analysis(beerschi$spectra.norm[as.character(beerschi$info$exp.G),],
                           beerschi$ppm,
                           beerschi$info$color)

# PCA
beer_H <- feat.PCA(beer_H, max.na = 0.05)
beer_L <- feat.PCA(beer_L, max.na = 0.05)
beer_G <- feat.PCA(beer_G, max.na = 0.05)


# visual inspection for PC selection (drop axes that explain only a signal outlier)
pdf("analysis/v1/plots/bierschinken_PCA_12345.pdf")
pairs(beer_L$PCA$scores[,2:6],
      main = "Bierschinken (L)", col = as.factor(beerschi$info$color2))
pairs(beer_H$PCA$scores[,2:6],
      main = "Bierschinken (H)", col = as.factor(beerschi$info$color2))
pairs(beer_G$PCA$scores[,2:6],
      main = "Bierschinken (G)", col = as.factor(beerschi$info$color2))
dev.off()

# # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# # Combined PCA axis for L+H #
# # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# beer_HL <- list()
# beer_HL$color <- cbind(beer_H$color[,.(sampleH = sample, color)],
#                        beer_L$color[,.(sampleL = sample)],
#                        beerschi$info[,.(newID,rack,color2)])
# beer_HL$filled <- rbind(merge(beer_H$filled, beer_HL$color[,.(sample = sampleH,
#                                                               newID)],
#                               by = "sample"),
#                         merge(beer_L$filled, beer_HL$color[,.(sample = sampleL,
#                                                               newID)],
#                               by = "sample"))
# beer_HL$filled$mode <- c(rep("H",nrow(beer_H$filled)),
#                          rep("L",nrow(beer_L$filled)))
# beer_HL$filled[, peakIndexMode := paste0(mode,"-", peakIndex)]
# tmpH <- beer_H$feat.tab; colnames(tmpH) <- paste0("H-",colnames(tmpH))
# tmpL <- beer_L$feat.tab; colnames(tmpL) <- paste0("L-",colnames(tmpL))
# beer_HL$feat.tab <- cbind(tmpH,tmpL)
# rm(tmpH); rm(tmpL)
# tmp_peaks.features <- dcast(beer_HL$filled, newID ~ peakIndexMode,
#                             value.var = "peakValue", fill = 0)
# setkey(tmp_peaks.features, "newID")
# tmp_peaks.features <- tmp_peaks.features[beer_HL$color$newID]
# tmp_peaks.features <- as.matrix(tmp_peaks.features[,-1])
# rownames(tmp_peaks.features) <- beer_HL$color$newID
# # hier evtl features rausfiltern, die selten >LOD sind
# tmp_peaks.features <- tmp_peaks.features[,apply(tmp_peaks.features,
#                                                 2,
#                                                 function(x) sum(x==0)) <= nrow(tmp_peaks.features) * 0.025]
# # imputation with kNN (only features, which occur in at least 97.5% of samples)
# tmp_peaks.features[which(tmp_peaks.features == 0, arr.ind = TRUE)] <- NA
# tmp_peaks.features <- t(bnstruct::knn.impute(t(tmp_peaks.features)))
#
# # exclude features within specific ranges
# keep_feat <- beer_HL$filled[!(ppm < 0.4 | (ppm > 4.41 & ppm < 6)), unique(peakIndexMode)]
# keep_feat <- keep_feat[keep_feat %in% colnames(tmp_peaks.features)]
# tmp_ft.tab.scaled <- speaq::SCANT(data.matrix = tmp_peaks.features[,keep_feat],
#                                   type = c("center","unit"))
# # the actual PCA
# ft.pca <- prcomp(tmp_ft.tab.scaled, scale. = FALSE, center = FALSE)
# dt.pca <- cbind(data.table(newID = rownames(tmp_ft.tab.scaled)),
#                 ft.pca$x)
# dt.pca <- merge(dt.pca, beer_HL$color)
# dt.pca.loading <- cbind(data.table(feat = rownames(ft.pca$rotation)),
#                         ft.pca$rotation)
# dt.pca.loading <- merge(dt.pca.loading, beer_HL$filled[,.(ppm = median(ppm),
#                                                           peakVal = median(peakValue)),
#                                                        by = peakIndexMode],
#                         by.x = "feat", by.y = "peakIndexMode")
# dt.pca.loading[, scale := ft.pca$scale[feat]]
# importance.pca <- summary(ft.pca)$importance
#
# pairs(dt.pca[,2:6],
#       col = as.factor(dt.pca$color), main = "Bierschinken (HL)")
# pairs(dt.pca[,2:6],
#       col = as.factor(dt.pca$rack), main = "Bierschinken (HL)")
#
# plot_pca_feat <- function(PC_OI, legend.position = "none") {
#   # identify top n features with max sqrt(axis1^2 + axis2^2) loadings
#   n <- 10
#   dt.pca.loading[order(abs(get(PC_OI[1])), decreasing = T)][1:n, feat]
#   # relfeats <- c(dt.pca.loading[order(abs(get(PC_OI[1])), decreasing = T)][1:n, feat],
#   #               dt.pca.loading[order(abs(get(PC_OI[2])), decreasing = T)][1:n, feat])
#
#   relfeats <- dt.pca.loading[order(sqrt((get(PC_OI[1])*scale)^2 + (get(PC_OI[2])*scale)^2), decreasing = T)][1:n, feat]
#
#   p_scores <- ggplot(dt.pca, aes(get(PC_OI[1]), get(PC_OI[2]), fill = color2)) +
#     geom_point(shape = 21) +
#     scale_fill_manual(values = beerschi_colors2) +
#     labs(x = paste0(PC_OI[1]," (",
#                     round(100*importance.pca["Proportion of Variance",PC_OI[1]]),"%)"),
#          y = paste0(PC_OI[2]," (",
#                     round(100*importance.pca["Proportion of Variance",PC_OI[2]]),"%)"),
#          fill = "Color of sample") +
#     theme_bw() +
#     theme(panel.grid = element_blank(),
#           axis.text = element_text(color = "black"),
#           legend.position = legend.position)
#
#   scaling_fac <- 1
#   tmpmat <- as.matrix(dt.pca.loading[feat %in% relfeats,.(get(PC_OI[1]),get(PC_OI[2]))])
#   for(isf in 2:50) {
#     tmptmpmat <- tmpmat * isf
#     if(min(tmptmpmat[,1]) < min(dt.pca[[PC_OI[1]]]) |
#        max(tmptmpmat[,1]) > max(dt.pca[[PC_OI[1]]]) |
#        min(tmptmpmat[,2]) < min(dt.pca[[PC_OI[2]]]) |
#        max(tmptmpmat[,2]) > max(dt.pca[[PC_OI[2]]]))
#       break
#     scaling_fac <- isf
#   }
#
#   p_loading <- ggplot(dt.pca.loading[feat %in% relfeats],
#                       aes(xend = get(PC_OI[1]) * scaling_fac,
#                           yend = get(PC_OI[2]) * scaling_fac,
#                           col = ifelse(grepl("H", feat),"H","L"))) +
#     geom_segment(x = 0, y = 0, alpha = 0.4)  +
#     geom_text_repel(aes(x=get(PC_OI[1]) * scaling_fac,
#                         y=get(PC_OI[2]) * scaling_fac,
#                         label = round(ppm, digits = 2)),
#                     max.overlaps = 100, show.legend = FALSE,
#                     min.segment.length = 100) +
#     scale_x_continuous(limits = c(min(dt.pca[[PC_OI[1]]]),max(dt.pca[[PC_OI[1]]]))) +
#     scale_y_continuous(limits = c(min(dt.pca[[PC_OI[2]]]),max(dt.pca[[PC_OI[2]]]))) +
#     scale_color_manual(values = c("#d01c8b","#4dac26")) +
#     labs(x = paste0(PC_OI[1]," (",
#                     round(100*importance.pca["Proportion of Variance",PC_OI[1]]),"%)"),
#          y = paste0(PC_OI[2]," (",
#                     round(100*importance.pca["Proportion of Variance",PC_OI[2]]),"%)"),
#          col = "Phase") +
#     theme_bw() +
#     theme(panel.grid = element_blank(),
#           axis.text = element_text(color = "black"),
#           legend.position = legend.position)
#
#   p_out <- ggpubr::ggarrange(p_scores, p_loading, nrow = 1)
#   return(p_out)
# }
#
# pcaHL_12 <- plot_pca_feat(c("PC1","PC2"), legend.position = "right")
# pcaHL_13 <- plot_pca_feat(c("PC1","PC3"), legend.position = "bottom")
# pcaHL_23 <- plot_pca_feat(c("PC2","PC3"), legend.position = "bottom")
#
#
# # ~ ~ ~ ~ ~ ~ #
# # Beer plot   #
# # ~ ~ ~ ~ ~ ~ #
# dt_scaled <- data.table(as.table(tmp_ft.tab.scaled))
# setnames(dt_scaled, c("newID","feat","signal"))
# dt_scaled <- merge(dt_scaled, beerschi$info, by = "newID")
# dt_scaled[, pos.y.new := pos.y]
# dt_scaled[rack == 2, pos.y.new := pos.y.new + 12]
# dt_scaled[, pos.x.new := pos.x]
# dt_scaled[rack == 3 & pos.x == "A", pos.x.new := "I"]
# dt_scaled[rack == 3 & pos.x == "B", pos.x.new := "J"]
#
# ggplot(dt_scaled[feat == "H-14100"], aes(x = pos.y.new, y = pos.x.new,
#                                          fill = signal)) +
#   geom_tile()
#
# dt_pca_beer <- copy(dt.pca[,1:10])
# dt_pca_beer <- merge(dt_pca_beer, beerschi$info, by = "newID")
# dt_pca_beer[, pos.y.new := pos.y]
# dt_pca_beer[rack == 2, pos.y.new := pos.y.new + 12]
# dt_pca_beer[, pos.x.new := pos.x]
# dt_pca_beer[rack == 3 & pos.x == "A", pos.x.new := "I"]
# dt_pca_beer[rack == 3 & pos.x == "B", pos.x.new := "J"]
#
#
# ggplot(dt_pca_beer, aes(x = pos.y.new, y = pos.x.new,
#                         fill = PC2)) +
#   geom_tile() +
#   scale_y_discrete(limits = rev)


# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# PCA Plot only for (L) #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
PC_OI <- c("PC1","PC2")

dt.pca.loading <- copy(beer_L$PCA$loadings)
tmp_filled <- beer_L$filled[,.(ppm = median(ppm),
                               peakVal = median(peakValue)),
                            by = peakIndex]
tmp_filled[, peakIndex := as.character(peakIndex)]
dt.pca.loading <- merge(dt.pca.loading, tmp_filled,
                        by.x = "feat", by.y = "peakIndex")
dt.pca.loading[, scale := beer_L$PCA$pca.object$scale[feat]]
n <- 10
#relfeats <-dt.pca.loading[order(sqrt((get(PC_OI[1])*scale)^2 + (get(PC_OI[2])*scale)^2), decreasing = T)][1:n, feat]
relfeats <-dt.pca.loading[order(sqrt((get(PC_OI[1]))^2 + (get(PC_OI[2]))^2), decreasing = T)][1:n, feat]
tmp_scores <- merge(beer_L$PCA$scores, beerschi$info[, .(sample = as.character(exp.L), color2)])

p_scores <- ggplot(tmp_scores, aes(get(PC_OI[1]), get(PC_OI[2]), fill = color2)) +
  geom_point(shape = 21) +
  scale_fill_manual(values = beerschi_colors2) +
  labs(x = paste0(PC_OI[1]," (",
                  round(100*beer_L$PCA$importance["Proportion of Variance",PC_OI[1]]),"%)"),
       y = paste0(PC_OI[2]," (",
                  round(100*beer_L$PCA$importance["Proportion of Variance",PC_OI[2]]),"%)"),
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
                  round(100*beer_L$PCA$importance["Proportion of Variance",PC_OI[1]]),"%)"),
       y = paste0(PC_OI[2]," (",
                  round(100*beer_L$PCA$importance["Proportion of Variance",PC_OI[2]]),"%)"),
       col = "Phase") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "right")

P_L_12 <- egg::ggarrange(p_scores, p_loading, nrow = 1)
ggsave("analysis/v1/plots/bierschinken_L_PCA.pdf", plot = P_L_12,
       width = 7.5, height = 3.2)

# Bear plot
dt_scaled <- data.table(as.table(beer_L$feat.tab))
dt_scaled <- cbind(dt_scaled[,.(feat = V2, signal = N)], beer_L$color[,.(exp.L = as.integer(sample))])
#setnames(dt_scaled, c("newID","feat","signal"))
dt_scaled <- merge(dt_scaled, beerschi$info, by = "exp.L")
dt_scaled[, pos.y.new := pos.y]
dt_scaled[rack == 2, pos.y.new := pos.y.new + 12]
dt_scaled[, pos.x.new := pos.x]
dt_scaled[rack == 3 & pos.x == "A", pos.x.new := "I"]
dt_scaled[rack == 3 & pos.x == "B", pos.x.new := "J"]
#dt_scaled[, signal_z := scale(signal), by = "feat"]
dt_scaled[, signal_z := (signal-min(signal))/(max(signal)-min(signal)), by = "feat"]
dt_scaled <- merge(dt_scaled, beer_L$filled[!duplicated(peakIndex),
                                            .(ppm = paste(round(ppm, digits = 2),"ppm"),
                                              feat = as.character(peakIndex))],
                   by = "feat")

p_relfeat_dist <- ggplot(dt_scaled[feat %in% relfeats], aes(x = pos.y.new, y = pos.x.new,
                                       fill = signal_z, col = signal_z)) +
  geom_tile() +
  facet_wrap("ppm", scales = "free", nrow = 2) +
  labs(x = "Index (x)", y = "Index (y)", fill = "range-scaled\nsignal",
       col = "range-scaled\nsignal") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(option = "magma") + scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
ggsave("analysis/v1/plots/bierschinken_L_relfeats.pdf", plot = p_relfeat_dist,
       width = 10.5, height = 3.35)

dt_pca_beer <- copy(tmp_scores[,1:10][, exp.L := as.integer(sample)])
dt_pca_beer <- merge(dt_pca_beer, beerschi$info, by = "exp.L")
dt_pca_beer[, pos.y.new := pos.y]
dt_pca_beer[rack == 2, pos.y.new := pos.y.new + 12]
dt_pca_beer[, pos.x.new := pos.x]
dt_pca_beer[rack == 3 & pos.x == "A", pos.x.new := "I"]
dt_pca_beer[rack == 3 & pos.x == "B", pos.x.new := "J"]


p_dist_pc1 <- ggplot(dt_pca_beer, aes(x = pos.y.new, y = pos.x.new,
                                      fill = get(PC_OI[1]), col = get(PC_OI[1]))) +
  geom_tile() +
  labs(x = "Index (x)", y = "Index (y)", title = PC_OI[1],
       fill = PC_OI[1], col = PC_OI[1]) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(option = "magma") +
  scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))

p_dist_pc2 <- ggplot(dt_pca_beer, aes(x = pos.y.new, y = pos.x.new,
                                      fill = get(PC_OI[2]), col = get(PC_OI[2]))) +
  geom_tile() +
  labs(x = "Index (x)", y = "Index (y)", title = PC_OI[2],
       fill = PC_OI[2], col = PC_OI[2]) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(option = "magma") +
  scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))

p_dist_pc12 <- egg::ggarrange(p_dist_pc1, p_dist_pc2, nrow = 1)
ggsave("analysis/v1/plots/bierschinken_L_PC1_PC2.pdf", plot = p_dist_pc12,
       width = 6.5, height = 2.05)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# PCA Plot only for (H) #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
PC_OI <- c("PC1","PC2")

dt.pca.loading <- copy(beer_H$PCA$loadings)
tmp_filled <- beer_H$filled[,.(ppm = median(ppm),
                               peakVal = median(peakValue)),
                            by = peakIndex]
tmp_filled[, peakIndex := as.character(peakIndex)]
dt.pca.loading <- merge(dt.pca.loading, tmp_filled,
                        by.x = "feat", by.y = "peakIndex")
dt.pca.loading[, scale := beer_H$PCA$pca.object$scale[feat]]
n <- 10
#relfeats <-dt.pca.loading[order(sqrt((get(PC_OI[1])*scale)^2 + (get(PC_OI[2])*scale)^2), decreasing = T)][1:n, feat]
relfeats <-dt.pca.loading[order(sqrt((get(PC_OI[1]))^2 + (get(PC_OI[2]))^2), decreasing = T)][1:n, feat]
tmp_scores <- merge(beer_H$PCA$scores, beerschi$info[, .(sample = as.character(exp.H), color2)])

p_scores <- ggplot(tmp_scores, aes(get(PC_OI[1]), get(PC_OI[2]), fill = color2)) +
  geom_point(shape = 21) +
  scale_fill_manual(values = beerschi_colors2) +
  labs(x = paste0(PC_OI[1]," (",
                  round(100*beer_H$PCA$importance["Proportion of Variance",PC_OI[1]]),"%)"),
       y = paste0(PC_OI[2]," (",
                  round(100*beer_H$PCA$importance["Proportion of Variance",PC_OI[2]]),"%)"),
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
                  round(100*beer_H$PCA$importance["Proportion of Variance",PC_OI[1]]),"%)"),
       y = paste0(PC_OI[2]," (",
                  round(100*beer_H$PCA$importance["Proportion of Variance",PC_OI[2]]),"%)"),
       col = "Phase") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "right")

P_L_12 <- egg::ggarrange(p_scores, p_loading, nrow = 1)
ggsave("analysis/v1/plots/bierschinken_H_PCA.pdf", plot = P_L_12,
       width = 7.5, height = 3.2)

# Bear plot
dt_scaled <- data.table(as.table(beer_H$feat.tab))
dt_scaled <- cbind(dt_scaled[,.(feat = V2, signal = N)], beer_H$color[,.(exp.H = as.integer(sample))])
#setnames(dt_scaled, c("newID","feat","signal"))
dt_scaled <- merge(dt_scaled, beerschi$info, by = "exp.H")
dt_scaled[, pos.y.new := pos.y]
dt_scaled[rack == 2, pos.y.new := pos.y.new + 12]
dt_scaled[, pos.x.new := pos.x]
dt_scaled[rack == 3 & pos.x == "A", pos.x.new := "I"]
dt_scaled[rack == 3 & pos.x == "B", pos.x.new := "J"]
#dt_scaled[, signal_z := scale(signal), by = "feat"]
dt_scaled[, signal_z := (signal-min(signal))/(max(signal)-min(signal)), by = "feat"]
dt_scaled <- merge(dt_scaled, beer_H$filled[!duplicated(peakIndex),
                                            .(ppm = paste(round(ppm, digits = 2),"ppm"),
                                              feat = as.character(peakIndex))],
                   by = "feat")

p_relfeat_dist <- ggplot(dt_scaled[feat %in% relfeats], aes(x = pos.y.new, y = pos.x.new,
                                                            fill = signal_z, col = signal_z)) +
  geom_tile() +
  facet_wrap("ppm", scales = "free", nrow = 2) +
  labs(x = "Index (x)", y = "Index (y)", fill = "range-scaled\nsignal",
       col = "range-scaled\nsignal") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(option = "magma") + scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
ggsave("analysis/v1/plots/bierschinken_H_relfeats.pdf", plot = p_relfeat_dist,
       width = 10.5, height = 3.35)

dt_pca_beer <- copy(tmp_scores[,1:10][, exp.H := as.integer(sample)])
dt_pca_beer <- merge(dt_pca_beer, beerschi$info, by = "exp.H")
dt_pca_beer[, pos.y.new := pos.y]
dt_pca_beer[rack == 2, pos.y.new := pos.y.new + 12]
dt_pca_beer[, pos.x.new := pos.x]
dt_pca_beer[rack == 3 & pos.x == "A", pos.x.new := "I"]
dt_pca_beer[rack == 3 & pos.x == "B", pos.x.new := "J"]


p_dist_pc1 <- ggplot(dt_pca_beer, aes(x = pos.y.new, y = pos.x.new,
                                      fill = get(PC_OI[1]), col = get(PC_OI[1]))) +
  geom_tile() +
  labs(x = "Index (x)", y = "Index (y)", title = PC_OI[1],
       fill = PC_OI[1], col = PC_OI[1]) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(option = "magma") +
  scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))

p_dist_pc2 <- ggplot(dt_pca_beer, aes(x = pos.y.new, y = pos.x.new,
                                      fill = get(PC_OI[2]), col = get(PC_OI[2]))) +
  geom_tile() +
  labs(x = "Index (x)", y = "Index (y)", title = PC_OI[2],
       fill = PC_OI[2], col = PC_OI[2]) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(option = "magma") +
  scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))

p_dist_pc12 <- egg::ggarrange(p_dist_pc1, p_dist_pc2, nrow = 1)
ggsave("analysis/v1/plots/bierschinken_H_PC1_PC2.pdf", plot = p_dist_pc12,
       width = 6.5, height = 2.05)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
# PCA Plot only for (G) #
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ #
PC_OI <- c("PC1","PC2")

dt.pca.loading <- copy(beer_G$PCA$loadings)
tmp_filled <- beer_G$filled[,.(ppm = median(ppm),
                               peakVal = median(peakValue)),
                            by = peakIndex]
tmp_filled[, peakIndex := as.character(peakIndex)]
dt.pca.loading <- merge(dt.pca.loading, tmp_filled,
                        by.x = "feat", by.y = "peakIndex")
dt.pca.loading[, scale := beer_G$PCA$pca.object$scale[feat]]
n <- 10
#relfeats <-dt.pca.loading[order(sqrt((get(PC_OI[1])*scale)^2 + (get(PC_OI[2])*scale)^2), decreasing = T)][1:n, feat]
relfeats <-dt.pca.loading[order(sqrt((get(PC_OI[1]))^2 + (get(PC_OI[2]))^2), decreasing = T)][1:n, feat]
tmp_scores <- merge(beer_G$PCA$scores, beerschi$info[, .(sample = as.character(exp.G), color2)])

p_scores <- ggplot(tmp_scores, aes(get(PC_OI[1]), get(PC_OI[2]), fill = color2)) +
  geom_point(shape = 21) +
  scale_fill_manual(values = beerschi_colors2) +
  labs(x = paste0(PC_OI[1]," (",
                  round(100*beer_G$PCA$importance["Proportion of Variance",PC_OI[1]]),"%)"),
       y = paste0(PC_OI[2]," (",
                  round(100*beer_G$PCA$importance["Proportion of Variance",PC_OI[2]]),"%)"),
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
                  round(100*beer_G$PCA$importance["Proportion of Variance",PC_OI[1]]),"%)"),
       y = paste0(PC_OI[2]," (",
                  round(100*beer_G$PCA$importance["Proportion of Variance",PC_OI[2]]),"%)"),
       col = "Phase") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "right")

P_L_12 <- egg::ggarrange(p_scores, p_loading, nrow = 1)
ggsave("analysis/v1/plots/bierschinken_G_PCA.pdf", plot = P_L_12,
       width = 7.5, height = 3.2)

# Bear plot
dt_scaled <- data.table(as.table(beer_G$feat.tab))
dt_scaled <- cbind(dt_scaled[,.(feat = V2, signal = N)], beer_G$color[,.(exp.G = as.integer(sample))])
#setnames(dt_scaled, c("newID","feat","signal"))
dt_scaled <- merge(dt_scaled, beerschi$info, by = "exp.G")
dt_scaled[, pos.y.new := pos.y]
dt_scaled[rack == 2, pos.y.new := pos.y.new + 12]
dt_scaled[, pos.x.new := pos.x]
dt_scaled[rack == 3 & pos.x == "A", pos.x.new := "I"]
dt_scaled[rack == 3 & pos.x == "B", pos.x.new := "J"]
#dt_scaled[, signal_z := scale(signal), by = "feat"]
dt_scaled[, signal_z := (signal-min(signal))/(max(signal)-min(signal)), by = "feat"]
dt_scaled <- merge(dt_scaled, beer_G$filled[!duplicated(peakIndex),
                                            .(ppm = paste(round(ppm, digits = 2),"ppm"),
                                              feat = as.character(peakIndex))],
                   by = "feat")

p_relfeat_dist <- ggplot(dt_scaled[feat %in% relfeats], aes(x = pos.y.new, y = pos.x.new,
                                                            fill = signal_z, col = signal_z)) +
  geom_tile() +
  facet_wrap("ppm", scales = "free", nrow = 2) +
  labs(x = "Index (x)", y = "Index (y)", fill = "range-scaled\nsignal",
       col = "range-scaled\nsignal") +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(option = "magma") + scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA))
ggsave("analysis/v1/plots/bierschinken_G_relfeats.pdf", plot = p_relfeat_dist,
       width = 10.5, height = 3.35)

dt_pca_beer <- copy(tmp_scores[,1:10][, exp.G := as.integer(sample)])
dt_pca_beer <- merge(dt_pca_beer, beerschi$info, by = "exp.G")
dt_pca_beer[, pos.y.new := pos.y]
dt_pca_beer[rack == 2, pos.y.new := pos.y.new + 12]
dt_pca_beer[, pos.x.new := pos.x]
dt_pca_beer[rack == 3 & pos.x == "A", pos.x.new := "I"]
dt_pca_beer[rack == 3 & pos.x == "B", pos.x.new := "J"]


p_dist_pc1 <- ggplot(dt_pca_beer, aes(x = pos.y.new, y = pos.x.new,
                                      fill = get(PC_OI[1]), col = get(PC_OI[1]))) +
  geom_tile() +
  labs(x = "Index (x)", y = "Index (y)", title = PC_OI[1],
       fill = PC_OI[1], col = PC_OI[1]) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(option = "magma") +
  scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))

p_dist_pc2 <- ggplot(dt_pca_beer, aes(x = pos.y.new, y = pos.x.new,
                                      fill = get(PC_OI[2]), col = get(PC_OI[2]))) +
  geom_tile() +
  labs(x = "Index (x)", y = "Index (y)", title = PC_OI[2],
       fill = PC_OI[2], col = PC_OI[2]) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0), limits = rev) +
  scale_fill_viridis_c(option = "magma") +
  scale_color_viridis_c(option = "magma") +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.title = element_text(hjust = 0.5))

p_dist_pc12 <- egg::ggarrange(p_dist_pc1, p_dist_pc2, nrow = 1)
ggsave("analysis/v1/plots/bierschinken_G_PC1_PC2.pdf", plot = p_dist_pc12,
       width = 6.5, height = 2.05)










pairs(beer_L$PCA$scores[,2:6],
      col = as.factor(beer_L$PCA$scores$color), main = "Bierschinken (L)")
pairs(beer_H$PCA$scores[,2:6],
      col = as.factor(beer_H$PCA$scores$color), main = "Bierschinken (H)")
pairs(beer_G$PCA$scores[,2:6],
      col = as.factor(beer_G$PCA$scores$color), main = "Bierschinken (G)")
dev.off()


ROI.ppm <- 6.1
roiWidth.ppm <- 0.1

ROIplot_SW(beer_H, ROI = ROI.ppm, ROI.width = roiWidth.ppm,
           color.code = MetBrewer::met.brewer("Isfahan2",7))
