spectra_analysis <- function(i_spectra,
                             i_ppm,
                             i_color,
                             PCA.minFeatNr = ceiling(nrow(i_spectra)/2)) {
  # exp
  # i_spectra <- baerchi$spectra.norm[as.character(baerchi$info$exp.H),]
  # i_ppm <- baerchi$ppm
  # i_color <- baerchi$info$color

  spq_summary <- list()
  spq_summary$color <- data.table(sample = rownames(i_spectra),
                                  color  = i_color)
  spq_summary$spectra <- melt(cbind(data.table(sample = rownames(i_spectra)),
                                    data.table(i_spectra)), id.vars = "sample",
                              variable.name = "ppm", value.name = "value")
  spq_summary$spectra[, ppm := as.numeric(as.character(ppm))]


  peaks <- speaq::getWaveletPeaks(Y.spec=i_spectra,
                                  X.ppm=i_ppm,
                                  baselineThresh = 10,
                                  SNR.Th = -1,
                                  nCPU = 6,
                                  include_nearbyPeaks = TRUE) # nCPU set to 2 for the


  spq_summary$peaks <- data.table(sample = rownames(i_spectra)[peaks$Sample],
                                  peak.index = peaks$peakIndex,
                                  ppm = peaks$peakPPM,
                                  value = peaks$peakValue,
                                  peak.SNR = peaks$peakSNR,
                                  peak.scale =peaks$peakScale)

  peaks.grouped <- speaq::PeakGrouper(Y.peaks = peaks,
                                      min.samp.grp = 5,
                                      grouping.window.width = 75)

  spq_summary$grouped <- data.table(peaks.grouped)
  spq_summary$grouped[, ppm := median(peakPPM), by = peakIndex]
  spq_summary$grouped[, sample := rownames(i_spectra)[Sample]]
  spq_summary$grouped[, Sample := NULL]

  peaks.filled <- speaq::PeakFilling(Y.grouped = peaks.grouped,
                                     Y.spec = i_spectra,
                                     max.index.shift = 5,
                                     nCPU = 6)

  spq_summary$filled <- data.table(peaks.filled)
  spq_summary$filled[, ppm := median(peakPPM), by = peakIndex]
  spq_summary$filled[, sample := rownames(i_spectra)[Sample]]
  spq_summary$filled[, Sample := NULL]

  # # do PCA
  # peaks.features <- speaq::BuildFeatureMatrix(peaks.filled)
  # peaks.features <- peaks.features[,apply(peaks.features,
  #                                         2,
  #                                         function(x) sum(x==0)) <= nrow(peaks.features) * 0.025]
  # message("PCA on ",ncol(peaks.features)," features.")
  # # imputation with kNN (only features, which occur in at least 97.5% of samples)
  # peaks.features[which(peaks.features == 0, arr.ind = TRUE)] <- NA
  # #peaks.features <- t(bnstruct::knn.impute(t(peaks.features), k = 5))
  # peaks.features <- missForest::missForest(peaks.features)$ximp
  #
  # spq_summary$feat.tab <- peaks.features
  # keep_feat <- spq_summary$filled[!(ppm < 0.4 | (ppm > 4.41 & ppm < 6)),
  #                                 as.character(unique(peakIndex))]
  # keep_feat <- keep_feat[keep_feat %in% colnames(peaks.features)]
  # ft.tab.scaled <- speaq::SCANT(data.matrix = peaks.features[,keep_feat],
  #                               type = c("center","unit"))
  # ft.pca <- prcomp(ft.tab.scaled, scale. = FALSE)
  # dt.pca <- cbind(data.table(sample = rownames(i_spectra)),
  #                 ft.pca$x)
  # dt.pca <- merge(dt.pca, spq_summary$color)
  # dt.pca.loading <- cbind(data.table(feat = rownames(ft.pca$rotation)),
  #                         ft.pca$rotation)
  # spq_summary$PCA <- list(scores = dt.pca,
  #                         loadings = dt.pca.loading,
  #                         importance = summary(ft.pca)$importance,
  #                         pca.object = ft.pca)


  # output
  return(spq_summary)
}

# Includes PCA and ICA
feat.PCA <- function(spq_summary, max.na = 0.025, min.SNR = 2) {
  require(missForest)
  #require(ica)
  #require(fastICA)
  
  # Build feature matrix
  peaks.features <- dcast(spq_summary$grouped[peakSNR >= min.SNR],
                          sample ~ peakIndex,
                          value.var = "peakValue")
  tmpnames <- rownames(peaks.features)
  peaks.features <- as.matrix(peaks.features[,-1])
  rownames(peaks.features) <- tmpnames
  peaks.features <- peaks.features[,apply(peaks.features,
                                          2,
                                          function(x) sum(is.na(x))) <= nrow(peaks.features) * max.na]

  # replace obvious outliers with NA
  peaks.features_ci <- apply(peaks.features, 2, quantile, na.rm = TRUE, prob = c(0.025, 0.975)) # 95% CI
  up_rm <- which(t(peaks.features) > peaks.features_ci[2,] + (peaks.features_ci[2,] - peaks.features_ci[1,]), arr.ind = TRUE)
  dw_rm <- which(t(peaks.features) < peaks.features_ci[1,] - (peaks.features_ci[2,] - peaks.features_ci[1,]), arr.ind = TRUE)
  if(nrow(up_rm) > 0)
    peaks.features[up_rm[,c(2,1)]] <- NA
  if(nrow(dw_rm) > 0)
    peaks.features[dw_rm[,c(2,1)]] <- NA

  peaks.features <- peaks.features[,apply(peaks.features,
                                          2,
                                          function(x) sum(is.na(x))) <= nrow(peaks.features) * max.na]


  # imputation with random forests
  peaks.features <- missForest::missForest(peaks.features)$ximp

  spq_summary$feat.tab <- peaks.features
  keep_feat <- spq_summary$filled[!(ppm < 0.4 | (ppm > 4.4 & ppm < 6.01)),
                                  as.character(unique(peakIndex))]
  keep_feat <- keep_feat[keep_feat %in% colnames(peaks.features)]
  # ft.tab.scaled <- speaq::SCANT(data.matrix = peaks.features[,keep_feat],
  #                               type = c("center","unit"))
  # ft.tab.scaled <- speaq::SCANT(data.matrix = peaks.features[,keep_feat],
  #                               type = c("center","pareto"))
  ft.tab.scaled <- apply(peaks.features[,keep_feat],2,
                         function(x) {
                           x <- (x - mean(x)) / sd(x)
                           return(x)
                         })

  message("PCA on ",ncol(ft.tab.scaled)," features.")

  ft.pca <- prcomp(ft.tab.scaled, scale. = FALSE)
  dt.pca <- cbind(data.table(sample = spq_summary$color$sample),
                  ft.pca$x)
  dt.pca <- merge(dt.pca, spq_summary$color)
  dt.pca.loading <- cbind(data.table(feat = rownames(ft.pca$rotation)),
                          ft.pca$rotation)
  spq_summary$PCA <- list(scores = dt.pca,
                          loadings = dt.pca.loading,
                          importance = summary(ft.pca)$importance,
                          pca.object = ft.pca)

  
  #––––––––––––––––––––––––––––––––––––––#
  # Independent component analysis (ICA) #
  #––––––––––––––––––––––––––––––––––––––#
  #nics <- ncol(ft.tab.scaled)
  # nics <- 2
  # ft.ica <- ica(ft.tab.scaled, nics, method = "fast", maxit = 5000)
  # dt.ica <- cbind(data.table(sample = spq_summary$color$sample),
  #                 ft.ica$S)
  # dt.ica <- merge(dt.ica, spq_summary$color)
  # dt.ica.loading <- cbind(data.table(feat = colnames(ft.tab.scaled)),
  #                         ft.ica$M)
  # spq_summary$ICA <- list(scores = dt.ica,
  #                         loadings = dt.ica.loading,
  #                         importance = ft.ica$vafs,
  #                         ica.object = ft.ica)
  # names(spq_summary$ICA$importance) <- paste0("V",1:nics)
  
  # 2 
  nics <- 3
  ft.ica <- ica::icafast(ft.tab.scaled, nics)
  dt.ica <- cbind(data.table(sample = spq_summary$color$sample),
                  ft.ica$S)
  setnames(dt.ica, paste0("V",1:nics),paste0("IC",1:nics))
  dt.ica <- merge(dt.ica, spq_summary$color)
  dt.ica.loading <- cbind(data.table(feat = colnames(ft.tab.scaled)),
                          t(ft.ica$W))
  # dt.ica.loading <- cbind(data.table(feat = colnames(ft.tab.scaled)),
  #                         ft.ica$M)
  setnames(dt.ica.loading, paste0("V",1:nics),paste0("IC",1:nics))
  spq_summary$ICA <- list(scores = dt.ica,
                          loadings = dt.ica.loading,
                          importance = ft.ica$vafs,
                          ica.object = ft.ica)
  names(spq_summary$ICA$importance) <- paste0("IC",1:nics)
  
  # output
  return(spq_summary)

}

