
baerchi_colors_tmp <- c(red = "#b30000", lightred = "#d619ba", white = "#ddc965")

p <- ROIplot_SW(baer_L, ROI = 2.92, ROI.width = 0.05, incl.plots = 1:3,
                color.code = baerchi_colors_tmp, alpha = 0.25)
ggsave("analysis/v1/plots/rev_baer_L_2.92.pdf", plot = p, width = 9, height = 7)

p <- ROIplot_SW(baer_L, ROI = 1.75, ROI.width = 0.075, incl.plots = 1:3,
                color.code = baerchi_colors_tmp, alpha = 0.25)
ggsave("analysis/v1/plots/rev_baer_L_1.75.pdf", plot = p, width = 9, height = 7)

p <- ROIplot_SW(baer_L, ROI = 2.16, ROI.width = 0.075, incl.plots = 1:3,
                color.code = baerchi_colors_tmp, alpha = 0.25)
ggsave("analysis/v1/plots/rev_baer_L_2.16.pdf", plot = p, width = 9, height = 7)

# p <- ROIplot_SW(baer_L, ROI = 4.38, ROI.width = 0.05, incl.plots = 1:3, 
#                 color.code = baerchi_colors_tmp, alpha = 0.25)
# ggsave("analysis/v1/plots/rev_baer_L_4.38.pdf", plot = p, width = 9, height = 7)
# 
# 
# ROIplot_SW(baer_L, ROI = 1.71, ROI.width = 0.075, incl.plots = 1:3,
#            color.code = baerchi_colors_tmp, alpha = 0.25)
# 
# ROIplot_SW(baer_H, ROI = 2.66, ROI.width = 0.025, incl.plots = 1:3,
#            color.code = baerchi_colors_tmp, alpha = 0.25)
