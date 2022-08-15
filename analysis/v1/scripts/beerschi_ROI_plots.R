
beerschi_colors2_spec <- c(w = "#ddc965", r = "#980043", g = "#df65b0")


p <- ROIplot_SW(beer_H, ROI = 0.925, ROI.width = 0.025, incl.plots = 1:3,
                color.code = beerschi_colors2_spec, alpha = 0.25)
ggsave("analysis/v1/plots/rev_beer_H_0.925.pdf", plot = p, width = 8, height = 8)

p <- ROIplot_SW(beer_L, ROI = 4.38, ROI.width = 0.025, incl.plots = 1:3,
                color.code = beerschi_colors2_spec, alpha = 0.25)
ggsave("analysis/v1/plots/rev_beer_L_4.38.pdf", plot = p, width = 8, height = 8)
