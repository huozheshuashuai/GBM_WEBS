setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
gm <- readRDS("results/geomx/geomx_bundle.rds")
res <- gm$roi_meta
pt <- gm$pt

# panel d: corrected title reflecting mixed-model result
p1 <- ggplot(pt, aes(stain, PC2, fill = stain)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = patient), alpha = 0.3, color = "grey60") +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5, size = 1) +
  scale_fill_manual(values = c("CD64" = "#7B4EA3", "GFAP" = "#3D9B8C")) +
  facet_wrap(~ timepoint) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(y = "Frozen KNCC PC2 axis score (ROI level)", x = NULL,
       title = "Axis detectable in GFAP-enriched ROIs, not confined to myeloid niches")
ggsave("figures_kncc/geomx_pc2_by_stain.pdf", p1, width = 7, height = 5)

# panel e: neutral title
p2 <- ggplot(res, aes(px_Microglia, PC2, color = stain)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("CD64" = "#7B4EA3", "GFAP" = "#3D9B8C")) +
  geom_smooth(method = "lm", se = FALSE, color = "grey30") +
  theme_bw(base_size = 12) +
  labs(x = "Microglia score (ROI)", y = "PC2 axis score (ROI)",
       title = "PC2 axis vs microglial signal across ROIs",
       caption = "microglia proxy validates stain-based ROI selection (CD64 > GFAP, p = 0.001)")
ggsave("figures_kncc/geomx_pc2_microglia.pdf", p2, width = 6, height = 5)
cat("panels fixed\n")
