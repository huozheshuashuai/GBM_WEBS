setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(dplyr); library(ggplot2) })
b <- readRDS("results/kncc/dual_omics_bundle.rds")
ph <- b$phos_delta
pairs_per_site <- rowSums(!is.na(ph))
df <- data.frame(pairs = pairs_per_site)
p1 <- ggplot(df, aes(pairs)) +
  geom_histogram(bins = 45, fill = "#7B4EA3") +
  geom_vline(xintercept = c(25, 40, 60), linetype = "dashed", color = "grey40") +
  annotate("text", x = c(25, 40, 60), y = Inf, label = c("≥25: 8,260", "≥40: 5,749", "≥60: 3,528"),
           vjust = 1.5, hjust = -0.1, size = 3.2) +
  theme_bw(base_size = 11) +
  labs(x = "Paired patients with non-missing delta per phosphosite",
       y = "Phosphosites",
       title = "Phosphoproteome effective coverage",
       subtitle = "median 12 pairs/site; sparse data limit site-level power")
sites_per_pt <- colSums(!is.na(ph))
df2 <- data.frame(sites = sites_per_pt)
p2 <- ggplot(df2, aes(sites)) +
  geom_histogram(bins = 40, fill = "#3D9B8C") +
  theme_bw(base_size = 11) +
  labs(x = "Phosphosites quantified per patient-pair", y = "Patient pairs",
       title = "Per-patient phosphosite coverage")
library(patchwork)
ggsave("figures_kncc/phospho_coverage_fig.pdf", p1 + p2, width = 11, height = 4.5)
cat("done\n")
