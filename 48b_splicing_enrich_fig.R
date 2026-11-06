setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(ggplot2))
ego <- read.csv("results/splicing_enrichment.csv")
ego <- ego[order(ego$p.adjust), ][1:8, ]
ego$Description <- factor(ego$Description, levels = rev(ego$Description))
ego$neglog <- -log10(ego$p.adjust)
p <- ggplot(ego, aes(neglog, Description, size = Count, color = p.adjust)) +
  geom_point() +
  scale_color_gradient(low = "#7B4EA3", high = "grey60") +
  theme_bw(base_size = 12) +
  labs(title = "GO enrichment: genes with WEBS-associated AS events",
       x = "-log10(adj. p)", y = NULL, size = "Genes", color = "adj. p")
ggsave("figures_final/fig_splicing_enrichment.pdf", p, width = 8, height = 4.5)
cat("saved\n")
