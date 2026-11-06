setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages(library(ggplot2))
d <- read.csv("results/wgcna_module_webs.csv")
d$sig <- ifelse(d$q < 0.05, "q<0.05", "n.s.")
d <- d[order(d$correlation), ]
d$module <- factor(d$module, levels = d$module)
p <- ggplot(d, aes(module, correlation, fill = sig)) +
  geom_bar(stat = "identity") + coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("q<0.05" = "#7B4EA3", "n.s." = "grey70")) +
  theme_bw(base_size = 12) + theme(legend.position = "none") +
  labs(title = "WGCNA module eigengenes vs WEBS (CPTAC)",
       subtitle = "Top module ME3 (r = 0.72) is neural-lineage, not spliceosome",
       x = NULL, y = "Pearson correlation with WEBS")
ggsave("figures_final/fig_wgcna_moduletrait.pdf", p, width = 6.5, height = 4.5)
cat("ok\n")
