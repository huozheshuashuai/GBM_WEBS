# Regenerate S8: acetylation heatmap of top spliceosome sites
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(dplyr); library(tidyr) })

tt <- read.csv("results/acetyl_sites_limma.csv")
tt$gene <- sub("-.*$", "", tt$site)
source("08_cptac_proteomics.R")
splice <- get_spliceosome_genes()
top <- tt %>% filter(gene %in% splice, adj.P.Val < 0.05) %>%
  arrange(desc(abs(t))) %>% head(40)

ac <- read.delim("E:/KIMI/data/HS_CPTAC_GBM_acetylome_mssm_site.cct", check.names = FALSE)
webs <- read.csv("results/cptac_webs.csv")
sub <- ac[match(top$site, ac[[1]]), ]
m <- as.matrix(sapply(sub[, -1], as.numeric)); rownames(m) <- top$site
common <- intersect(colnames(m), webs$sample)
m <- m[, common]
webs <- webs[match(common, webs$sample), ]
ord <- order(webs$WEBS)
m <- m[, ord]
grp <- webs$WEBS_group[ord]

df <- as.data.frame(m)
df$site <- rownames(m)
long <- pivot_longer(df, -site, names_to = "sample", values_to = "val")
long$site <- factor(long$site, levels = rev(top$site))
long$sample <- factor(long$sample, levels = colnames(m))
p <- ggplot(long, aes(sample, site, fill = val)) +
  geom_tile() +
  scale_fill_gradient2(low = "#3D9B8C", mid = "white", high = "#7B4EA3",
                       midpoint = 0, na.value = "grey90") +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 6)) +
  labs(title = "Top 40 differential spliceosome acetylation sites (FDR<0.05)",
       subtitle = "Samples ordered by WEBS (left=LLS, right=LHS)",
       x = paste0("99 CPTAC-GBM samples (LHS: ", sum(grp == "LHS"), " right)"),
       y = NULL, fill = "Acetylation")
ggsave("figures_final/figS_acetyl_heatmap.pdf", p, width = 8, height = 8)
cat("ok\n")
