# B/W-readability patches for BJC submission
# 1. 31_fig_landscape.R: F2D add shape encoding for WEBS groups
s = open("31_fig_landscape.R", encoding="utf-8").read()
s = s.replace('pD <- DimPlot(obj, group.by = "WEBS_group", cols = COL_WEBS) +',
              'pD <- DimPlot(obj, group.by = "WEBS_group", cols = COL_WEBS,\n               shape.by = "WEBS_group") +\n  scale_shape_manual(values = c("LLS" = 16, "LDTS" = 17, "LHS" = 15), name = "WEBS group") +')
open("31_fig_landscape.R", "w", encoding="utf-8").write(s)
print("F2D ok")

# 2. 67_gnd_figures.R: F3a spliceosome-up points -> triangle shape
s = open("67_gnd_figures.R", encoding="utf-8").read()
s = s.replace('lim$sz <- ifelse(lim$cat == "Spliceosome up", 1.8, 0.4)',
              'lim$sz <- ifelse(lim$cat == "Spliceosome up", 2.0, 0.4)\nlim$shp <- ifelse(lim$cat == "Spliceosome up", 17, 16)')
s = s.replace('geom_point(aes(color = cat, size = sz), alpha = 0.6) +',
              'geom_point(aes(color = cat, size = sz, shape = shp), alpha = 0.65) +\n  scale_shape_identity() +')
open("67_gnd_figures.R", "w", encoding="utf-8").write(s)
print("F3a ok")

# 3. 69_figure_fixes.R: F4b sig points -> triangles; F6 dashed lines; F7d sig labels
s = open("69_figure_fixes.R", encoding="utf-8").read()
s = s.replace('ph$sz <- ifelse(grepl("Spliceosome sig", ph$cat), 1.6, 0.4)',
              'ph$sz <- ifelse(grepl("Spliceosome sig", ph$cat), 1.8, 0.4)\nph$shp <- ifelse(grepl("Spliceosome sig", ph$cat), 17, 16)')
s = s.replace('geom_point(aes(color = cat, size = sz), alpha = 0.55) +',
              'geom_point(aes(color = cat, size = sz, shape = shp), alpha = 0.55) +\n  scale_shape_identity() +')
# KM: add linetype by group
s = s.replace('ggsurvplot(fit, data = d, pval = TRUE, pval.method = TRUE,',
              'ggsurvplot(fit, data = d, pval = TRUE, pval.method = TRUE,\n             linetype = "strata",')
# F7d: significance text labels on bars
s = s.replace('''p7d <- ggplot(key, aes(gene, spearman_r, fill = sig)) +
  geom_bar(stat = "identity") + coord_flip() +''',
'''p7d <- ggplot(key, aes(gene, spearman_r, fill = sig)) +
  geom_bar(stat = "identity") + coord_flip() +
  geom_text(aes(label = ifelse(sig == "n.s.", "n.s.", ifelse(sig == "nominal p < 0.05", "nom.", "FDR")),
                hjust = ifelse(spearman_r < 0, 1.15, -0.15)), size = 2.6) +''')
open("69_figure_fixes.R", "w", encoding="utf-8").write(s)
print("F4b/F6/F7d ok")
