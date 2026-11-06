# =============================================================================
# 70. Figure 1 rebuild: clean vector study-design flowchart (replaces AI image)
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(grid) })

PURPLE <- "#7B4EA3"; TEAL <- "#3D9B8C"; RED <- "#d62728"; BLUE <- "#1f77b4"; GOLD <- "#D4A017"

box <- function(x, y, w, h, label, fill, size = 2.6, color = "grey25") {
  list(
    annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2,
             fill = fill, color = color, linewidth = 0.4, alpha = 0.92),
    annotate("text", x = x, y = y, label = label, size = size, lineheight = 0.95)
  )
}
arr <- function(x1, y1, x2, y2) {
  annotate("segment", x = x1, y = y1, xend = x2, yend = y2,
           arrow = grid::arrow(length = unit(0.15, "cm")), color = "grey40", linewidth = 0.5)
}

p <- ggplot() +
  # ---- Data layer (left column) ----
  box(1.5, 8.6, 2.5, 1.3, "Single-cell & spatial\nGSE138794 (35,823 cells)\nGSE131928, GSE194329", "#E8E4F3") +
  box(1.5, 6.9, 2.5, 1.3, "CPTAC-GBM (n = 99)\nProteome, acetylome,\nphosphoproteome, RNA-seq", "#E4F0EC") +
  box(1.5, 5.2, 2.5, 1.3, "Bulk cohorts\nTCGA, CGGA-693/325, GLASS,\nGSE7696, GSE83300", "#FDF3E3") +
  box(1.5, 3.5, 2.5, 1.3, "Functional resources\nTCGASpliceSeq (PSI)\nDepMap, GDSC2", "#F3E4E4") +

  # ---- WEBS definition (center) ----
  box(4.8, 6.5, 2.9, 1.9, "WEBS\nUCell(writers) \u2212 UCell(erasers)\nLLS / LDTS / LHS =\nlow / intermediate / high-WEBS state",
      "#DDD3EA", size = 2.8, color = PURPLE) +

  # ---- Analysis layer (right column) ----
  box(8.3, 8.6, 2.9, 1.3, "Malignant-state characterization\nNeftel states, pseudotime,\nspatial niches", "#E8E4F3") +
  box(8.3, 6.9, 2.9, 1.3, "Spliceosome accumulation\n107/111 proteins up (CPTAC)\nacetylation + phosphorylation", "#E4F0EC") +
  box(8.3, 5.2, 2.9, 1.3, "Splicing output\n4,387 WEBS-associated AS events\n(replicated in TCGA-LGG)", "#FDF3E3") +
  box(8.3, 3.5, 2.9, 1.3, "Clinical association &\ncandidate therapeutic\nvulnerabilities", "#F3E4E4") +

  # ---- arrows ----
  arr(2.75, 8.6, 3.35, 7.1) + arr(2.75, 6.9, 3.35, 6.7) +
  arr(2.75, 5.2, 3.35, 6.3) + arr(2.75, 3.5, 3.35, 5.9) +
  arr(6.25, 7.1, 6.85, 8.5) + arr(6.25, 6.9, 6.85, 6.9) +
  arr(6.25, 6.6, 6.85, 5.3) + arr(6.25, 6.2, 6.85, 3.7) +

  # ---- headers ----
  annotate("text", x = 1.5, y = 9.6, label = "Data", fontface = "bold", size = 3.4) +
  annotate("text", x = 4.8, y = 9.6, label = "Score", fontface = "bold", size = 3.4) +
  annotate("text", x = 8.3, y = 9.6, label = "Findings", fontface = "bold", size = 3.4) +

  # ---- footer ----
  annotate("text", x = 4.9, y = 2.75, size = 2.3, color = "grey35",
           label = "Code and the WEBS implementation are publicly available (GitHub / Zenodo); datasets are accessible\nthrough the original repositories subject to their respective access requirements.") +

  coord_cartesian(xlim = c(0, 9.8), ylim = c(1.9, 10), expand = FALSE) +
  theme_void()

ggsave("figures_GND/Figure_1.pdf", p, width = 10, height = 7)
ggsave("figures_GND/Figure_1.png", p, width = 10, height = 7, dpi = 300)
cat("Figure 1 rebuilt.\n")
