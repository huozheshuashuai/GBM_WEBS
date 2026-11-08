# =============================================================================
# 76. Figure 1 BEAUTIFIED: navy-header boxes + central balance scale + footer bar
# =============================================================================
setwd("E:/KIMI/GBM_WEBS_Project")
suppressMessages({ library(ggplot2); library(grid) })

NAVY  <- "#1F3B5C"; PURPLE <- "#7B4EA3"; PURPLE_BG <- "#EDE7F6"
TEAL  <- "#3D9B8C"; GOLD <- "#D4A017"; GREEN_BG <- "#E8F5E9"; GREEN_TXT <- "#2E5E32"
ARR   <- "grey45"

# helper: box with navy header + pastel body
box2 <- function(x, y, w, h, title, body, fill) {
  hh <- h * 0.30
  list(
    annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2,
             fill = fill, color = NAVY, linewidth = 0.8),
    annotate("rect", xmin = x - w/2, xmax = x + w/2, ymin = y + h/2 - hh, ymax = y + h/2,
             fill = NAVY, color = NAVY),
    annotate("text", x = x, y = y + h/2 - hh/2, label = title,
             color = "white", fontface = "bold", size = 3.0),
    annotate("text", x = x, y = y - h/2 + (h - hh)/2, label = body,
             color = "grey15", size = 2.5, lineheight = 1.05)
  )
}
arr <- function(x1, y1, x2, y2, col = ARR, lw = 0.6) {
  annotate("segment", x = x1, y = y1, xend = x2, yend = y2,
           arrow = grid::arrow(length = unit(0.16, "cm"), type = "closed"),
           color = col, linewidth = lw)
}

# balance scale elements (centered at cx, cy)
balance <- function(cx, cy, s = 1) {
  post_h <- 0.85 * s; beam_w <- 1.9 * s; drop <- 0.42 * s; pan_w <- 0.62 * s; pan_h <- 0.30 * s
  list(
    # post + base
    annotate("segment", x = cx, y = cy + post_h/2, xend = cx, yend = cy - post_h/2,
             color = NAVY, linewidth = 1.1),
    annotate("segment", x = cx - 0.35*s, y = cy - post_h/2, xend = cx + 0.35*s, yend = cy - post_h/2,
             color = NAVY, linewidth = 1.6),
    # beam
    annotate("segment", x = cx - beam_w/2, y = cy + post_h/2, xend = cx + beam_w/2, yend = cy + post_h/2,
             color = NAVY, linewidth = 1.1),
    annotate("point", x = cx, y = cy + post_h/2, color = NAVY, size = 2.5),
    # left hanger + pan (writers, gold)
    annotate("segment", x = cx - beam_w/2, y = cy + post_h/2, xend = cx - beam_w/2 - pan_w/2,
             yend = cy + post_h/2 - drop, color = NAVY, linewidth = 0.6),
    annotate("segment", x = cx - beam_w/2, y = cy + post_h/2, xend = cx - beam_w/2 + pan_w/2,
             yend = cy + post_h/2 - drop, color = NAVY, linewidth = 0.6),
    annotate("polygon", x = c(cx - beam_w/2 - pan_w/2, cx - beam_w/2 + pan_w/2, cx - beam_w/2),
             y = c(cy + post_h/2 - drop, cy + post_h/2 - drop, cy + post_h/2 - drop - pan_h),
             fill = GOLD, color = NAVY, linewidth = 0.5),
    annotate("text", x = cx - beam_w/2, y = cy + post_h/2 - drop - pan_h - 0.12*s,
             label = "Writers", color = "#B8860B", fontface = "bold", size = 2.7),
    # right hanger + pan (erasers, teal)
    annotate("segment", x = cx + beam_w/2, y = cy + post_h/2, xend = cx + beam_w/2 - pan_w/2,
             yend = cy + post_h/2 - drop, color = NAVY, linewidth = 0.6),
    annotate("segment", x = cx + beam_w/2, y = cy + post_h/2, xend = cx + beam_w/2 + pan_w/2,
             yend = cy + post_h/2 - drop, color = NAVY, linewidth = 0.6),
    annotate("polygon", x = c(cx + beam_w/2 - pan_w/2, cx + beam_w/2 + pan_w/2, cx + beam_w/2),
             y = c(cy + post_h/2 - drop, cy + post_h/2 - drop, cy + post_h/2 - drop - pan_h),
             fill = TEAL, color = NAVY, linewidth = 0.5),
    annotate("text", x = cx + beam_w/2, y = cy + post_h/2 - drop - pan_h - 0.12*s,
             label = "Erasers", color = "#2A7A6C", fontface = "bold", size = 2.7)
  )
}

W <- 12; H <- 7.4
p <- ggplot() +
  # ---- left data boxes ----
  box2(1.9, 6.35, 3.2, 1.55, "Single-cell & spatial",
       "GSE138794 (35,823 cells)\nGSE131928, GSE194329", "#E8E4F3") +
  box2(1.9, 4.55, 3.2, 1.55, "CPTAC-GBM (n = 99)",
       "Proteome, acetylome,\nphosphoproteome, RNA-seq", "#E4F0EC") +
  box2(1.9, 2.75, 3.2, 1.55, "Bulk cohorts",
       "TCGA, CGGA-693/325, GLASS,\nGSE7696, GSE83300", "#FDF3E3") +
  box2(1.9, 0.95, 3.2, 1.55, "Functional resources",
       "TCGASpliceSeq (PSI),\nDepMap, GDSC2", "#F3E4E4") +

  # ---- center WEBS box ----
  annotate("rect", xmin = 4.4, xmax = 7.6, ymin = 2.3, ymax = 5.9,
           fill = PURPLE_BG, color = PURPLE, linewidth = 1.2) +
  annotate("text", x = 6.0, y = 5.55, label = "WEBS", fontface = "bold", size = 5.2, color = NAVY) +
  annotate("text", x = 6.0, y = 5.05, label = "UCell(writers) \u2212 UCell(erasers)",
           size = 2.9, color = "grey15") +
  annotate("text", x = 6.0, y = 4.72, label = "LLS / LDTS / LHS =",
           size = 2.5, color = "grey30") +
  annotate("text", x = 6.0, y = 4.44, label = "low / intermediate / high-WEBS state",
           size = 2.5, color = "grey30") +
  balance(6.0, 3.35, 0.95) +

  # ---- right findings boxes ----
  box2(10.1, 6.35, 3.4, 1.55, "Malignant-state characterization",
       "Neftel states, pseudotime,\nspatial niches", "#E8E4F3") +
  box2(10.1, 4.55, 3.4, 1.55, "Spliceosome accumulation",
       "107/111 spliceosome proteins\nup-regulated (CPTAC)", "#E4F0EC") +
  box2(10.1, 2.75, 3.4, 1.55, "Splicing output",
       "4,387 WEBS-associated AS events\n(replicated in TCGA-LGG)", "#FDF3E3") +
  box2(10.1, 0.95, 3.4, 1.55, "Clinical association",
       "Clinical association &\ncandidate therapeutic vulnerabilities", "#F3E4E4") +

  # ---- arrows ----
  arr(3.5, 6.35, 4.4, 4.9) + arr(3.5, 4.55, 4.4, 4.4) +
  arr(3.5, 2.75, 4.4, 3.9) + arr(3.5, 0.95, 4.4, 3.4) +
  arr(7.6, 4.9, 8.4, 6.35) + arr(7.6, 4.4, 8.4, 4.55) +
  arr(7.6, 3.9, 8.4, 2.75) + arr(7.6, 3.4, 8.4, 0.95) +

  # ---- footer bar ----
  annotate("rect", xmin = 0.25, xmax = 11.75, ymin = -0.35, ymax = -1.05,
           fill = GREEN_BG, color = "#A5D6A7", linewidth = 0.8) +
  annotate("text", x = 6.0, y = -0.70, color = GREEN_TXT, size = 2.5,
           label = "Code and the WEBS implementation are publicly available (GitHub); datasets are accessible\nthrough the original repositories subject to their respective access requirements.") +

  coord_cartesian(xlim = c(0, W), ylim = c(-1.25, 7.4), expand = FALSE) +
  theme_void()

ggsave("figures_GND/Figure_1.pdf", p, width = 12, height = 7.4)
ggsave("figures_GND/Figure_1.png", p, width = 12, height = 7.4, dpi = 300)
cat("beautified Figure 1 saved\n")
