# =====================================================================
# MC-MED — Figure 2: distribution of the bolus effect on MAP at 60 min by
# phenotype and between-phenotype contrasts (half-eye plots of the
# cluster-bootstrap draws; region of practical equivalence +/- 2 mm Hg).
# Reads 05_effects.csv and 08_bootstrap_draws_primary.csv. Saves TIFF (LZW,
# 600 dpi), PNG and PDF. No title inside the figure; legend is separate.
# Font: Lato if installed (free, Google Fonts), otherwise Arial/Helvetica.
# Requires: data.table, ggplot2, ragg, systemfonts
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
home    <- if (.Platform$OS.type == "windows") Sys.getenv("USERPROFILE") else path.expand("~")
out_dir <- Sys.getenv("MCMED_OUT_DIR", unset = file.path(home, "Desktop", "mimiciv", "analysis_dap_pp_bolus_seqtrials"))
if (!dir.exists(out_dir)) out_dir <- getwd()
mcid <- 2
fams <- unique(systemfonts::system_fonts()$family)
fam  <- c("Lato", "Arial", "Helvetica", "Liberation Sans")[c("Lato", "Arial", "Helvetica", "Liberation Sans") %in% fams][1]
if (is.na(fam)) fam <- "sans"
minus <- function(x) gsub("-", "−", x, fixed = TRUE)

dr  <- fread(file.path(out_dir, "08_bootstrap_draws_primary.csv"))
eff <- fread(file.path(out_dir, "05_effects.csv"))[analysis == "Primary (sequential trials, IPW + outcome model)" & grepl("^MAP", outcome)]
grp <- data.table(
  term  = c("map_Overall", "map_Low pulse pressure", "map_Vasodilatory", "map_Mixed",
            "map_contrast_Vasodilatory vs Low pulse pressure", "map_contrast_Mixed vs Low pulse pressure"),
  group = c("Overall", "Low pulse pressure", "Vasodilatory", "Mixed",
            "contrast_Vasodilatory vs Low pulse pressure", "contrast_Mixed vs Low pulse pressure"),
  label = c("All trial-visits", "Low pulse pressure", "Vasodilatory", "Mixed",
            "Vasodilatory vs low pulse pressure", "Mixed vs low pulse pressure"),
  block = c(rep("bolus", 4), rep("contrast", 2)),
  y     = c(6.3, 5.2, 4.2, 3.2, 1.25, 0.25))
pal <- c("Overall" = "#3B3B3B", "Low pulse pressure" = "#0072B2", "Vasodilatory" = "#E69F00", "Mixed" = "#D55E00",
         "contrast_Vasodilatory vs Low pulse pressure" = "#8C8C8C", "contrast_Mixed vs Low pulse pressure" = "#8C8C8C")
ink <- "#222222"; mid <- "#666666"; light <- "#9A9A9A"

dens <- rbindlist(lapply(seq_len(nrow(grp)), function(i) {
  x <- dr[[grp$term[i]]]; x <- x[is.finite(x)]
  d <- density(x, adjust = 1.2, n = 1024)
  data.table(group = grp$group[i], x = d$x, h = d$y / max(d$y) * 0.72, y0 = grp$y[i])
}))
dens <- dens[x >= -4.4 & x <= 4.4]
dens[, region := fifelse(abs(x) >= mcid, "out", "in")]
dens[, seg := paste(group, rleid(region)), by = group]
pt <- merge(grp, eff[, .(group, est, lo, hi)], by = "group")
pt[, p_gt := sapply(term, function(t) mean(dr[[t]] > mcid, na.rm = TRUE))]
pt[, ci  := minus(sprintf("%.1f (%.1f to %.1f)", est, lo, hi))]
pt[, pg  := fifelse(block == "bolus", fifelse(p_gt < 0.01, "<1%", sprintf("%.0f%%", 100 * p_gt)), "")]

xl <- c(-4.4, 4.4); xc1 <- 5.0; xc2 <- 11.0        # text columns (data units)
p <- ggplot() +
  annotate("rect", xmin = -mcid, xmax = mcid, ymin = -Inf, ymax = 6.95, fill = "#EFEFEF") +
  annotate("segment", x = 0, xend = 0, y = -0.05, yend = 6.95, linewidth = 0.35, colour = mid) +
  geom_ribbon(data = dens, aes(x = x, ymin = y0, ymax = y0 + h, group = seg,
                               fill = group, alpha = region), colour = NA) +
  geom_line(data = dens, aes(x = x, y = y0 + h, group = group, colour = group), linewidth = 0.45) +
  geom_segment(data = pt, aes(x = lo, xend = hi, y = y, yend = y), linewidth = 0.9, colour = ink, lineend = "round") +
  geom_point(data = pt, aes(x = est, y = y), size = 2.4, shape = 21, fill = "white", colour = ink, stroke = 0.9) +
  # right-hand columns
  geom_text(data = pt, aes(x = xc1, y = y, label = ci), hjust = 0, size = 3.3, family = fam, colour = ink) +
  geom_text(data = pt, aes(x = xc2, y = y, label = pg), hjust = 1, size = 3.3, family = fam, colour = ink) +
  annotate("text", x = xc1, y = 7.25, label = "Effect, mm Hg (95% CI)", hjust = 0, size = 3.1, family = fam, colour = mid) +
  annotate("text", x = xc2, y = 7.25, label = "Estimates > 2 mm Hg", hjust = 1, size = 3.1, family = fam, colour = mid) +
  annotate("segment", x = xc1, xend = xc2, y = 6.98, yend = 6.98, linewidth = 0.3, colour = light) +
  # block headings and region label
  annotate("text", x = -9.6, y = 7.25, label = "EFFECT OF THE BOLUS", hjust = 0, size = 2.8, family = fam, colour = mid, fontface = "bold") +
  annotate("text", x = -9.6, y = 2.25, label = "DIFFERENCE BETWEEN PHENOTYPES", hjust = 0, size = 2.8, family = fam, colour = mid, fontface = "bold") +
  annotate("text", x = 0, y = 7.25, label = "Region of practical equivalence", size = 3.0, family = fam, colour = mid, fontface = "italic") +
  # row labels drawn as text so they align with the headings
  geom_text(data = grp, aes(x = -9.6, y = y, label = label), hjust = 0, size = 3.6, family = fam, colour = ink) +
  scale_y_continuous(expand = expansion(add = c(0.25, 0.1))) +
  scale_x_continuous(breaks = seq(-4, 4, 1), labels = function(z) minus(format(z))) +
  scale_fill_manual(values = pal, guide = "none") + scale_colour_manual(values = pal, guide = "none") +
  scale_alpha_manual(values = c(`in` = 0.30, out = 0.85), guide = "none") +
  coord_cartesian(xlim = c(-9.6, xc2), ylim = c(0, 7.4), clip = "off") +
  labs(x = "Effect of a crystalloid bolus on MAP at 60 min (mm Hg)", y = NULL) +
  theme_minimal(base_size = 10.5, base_family = fam) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = mid, size = 9),
        axis.title.x = element_text(colour = ink, size = 10, margin = margin(t = 6), hjust = 0.48),
        plot.margin = margin(10, 14, 6, 8))
# x axis drawn only under the effect scale
p <- p + annotate("segment", x = -4.4, xend = 4.4, y = -0.05, yend = -0.05, linewidth = 0.4, colour = mid) +
  annotate("segment", x = seq(-4, 4, 1), xend = seq(-4, 4, 1), y = -0.05, yend = -0.16, linewidth = 0.4, colour = mid)

w <- 10.5; h <- 5.2
ggsave(file.path(out_dir, "fig2_effect_distributions_map.tiff"), p, width = w, height = h, dpi = 600,
       device = ragg::agg_tiff, compression = "lzw", bg = "white")
ggsave(file.path(out_dir, "fig2_effect_distributions_map.png"),  p, width = w, height = h, dpi = 300, device = ragg::agg_png, bg = "white")
ggsave(file.path(out_dir, "fig2_effect_distributions_map.pdf"),  p, width = w, height = h, device = cairo_pdf)
message("Figure 2 saved (font: ", fam, ") in ", normalizePath(out_dir))
