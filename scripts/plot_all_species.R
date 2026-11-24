#!/usr/bin/env Rscript
# Draw complete and truncated CDS F1 from the frozen species summary.
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
d <- read.delim(file.path(args[1], "mean_sd.tsv"), na.strings = "NA")
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- subset(d, metric == "f1")
d$group <- factor(d$group, c("genome", "10-100kb", "5-10kb", "2-5kb", "1-2kb"),
                  c("Original", "10–100 kb", "5–10 kb", "2–5 kb", "1–2 kb"))
d$group <- droplevels(d$group)
d$level <- factor(d$level, c("nucleotide", "exon", "gene"),
                  c("Nucleotide", "Exon", "Gene / CDS"))
d$lower <- d$mean - d$sd
d$upper <- d$mean + ifelse(is.na(d$sd), 0, d$sd)
tool_colors <- c(ANNEVO = "#738B91", Tiberius = "#BDA982", `Carbon-A` = "#B76E60", Helixer = "#92919F")
type_colors <- c("#514A46", unname(tool_colors[1:3]))

draw <- function(x, title, colors) {
  dodge <- position_dodge(width = .80)
  ggplot(x, aes(group, mean, fill = category)) +
    geom_col(position = dodge, width = .72, alpha = .88, na.rm = TRUE) +
    geom_errorbar(aes(ymin = lower, ymax = mean + sd), position = dodge,
                  width = .13, linewidth = .40, colour = "#444A50", na.rm = TRUE) +
    geom_text(aes(y = ifelse(is.na(mean), 0, upper),
                  label = ifelse(is.na(mean), "NA", sprintf("%.1f", mean))),
              position = dodge, vjust = -.65, size = if (nlevels(d$group) > 4) 2.2 else 2.7, colour = "#30363B") +
    facet_wrap(~level, nrow = 1) +
    scale_fill_manual(values = if (is.null(names(colors))) setNames(colors, levels(x$category)) else colors, drop = FALSE) +
    scale_y_continuous(breaks = seq(0, 125, 25), expand = expansion(mult = c(.01, .14))) +
    labs(title = title, x = NULL, y = "F1 (%)", fill = NULL) +
    theme_minimal(base_size = 11, base_family = "sans") +
    theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(colour = "#E8EBED", linewidth = .35),
          panel.spacing = grid::unit(1.2, "lines"),
          strip.text = element_text(face = "bold", size = 11, margin = margin(b = 10)),
          legend.position = "top", legend.key.height = grid::unit(.32, "cm"),
          legend.key.width = grid::unit(.55, "cm"),
          legend.margin = margin(t = 2, b = 6),
          plot.title = element_text(face = "bold", size = 15, margin = margin(b = 8)),
          plot.margin = margin(12, 14, 10, 12),
          axis.text = element_text(colour = "#454B50"))
}

full <- subset(d, stratum == "complete")
full$category <- droplevels(factor(full$tool, c("helixer", "annevo", "tiberius", "carbon"),
                                  c("Helixer", "ANNEVO", "Tiberius", "Carbon-A")))
stopifnot(nrow(full) == 3 * nlevels(full$category) * nlevels(d$group),
          all(is.na(full$mean) == (full$n_species == 0)))
plots <- list(complete_cds = draw(full, "Complete reference CDS", tool_colors))
for (tool in intersect(c("helixer", "carbon"), unique(d$tool))) {
  partial <- d[d$tool == tool & d$group != "Original" &
                 d$stratum %in% c("partial", "missing_5prime", "missing_3prime", "missing_both"), ]
  partial$group <- droplevels(partial$group)
  partial$category <- factor(partial$stratum,
                             c("partial", "missing_5prime", "missing_3prime", "missing_both"),
                             c("Overall", "Missing 5′ end", "Missing 3′ end", "Missing both ends"))
  stopifnot(nrow(partial) == 12 * (nlevels(d$group) - 1),
            all(is.na(partial$mean) == (partial$n_species == 0)))
  label <- if (tool == "carbon") "Carbon-A" else "Helixer"
  plots[[paste0(tool, "_partial_cds")]] <- draw(partial, paste0(label, ": truncated reference CDS"), type_colors)
}
if (all(c("carbon", "helixer") %in% unique(d$tool))) {
  paired <- subset(d, stratum == "partial" & group != "Original")
  paired$group <- droplevels(paired$group)
  paired$category <- factor(paired$tool, c("helixer", "carbon"), c("Helixer", "Carbon-A"))
  plots$truncated_cds <- draw(paired, "Truncated reference CDS", tool_colors)
}
for (name in names(plots)) {
  ggsave(file.path(out, paste0(name, ".png")), plots[[name]],
         width = 18, height = 4.4, dpi = 220, bg = "white")
  ggsave(file.path(out, paste0(name, ".pdf")), plots[[name]],
         width = 18, height = 4.4, device = cairo_pdf, bg = "white")
}
