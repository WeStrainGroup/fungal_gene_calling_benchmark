#!/usr/bin/env Rscript
# Plot equal-species means after taking each input's replicate median.
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
d <- read.delim(args[1])
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
stopifnot(nrow(d) == 20, all(d$species == 17), length(unique(d$replicates)) == 1)
rounds <- unique(d$replicates)
stopifnot(rounds %in% c(1, 3))
d$group <- factor(d$group, c("genome", "10-100kb", "5-10kb", "2-5kb", "1-2kb"),
                  c("Original", "10–100 kb", "5–10 kb", "2–5 kb", "1–2 kb"))
d$tool <- factor(d$tool, c("annevo", "tiberius", "helixer", "carbon"),
                 c("ANNEVO", "Tiberius", "Helixer", "Carbon-A"))
colors <- c(ANNEVO = "#738B91", Tiberius = "#BDA982", Helixer = "#92919F", `Carbon-A` = "#B76E60")
subtitle <- if (rounds == 3) "Mean ± SD across 17 species; each input represented by its median of three runs" else
  "Round 1 preview: mean ± SD across 17 species; three-run results pending"
format_value <- function(x) ifelse(x < 1, sprintf("%.2f", x), sprintf("%.1f", x))
base <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#E8EBED", linewidth = .35),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text = element_text(colour = "#454B50"),
        strip.text = element_text(face = "bold", margin = margin(b = 9)),
        legend.position = "top", plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(size = 10, colour = "#646A70"),
        plot.margin = margin(12, 15, 10, 12), panel.spacing = grid::unit(1.2, "lines"))
p <- ggplot(d, aes(tool, mean_mbp_per_minute, fill = tool)) +
  geom_col(width = .68, alpha = .85) +
  geom_errorbar(aes(ymin = mean_mbp_per_minute - sd_mbp_per_minute,
                    ymax = mean_mbp_per_minute + sd_mbp_per_minute),
                width = .16, linewidth = .45, colour = "#555B60") +
  geom_text(aes(y = mean_mbp_per_minute + sd_mbp_per_minute,
                label = format_value(mean_mbp_per_minute)), vjust = -.5, size = 3, colour = "#42474B") +
  facet_wrap(~group, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = colors) +
  scale_y_continuous(expand = expansion(mult = c(0, .20))) +
  labs(title = "Gene prediction throughput", subtitle = subtitle, x = NULL,
       y = "Throughput (Mb/min)", fill = NULL) + base
ggsave(file.path(out, "throughput.png"), p, width = 12.5, height = 4.1, dpi = 200, bg = "white")
ggsave(file.path(out, "throughput.pdf"), p, width = 12.5, height = 4.1, device = cairo_pdf, bg = "white")
