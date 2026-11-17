#!/usr/bin/env Rscript
# Three evaluation levels, with consistent colors across panels.
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[1] else 'results/benchmark.tsv'
output <- if (length(args) >= 2) args[2] else 'figures'
dir.create(output, recursive = TRUE, showWarnings = FALSE)
data <- read.delim(input)
data$group <- factor(data$group, c('genome', '50-100kb', '5-10kb', '1-2kb'),
                     c('Original', '50-100 kb', '5-10 kb', '1-2 kb'))
data$tool <- factor(data$tool, c('annevo', 'tiberius', 'carbon'), c('ANNEVO', 'Tiberius', 'Carbon-A'))
data$species <- factor(data$species, c('Saccharomyces_cerevisiae', 'Candida_albicans'),
                       c('S. cerevisiae', 'C. albicans'))
data$level <- factor(data$level, c('nucleotide', 'exon', 'gene'), c('Nucleotide', 'Exon', 'Gene'))
colors <- c('ANNEVO' = '#6E89A6', 'Tiberius' = '#91AA9D', 'Carbon-A' = '#C78180')
p <- ggplot(data, aes(group, f1, fill = tool)) +
  geom_col(position = position_dodge(width = .82), width = .74, alpha = .90) +
  geom_text(aes(label = sprintf('%.1f', f1)), position = position_dodge(width = .82),
            vjust = -.4, size = 2.55) +
  facet_grid(level ~ species) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0, 108), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
  labs(x = NULL, y = 'F1 (%)', fill = NULL,
       caption = 'Full nuclear genomes; one fixed partition per length range.\nOfficial Carbon-A metric definitions; no external confidence or CDS-length filter.') +
  theme_minimal(base_size = 11) +
  theme(legend.position = 'top', panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), strip.text.y = element_text(angle = 0),
        plot.caption = element_text(size = 9, hjust = 0), panel.spacing = unit(1.1, 'lines'))
ggsave(file.path(output, 'benchmark_f1.pdf'), p, width = 12, height = 8)
ggsave(file.path(output, 'benchmark_f1.png'), p, width = 12, height = 8, dpi = 180)
