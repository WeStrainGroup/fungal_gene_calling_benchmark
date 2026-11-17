#!/usr/bin/env Rscript
# Means across equally weighted species.
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[1] else 'results/species_summary.tsv'
output <- if (length(args) >= 2) args[2] else 'figures'
dir.create(output, recursive = TRUE, showWarnings = FALSE)
data <- read.delim(input)
data <- subset(data, metric == 'f1')
stopifnot(length(unique(data$n_species)) == 1)
data$group <- factor(data$group, c('genome','50-100kb','5-10kb','1-2kb'), c('Original','50-100 kb','5-10 kb','1-2 kb'))
data$tool <- factor(data$tool, c('annevo','tiberius','carbon'), c('ANNEVO','Tiberius','Carbon-A'))
data$level <- factor(data$level, c('nucleotide','exon','gene'), c('Nucleotide','Exon','Gene'))
colors <- c('ANNEVO'='#6E89A6','Tiberius'='#91AA9D','Carbon-A'='#C78180')
p <- ggplot(data, aes(group, mean, fill=tool)) +
  geom_col(position=position_dodge(.82), width=.74, alpha=.9) +
  geom_text(aes(label=sprintf('%.1f',mean)), position=position_dodge(.82), vjust=-.4, size=2.3) +
  facet_wrap(~level, nrow=1) +
  scale_fill_manual(values=colors) +
  scale_y_continuous(breaks=seq(0,100,25), expand=expansion(mult=c(.02,.07))) +
  labs(x=NULL,y='F1 (%)',fill=NULL,
       caption=paste0('Mean across ',unique(data$n_species),' species; each species has equal weight. No pooling across benchmark sets.')) +
  theme_minimal(base_size=11) +
  theme(legend.position='top', panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        plot.caption=element_text(size=9,hjust=0), panel.spacing=unit(1.2,'lines'))
ggsave(file.path(output,'benchmark_mean.pdf'),p,width=12,height=4.7)
ggsave(file.path(output,'benchmark_mean.png'),p,width=12,height=4.7,dpi=180)
