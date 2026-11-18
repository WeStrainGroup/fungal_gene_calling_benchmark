#!/usr/bin/env Rscript
# Species-equal means; no SD or error bars.
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]
out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
stopifnot(jsonlite::fromJSON(file.path(root, 'STATUS.json'))$status == 'COMPLETE')
d <- read.delim(file.path(root, 'species_mean.tsv'), na.strings='NA')
d$group <- factor(d$group, c('genome','50-100kb','5-10kb','1-2kb'),
                  c('Original','50-100 kb','5-10 kb','1-2 kb'))
d$tool <- factor(d$tool, c('annevo','tiberius','carbon'), c('ANNEVO','Tiberius','Carbon-A'))
d$level <- factor(d$level, c('nucleotide','exon','gene'), c('Nucleotide','Exon','Gene / CDS'))
colors <- c('ANNEVO'='#6E89A6','Tiberius'='#91AA9D','Carbon-A'='#C78180')
base <- theme_minimal(base_size=11) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        legend.position='top', plot.caption=element_text(size=9,hjust=0))
full <- subset(d, stratum=='complete' & metric=='f1')
p <- ggplot(full,aes(group,mean,fill=tool)) +
  geom_col(position=position_dodge(.82),width=.74,alpha=.9,na.rm=TRUE) +
  geom_text(aes(label=ifelse(is.na(mean),'NA',sprintf('%.1f',mean)),y=ifelse(is.na(mean),0,mean)),
            position=position_dodge(.82),vjust=-.5,size=2.3) +
  facet_wrap(~level,nrow=1) + scale_fill_manual(values=colors) +
  scale_y_continuous(limits=c(0,108),breaks=seq(0,100,25),expand=expansion(mult=c(0,.01))) +
  labs(x=NULL,y='F1 (%)',fill=NULL,title='Complete reference CDS',
       subtitle='Common input >=1 kb; CDS >=200 bp; Carbon-A confidence >=0.1 or NA',
       caption='Species-equal means. Strand-aware scoring regions include background; cross-scope predictions are reported separately.') + base
p2data <- subset(d, stratum!='complete' & metric=='f1' & group!='Original')
p2data$stratum <- factor(p2data$stratum,c('partial','missing_5prime','missing_3prime','missing_both'),
                         c('All truncated CDS',"Missing 5' end","Missing 3' end",'Missing both ends'))
p2 <- ggplot(p2data,aes(group,mean)) + geom_col(fill=colors['Carbon-A'],alpha=.9,width=.65,na.rm=TRUE) +
  geom_text(aes(label=ifelse(is.na(mean),'NA',sprintf('%.1f',mean)),y=ifelse(is.na(mean),0,mean)),vjust=-.5,size=2.5) +
  facet_grid(stratum~level) +
  scale_y_continuous(limits=c(0,108),breaks=seq(0,100,25),expand=expansion(mult=c(0,.01))) +
  labs(x=NULL,y='F1 (%)',title='Carbon-A: truncated reference CDS',
       subtitle='Conditional scores within truncated-reference regions; same confidence and length filters',
       caption='Species-equal means; NA means no eligible reference. Original assemblies have no clipped reference CDS.\nMixed-end isoforms and unassignable false positives remain in the aggregate diagnostics.') + base
for (name in c('complete_cds','carbon_partial_cds')) {
  plot <- if (name=='complete_cds') p else p2
  height <- if (name=='complete_cds') 4.8 else 9
  ggsave(file.path(out,paste0(name,'.pdf')),plot,width=12,height=height)
  ggsave(file.path(out,paste0(name,'.png')),plot,width=12,height=height,dpi=180)
}
