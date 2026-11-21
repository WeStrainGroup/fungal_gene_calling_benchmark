#!/usr/bin/env Rscript
# Completed species; bars show means and error bars show sample SD.
library(ggplot2)
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]
out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
status <- jsonlite::fromJSON(file.path(root, 'STATUS.json'))
stopifnot(status$status %in% c('COMPLETE','INTERIM'))
n_species <- if (is.null(status$included_species)) status$expected_species else status$included_species
scope <- if (status$status=='INTERIM') 'Interim: only species-condition pairs with all three tools validated' else sprintf('%d species; all tools and conditions complete',n_species)
groups <- c('genome','50-100kb','5-10kb','1-2kb')
group_labels <- c('Original','50-100 kb','5-10 kb','1-2 kb')
if (!is.null(status$species_by_condition)) {
  counts <- vapply(groups,function(g) length(status$species_by_condition[[g]]),integer(1))
  group_labels <- paste0(group_labels,'\n(n=',counts,')')
}
d <- read.delim(file.path(root, 'mean_sd.tsv'), na.strings='NA')
d$upper <- d$mean + ifelse(is.na(d$sd), 0, d$sd)
d$lower <- d$mean - d$sd
d$original <- d$group=='genome'
d$group <- factor(d$group, groups, group_labels)
d$tool <- factor(d$tool, c('annevo','tiberius','carbon'), c('ANNEVO','Tiberius','Carbon-A'))
d$level <- factor(d$level, c('nucleotide','exon','gene'), c('Nucleotide','Exon','Gene / CDS'))
colors <- c('ANNEVO'='#6E89A6','Tiberius'='#91AA9D','Carbon-A'='#C78180')
base <- theme_minimal(base_size=11) +
  theme(panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
        legend.position='top', plot.caption=element_text(size=9,hjust=0))
full <- subset(d, stratum=='complete' & metric=='f1')
p <- ggplot(full,aes(group,mean,fill=tool)) +
  geom_col(position=position_dodge(.82),width=.74,alpha=.9,na.rm=TRUE) +
  geom_errorbar(aes(ymin=lower,ymax=mean+sd),position=position_dodge(.82),width=.15,na.rm=TRUE) +
  geom_text(aes(label=ifelse(is.na(mean),'NA',sprintf('%.1f',mean)),y=ifelse(is.na(mean),0,upper)),
            position=position_dodge(.82),vjust=-.5,size=2.3) +
  facet_wrap(~level,nrow=1) + scale_fill_manual(values=colors) +
  scale_y_continuous(expand=expansion(mult=c(.03,.12))) +
  labs(x=NULL,y='F1 (%)',fill=NULL,title='Complete reference CDS',
       subtitle=scope,
       caption='Species-equal mean +/- sample SD. n: species with all three tools validated in that condition.\nMetric-specific defined counts are reported in the table.') + base
p2data <- subset(d, stratum!='complete' & metric=='f1' & !original)
p2data$stratum <- factor(p2data$stratum,c('partial','missing_5prime','missing_3prime','missing_both'),
                         c('All truncated CDS',"Missing 5' end","Missing 3' end",'Missing both ends'))
p2 <- ggplot(p2data,aes(group,mean)) + geom_col(fill=colors['Carbon-A'],alpha=.9,width=.65,na.rm=TRUE) +
  geom_errorbar(aes(ymin=lower,ymax=mean+sd),width=.15,na.rm=TRUE) +
  geom_text(aes(label=ifelse(is.na(mean),'NA',sprintf('%.1f',mean)),y=ifelse(is.na(mean),0,upper)),vjust=-.5,size=2.5) +
  facet_grid(stratum~level) +
  scale_y_continuous(expand=expansion(mult=c(.03,.12))) +
  labs(x=NULL,y='F1 (%)',title='Carbon-A: truncated reference CDS',
       subtitle=paste(scope,'\nConditional scores within truncated-reference regions'),
       caption='Species-equal mean +/- sample SD; conditional truncated-region scores. SD is undefined with fewer than two defined species.\nAxis n: species completing all three tools in that condition; per-metric defined counts are reported in the table.') + base
for (name in c('complete_cds','carbon_partial_cds')) {
  plot <- if (name=='complete_cds') p else p2
  height <- if (name=='complete_cds') 4.8 else 9
  ggsave(file.path(out,paste0(name,'.pdf')),plot,width=12,height=height)
  ggsave(file.path(out,paste0(name,'.png')),plot,width=12,height=height,dpi=180)
}
