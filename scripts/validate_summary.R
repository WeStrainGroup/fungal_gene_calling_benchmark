args <- commandArgs(trailingOnly=TRUE)
root <- args[1]
m <- read.delim(file.path(root, 'metrics.tsv'), na.strings='NA')
s <- read.delim(file.path(root, 'mean_sd.tsv'), na.strings='NA')
status <- jsonlite::fromJSON(file.path(root,'STATUS.json'))
key <- c('species','group','tool','stratum','level')
stopifnot(!anyDuplicated(m[key]))
close <- function(a,b) {if (is.na(a) || is.na(b)) return(is.na(a) && is.na(b)); abs(a-b)<1e-9}
for (i in seq_len(nrow(s))) {
 r <- s[i,]
 rows <- m[m$group==r$group & m$tool==r$tool & m$stratum==r$stratum & m$level==r$level,]
 stopifnot(setequal(rows$species,status$species_by_condition[[r$group]]))
 v <- rows[[r$metric]]
 v <- v[!is.na(v)]
 stopifnot(length(v)==r$n_species,nrow(rows)==r$total_species)
 stopifnot(close(if(length(v)) mean(v) else NA_real_,r$mean))
 stopifnot(close(if(length(v)>1) sd(v) else NA_real_,r$sd))
}
cases <- unique(m[c('species','group','tool')])
units <- unique(m[c('species','group')])
stopifnot(nrow(cases)==length(unique(m$tool))*nrow(units), nrow(cases)==status$completed_cases)
receipt <- list(status='PASS',summary_rows_independently_checked=nrow(s),mean_sd_engine='R base mean and sample sd',paired_cases=nrow(cases),paired_species_conditions=nrow(units),scoring_version=status$scoring_version,validated_at=format(Sys.time(),tz='Asia/Shanghai',usetz=TRUE))
jsonlite::write_json(receipt,file.path(root,'validation.json'),pretty=TRUE,auto_unbox=TRUE)
print(receipt)
