### TOOLS
library(tidyverse)
source('functions.r')
pops <- c('CH','PT','GB','CH13','CH13_2')
sizes <- c('10kb','100kb','1Mb','5Mb')
  
### PYRHO FILES
for(size in sizes){

for(pop in pops){
  dt <- read_delim(paste0('../1.data/2.pyrho/scaled_datasets/',pop,'_',size,'_scaled.txt'))[,c(1,2,3,6)]
  names(dt)[4] <- paste0(pop,'_cM')
  if(pop != pops[1]){
    d <- left_join(d,dt,by=c('ss','start','end'))
  } else{
    d <- dt 
  }
}

if(size != '5Mb'){
### PI 
for(pop in pops){
pi <- read_delim(paste0('../1.data/3.windows/pi/',size,'/',pop,'_merged_corrected_',size,'_pi.table'))[,seq(1,5)]
names(pi) <- c('ss','start','end',paste0(pop,'_snps'),paste0(pop,'_pi'))
pi$start <- pi$start - 1 
# merge 
d <- dplyr::left_join(d,pi,by=c('ss','start','end'))
}

### GC
gc <- read_delim(paste0('../1.data/3.windows/GC_windows/GC_',size))[,c(1,2,3,5)]
names(gc) <- c('ss','start','end','gc')
# merge
d <- dplyr::left_join(d,gc,by=c('ss','start','end'))

### CpGIs
cpgi <- read_delim((paste0('../1.data/3.windows/GC_windows/cpgi_',size)),
                   col_names = c('ss','start','end','tmp1','tmp2','tmp3','overlap'))
cpgi$tag <- paste0(cpgi$ss,cpgi$end)
tmp <- aggregate(cpgi$overlap,by=list(cpgi$tag), FUN=sum, na.rm=T)
tmp$ss <- cpgi$ss[match(tmp$Group.1,cpgi$tag)]
tmp$start <- cpgi$start[match(tmp$Group.1,cpgi$tag)]
tmp$end <- cpgi$end[match(tmp$Group.1,cpgi$tag)]
tmp <- tmp[,-1]
names(tmp)[1] <- 'cpgi'
# merge
d <- dplyr::left_join(d, tmp, by=c('ss','start','end'))

### Genes 
tss <- read_delim(paste0('../1.data/3.windows/Gene_windows/',size,'_tss.table'),
                  col_names=c('ss','start','end','tss'))
d <- dplyr::left_join(d, tss, by=c('ss','start','end'))
tes <- read_delim(paste0('../1.data/3.windows/Gene_windows/',size,'_tes.table'),
                  col_names=c('ss','start','end','tes'))
d <- dplyr::left_join(d, tes, by=c('ss','start','end'))
g <- read_delim(paste0('../1.data/3.windows/Gene_windows/',size,'_genes.table'),
                col_names=c('ss','start','end','tmp1','tmp2','tmp3','overlap'))
tmp <- g %>% group_by(ss,start,end) %>% summarize(genes=sum(overlap))
d <- dplyr::left_join(d, tmp, by=c('ss','start','end'))


### PED for 1 Mb
if(size == '1Mb'){
  lm <- read_delim('../1.data/1.LepMAP/orders/fitted_values_gams.table')
  win <- read_delim('../1.data/3.windows/genome_1mb_windows.bed', 
                    col_names = c('ss','start','end'))
  win$ss <- rename.ss(win$ss)
  win$ped <- window.ped(win,ped = lm)
  win$ss <- gsub('ss','Super-Scaffold_',win$ss)
  d <- left_join(d,win,by=c('ss','start','end'))
}
} else {
  lm <- read_delim('../1.data/1.LepMAP/orders/fitted_values_gams.table')
  win <- read_delim('../1.data/3.windows/genome_5mb_windows.bed', 
                    col_names = c('ss','start','end'))
  win$ss <- rename.ss(win$ss)
  win$ped <- window.ped(win,ped = lm)
  win$ss <- gsub('ss','Super-Scaffold_',win$ss)
  d <- left_join(d,win,by=c('ss','start','end'))
}
  
  write.table(d,file=paste0('../1.data/3.windows/total_',size,'.table'),quote=F,row.names=F)
  
}
