### TOOLS
library(tidyverse)
source('functions.r')
pops <- c('CH','PT','GB','CH13','CH13_2','CH13_3','CH13_4','CH13_5')
sizes <- c('1kb','10kb','100kb','1Mb','5Mb')
size='1kb'
### PYRHO FILES
for(size in sizes){

  
for(pop in pops){
  
  if(size != '1kb'){
  dt <- read_delim(paste0('../1.data/2.pyrho/scaled_datasets/',pop,'_',size,'_scaled.txt'))
  names(dt)[4] <- paste0(pop,'_cM')
  } else { 
  dt <- read_delim(paste0('../1.data/2.pyrho/scaled_datasets/',pop,'_total.theta.intersect.bed_1000.thetaxN.scaled'))
  }
  
  if(pop != pops[1]){
    d <- left_join(d,dt,by=c('ss','start','end'))
  } else{
    d <- dt 
  }
}
if(size != '5Mb'){
### PI 
if(size != '1kb'){
  for(pop in pops){
    pi <- read_delim(paste0('../1.data/3.windows/pi/',size,'/',pop,'_merged_corrected_',size,'_pi.table'))[,seq(1,5)]
    names(pi) <- c('ss','start','end',paste0(pop,'_snps'),paste0(pop,'_pi'))
    pi$ss[pi$ss =='Super-Scaffold_2' & pi$end < len2.2] <- 'Super-Scaffold_2.2'
    pi$start <- pi$start - 1 
    pi <- pi[,c(1,2,4,5)]
    # merge 
    d <- dplyr::left_join(d,pi,by=c('ss','start'))
    }
}

### GC
gc <- read_delim(paste0('../1.data/3.windows/GC_windows/GC_',size))[,c(1,2,3,5)]
names(gc) <- c('ss','start','end','gc')
gc$ss[gc$ss =='Super-Scaffold_2' & gc$end < len2.2] <- 'Super-Scaffold_2.2'
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
tmp$ss[tmp$ss =='Super-Scaffold_2' & tmp$end < len2.2] <- 'Super-Scaffold_2.2'
tmp <- tmp[,-1]
names(tmp)[1] <- 'cpgi'
# merge
d <- dplyr::left_join(d, tmp, by=c('ss','start','end'))

### Genes 
tss <- read_delim(paste0('../1.data/3.windows/Gene_windows/',size,'_tss.table'),
                  col_names=c('ss','start','end','tss'))
tss$ss[tss$ss=='Super-Scaffold_2' & tss$end < len2.2] <- 'Super-Scaffold_2.2'
d <- dplyr::left_join(d, tss, by=c('ss','start','end'))
g <- read_delim(paste0('../1.data/3.windows/Gene_windows/',size,'_genes.table'),
                col_names=c('ss','start','end','tmp1','tmp2','tmp3','overlap'))
tmp <- g %>% group_by(ss,start,end) %>% summarize(Lgenes=sum(overlap),ngenes=n())
tmp$ss[tmp$ss =='Super-Scaffold_2' & tmp$end < len2.2] <- 'Super-Scaffold_2.2'

d <- dplyr::left_join(d, tmp, by=c('ss','start','end'))

# Repeats 
r <- read_delim(paste0('../1.data/3.windows/TE_windows/',size,'_repeats.table'),
                col_names=c('ss','start','end','tmp1','tmp2','tmp3','overlap'))
tmp <- r %>% group_by(ss,start,end) %>% summarize(Lrepeats=sum(overlap),nrepeats=n())
tmp$ss[tmp$ss =='Super-Scaffold_2' & tmp$end < len2.2] <- 'Super-Scaffold_2.2'

d <- dplyr::left_join(d, tmp, by=c('ss','start','end'))

### PED for 1 Mb
if(size == '1Mb'){
  lm <- read_delim('../1.data/1.LepMAP/R.AllSNPs_physical/physical_orders_full.txt')
  lm$ss[lm$lg==24] <- 'Super-Scaffold_2.2'
  win <- read_delim('../1.data/3.windows/genome_1mb_windows.bed', 
                    col_names = c('ss','start','end'))
  win$ss[win$end < len2.2 & win$ss=='Super-Scaffold_2'] <- 'Super-Scaffold_2.2'
  win <- cbind.data.frame(win,window.ped(win,lm,sex.based=T))
  win$av <- window.ped(win,ped = lm)
  d <- left_join(d,win,by=c('ss','start','end'))
}
} else {
  lm <- read_delim('../1.data/1.LepMAP/R.AllSNPs_physical/physical_orders_full.txt')
  lm$ss[lm$lg==24] <- 'Super-Scaffold_2.2'
  win <- read_delim('../1.data/3.windows/genome_5mb_windows.bed', 
                    col_names = c('ss','start','end'))
  win$ss[win$end < len2.2 & win$ss=='Super-Scaffold_2'] <- 'Super-Scaffold_2.2'
  win <- cbind.data.frame(win,window.ped(win,lm,sex.based=T))
  win$av <- window.ped(win,ped = lm)
  d <- left_join(d,win,by=c('ss','start','end'))
}
  
  write.table(d,file=paste0('../1.data/3.windows/total_',size,'.table'),quote=F,row.names=F)
  
}
