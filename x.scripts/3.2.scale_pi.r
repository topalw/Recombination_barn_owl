### TOOLS
library(tidyverse)
source('functions.r')
sizes <- c('10kb','100kb','1Mb','5Mb')
pops <- c('CH','PT','GB','CH13','CH13_2','CH13_3','CH13_4','CH13_5')
size.translator <- c(10000,100000,1000000,5000000)
options(scipen=100000)

for(size in sizes){
  tmps <- size.translator[sizes==size]
  path=paste0('../1.data/3.windows/pi/',size,'/')
  mask <- read_delim(paste0('../1.data/3.windows/mask/mask_overlap_',size,'.bed'),
                     col_names = c('ss1','start1','end1',
                                   'ss2','start2','end2','overlap'))
  mask$tag <- paste0(mask$ss1,':',mask$end1)
  w <- mask$end1[1] - mask$start1[1]
  mask <- aggregate(mask$overlap,by=list(mask$tag),FUN=sum)
  # should multiply all estimates by win/(win-overlap) but make sure to check if mask=w
  # in that case set pi to  NA since we cannot tell how many snps are there 
  mask$scale <- ifelse(w/(w-mask$x)!=Inf,w/(w-mask$x),NA)
  for(pop in pops){
    tmp <- list.files(path,pattern=paste0(pop,'_Super'))
    d <- read_delim(paste0(path,tmp))
    d$tag <- paste0(d$CHROM,':',d$BIN_END)
    d$pi.corrected <- d$PI*mask$scale[match(d$tag,mask$Group.1)]
    write.table(d,paste0(path,pop,'_merged_corrected_',size,'_pi.table'),quote=F,row.names=F)
  }
}
