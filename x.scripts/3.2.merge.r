### TOOLS
library(tidyverse)
source('functions.r')

### PYRHO FILES
ch <- read_delim('../1.data/2.pyrho/scaled_datasets/ch_scaled_all.txt')[,c(1,2,3,6)]
gb <- read_delim('../1.data/2.pyrho/scaled_datasets/gb_scaled_all.txt')[,c(1,2,3,6)]
pt <- read_delim('../1.data/2.pyrho/scaled_datasets/pt_scaled_all.txt')[,c(1,2,3,6)]
ch13_1 <- read_delim('../1.data/2.pyrho/scaled_datasets/ch13_1_scaled_all.txt')[,c(1,2,3,6)]
ch13_2 <- read_delim('../1.data/2.pyrho/scaled_datasets/ch13_2_scaled_all.txt')[,c(1,2,3,6)]
d <- full_join(ch,gb,by=c('ss','a','b'))
d <- full_join(d,pt,by=c('ss','a','b'))
d <- full_join(d,ch13_1,by=c('ss','a','b'))
d <- full_join(d,ch13_2,by=c('ss','a','b'))
d <- d[order(d$ss,d$a,d$b),]
names(d) <- c('ss','a','b','ch','gb','pt','ch13_1','ch13_2')
rm(ch,gb,pt,ch13_1,ch13_2)
d$ss <- rename.ss(d$ss)

### homogenize in windows 
k <- read_delim('../1.data/3.windows/genome_1kb_windows.bed',
                col_names=c('ss','start','end'))

### go from cM to cM/Mb 
drate <- d
for(i in 4:8){
tmp <- d[,c(2,3,i)]
indx <- which(!is.na(tmp[,3]))
drate[indx,i]<- apply(tmp[indx,],1,FUN=cmmb)
}
