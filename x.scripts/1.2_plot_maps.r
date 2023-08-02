library(readr)
library(tidyverse)
library(viridis)
source('functions.r')

##################################################

# read real distribution
real <- read_delim('../1.data/1.LepMAP/maps/Call2_tf1_maf05_miss05_1kb_pruned_f.call_dist.csv',
                   col_names = c('nvar','ss'))
real$ss <- rename.ss(real$ss)
# order
real <- real[order(real$nvar,decreasing=T),]
# remove some small shits 
real <- real[real$nvar > 20,]
# add colour pallete 
real$col <- colorRampPalette(c('#472183','#82C3EC','#F1F6F5'))(nrow(real))
pdf('plots/distr_lod_all.pdf',width=12,height=12)
par(mfrow=c(4,3))
# true distribution
barplot(real$nvar, ylab='# of SNPs', ylim=c(0,14000),las=2,
        col=real$col, main='real distribution of SNPs',space=0,border = NA)
axis(1,at=seq(0,nrow(real),10),labels = seq(0,nrow(real),10),las=2)
# loop through lods 
for(lod in seq(11,21)){
  # read file for lod
  file <- list.files(path='../1.data/1.LepMAP/maps/',pattern=paste0('*lod',lod,'_theta0.03_map.txt.positions'))
  tmp <- read_delim(paste0('../1.data/1.LepMAP/maps/',file),col_names=c('ss','pos','lg'))
  # rename 
  tmp$ss <- rename.ss(tmp$ss)
  # rm small shits 
  tmp <- tmp[tmp$ss %in% real$ss,]
  # summarize 
  tmp2 <- tmp %>% group_by(lg,ss) %>% summarize('nvar'=length(pos))
  tmp2 <- tmp2[order(tmp2$ss),]
  maxx <- nrow(real) #ifelse(length(unique(tmp2$lg)) < nrow(real), nrow(real), length(unique(tmp2$lg)))
  maxy <- max(aggregate(tmp2$nvar, by=list(tmp2$lg), FUN=sum))
  # initiate plot with right dimensions 
  plot(0, pch='', 
       xlim=c(0, maxx),
       ylim=c(0, maxy), 
       ylab='',
       xlab='',
       main=paste('LOD',lod),
       axes=F)
  # start stacking / lg
  for(i in seq(1,maxx)){
    xlims <- c(i - .5, i + .5)
    lgd <- tmp2[tmp2$lg == unique(tmp2$lg)[order(unique(tmp2$lg))][i], ]
    bot = 0
    lgd <- lgd[order(lgd$nvar, decreasing=T), ]
    lgd$cum <- cumsum(lgd$nvar)
    for(j in 1:nrow(lgd)){
      top = lgd$cum[j]
      rect(xlims[1], bot, xlims[2], top, col=real$col[match(lgd$ss[j], real$ss)],
           border=NA)
      bot = top
    }
  }
  axis(1, at=seq(0, nrow(real), 10), labels = seq(0,nrow(real), 10), las=2)
  axis(2, at=seq(0, maxy, maxy/4), labels = round(seq(0, maxy, maxy/4),-1), las=2 ) 
}
dev.off()

##################################################

# Plot on lod 15
lod = 15
pdf('plots/lod15.pdf',width=6,height=6)
file <- list.files(path='../1.data/1.LepMAP/maps/',pattern=paste0('*lod',lod,'_theta0.03_map.txt.positions'))
tmp <- read_delim(paste0('../1.data/1.LepMAP/maps/',file),col_names=c('ss','pos','lg'))
# rename 
tmp$ss <- rename.ss(tmp$ss)
# rm small shits 
tmp <- tmp[tmp$ss %in% real$ss,]
# summarize 
tmp2 <- tmp %>% group_by(lg,ss) %>% summarize('nvar'=length(pos))
tmp2 <- tmp2[order(tmp2$ss),]
maxx <- nrow(real) 
maxy <- max(aggregate(tmp2$nvar, by=list(tmp2$lg), FUN=sum))
# initiate plot with right dimensions 
plot(0, pch='', 
     xlim=c(0, maxx),
     ylim=c(0, maxy), 
     ylab='',
     xlab='',
     main=paste('LOD',lod),
     axes=F)
# start stacking / lg
for(i in seq(1,maxx)){
  xlims <- c(i - .5, i + .5)
  lgd <- tmp2[tmp2$lg == unique(tmp2$lg)[order(unique(tmp2$lg))][i], ]
  bot = 0
  lgd <- lgd[order(lgd$nvar, decreasing=T), ]
  lgd$cum <- cumsum(lgd$nvar)
  for(j in 1:nrow(lgd)){
    top = lgd$cum[j]
    rect(xlims[1], bot, xlims[2], top, col=real$col[match(lgd$ss[j], real$ss)],
         border=NA)
    bot = top
  }
}
axis(1, at=seq(0, nrow(real), 10), labels = seq(0,nrow(real), 10), las=2)
axis(2, at=seq(0, maxy, maxy/4), labels = round(seq(0, maxy, maxy/4),-1), las=2 ) 

dev.off()

##################################################

### Plot unique ss gain / LG

count.uniq <- function(df){
  uniq.vec <- rep(0,maxLGs)
  tmp.book <- ''
  for(i in 1:maxLGs){
    lgs <- as.character(df$ss[df$lg==i])
    new.lgs <- lgs[! lgs %in% tmp.book]
    tmp.book <- append(tmp.book, new.lgs)
    uniq.vec[i] <- length(new.lgs)
  }
  return(cumsum(as.numeric(uniq.vec)))
}
pdf('plots/plot_cumulative_lines.pdf',width=6,height=6)
plot(0,pch='', 
     ylab='# of unique scaffolds added', ylim=c(0,50), 
     xlim=c(0,50), xlab='LG',axes=F)
axis(1,at=seq(0,50,10),labels = seq(0,50,10))
axis(2,at=seq(0,50,10),labels = seq(0,50,10))
lods <- seq(11,21)
p <- c('grey85','grey60','grey30','grey1', viridis(7))
for(lod in lods){
  file <- list.files(path='../1.data/1.LepMAP/maps/',pattern=paste0('*lod',lod,'_theta0.03_map.txt.positions'))
  tmp <- read_delim(paste0('../1.data/1.LepMAP/maps/',file),col_names=c('ss','pos','lg'))
  # rename 
  tmp$ss <- rename.ss(tmp$ss)
  # rm small shits 
  tmp <- tmp[tmp$ss %in% real$ss,]
  # summarize 
  tmp2 <- tmp %>% group_by(lg,ss) %>% summarize('nvar'=length(pos))
  tmp2 <- tmp2[order(tmp2$ss),]
  maxLGs <- ifelse(length(unique(tmp2$lg)) > nrow(real), nrow(real), length(unique(tmp2$lg)))
  lines(count.uniq(tmp2)~seq(1:maxLGs),col=p[which(lods==lod)],lwd=2.5)
}
abline(0,1,lty=2)
legend('bottomright',lwd=3,col=p,legend=lods, title = 'LOD')
dev.off()
