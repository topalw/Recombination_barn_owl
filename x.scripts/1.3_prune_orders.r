args <- commandArgs(trailingOnly=T)
library(tidyverse)
source('functions.r')
shit <- c('ss24','ss25','ss100000100078','ss200000178', 'ss100000100064', 'ss4')
p <- c('orange','blue')
focus <- args[1]
files <- fs::dir_ls(path= paste0('../1.data/1.LepMAP/orders/',focus), glob='*.perf')
# this will be the unpruned dataset 
unpruned  <- data.frame('ss'=character(),'pos'=integer(),
                        'm'=double(),'f'=double(),'lg'=integer())
# this will be the pruned + oriented dataset before sd prunning
d  <- data.frame('ss'=character(),'pos'=integer(),
                 'm'=double(),'f'=double(),'lg'=integer())
# while reading prune male and females on edges 
# unpruned goes to unpruned, pruned to d 
for(f in files){
  d1 <- read_delim(f, col_names = c('ss','pos','m','f'))
  d1$pos <- as.integer(d1$pos)
  lg <- strsplit(basename(f),'_')[[1]][9]
  d1$ss <- rename.ss(d1$ss)
  d1$lg <- rep(lg, nrow(d1))
  unpruned <- rbind(unpruned,d1)
  d1$m <- prune.ends(d1$m,100,2)
  d1$f <- prune.ends(d1$f,100,2)
  d <- rbind(d,d1)
}
d <- d[! d$ss %in% shit,]
unpruned <- unpruned[! unpruned$ss %in% shit,]
# make sex averaged 
unpruned$av <- (unpruned$m + unpruned$f) / 2
d$av <- (d$m + d$f)/2
# make merged positions for lg20
newpos <- unpruned[unpruned$lg==20,]
newpos$pos[newpos$ss=='ss49'] <- newpos$pos[newpos$ss=='ss49']
newpos$pos[newpos$ss=='ss3'] <- newpos$pos[newpos$ss == 'ss3']  + 
  max(newpos$pos[newpos$ss=='ss49'])
unpruned[unpruned$lg==20,] <- newpos
newpos <- d[d$lg==20,]
newpos$pos[newpos$ss=='ss49'] <- newpos$pos[newpos$ss=='ss49']
newpos$pos[newpos$ss=='ss3'] <- newpos$pos[newpos$ss == 'ss3']  + 
  max(newpos$pos[newpos$ss=='ss49'])
d[d$lg==20,] <- newpos
# in the loop that follows d is being modified but the sd prunning happens in the 
# temporary dataframe tmp, thus we need to make a new df to save our pruned results 
# thats sdpruned. Care that this code prunes the sd on on sex-averaged and not on male -female specifici
sdpruned <-  data.frame('ss'=character(),'pos'=integer(),
                        'm'=double(),'f'=double(),'lg'=integer(),
                        'av' = double())
# pdf 1 to compare prunning
pdf(paste0('plots/comparing_LM_pruning_',focus,'.pdf'),12,12)
par(mfrow=c(3,3))
for(lg in unique(d$lg)){
  # make sure it starts from 0 
  d$av[d$lg==lg] <- d$av[d$lg==lg] - min(d$av[d$lg==lg],na.rm=T)
  # fix maps so that 0 cm is at the beggining ! 
  lm1 <- lm(d$av[d$lg==lg]~d$pos[d$lg==lg])
  if(lm1$coefficients[2] < 0){
    d$av[d$lg==lg] <-abs( d$av[d$lg==lg] - max(d$av[d$lg==lg], na.rm=T) )
    d$m[d$lg==lg] <-abs( d$m[d$lg==lg] - max(d$m[d$lg==lg], na.rm=T) )
    d$f[d$lg==lg] <-abs( d$f[d$lg==lg] - max(d$f[d$lg==lg], na.rm=T) )
  }
  # Prune based on residuals of order ~ physical 
  tmp <- d[d$lg == lg & ! is.na(d$av), ] #only check av for NA since it will include either or 
  # make physical order
  tmp <- tmp[order(tmp$pos,decreasing = F),]
  tmp$ph <- 1:length(tmp$pos)
  # make genetic order 
  tmp <- tmp[order(tmp$m, decreasing = F),]
  tmp$g <- 1:length(tmp$pos)
  lmor <- lm(tmp$g~tmp$ph)
  # get "bad" points with residuals > half an SD from the mean
  bad <- which(abs(lmor$residuals) > (mean(lmor$residuals) + sd(lmor$residuals)*.5))
  # set to NA
  tmp[bad,c(3,4,6)] <- NA
  # make 'em start at 0 
  tmp$av <- tmp$av - min(tmp$av,na.rm=T)
  tmp$m <- tmp$m - min(tmp$m,na.rm=T)
  tmp$f <- tmp$f - min(tmp$f,na.rm=T)
  sdpruned <- rbind(sdpruned,tmp)
  # PLOT ALL to compare prunning with a certain span  # 
  # unpruned
  plot(unpruned$m[unpruned$lg==lg] ~ unpruned$pos[unpruned$lg==lg],pch=20,main='unpruned males',xlab='',ylab='cM')
  lmm1 <- loess(unpruned$m[unpruned$lg==lg] ~ unpruned$pos[unpruned$lg==lg],span=0.1)
  lines(lmm1$fitted[order(lmm1$x)]~lmm1$x[order(lmm1$x)], col=p[1])
  plot(unpruned$f[unpruned$lg==lg] ~ unpruned$pos[unpruned$lg==lg],pch=20,main='unpruned females',xlab='',ylab='')
  lmm1 <- loess(unpruned$f[unpruned$lg==lg] ~ unpruned$pos[unpruned$lg==lg],span=0.1)
  lines(lmm1$fitted[order(lmm1$x)]~lmm1$x[order(lmm1$x)], col=p[1])
  plot(unpruned$av[unpruned$lg==lg] ~ unpruned$pos[unpruned$lg==lg],pch=20,main='unpruned averaged',xlab='',ylab='')
  lmm1 <- loess(unpruned$av[unpruned$lg==lg] ~ unpruned$pos[unpruned$lg==lg],span=0.1)
  lines(lmm1$fitted[order(lmm1$x)]~lmm1$x[order(lmm1$x)], col=p[1])
  # pruned before sd
  plot(d$m[d$lg==lg] ~ d$pos[d$lg==lg],pch=20,main='edges males',xlab='',ylab='cM')
  lmm1 <- loess(d$m[d$lg==lg] ~ d$pos[d$lg==lg],span=0.1)
  lines(lmm1$fitted[order(lmm1$x)]~lmm1$x[order(lmm1$x)], col=p[1])
  plot(d$f[d$lg==lg] ~ d$pos[d$lg==lg],pch=20,main='edges females',xlab='',ylab='')
  lmm1 <- loess(d$f[d$lg==lg] ~ d$pos[d$lg==lg],span=0.1)
  lines(lmm1$fitted[order(lmm1$x)]~lmm1$x[order(lmm1$x)], col=p[1])
  plot(d$av[d$lg==lg] ~ d$pos[d$lg==lg],pch=20,main='edges averaged',xlab='',ylab='')
  lmm1 <- loess(d$av[d$lg==lg] ~ d$pos[d$lg==lg],span=0.1)
  lines(lmm1$fitted[order(lmm1$x)]~lmm1$x[order(lmm1$x)], col=p[1])
  # pruned after sd
  plot(tmp$m ~ tmp$pos,pch=20,main='sd males',xlab='',ylab='cM')
  lmm1 <- loess(tmp$m ~ tmp$pos,span=0.1)
  lines(lmm1$fitted[order(lmm1$x)]~lmm1$x[order(lmm1$x)], col=p[1])
  plot(tmp$f ~ tmp$pos,pch=20,main='sd females',xlab='',ylab='')
  lmm1 <- loess(tmp$f ~ tmp$pos,span=0.1)
  lines(lmm1$fitted[order(lmm1$x)]~lmm1$x[order(lmm1$x)], col=p[1])
  plot(tmp$av ~ tmp$pos,pch=20,main='sd averaged',xlab='',ylab='')
  lmm1 <- loess(tmp$av ~ tmp$pos,span=0.1)
  lines(lmm1$fitted[order(lmm1$x)]~lmm1$x[order(lmm1$x)], col=p[1])
  mtext(paste('lg',lg,unique(d$ss[d$lg==lg])),line = -1.5,outer = T)
}
dev.off()

write.table(d,paste0('../1.data/1.LepMAP/orders/',focus,'/full_orders_fixed349_nosmall_ascendingAll_mor.txt'),
            quote=F,row.names=F)

