lm <- read.table('../1.data/1.LepMAP/orders/fitted_values_gams.table',h=T)
source('functions.r')
d1m <- read.table('../1.data/3.windows/genome_1mb_windows.bed')
names(d1m) <- c('ss','start','end')
d1m$ss <- rename.ss(d1m$ss)
d1m$ss[d1m$ss == 'ss2' & d1m$end <= 28692481] <- 'ss2.2'
sbped <- window.ped(windows = d1m, ped=lm,sex.based = T)
d1m$m <- sbped$m
d1m$f <- sbped$f
# pdf('plots/Comparing_sexes_windowed.pdf')
# for(ss in unique(lm$ss)){
# plot(d1m$m[d1m$ss==ss]~d1m$start[d1m$ss==ss],type='l',col='navy',lwd=2,
#      main=paste0(ss),ylim=c(0,max(c(d1m$m[d1m$ss == ss],d1m$f[d1m$ss == ss]),na.rm=T)))
# lines(d1m$f[d1m$ss==ss]~d1m$start[d1m$ss==ss],col='magenta',lwd=2)
#}
#dev.off()

d1m <- d1m[!is.na(d1m$m) & ! is.na(d1m$f),]
d1m <- d1m[d1m$m != -Inf & d1m$f != -Inf, ]

plot(d1m$m~d1m$f,log='xy')
hist(d1m$m,101,col=add.alpha('navy',.4)) -> h 
hist(d1m$f,h$breaks,col=add.alpha('red',.4),add=T)
hist(d1m$m - d1m$f,41)
abline(v=5,col='navy')
abline(v=-5,col='red')
d1m$diff <- d1m$m - d1m$f
d1m$dist <- NA
for(ss in unique(d1m$ss)){
  indx <- which(d1m$ss == ss)
  if(length(indx > 1)){
  d1m$dist[indx] <- seq(0,length(indx)-1,1)/(length(indx)-1)
  }
}

d1m$dist[which(d1m$dist > 0.5)] <- abs(d1m$dist[which(d1m$dist > 0.5)] - 1)

mdiff <- d1m$dist[d1m$diff > 2 & ! is.na(d1m$diff)]
fdiff <- d1m$dist[d1m$diff < -2 & ! is.na(d1m$diff)]
boxplot(mdiff,fdiff,pch='')
points(mdiff~jitter(rep(1,length(mdiff))),pch=20,col='blue')
points(fdiff~jitter(rep(2,length(fdiff))),pch=20,col='red')
t.test(mdiff,fdiff)

hist(d1m$dist[d1m$diff > 2 & ! is.na(d1m$diff)],41, 
     col=add.alpha('blue',0.4)) -> h
hist(d1m$dist[d1m$diff < -2 & !is.na(d1m$diff)],h$breaks, add=T,
     col=add.alpha('red4',.4))
tmp1 <- d1m$dist[d1m$diff > 2 & !is.na(d1m$diff)]
tmp2 <- d1m$dist[d1m$diff < -2 & !is.na(d1m$diff)]
t.test(tmp1,tmp2)

# differences are not significantly different in location on chr 
# but usually concentrated on ends which could be due to centro / telo 
# effect  but also due to differential prunning of edge markers !!! 

