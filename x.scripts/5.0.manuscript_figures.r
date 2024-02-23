# this script makes the main text figures 
# CONTENTS 
  # TOOLS NEEDED - lines 12 - 24
  # DATA ACQUISITION AND CURATION - lines 25 - 152
  # FIGURE 1  - lines 154 - 231
  # FIGURE 2  - lines 232 - 437
  # FIGURE 3  - lines 442 - 574
  # FIGURE 4  - lines 575 - 752 


{ # TOOLS NEEDED 
  # only works in Rstudio
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source('functions.r')
library(tidyverse)
library(corrplot)
library(vioplot)
library(plotrix)
library(car)
library(grid)
library(DescTools)
library(recolorize)
library(rworldxtra)
library(tmap)
library(sf)
}
{  # ACQUIRE DATA AND CURATE #
# LOAD PEDIGREE data 
{
  lm <- read.table('../1.data/1.LepMAP/orders/fitted_values_gams.table',h=T)
  co <- read.table('../1.data/1.LepMAP/orders/COs_male_female.table',h=T)    
  co$distunf <- NA
  for(lg in unique(co$lg)){
    co$distunf[co$lg==lg] <-  co$pos[co$lg==lg] / max(co$pos[co$lg==lg])
  }
  co$distunf[co$distunf > .5 ] <- abs(co$distunf[co$distunf > .5 ] - 1)
  lm$ss[lm$lg == 24] <- 'ss2.2'
  lm <- left_join(lm,co,by=c('lg','pos'))
  lms <- lm %>% group_by(lg) %>% summarize(lg=unique(lg),
                                           ss=unique(ss),
                                           m=max(m,na.rm=T),
                                           f=max(f,na.rm=T),
                                           av=max(av,na.rm=T),
                                           len=(max(pos)-min(pos))/1e6,
                                           loci = n(),
                                           CO_male = sum(maleCO,na.rm=T),
                                           CO_female = sum(femaleCO,na.rm=T))
  sum(lms$CO_female,na.rm=T)
  ns <- read.table('../1.data/3.windows/Ns/reference_Ns.bed')
  ns$V1[ns$V1 == 'Super-Scaffold_2' & ns$V3 <= 28692481] <- 'Super-Scaffold_2.2'
  ag <- ns %>% group_by(V1) %>% summarize('N' = sum(V3 - V2))
  lms$Ns <- ag$N[match(lms$ss,rename.ss(ag$V1))]
  lms$lenN <- lms$len - (lms$Ns)/1e6
  rintra <- read.table('../1.data/1.LepMAP/orders/rintra.table',h=T)
  lms <- left_join(lms,rintra[,c(1,3,4,5,6)],by='lg')
}
# load LD data 
d100 <- read.table('../1.data/3.windows/total_100kb.table',h=T)
d10 <- read.table('../1.data/3.windows/total_10kb.table',h=T)
d1 <- read.table('../1.data/3.windows/total_1kb.table',h=T)
d1m <- read.table('../1.data/3.windows/total_1Mb.table',h=T)
d1m <- d1m[d1m$ped  != -Inf,]

# MAKE CM TO CM/MB
for(i in seq(4,8)){
  d100[,i] <- d100[,i] *10
  d10[,i] <- d10[,i] * 100
  d1[,i] <- d1[,i] * 1000
}

# NEED TO FIX CPGIS 
d10$cpgi[d10$cpgi >1e4 & !is.na(d10$cpgi)] <- NA
d1$cpgi[d1$cpgi >1e3 & !is.na(d1$cpgi)] <- NA

# MAKE DISTANCE TO END
d1m$dist <- distance(d1m)
d100$dist <- distance(d100)
d10$dist <- distance(d10)
d1$dist <- distance(d1)

# REMOVE N WINDOWS 
# 100kb
nr <- nrow(d100)
tmp <- read.table('../1.data/3.windows/Ns/100kb_N.overlap')
tmp <- tmp %>% group_by(V1,V2,V3) %>% summarize('N' = sum(V7))
names(tmp) <- c('ss','start','end','N')
d100 <- left_join(d100,tmp,by=c('ss','start','end'))
d100 <- d100[is.na(d100$N),]
print(paste('100kb lost',100*(nr-nrow(d100))/nr,'% rows'))
# 10kb
nr <- nrow(d10)
tmp <- read.table('../1.data/3.windows/Ns/10kb_N.overlap')
tmp <- tmp %>% group_by(V1,V2,V3) %>% summarize('N' = sum(V7))
names(tmp) <- c('ss','start','end','N')
d10 <- left_join(d10,tmp,by=c('ss','start','end'))
d10 <- d10[is.na(d10$N),]
print(paste('10kb lost',100*(nr-nrow(d10))/nr,'% rows'))
# 1kb
nr <- nrow(d1)
tmp <- read.table('../1.data/3.windows/Ns/1kb_N.overlap')
tmp <- tmp %>% group_by(V1,V2,V3) %>% summarize('N' = sum(V7))
names(tmp) <- c('ss','start','end','N')
d1 <- left_join(d1,tmp,by=c('ss','start','end'))
d1 <- d1[is.na(d1$N),]
print(paste('1kb lost',100*(nr-nrow(d1))/nr,'% rows'))
# 1Mb
nr <- nrow(d1m)
tmp <- read.table('../1.data/3.windows/Ns/1mb_N.overlap')
tmp <- tmp %>% group_by(V1,V2,V3) %>% summarize('N' = sum(V7))
names(tmp) <- c('ss','start','end','N')
d1m <- left_join(d1m,tmp,by=c('ss','start','end'))
d1m <- d1m[is.na(d1m$N),]
print(paste('1mb lost',100*(nr-nrow(d1m))/nr,'% rows'))

# REMOVE ROWS WITH MISSING RECOMBINATION DATA 
d100 <- d100[! is.na(d100$CH_cM) & ! is.na(d100$PT_cM) & ! is.na(d100$GB_cM),]
d10 <- d10[! is.na(d10$CH_cM) & ! is.na(d10$PT_cM) & ! is.na(d10$GB_cM),]
d1 <- d1[! is.na(d1$CH_cM) & ! is.na(d1$PT_cM) & ! is.na(d1$GB_cM),]
d1m <- d1m[! is.na(d1m$CH_cM) & ! is.na(d1m$PT_cM) & ! is.na(d1m$GB_cM),]

# REMOVE CH OUTLIERS 
d100 <- d100[ d100$CH_cM < 50,]
d10 <- d10[ d10$CH_cM < 50,]
d1 <- d1[ d1$CH_cM < 50,]
d1m <- d1m[ d1m$CH13_cM < 50,]

# REMOVE Z SCAFFOLDS
z <- c('Super-Scaffold_13','Super-Scaffold_42')
d100z <- d100[ d100$ss %in% z, ]
d100 <- d100[! d100$ss %in% z, ]
d10z <- d10[ d10$ss %in% z, ]
d10 <- d10[! d10$ss %in% z, ]
d1z <- d1[ d1$ss %in% z, ]
d1 <- d1[! d1$ss %in% z, ]
d1mz <- d1m[d1m$ss %in% z,]
d1m <- d1m[! d1m$ss %in% z,]

# HOTSPOTS AS MORE THAN 5 TIMES THE MEAN 
hot1 <- which(d1$CH_hi >= 5)
pthot <- which(d1$PT_hi >= 5)
gbhot <- which(d1$GB_hi >= 5)
ch13hot <- which(d1$CH13_hi >= 5)
ch13_2hot <- which(d1$CH13_2_hi >= 5)

# GLOBAL HOTSPOTS AS 10 TIMES AVERAGE OVERAL RATES 
hot10 <- which(d1$CH_cM /mean(d1$CH_cM,na.rm=T)  >= 10) # CH
pthot10 <- which(d1$PT_cM /mean(d1$PT_cM,na.rm=T)  >= 10) # PT
gbhot10 <- which(d1$GB_cM /mean(d1$GB_cM,na.rm=T)  >= 10) # GB
ch13hot10 <- which(d1$CH13_cM /mean(d1$CH13_cM,na.rm=T)  >= 10) # CH
ch13_2hot10 <- which(d1$CH13_2_cM /mean(d1$CH13_2_cM,na.rm=T)  >= 10) # CH
overlaps <- read.table('overlaps.table')
} # END DATA ACQUISITION AND CURATION


###################################
###################################
###################################
###################################
#### FIGURE 1 #####################
###################################
###################################
###################################
###################################
{
# general 
  plotit = T
axis.size=1.4
mtext.size=1.3
# colours
crp <- viridis::mako(nrow(lms),alpha=1,direction=1)
lms$col[order(lms$lenN,decreasing = T)] <- crp
malecolor=viridis::mako(4)[3]
femalecolor=viridis::mako(4)[2]
# file 
if(plotit){
pdf('plots/v1_figure1.pdf',height=12,width=10)
}
# layout of image 
layout(matrix(c(1,5,
                1,4,
                2,4,
                3,4),byrow = T,nrow=4,ncol=2),
       heights=c(.5,.5,1,1),
       widths=c(1,1))
# PLOT 1 - cM~Mb
{ 
  # regression 
  lm1 <- lm(av~lenN,data=lms)
  # plot 
  plot(y=lms$av,x=lms$lenN, col=lms$col, xlab='', 
       ylab='', main='', axes =F,pch=20, 
       cex=6, xlim=c(0,70),ylim=c(0,100))
  # peripherals 
  axis(1, at = c(0,20,40,60,70), labels=c(0,20,40,60,70),cex.axis=axis.size)
  axis(2, at = c(0,20,41,60,80,100), labels= c(0,20,41,60,80,100),
       cex.axis=axis.size,las=2) 
  mtext('Genetic length (cM)',2,2.5,cex=mtext.size)
  mtext('Physical length (Mb)',1,2.5,cex=mtext.size)
  lnew <- data.frame(lenN=seq(0,80,length.out=100))
  p <- predict(lm1,newdata=lnew,type = 'response',se.fit = T)
  lines(p$fit~lnew$len,lwd=2,col='black')
  lines(y=p$fit-1.96*p$se.fit,x=lnew$len,lty=2,lwd=2,col='black')
  lines(y=p$fit+1.96*p$se.fit,x=lnew$len,lty=2,lwd=2,col='black')
  mtext('A',3,at=0,cex=2)
  # legend 
  crp <- viridis::mako(nrow(lms),alpha=1,direction=-1)
  for(i in seq(1,39,3)){
    points(68,11+i/3,pch=20,cex=3,col=crp[i])
  }
  text(68, 31, "Mbp" , cex = .8)
  text(66, 10, "5" , cex = .8)
  text(66, 26, "61" , cex = .8)
  segments(66, 12, 66,23 )

}
### plot 2 - SEXES 
{
  plot(lms$f~lms$m,xlab='', ylab='',
       col=lms$col, pch=16, cex=4,
       xlim=c(18,82),ylim=c(18,82), axes=F, main ='') 
  axis(1,at=seq(20,80,20),labels=seq(20,80,20),cex.axis=axis.size)
  axis(2,at=seq(20,80,20),labels=seq(20,80,20),cex.axis=axis.size,las=2)
  mtext('Male genetic length (cM)',1,line=2.5,cex=mtext.size)
  mtext('Female genetic length (cM)',2,line=2.5,cex=mtext.size)
  abline(0,1,lty=2,lwd=2)
  mtext('B',3,at=16,cex=2)
  # legend 
  for(i in seq(1,39,3)){
    points(81,21+i/4,pch=20,cex=3,col=crp[i])
  }
  text(81, 35, "Mbp" , cex = .8)
  text(79, 20, "5" , cex = .8)
  text(79, 32, "61" , cex = .8)
  segments(79, 22.5, 79, 28.5 )
}
### plot 3 - CO location density 
{
  maley <- density(co$distunf[co$maleCO!=0],from=0,to=.5)$y
  malex <- density(co$distunf[co$maleCO!=0],from=0,to=.5)$x
  plot(x=malex,y=maley/sum(maley),col=malecolor,lwd=4,type='l',axes=F,xlab='',ylab='')
  axis(1,at=seq(0,0.5,0.1),labels = seq(0,50,10),cex.axis=1.2)
  axis(2,at=seq(0,.004,.001),labels=seq(0,4,1),cex.axis=1.2,las=2)
  femaley <- density(co$distunf[co$femaleCO!=0],from=0,to=.5)$y
  femalex <- density(co$distunf[co$femaleCO!=0],from=0,to=.5)$x
  lines(x=femalex,y=femaley/sum(femaley),col=femalecolor,lwd=4,lty=1)
  mtext('Density (x 1e-3)',2,2.5,cex=mtext.size)
  mtext('Distance from LG end (% of LG length)',1,2.5,cex=mtext.size)
  mtext('C',3,at=-.02,cex=2)
  legend(x=0.32,y=.004,lty=c(1,1),col=c(malecolor,femalecolor),lwd=4,
         legend=c('Male','Female'),cex=1.5,bty='n')
  # test significance 
  ks.test(maley,femaley)
}
### plot 4 - Barplot of r intra
{
  bp <- as.matrix(lms[order(lms$lenN,decreasing=T),c(14,15)])
  lms$difr <- rep(NA,nrow(bp))
  lms$diff <- rep(NA,nrow(bp))
  for(i in 1:nrow(lms)){
    lms$difr[i] <- lms$rintra_M[i]/lms$rintra_F[i]
    lms$diff[i] <- ( lms$m[i] - lms$f[i] ) / lms$m[i]
  }
  par(mar=c(5,4,4,4))
  barplot(lms$difr[order(lms$len,decreasing=F)] -1,
          xlab='', width=2,
          col=c(malecolor,femalecolor)[
            ifelse(lms$difr[order(lms$len,decreasing=F)] > 1,1,2)],
          axes=F, space=1, xlim=c(-.6,.6), ylim=c(0,150),
          horiz=T, border=NA)
  axis(1, at = c(-.5,0,.5),labels=c(0.5,1,1.5),cex.axis=1.3,las=1)
  mtext('Intrachromosomal shuffling (Male / Female)',1,2.5,cex=1)
  axis(2, at=seq(3,155,4), lms$lg[order(lms$lenN,decreasing=F)],
       cex.axis=1.1, las=2, adj=1)
  mtext('Linkage Group',2,2.5,cex=mtext.size)
  abline(v=0,lwd=2)
  legend('topleft', pch=c(15,15), pt.cex=c(2,2), col=c(malecolor,femalecolor),
         legend=c('Male > Female','Female > Male'), cex=1.5, bty='n')
  mtext('D',3,at=-.6,cex=2)
  axis(4,at=seq(3,155,4), labels=round(lms$len[order(lms$len)],0),
       cex.axis=1,las=2,adj=1)
  mtext('Physical Length (Mb)',4,2.5,cex=mtext.size)
}
if(plotit){
dev.off()
}
}


###################################
###################################
###################################
###################################
#### FIGURE 2 #####################
###################################
###################################
###################################
###################################
{
  plotit=T
  # global figure options 
  label.size = 1.5
  axis.size = 1.5
  # colours
  v1 <- viridis::mako(5,alpha = 1)
  v1a <- viridis::mako(5,alpha = .6)
  v2a <- viridis::mako(5,alpha = .2)
  # layout 
  if(plotit){
    pdf('plots/v1_figure2.pdf',width=12,height=14)
  }
  layout( matrix(c(8,1,1,9,10,
                   2,3,4,4,4,
                   5,5,6,6,7),
                 nrow=3, ncol=5, byrow = T), 
          widths = c(10,15,5,10,15),  
          heights = c(1,1,1) )
  # three scaffold examples 
  ss.ex <- c('Super-Scaffold_30',
             'Super-Scaffold_39',
             'Super-Scaffold_8')
  lg.ex <- c(36,32,22)
  {
    # PLOT 1 - CORRELATION
    # remove bad pedigree window 
    d1m <- d1m[d1m$ped != -Inf,]
    # make ln (x+1)
    d1m$x <- log1p(d1m$CH_cM)
    d1m$y <- log1p(d1m$ped)
    plot(log1p(d1m$ped)~log1p(d1m$CH_cM),pch=20,col=v2a[2],xlim=c(0,3),ylim=c(0,3),
         axes=F,xlab='',ylab='',cex=3)
    lm2 <- lm(y~x, data=d1m)
    p2 <- predict(lm2, newdata=data.frame('x'=seq(0,3,length.out=1000)),type='response',se.fit = T)
    lines(y=p2$fit,x=seq(0,3,length.out=1000),lwd=2,col=v1[1])
    lines(y=p2$fit-1.96*p2$se.fit,x=seq(0,3,length.out=1000),lty=2,lwd=2,col=v1[1])
    lines(y=p2$fit+1.96*p2$se.fit,x=seq(0,3,length.out=1000),lty=2,lwd=2,col=v1[1])
    axis(1,at=seq(0,3,1),labels=round(expm1(seq(0,3,1)),1),cex.axis=1.2)
    axis(2,at=seq(0,3,1),labels=round(expm1(seq(0,3,1)),1),las=2,cex.axis=1.2)
    mtext('pyrho (log(x+1))',1,2.4,cex=1.2)
    mtext('LepMap3 (log(x+1))',2,2.8,cex=1.2)
    mtext('A',3,at=0,cex=2)
  }
  {
    # PLOT 2 
    ## PLOT SCAFFOLDS
    plot(d10$CH_cM[d10$ss=='Super-Scaffold_30']~
           d10$start[d10$ss=='Super-Scaffold_30'],type='l',
         axes=F,xlim=c(0,8e6),ylim=c(0,25),
         xlab='',ylab='',main='',col=v1[2],lwd=2)
    abline(h=mean(d10$CH_cM,na.rm=T),lty=2,col='black')
    axis(1,at=c(0,8e6),labels=c(0,8),cex.axis=axis.size)
    mtext('Mb',1,2.5,cex = label.size)
    points(max(d10$start[d10$ss=='Super-Scaffold_30'],na.rm=T)/2,22,pch=15,cex=4,col=v1[2])
    axis(2,at=seq(0,25,10),labels=seq(0,25,10),cex.axis=axis.size,las=2)
    mtext('cM/Mb',2,2.5,cex = label.size)
    mtext('LG 36',3,0,cex = label.size)
    mtext('B',3,at=1,cex=2)
    #
    plot(d10$CH_cM[d10$ss=='Super-Scaffold_39']~
           d10$start[d10$ss=='Super-Scaffold_39'],type='l',
         axes=F,xlim=c(0,15e6),ylim=c(0,25),
         xlab='',ylab='',main='',col=v1[3],lwd=2)
    axis(1,at=c(0,8e6,15e6),labels=c(0,8,15),cex.axis=axis.size)
    abline(h=mean(d10$CH_cM,na.rm=T),lty=2,col='black')
    mtext('Mb',1,2.5,cex = label.size)
    mtext('LG 32',3,0,cex = label.size)
    points(max(d10$start[d10$ss=='Super-Scaffold_39'],na.rm=T)/2,22,pch=18,cex=4,col=v1[3])
    axis(2,at=seq(0,25,10),labels=seq(0,25,10),cex.axis=axis.size,las=2)
    plot(d10$CH_cM[d10$ss=='Super-Scaffold_8']~
           d10$start[d10$ss=='Super-Scaffold_8'],type='l',
         axes=F,xlim=c(0,25e6),ylim=c(0,25),
         xlab='',ylab='', main='',col=v1[4],lwd=2)
    points(max(d10$start[d10$ss=='Super-Scaffold_8'],na.rm=T)/2,22,pch=17,cex=4,col=v1[4])
    abline(h=mean(d10$CH_cM,na.rm=T),lty=2,col='black')
    axis(1,at=c(0,8e6,15e6,25e6),labels=c(0,8,15,25),
         cex.axis=axis.size)
    mtext('Mb',1,2.5,cex = label.size)
    mtext('LG 22',3,0,cex = label.size)
    axis(2,at=seq(0,25,10),labels=seq(0,25,10),cex.axis=axis.size,las=2)
  }
  {
    # PLOT 3 - HOTSPOT PLOT 
    # use only linkage map scaffolds 
    d10 <- d10[d10$ss %in% gsub('ss','Super-Scaffold_',lm$ss),]
    # dataset with propseq and late gini 
    dfd <- data.frame(ss=unique(d10$ss),propsec=NA)
    dfd$gini <- NA 
    for(ss in dfd$ss){
      dfd$gini[dfd$ss==ss] <- Gini(d10$CH_cM[d10$ss==ss])
    }
    # plot 
    plot(0,pch='',xlab='', xlim=c(0,1),
         ylab= '', ylim=c(0,1), axes=F)
    # add lines per scaffold 
    for(ss in unique(d10$ss)){
      tmp <- d10[d10$ss == ss,] # subset
      tmp$CH_cM <- tmp$CH_cM / 100 
      # make x and y for lines based on ordering
      x1 <- cumsum(tmp$CH_cM[order(tmp$CH_cM,decreasing = T)])/sum(tmp$CH_cM,na.rm=T)
      y1 <- seq(0,1,length.out=nrow(tmp))
      # set colour for example scaffolds 
      lc1 <- ifelse(ss %in% ss.ex, 
                    v1[which(ss.ex == ss) + 1] ,
                    add.alpha('black',.3))
      lines(x= x1, y=y1, lwd=3, col=lc1)
      dfd$propsec[dfd$ss==ss] <- y1[which(abs(x1-0.8) == min(abs(x1-0.8)))]
    }
    # make overall lines 
    x1 <- cumsum(d10$CH_cM[order(d10$CH_cM,decreasing = T)])/sum(d10$CH_cM,na.rm=T)
    y1 <- seq(0,1,length.out=nrow(d10))
    lines(x=x1, y=y1, col='grey',lwd=4,lty=2)
    # add example points - maybe remove 
    points(x=.65,y=.40,pch=15,cex=3,col=v1[2])
    points(x=.75,y=.32,pch=18,cex=4,col=v1[3])
    points(x=.85,y=.28,pch=17,cex=3,col=v1[4])
    # peripherals 
    abline(0,1,lty=3)
    axis(1,at=seq(0,1,.2),labels=seq(0,1,.20),cex.axis=axis.size)
    mtext('cumulative recombination',1,2.6,cex = label.size)
    axis(2,at=seq(0,1,.2),labels=seq(0,1,.20),cex.axis=axis.size,las=2)
    mtext('cumulative sequence',2,2.6,cex = label.size)
    mtext('C',3,at=0,cex=2)
  }
  {
    # PLOT 4 & 5 GINIS 
    # make dataset 
    lm <- read.table('../1.data/1.LepMAP/orders/fitted_values_gams.table',h=T)
    lm$ss[lm$lg == 24] <- 'ss2.2'
    d10.tmp <- d10[d10$ss %in% gsub('ss','Super-Scaffold_',unique(lm$ss)),]
    ss10 <- d10.tmp %>% group_by(ss) %>% summarize(cM = sum(CH_cM) / 100, 
                                                   pi = mean(CH_pi,na.rm=T),
                                                   gc = mean(gc,na.rm=T),
                                                   cpgi = sum(cpgi,na.rm=T),
                                                   len = max(end)-min(start),
                                                   ts = sum(tes) + sum(tss),
                                                   se.pi = plotrix::std.error(CH_pi,na.rm=T),
                                                   loci = sum(CH_snps,na.rm=T))
    ss10 <- left_join(ss10,dfd,by='ss')
    # set pch, cex and col for examples 
    pch1 <- rep(16,nrow(ss10))
    cex1 <- rep(4,nrow(ss10))
    col1 <- rep(add.alpha('black',.4),nrow(ss10))
    pch1[ss10$ss == 'Super-Scaffold_30'] <- 15
    cex1[ss10$ss == 'Super-Scaffold_30'] <- 4
    col1[ss10$ss == 'Super-Scaffold_30'] <- v1[2]
    pch1[ss10$ss == 'Super-Scaffold_8'] <- 17
    cex1[ss10$ss == 'Super-Scaffold_8'] <- 4
    col1[ss10$ss == 'Super-Scaffold_8'] <- v1[4]
    pch1[ss10$ss == 'Super-Scaffold_39'] <- 18
    cex1[ss10$ss == 'Super-Scaffold_39'] <- 5
    col1[ss10$ss == 'Super-Scaffold_39'] <- v1[3]
    plot(y=ss10$gini,x=ss10$len,pch=pch1,cex=cex1,xlab='',
         ylab='',ylim=c(.3,.75),xlim=c(0,65e6),axes=F,col=col1)
    axis(2,at=seq(0.3,.7,.1),labels = seq(0.3,.7,.1),cex.axis=axis.size,las=2)
    abline(h=0.6,lty=2)
    mtext('Gini Coefficient',2,2.6,cex = label.size)
    axis(1,at=c(0,25e6,50e6,65e6),labels=c(0,25,50,65),cex.axis=axis.size)
    mtext('Length (Mb)',1,2.6,cex = label.size)
    mtext('D',3,cex=2,at=.39)
    plot(y=ss10$gini,x=ss10$pi,pch=pch1,cex=cex1,axes=F,xlab='',
         ylab='',ylim=c(.3,.75),xlim=c(0.0013,0.0027),col=col1)
    axis(2,at=seq(0.3,.7,.1),labels = seq(0.3,.7,.1),cex.axis=axis.size,las=2)
    abline(h=0.6,lty=2)
    axis(1,at=c(0.0015,0.002,0.0025,0.0027),labels=c(1.5,2,2.5,2.7),cex.axis=axis.size)
    mtext('pi (x1000)',1,2.6,cex = label.size)
    mtext('E',3,at=1.5e-3,cex=2)
  }
  if(plotit){
  dev.off()
  }
}


###################################
###################################
###################################
###################################
#### FIGURE 3 #####################
###################################
###################################
###################################
###################################

{
  # PREFIGURE DATA
# matrix of TES and TSS 
d1$ts <- d1$tes + d1$tss
ts <- which(d1$ts > 0)
ts <- ts[ts > 41]
tsmat <- matrix(data=NA,nrow=length(ts),ncol=81)
for(i in 1:length(ts)){
  r1 <- (ts[i] - 40) : (ts[i] + 40)
  tsmat[i,] <- d1$CH_cM[r1] / mean(d1$CH_cM[r1],na.rm=T)
}
# matrix of CPGIs
cp <- which(d1$cpgi == 1000)
cp <- cp[cp > 41]
cpmat <- matrix(data=NA,nrow=length(cp),ncol=81)
for(i in 1:length(cp)){
  r1 <- (cp[i] - 40) : (cp[i] + 40)
  cpmat[i,] <- d1$CH_cM[r1] / mean(d1$CH_cM[r1],na.rm=T)
}
# matrix of hotspots for example plot 
hot1 <- hot1[hot1 > 41]
hotmat <- matrix(data=NA,nrow=length(hot1),ncol=81)
for(i in 1:length(hot1)){
  r1 <- (hot1[i]-40):(hot1[i]+40)
  hotmat[i,] <- d1$CH_cM[r1]/mean(d1$CH_cM[r1],na.rm=T)
}
}
{
  #general 
  plotit = T
  if(plotit){
    pdf('plots/v1_figure3.pdf',width=12,height=12)
  }
  #jpeg('plots/v1_figure3.jpg',width=1440,height=1200,res=150)
  #
  layout(matrix(c(1,4,4,2,5,5,3,6,6),byrow = T,nrow=3,ncol=3),
         heights = c(1,1,1),widths=c(1,1,1))
  axis.size = 1.2
  mtext.size = 1 
  viocol = viridis::mako(4)[c(4,3,2)]
  # TSS+TES / CPGIS
  elcol <- viridis::mako(4)[c(2,3)]
  plot(colMeans(tsmat,na.rm=T),type='l',axes=F,xlab='',ylab='',lwd=4,
       xlim=c(1,81),ylim=c(.9,1.6),main='',
       col=elcol[1],lty=1)
  lines(colMeans(cpmat,na.rm=T),col=elcol[2],lwd=4,lty=2)
  axis(1,at=c(1,21,41,61,81),labels=c(-40,-20,0,20,40),cex.axis=axis.size)
  axis(2,at=c(.9,1.2,1.6),labels=c(.9,1.2,1.6),cex.axis=axis.size,las=2)
  mtext('Distance from element (kb)',1,2.5,cex=mtext.size)
  mtext('RRR80',side=2,line=2.5,cex=mtext.size)
  legend('topleft',lty=c(1,3),col=elcol,legend=c('TSS/TES','CGIs'),lwd=4,bty='n',
         cex=1.2)
  abline(v=41,lwd=1,lty=1)
  mtext('A',3,at=1,cex=2)
  # LOCAL HOTSPOTS
  plot(colMeans(hotmat,na.rm=T),axes=F,xlab='',ylab='',xlim=c(1,80),
       ylim=c(0,8),col=viocol[1],pch=16,cex=2,type='b',main='local hotspots')
  points(x=41,y=colMeans(hotmat,na.rm=T)[41],pch=16,col=viocol[2],cex=2)
  axis(1,at=c(1,21,41,61,81),labels=c(-40,-20,0,20,40),cex.axis=axis.size)
  mtext('Position (kb)',1,2.5,cex=mtext.size)
  axis(2,at=seq(0,8,2),labels = seq(0,8,2),las=2,cex=axis.size)
  mtext('RRR80',2, 2.5,cex=mtext.size)
  abline(h=5)
  abline(v=41,lty=2,lwd=2)
  mtext('B',3,at=1,cex=2)
  # GLOBAL HOTSPOTS
  hist(log(d1$CH_cM/mean(d1$CH_cM)),101,xlab='',ylab='',xlim=c(-5,5),
       axes=F,main='global hotspots',col=viocol[1],border = NA) -> h
  hist(log(d1$CH_cM[hot10]/mean(d1$CH_cM)),h$breaks,add=T,col=viocol[3],
       border=NA)
  abline(v=log(10),col=viocol[3],lwd=2,lty=2)
  axis(1,at=seq(-5,5,5),labels = seq(-5,5,5),cex.axis=axis.size)
  mtext('RRR (log)',1,2.5,cex=mtext.size)
  axis(2,at=seq(0,3e4,1e4),labels=seq(0,30,10),cex.axis=axis.size,las=2)
  mtext('Counts (x1000)',2,2.5,cex=mtext.size)
  mtext('C',3,at=-5,cex=2)
  # Example chr 
  tmp1 <- d1[hot1,]
  tmp2 <- d1[hot10,]
  plot(y=d1$CH_cM[d1$ss=='Super-Scaffold_39']/mean(d1$CH_cM),
       x= d1$start[d1$ss=='Super-Scaffold_39'],pch=20,
       col=add.alpha('grey70',.5),axes=F,xlab='',ylab='',
       ylim=c(0,18),cex=.8,xlim=c(0,15e6),main='')
  axis(1,at=seq(0,15e6,5e6),labels=seq(0,15,5),cex.axis=axis.size)
  mtext('Position (Mb)',1,2.5,cex=mtext.size)
  mtext('LG 32',3,0,cex=mtext.size)
  axis(2,at=c(0,5,10,15),labels=c(0,5,10,15),las=2,cex.axis=axis.size)
  mtext('RRR',2,2.5,cex=mtext.size)
  points(y=tmp1$CH_cM[tmp1$ss=='Super-Scaffold_39']/mean(d1$CH_cM),
         x=tmp1$start[tmp1$ss=='Super-Scaffold_39'], col=viocol[2],pch=20,
         cex=1.5)
  points(y=tmp2$CH_cM[tmp2$ss=='Super-Scaffold_39']/mean(d1$CH_cM),
         x=tmp2$start[tmp2$ss=='Super-Scaffold_39'], col=add.alpha(viocol[3],.6),
         cex=1.3)
  abline(h=10,col=viocol[3],lwd=2)
  mtext('D',3,at=0,cex=2)
  
  # VIOLIN PLOTS
  vioplot::vioplot(d1$dist,d1$dist[hot1],d1$dist[hot10],col=viocol,
                   ylim=c(0,.55), axes=FALSE, xaxt='n',yaxt='n',
                   main='')
  axis(1, at = c(1,2,3), labels=c('all','local','global'),cex.axis=axis.size)
  axis(2,at = seq(0,0.5,0.1),labels=seq(0,50,10),cex.axis=axis.size,las=2)
  mtext('***',line = -.8,at=2,cex=mtext.size)
  mtext('Distance from LG end (%)',2,2.5,cex=1.1)
  wilcox.test(d1$dist,d1$dist[hot1])
  #wilcox.test(d1$dist,d1$dist[hot5])
  mtext('E',3,at=.4,cex=2)
  ###
  vioplot::vioplot(d1$gc*100,100*d1$gc[hot1],100*d1$gc[hot10],col=viocol,main='',
                   ylim=c(20,88), axes=FALSE, xaxt='n',yaxt='n')
  axis(1, at = c(1,2,3), labels=c('all', 'local','global'),
       cex.axis=axis.size)
  axis(2,at = seq(20,80,20),labels=seq(20,80,20),cex.axis=axis.size,las=2)
  mtext('GC %',2,2.5,cex=mtext.size)
  mtext('windows',1,2.5,cex=mtext.size)
  wilcox.test(d1$gc,d1$gc[hot1])
  #wilcox.test(d1$gc,d1$gc[hot5])
  #wilcox.test(d1$gc[hot1],d1$gc[hot5])
  mtext('F',3,at=.4,cex=2)
  ######
  if(plotit){
    dev.off()
  }
}


###################################
###################################
###################################
###################################
#### FIGURE 4 #####################
###################################
###################################
###################################
###################################
{# MAP save to pdf and merge in illustrator
  metadata <- read.table('../2.metadata/metadata_Bioproject_v4.tsv',
                         h=T,sep='\t')
  popcol = viridis::mako(6,direction=1)[c(4,2,3,1,5)]
  metadata <- metadata[metadata$REF_Panel==1,]
  CH <- which(metadata$Population == 'CH')
  PT <- which(metadata$Population %in% c('PT','MA'))
  GB <- which(metadata$Population == 'GB')
  CH.d <- cbind.data.frame(metadata$lon[CH], metadata$lat[CH])
  CH.spdf <- SpatialPointsDataFrame(CH.d, metadata[CH,c(2,4,6,7,8,12,13,14,15)])
  CHsf <- st_as_sf(CH.spdf, crs = st_crs(3035) )
  PT.d <- cbind.data.frame(metadata$lon[PT], metadata$lat[PT])
  PT.spdf <- SpatialPointsDataFrame(PT.d, metadata[PT,c(2,4,6,7,8,12,13,14,15)])
  PTsf <- st_as_sf(PT.spdf, crs = st_crs(3035) )
  GB.d <- cbind.data.frame(metadata$lon[GB], metadata$lat[GB])
  GB.spdf <- SpatialPointsDataFrame(GB.d, metadata[GB,c(2,4,6,7,8,12,13,14,15)])
  GBsf <- st_as_sf(GB.spdf, crs = st_crs(3035) )
  data(countriesHigh)
  world <- countriesHigh %>%  st_as_sf
  world<- st_transform(world, crs = st_crs(3035))
  
  europe_bbox <- st_bbox(c(xmin = -10, xmax = 30, ymax =65, ymin = 30), crs = st_crs(4326))
  europe_bbox = st_bbox(
    st_transform(
      st_as_sfc(europe_bbox), 
      3035
    )
  )
  europe_map_cropped <-world %>%  st_as_sf %>%  sf::st_crop(europe_bbox)
  t <- tm_shape(europe_map_cropped)  + 
    tm_polygons() + 
    tm_graticules() + 
    tm_layout(title = "",
              frame=F,
              bg.color = "white",
              space.color = "white",
              earth.boundary.lwd=1.5)+ 
    tm_shape(CHsf)  + 
    tm_dots(col=popcol[3],shape=20,size=1.5,alpha=1) + 
    tm_shape(PTsf)  + 
    tm_dots(col=popcol[2],shape=20,size=1.5,alpha=1) + 
    tm_shape(GBsf)  + 
    tm_dots(col=popcol[1],shape=20,size=1.5,alpha=1) 
  
  pdf('plots/tmpmap.pdf',10,10)
  print(t)
  dev.off()
}
# Fig 4 part 2 
# make pdf and merge with map in illustrator
{
  # global parameters 
  mtext.size=1.5
  axis.size=1.2
  popcol = viridis::mako(6,direction=1)[c(4,2,3,1,5)]
  plotit = T
  if(plotit){
    pdf('plots/v1_figure4.pdf',12,12)
  }
  layout(matrix(c(1,3,2,4),nrow=2,ncol=2),heights = c(.6,.4),widths = c(.5,.5))
  
  # 1 empty to fill with map 
  plot(0,pch='',xlab='',axes=F,ylab='')
  mtext('A',3,cex=2,adj=-.1)
  
  ## 2 Example LG32 
  {
    par(mar=c(4,3,2,1)+0.1)
    ## GB
    plot(d1$GB_cM[d1$ss=='Super-Scaffold_39'& d1$end < 5e6]~
           d1$end[d1$ss=='Super-Scaffold_39'& d1$end < 5e6],
         type='l',lwd=2,ylim=c(0,200),col=popcol[1],
         xlab='',ylab='',axes=F,pch=20,xlim=c(0,5.4e6),main='LG 32')
    segments(x=0,x1=5e6,y0=10*mean(d1$CH13_cM),y1=10*mean(d1$CH13_cM),lty=2)
    tmp3 <- d1[gbhot,]
    points(y=rep(30,length(tmp3$end[tmp3$ss=='Super-Scaffold_39'& tmp3$end < 5e6])),
           x=tmp3$end[tmp3$ss=='Super-Scaffold_39'& tmp3$end < 5e6], col=popcol[1],
           pch=16,cex=1)
    ## PT
    lines(40+d1$PT_cM[d1$ss=='Super-Scaffold_39'& d1$end < 5e6]~
            d1$end[d1$ss=='Super-Scaffold_39'& d1$end < 5e6],
          type='l',lwd=2,col=popcol[2])
    segments(x=0,x1=5e6,y0=40+10*mean(d1$CH13_cM),y1=40+10*mean(d1$CH13_cM),lty=2)
    tmp2 <- d1[pthot,]
    points(y=rep(70,length(tmp2$end[tmp2$ss=='Super-Scaffold_39'& tmp2$end < 5e6])),
           x=tmp2$end[tmp2$ss=='Super-Scaffold_39'& tmp2$end < 5e6], col=popcol[2],
           pch=16, cex=1)
    # CH
    lines(80+d1$CH_cM[d1$ss=='Super-Scaffold_39' & d1$end < 5e6 ]~
            d1$end[d1$ss=='Super-Scaffold_39'& d1$end < 5e6],
          type='l',lwd=2,col=popcol[3])
    segments(x=0,x1=5e6,y0=80+10*mean(d1$CH13_cM),y1=80+10*mean(d1$CH13_cM),lty=2)
    tmp1 <- d1[hot1,]
    points(y=rep(110,length(tmp1$end[tmp1$ss=='Super-Scaffold_39'& tmp1$end < 5e6])),
           x=tmp1$end[tmp1$ss=='Super-Scaffold_39'& tmp1$end < 5e6], col=popcol[3],
           pch=16, cex=1)
    ## CH13 
    lines(120+d1$CH13_cM[d1$ss=='Super-Scaffold_39'& d1$end < 5e6]~
            d1$end[d1$ss=='Super-Scaffold_39'& d1$end < 5e6],
          type='l',lwd=2,col=popcol[4],lty=1)
    segments(x=0,x1=5e6,y0=120+10*mean(d1$CH13_cM),y1=120+10*mean(d1$CH13_cM),lty=2)
    tmp3 <- d1[ch13hot,]
    points(y=rep(150,length(tmp3$end[tmp3$ss=='Super-Scaffold_39'& tmp3$end < 5e6])),
           x=tmp3$end[tmp3$ss=='Super-Scaffold_39'& tmp3$end < 5e6], col=popcol[4],
           pch=16,cex=1)
    ## CH13_2
    lines(160+d1$CH13_2_cM[d1$ss=='Super-Scaffold_39'& d1$end < 5e6]~
            d1$end[d1$ss=='Super-Scaffold_39'& d1$end < 5e6],
          type='l',lwd=2,col=popcol[5],lty=1)
    segments(x=0,x1=5e6,y0=160+10*mean(d1$CH13_cM),y1=160+10*mean(d1$CH13_cM),lty=2)
    tmp3 <- d1[ch13_2hot,]
    points(y=rep(190,length(tmp3$end[tmp3$ss=='Super-Scaffold_39'& tmp3$end < 5e6])),
           x=tmp3$end[tmp3$ss=='Super-Scaffold_39'& tmp3$end < 5e6], col=popcol[5],
           pch=16,cex=1)
    # for all fig2 
    axis(1,at=seq(0,5e6,1e6),labels=seq(0,5,1),cex.axis=axis.size)
    mtext('Mb',1,2.2,cex=mtext.size)
    mtext('cM/Mb',2,2.2,cex=mtext.size)
    axis(2,at=c(0,20),labels=c(0,20),cex.axis=axis.size,las=2)
    axis(2,at=c(40,60),labels=c(0,20),cex.axis=axis.size,las=2)
    axis(2,at=c(80,100),labels=c(0,20),cex.axis=axis.size,las=2)
    axis(2,at=c(120,140),labels=c(0,20),cex.axis=axis.size,las=2)
    axis(2,at=c(160,180),labels=c(0,20),cex.axis=axis.size,las=2)
    text('GB',x=5.1e6,y=10)
    text('PT',x=5.1e6,y=50)
    text('CH',x=5.1e6,y=90)
    text('CH13',x=5.1e6,y=130,cex=.8)
    text('CH13_2',x=5.1e6,y=170,cex=.8)
    mtext('B',3,-2.2,cex=2,adj=-.1)
  }
  
  ## 3 Corrplot at 2 scales 
  {
    call <- cor(d1[,c(8,7,4,5,6)],use='complete.obs')
    tmp <- cor(d100[,c(8,7,4,5,6)],use='complete.obs')
    call[lower.tri(call)] <- tmp[lower.tri(tmp)]
    row.names(call) <- c('CH13_2', 'CH13','CH', 'PT', 'GB')
    colnames(call) <- c('CH13_2', 'CH13','CH', 'PT', 'GB')
    #layout(matrix(c(4,1,2,3),nrow=2,ncol=2,byrow = F))
    ax.cex = 1.2
    lab.cex = 1.2
    par(mar=c(1,1,1,1))
    tmpx <- viridis::mako(49,direction = -1)
    corrplot(call,method='color',is.corr=F,
             col = tmpx,#COL1(sequential = 'YlGn'),
             col.lim = c(.6,1),
             type='full', number.digits=2,
             addCoef.col = 'white',cl.pos = 'b',
             diag=F, tl.cex = 1.7, number.cex = 1.5,
             tl.col=popcol[c(5,4,3,2,1)], tl.srt=45,na.label = ' ', tl.pos = 'd')
    #mtext('Correlation',4,-2.2,cex=mtext.size)
    mtext('C',3,cex=2,adj=0.05)
  }
  
  ## 4. Venn diagram
  {
    venn <- VennDiagram::venn.diagram(list(chhot2,pthot2,gbhot2),
                                      category.names=c('CH','PT','GB'),
                                      filename=NULL,
                                      disable.logging = T,
                                      label.col = 'white',
                                      cex=1.5,
                                      main='Local Hotspot overlap',
                                      fill=popcol[c(3,2,1)],
                                      col=rep('black',3),
                                      lwd=rep(2,3))
    vp <- viewport(x=.52,y=.4,width=.4,height=.4,just=c('left','top'))
    pushViewport(vp)
    grid.draw(venn)
    grid.draw(venn)
    popViewport()
    plot(0,pch='',xlab='',axes=F,ylab='')
    mtext('D',3,cex=2,adj=0.05)
    par(mar=c(5.1, 4.1, 4.1, 2.1))
    
  }
  if(plotit){
    dev.off()
  }
}


##### ENO OF MANUSCRIPT FIGURES ######


