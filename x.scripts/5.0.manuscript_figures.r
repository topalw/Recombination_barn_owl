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
  
# ACQUIRE DATA AND CURATE #
# LOAD PEDIGREE data 
{
  lm <- read_delim('../1.data/1.LepMAP/R.AllSNPs_physical/physical_orders_full_COs.txt')
  lm$dist <- distance(lm)
  # read lengths (removing Ns from Bionano)
  lengths <- read.table('../1.data/lengths.txt',h=T)
  lengths <- rbind.data.frame(lengths,data.frame('ss'='Super-Scaffold_2.2',
                                      'len'=len2.2,
                                      'lg'=24))
  # read info / LG
  lms <- read.table('../1.data/1.LepMAP/R.AllSNPs_physical/linkage_map_perLG.txt',h=T)
  lms$ss[lms$lg==24] <- 'Super-Scaffold_2.2'
  lms$len[lms$lg==24] <- len2.2
  # read centromeres 
  centro <- read.csv('../3.Centromeres/centromeres.csv',h=T)
  lms$centro <- centro$centro[match(lms$lg,centro$lg)]
  # add crossovers / LG
  tmp <- aggregate(lm$maleCOs,by=list(lm$lg),FUN=sum,na.rm=T)
  lms$maleCOs <- tmp$x[match(lms$lg,tmp$Group.1)]
  tmp <- aggregate(lm$femaleCOs,by=list(lm$lg),FUN=sum,na.rm=T)
  lms$femaleCOs <- tmp$x[match(lms$lg,tmp$Group.1)]
  # make dataset with unique lg rows 
  lms2 <- lms %>% group_by(lg) %>% summarize(m=max(m),
                                             f=max(f),
                                             av=max(av),
                                             len=sum(len)/1e6,
                                             loci=sum(loci),
                                             l2 = sum(l2),
                                             rintra_m=sum(rintra_m),
                                             rintra_f=sum(rintra_f),
                                             mCOs=unique(maleCOs),
                                             fCOs=unique(femaleCOs),
                                             centro=unique(centro))
  lms2$centro_start <- ifelse(lms2$centro - lms2$len /2  > 0, lms2$len*1e6, 1)
}

# load LD data 
{
d100 <- read.table('../1.data/3.windows/total_100kb.table',h=T)
d10 <- read.table('../1.data/3.windows/total_10kb.table',h=T)
d1 <- read.table('../1.data/3.windows/total_1kb.table',h=T)
d1m <- read.table('../1.data/3.windows/total_1Mb.table',h=T)
# d1m <- d1m[d1m$ped  != -Inf,]

# MAKE CM TO CM/MB
for(i in seq(4,11)){
  d100[,i] <- d100[,i] *10
  d10[,i] <- d10[,i] * 100
}
for(i in c(5,7,9,11,13,15,17,19)){
  d1[,i] <- d1[,i] * 1000
}

# NEED TO FIX CPGIS 
d10$cpgi[d10$cpgi >1e4 & !is.na(d10$cpgi)] <- 1e4
d1$cpgi[d1$cpgi >1e3 & !is.na(d1$cpgi)] <- 1e3

# MAKE DISTANCE TO END
d1m$dist <- distance(d1m)
d100$dist <- distance(d100)
d10$dist <- distance(d10)
d1$dist <- distance(d1)

# REMOVE N WINDOWS with > Lwin/10
# 100kb
nr <- nrow(d100)
tmp <- read.table('../1.data/3.windows/Ns/100kb_N.overlap')
tmp <- tmp %>% group_by(V1,V2,V3) %>% summarize('N' = sum(V7))
names(tmp) <- c('ss','start','end','N')
d100 <- left_join(d100,tmp,by=c('ss','start','end'))
d100 <- d100[d100$N < 10000 | is.na(d100$N) ,]
print(paste('100kb lost',round(100*(nr-nrow(d100))/nr,2),'% rows'))
# 10kb
nr <- nrow(d10)
tmp <- read.table('../1.data/3.windows/Ns/10kb_N.overlap')
tmp <- tmp %>% group_by(V1,V2,V3) %>% summarize('N' = sum(V7))
names(tmp) <- c('ss','start','end','N')
d10 <- left_join(d10,tmp,by=c('ss','start','end'))
d10 <- d10[d10$N < 1000 | is.na(d10$N),]
print(paste('10kb lost',round(100*(nr-nrow(d10))/nr,2),'% rows'))
# 1kb
nr <- nrow(d1)
tmp <- read.table('../1.data/3.windows/Ns/1kb_N.overlap')
tmp <- tmp %>% group_by(V1,V2,V3) %>% summarize('N' = sum(V7))
names(tmp) <- c('ss','start','end','N')
d1 <- left_join(d1,tmp,by=c('ss','start','end'))
d1 <- d1[d1$N < 100 | is.na(d1$N),]
print(paste('1kb lost',round(100*(nr-nrow(d1))/nr,2),'% rows'))
# 1Mb
nr <- nrow(d1m)
tmp <- read.table('../1.data/3.windows/Ns/1mb_N.overlap')
tmp <- tmp %>% group_by(V1,V2,V3) %>% summarize('N' = sum(V7))
names(tmp) <- c('ss','start','end','N')
d1m <- left_join(d1m,tmp,by=c('ss','start','end'))
d1m <- d1m[d1m$N < 10000 | is.na(d1m$N),]
print(paste('1mb lost',round(100*(nr-nrow(d1m))/nr,2),'% rows'))

# REMOVE ROWS WITH MISSING RECOMBINATION DATA 
#d100 <- d100[! is.na(d100$CH_cM) & ! is.na(d100$PT_cM) & ! is.na(d100$GB_cM),]
#d10 <- d10[! is.na(d10$CH_cM) & ! is.na(d10$PT_cM) & ! is.na(d10$GB_cM),]
#d1 <- d1[! is.na(d1$CH_cM) & ! is.na(d1$PT_cM) & ! is.na(d1$GB_cM),]
#d1m <- d1m[! is.na(d1m$CH_cM) & ! is.na(d1m$PT_cM) & ! is.na(d1m$GB_cM),]

# fix genes / repeats so that NA = 0 
d1m$ngenes[is.na(d1m$ngenes)] <- 0 
d1m$nrepeats[is.na(d1m$nrepeats)] <- 0 
d1m$cpgi[is.na(d1m$cpgi)] <- 0
d100$ngenes[is.na(d100$ngenes)] <- 0 
d100$nrepeats[is.na(d100$nrepeats)] <- 0 
d100$cpgi[is.na(d100$cpgi)] <- 0
d10$ngenes[is.na(d10$ngenes)] <- 0 
d10$nrepeats[is.na(d10$nrepeats)] <- 0 
d10$cpgi[is.na(d10$cpgi)] <- 0
d1$ngenes[is.na(d1$ngenes)] <- 0 
d1$Lgenes[is.na(d1$Lgenes)] <- 0 
d1$nrepeats[is.na(d1$nrepeats)] <- 0 
d1$Lrepeats[is.na(d1$Lrepeats)] <- 0 
d1$cpgi[is.na(d1$cpgi)] <- 0



# REMOVE CH OUTLIERS 
d100 <- d100[ d100$CH_cM < 50,]
d10 <- d10[ d10$CH_cM < 50,]
d1 <- d1[ d1$CH_cM < 50 ,]
d1m <- d1m[ d1m$CH13_cM < 50,]
d1 <- d1[! is.na(d1$ss),]

# # REMOVE Z SCAFFOLDS
z <- c('Super-Scaffold_13','Super-Scaffold_42')
d100z <- d100[ d100$ss %in% z, ]
d100 <- d100[! d100$ss %in% z, ]
d10z <- d10[ d10$ss %in% z, ]
d10 <- d10[! d10$ss %in% z, ]
d1z <- d1[ d1$ss %in% z, ]
d1 <- d1[! d1$ss %in% z, ]
d1mz <- d1m[d1m$ss %in% z,]
d1m <- d1m[! d1m$ss %in% z,]

# HOTSPOTS AS MORE THAN 5 TIMES THE MEAN + at least 2 windows  
hot1 <- which(d1$CH_hi >= 5)
hot1 <- hot1[which(d1$CH_cM[hot1] > 1 & d1$CH_cM[hot1] < 10)]
hot1 <- hot1[which(d1$Lrepeats[hot1] <= 500)]
hot1 <- hot1[unique(c(which(diff(hot1)==1), which(diff(hot1)==1)+1))]
pthot <- which(d1$PT_hi >= 10)
pthot <- pthot[unique(c(which(diff(pthot)==1), which(diff(pthot)==1)+1))]
gbhot <- which(d1$GB_hi >= 5)
gbhot <- gbhot[unique(c(which(diff(gbhot)==1), which(diff(gbhot)==1)+1))]
ch13hot <- which(d1$CH13_hi >= 5)
ch13hot <- ch13hot[unique(c(which(diff(ch13hot)==1), which(diff(ch13hot)==1)+1))]
ch13_2hot <- which(d1$CH13_2_hi >= 5)
ch13_2hot <- ch13_2hot[unique(c(which(diff(ch13_2hot)==1), which(diff(ch13_2hot)==1)+1))]
ch13_3hot <- which(d1$CH13_3_hi >= 5)
ch13_3hot <- ch13_3hot[unique(c(which(diff(ch13_3hot)==1), which(diff(ch13_3hot)==1)+1))]
ch13_4hot <- which(d1$CH13_4_hi >= 5)
ch13_4hot <- ch13_4hot[unique(c(which(diff(ch13_4hot)==1), which(diff(ch13_4hot)==1)+1))]
ch13_5hot <- which(d1$CH13_5_hi >= 5)
ch13_5hot <- ch13_5hot[unique(c(which(diff(ch13_5hot)==1), which(diff(ch13_5hot)==1)+1))]

# overlap tables for Figure 4 

hotmat <- matrix(data=NA,nrow=length(hot1),ncol=81)
pthotmat <- matrix(data=NA,nrow=length(hot1),ncol=81)
gbhotmat <- matrix(data=NA,nrow=length(hot1),ncol=81)
ch13hotmat <- matrix(data=NA,nrow=length(hot1),ncol=81)
ch13_2hotmat <- matrix(data=NA,nrow=length(hot1),ncol=81)
ch13_3hotmat <- matrix(data=NA,nrow=length(hot1),ncol=81)
ch13_4hotmat <- matrix(data=NA,nrow=length(hot1),ncol=81)
ch13_5hotmat <- matrix(data=NA,nrow=length(hot1),ncol=81)

for(i in 1:length(hot1)){
  r1 <- (hot1[i]-40):(hot1[i]+40)
  hotmat[i,] <- d1$CH_cM[r1]/mean(d1$CH_cM[r1],na.rm=T)
  pthotmat[i,] <- d1$PT_cM[r1]/mean(d1$PT_cM[r1],na.rm=T)
  gbhotmat[i,] <- d1$GB_cM[r1]/mean(d1$GB_cM[r1],na.rm=T)
  ch13hotmat[i,] <- d1$CH13_cM[r1]/mean(d1$CH13_cM[r1],na.rm=T)
  ch13_2hotmat[i,] <- d1$CH13_2_cM[r1]/mean(d1$CH13_2_cM[r1],na.rm=T)
  ch13_3hotmat[i,] <- d1$CH13_3_cM[r1]/mean(d1$CH13_3_cM[r1],na.rm=T)
  ch13_4hotmat[i,] <- d1$CH13_4_cM[r1]/mean(d1$CH13_4_cM[r1],na.rm=T)
  ch13_5hotmat[i,] <- d1$CH13_5_cM[r1]/mean(d1$CH13_5_cM[r1],na.rm=T)
}

#overlaps <- read.table('overlaps.table')

d10 <- d10[!is.na(d10$ss),]
# dataset with propseq and late gini 
dfd <- data.frame(ss=unique(d10$ss),propseq = NA)
dfd$gini <- NA 
for(ss in dfd$ss){
  dfd$gini[dfd$ss==ss] <- Gini(d10$CH_cM[d10$ss==ss],na.rm=T)
}
# make dataset 
d10$lg <- lms$lg[match(d10$ss,lms$ss)]
ss10 <- d10 %>% group_by(lg) %>% summarize(ss= unique(ss),
                                           cM = sum(CH_cM,na.rm=T) / 100,
                                           cMMb = median(CH_cM,na.rm=T),
                                           pi = median(CH_pi,na.rm=T),
                                           gc = median(gc,na.rm=T),
                                           cpgi = sum(cpgi,na.rm=T),
                                           len = max(end)-min(start),
                                           cpgi_mb = sum(cpgi,na.rm=T) / (len ),
                                           loci = sum(CH_snps,na.rm=T))
ss10 <- left_join(ss10,dfd[,c(1,3)],by='ss')
# 
genes <- read.table('../1.data/3.windows/Gene_windows/genes.bed')
genes$dif <- genes$V3-genes$V2
gd <- genes %>% group_by(ss=V1) %>% summarize('genes' = n(),
                                              'lgenes' = sum(dif))
ss10$genes <- gd$genes[match(ss10$ss,gd$ss)]
ss10$lgenes <- gd$genes[match(ss10$ss,gd$ss)]
ss10$lgene_Mb <- ss10$lgenes / (ss10$len / 1e6)
ss10$gene_Mb <- ss10$genes / (ss10$len / 1e6)

d10 <- d10[d10$CH_snps < 200,]
d10 <- d10[d10$ss != 'Super-Scaffold_49',]
cold10 <- which(d10$CH_cM / mean(d10$CH_cM,na.rm=T)  <= .04) # CH
hot10 <- which(d10$CH_cM / mean(d10$CH_cM,na.rm=T)  >= 5) # CH

random10 <- which(d10$CH_cM / mean(d10$CH_cM,na.rm=T) >= .02 & 
                    d10$CH_cM / mean(d10$CH_cM,na.rm=T) <= 5) 
random10 <- sample(random10,4000)

# PREFIGURE DATA for Figure 3 
# matrix of TES and TSS 
d1$ts <- d1$tss
ts <- which(d1$tss > 0)
ts <- ts[ts > 41]
tsmat <- matrix(data=NA,nrow=length(ts),ncol=81)
for(i in 1:length(ts)){
  r1 <- (ts[i] - 40) : (ts[i] + 40)
  tsmat[i,] <- d1$CH_cM[r1] / mean(d1$CH_cM[r1],na.rm=T)
}
# matrix of CPGIs
cp <- which(d1$cpgi > 500)
cp <- cp[cp > 41]
cpmat <- matrix(data=NA,nrow=length(cp),ncol=81)
for(i in 1:length(cp)){
  r1 <- (cp[i] - 40) : (cp[i] + 40)
  cpmat[i,] <- d1$CH_cM[r1] / mean(d1$CH_cM[r1],na.rm=T)
}

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
# make sure lm is / lg 

# colours
crp <- viridis::mako(nrow(lms2),alpha=1,direction=1)
lms2$col <-''
lms2$col[order(lms2$len,decreasing = T)] <- crp
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
  lm1 <- lm(av~len,data=lms2)
  # plot 
  plot(y=lms2$av,x=lms2$len, col=lms2$col, xlab='', 
       ylab='', main='', axes=F,pch=16, 
       cex=4, xlim=c(0,100),ylim=c(0,100))
  # peripherals 
  axis(1, at = seq(0,100,25), labels=seq(0,100,25),cex.axis=axis.size)
  axis(2, at = c(0,20,41,60,80,100), labels= c(0,20,41,60,80,100),
       cex.axis=axis.size,las=2) 
  mtext('Genetic length (cM)',2,2.5,cex=mtext.size)
  mtext('Physical length (Mb)',1,2.5,cex=mtext.size)
  lnew <- data.frame(len=seq(0,100,length.out=1000))
  p <- predict(lm1,newdata=lnew,type = 'response',se.fit = T)
  lines(p$fit~lnew$len,lwd=2,col='black')
  # now adds white dot for less than 50cM
  points(y=lms2$av[lms2$av < 50],x=lms2$len[lms2$av < 50],cex=1,
         pch='X',col='white')
  points(y=lms2$av[lms2$lg==40],x=lms2$len[lms2$lg==40],cex=2,pch='z',col='white')
  lines(y=p$fit-1.96*p$se.fit,x=lnew$len,lty=2,lwd=2,col='black')
  lines(y=p$fit+1.96*p$se.fit,x=lnew$len,lty=2,lwd=2,col='black')
  mtext('A',3,at=0,cex=2)
  # legend 
  crp <- viridis::mako(nrow(lms2),alpha=1,direction=-1)
  for(i in seq(1,39,3)){
    points(98,11+i/2,pch=20,cex=5,col=crp[i])
  }
  text(98, 36, "Mb" , cex = 1.2)
  text(94, 10, "6" , cex = 1)
  text(94, 32, "90" , cex = 1)
  segments(94, 12, 94,30 )

}
### plot 2 - SEXES 
# remove Z 
#lms2 <- lms2[lms2$lg!=40,]

{
  plot(lms2$f~lms2$m,xlab='', ylab='',
       col=lms2$col, pch=16, cex=4,
       xlim=c(20,100),ylim=c(20,100), axes=F, main ='') 
  axis(1,at=seq(20,100,20),labels=seq(20,100,20),cex.axis=axis.size)
  axis(2,at=seq(20,100,20),labels=seq(20,100,20),cex.axis=axis.size,las=2)
  mtext('Male genetic length (cM)',1,line=2.5,cex=mtext.size)
  mtext('Female genetic length (cM)',2,line=2.5,cex=mtext.size)
  # Z 
  points(y=lms2$f[lms2$lg==40],x=lms2$m[lms2$lg==40],cex=2,pch='z',col='white')
  points(y=lms2$f[lms2$av < 50],x=lms2$m[lms2$av < 50],cex=1,
         pch='X',col='white')
  abline(0,1,lty=2,lwd=2)
  mtext('B',3,at=18,cex=2)
  # legend 
  crp <- viridis::mako(nrow(lms2),alpha=1,direction=-1)
  for(i in seq(1,39,3)){
    points(98,21+i/2,pch=20,cex=5,col=crp[i])
  }
  text(98, 45, "Mb" , cex = 1.2)
  text(94, 20, "6" , cex = 1)
  text(94, 42, "90" , cex = 1)
  segments(94, 22, 94,40 )
  
  
}
### plot 3 - CO location density
{
  #????#
  # lm$prop <- 0
  # for(lg in lgs){
  #   lm$prop[lm$lg==lg] <- lm$newpos[lm$lg==lg] / max(lm$newpos[lm$lg==lg])
  #   if(lg==1){
  #     lm$prop[lm$lg==lg] <- (lm$newpos[lm$lg==lg] - len2.2)/ (max(lm$newpos[lm$lg==lg]) - len2.2)
  #   }
  # }
  male <- density(lm$dist[lm$maleCOs>0 & ! is.na(lm$maleCOs) & lm$lg!=40],from=0,to=.5)
  female <- density(lm$dist[lm$femaleCOs>0& ! is.na(lm$femaleCOs)& lm$lg!=40],from=0,to=.5)
  plot(x=male$x,y=male$y/sum(male$y),
       col=malecolor,lwd=4,type='l',axes=F,xlab='',ylab='')
  axis(1,at=seq(0,.5,0.1),labels = seq(0,50,10),cex.axis=1.2)
  axis(2,at=seq(0,.006,.002),labels=seq(0,6,2),cex.axis=1.2,las=2)
  lines(x=female$x,y=female$y/sum(female$y),col=femalecolor,lwd=4,lty=1)
  mtext('Density (x 1e-3)',2,2.5,cex=mtext.size)
  mtext('Distance from LG end (% of LG length)',1,2.5,cex=mtext.size)
  mtext('C',3,at=-0.01,cex=2)
  legend('topright',lty=c(1,1),col=c(malecolor,femalecolor),lwd=4,
         legend=c('Male','Female'),cex=1.5,bty='n')
  # test significance 
  ks.test(male$y,female$y)
}
### plot 4 - Barplot of r intra
{
  lms2$difr <- NA
  lms2$diff <- NA
  for(i in 1:nrow(lms2)){
    lms2$difr[i] <- (lms2$rintra_m[i]-lms2$rintra_f[i]) / lms2$rintra_f[i]
    lms2$diff[i] <- ( lms2$m[i] - lms2$f[i] ) / lms2$f[i]
  }
  lms2$difr[lms2$lg==40] <- NA
  par(mar=c(5,4,4,4))
  barplot(lms2$difr[order(lms2$len,decreasing=F)],
          xlab='', width=2,
          col=c(malecolor,femalecolor)[
            ifelse(lms2$difr[order(lms2$len,decreasing=F)] > 0,1,2)],
          axes=F, space=1, xlim=c(-.8,.8), ylim=c(5,155),
          horiz=T, border=NA)
  axis(1, at = seq(-.8,.8,.4),labels=c(80,40,0,40,80),cex.axis=1.3,las=1)
  mtext('Proportional increase in intrachromosomal shuffling (%)',1,2.5,cex=1)
  axis(2, at=seq(3,159,4), lms2$lg[order(lms2$len,decreasing=F)],
       cex.axis=1.1, las=2, adj=1)
  mtext('Linkage Group',2,2.5,cex=mtext.size)
  abline(v=0,lwd=2)
  legend('topleft', pch=c(15,15), pt.cex=c(2,2), col=c(malecolor,femalecolor),
         legend=c('Male > Female','Female > Male'), cex=1.5, bty='n')
  mtext('D',3,at=-.9,cex=2)
  axis(4,at=seq(3,159,4), labels=round(lms2$len[order(lms2$len,decreasing = F)],0),
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
    # make ln (x+1)
    d1m$x <- log1p(d1m$CH_cM)
    d1m$y <- log1p(d1m$av)
    plot(log1p(d1m$av)~log1p(d1m$CH_cM),pch=20,col=v2a[2],xlim=c(0,3),ylim=c(0,3),
         axes=F,xlab='',ylab='',cex=3)
    lm2 <- lm(y~x, data=d1m)
    p2 <- predict(lm2, newdata=data.frame('x'=seq(0,3,length.out=1000)),type='response',se.fit = T)
    lines(y=p2$fit,x=seq(0,3,length.out=1000),lwd=2,col=v1[1])
    lines(y=p2$fit-1.96*p2$se.fit,x=seq(0,3,length.out=1000),lty=2,lwd=2,col=v1[1])
    lines(y=p2$fit+1.96*p2$se.fit,x=seq(0,3,length.out=1000),lty=2,lwd=2,col=v1[1])
    axis(1,at=seq(0,3,1),labels=round(expm1(seq(0,3,1)),1),cex.axis=1.2)
    axis(2,at=seq(0,3,1),labels=round(expm1(seq(0,3,1)),1),las=2,cex.axis=1.2)
    mtext('LD-based estimate - log(x+1)',1,2.4,cex=1.2)
    mtext('Family estimate - log(x+1)',2,2.8,cex=1.2)
    mtext('A',3,at=0,cex=2)
  }
  {
    # PLOT 2 
    ## PLOT SCAFFOLDS
    ss='Super-Scaffold_30'
    plot(d10$CH_cM[d10$ss==ss]~
           d10$start[d10$ss==ss],type='l',
         axes=F,xlim=c(0,8e6),ylim=c(-2,40),
         xlab='',ylab='',main='',col=v1[2],lwd=2)
    abline(h=mean(d10$CH_cM,na.rm=T),lty=2,col='black',lwd=2)
    abline(h=mean(d10$CH_cM[d10$ss==ss],na.rm=T),col=v1[2],lwd=2)
    axis(1,at=c(0,8e6),labels=c(0,8),cex.axis=axis.size)
    mtext('Mb',1,2.5,cex = label.size)
    points(max(d10$start[d10$ss==ss],na.rm=T)/2,35,pch=15,cex=4,col=v1[2])
    axis(2,at=seq(0,40,10),labels=seq(0,40,10),cex.axis=axis.size,las=2)
    segments(x0=lm$pos[lm$ss==ss & lm$femaleCOs >0 ],
             x1=lm$pos[lm$ss==ss & lm$femaleCOs >0 ],
             y0=-2,y1=-1.1,lwd=2,col=femalecolor)
    segments(x0=lm$pos[lm$ss==ss & lm$maleCOs >0 ],
             x1=lm$pos[lm$ss==ss & lm$maleCOs >0 ],
             y0=-1,y1=-.2,lwd=2,col=malecolor)
    mtext('cM/Mb',2,2.5,cex = label.size)
    mtext('LG 36',3,0,cex = label.size)
    mtext('B',3,at=1,cex=2)
    #
    ss = 'Super-Scaffold_39'
    plot(d10$CH_cM[d10$ss==ss]~
           d10$start[d10$ss==ss],type='l',
         axes=F,xlim=c(0,15e6),ylim=c(-2,40),
         xlab='',ylab='',main='',col=v1[3],lwd=2)
    axis(1,at=c(0,8e6,15e6),labels=c(0,8,15),cex.axis=axis.size)
    abline(h=mean(d10$CH_cM,na.rm=T),lty=2,col='black',lwd=2)
    abline(h=mean(d10$CH_cM[d10$ss==ss],na.rm=T),col=v1[3],lwd=2)
    segments(x0=lm$pos[lm$ss==ss & lm$femaleCOs >0 ],
             x1=lm$pos[lm$ss==ss & lm$femaleCOs >0 ],
             y0=-2,y1=-1.1,lwd=2,col=femalecolor)
    segments(x0=lm$pos[lm$ss==ss & lm$maleCOs >0 ],
             x1=lm$pos[lm$ss==ss & lm$maleCOs >0 ],
             y0=-1,y1=-.2,lwd=2,col=malecolor)
    mtext('Mb',1,2.5,cex = label.size)
    mtext('LG 32',3,0,cex = label.size)
    points(max(d10$start[d10$ss==ss],na.rm=T)/2,35,pch=18,cex=4,col=v1[3])
    axis(2,at=seq(0,40,10),labels=seq(0,40,10),cex.axis=axis.size,las=2)
    # 
    ss = 'Super-Scaffold_8'
    plot(d10$CH_cM[d10$ss==ss]~
           d10$start[d10$ss==ss],type='l',
         axes=F,xlim=c(0,25e6),ylim=c(-2,40),
         xlab='',ylab='', main='',col=v1[4],lwd=2)
    points(max(d10$start[d10$ss==ss],na.rm=T)/2,35,
           pch=17,cex=4,col=v1[4])
    abline(h=mean(d10$CH_cM,na.rm=T),lty=2,col='black',lwd=2)
    abline(h=mean(d10$CH_cM[d10$ss==ss],na.rm=T),col=v1[4],lwd=2)
    axis(1,at=c(0,8e6,15e6,25e6),labels=c(0,8,15,25),
         cex.axis=axis.size)
    mtext('Mb',1,2.5,cex = label.size)
    mtext('LG 22',3,0,cex = label.size)
    segments(x0=lm$pos[lm$ss==ss & lm$femaleCOs >0 ],
             x1=lm$pos[lm$ss==ss & lm$femaleCOs >0 ],
             y0=-2,y1=-1.1,lwd=2,col=femalecolor)
    segments(x0=lm$pos[lm$ss==ss & lm$maleCOs >0 ],
             x1=lm$pos[lm$ss==ss & lm$maleCOs >0 ],
             y0=-1,y1=-.2,lwd=2,col=malecolor)
    axis(2,at=seq(0,40,10),labels=seq(0,40,10),cex.axis=axis.size,las=2)
  }
  {
    # PLOT 3 - HOTSPOT PLOT 
    # use only linkage map scaffolds 
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
      dfd$propsec[dfd$ss==ss] <- y1[which(abs(x1-0.8) == min(abs(x1-0.8),na.rm=T))]
    }
    # make overall lines 
    x1 <- cumsum(d10$CH_cM[order(d10$CH_cM,decreasing = T)])/sum(d10$CH_cM,na.rm=T)
    y1 <- seq(0,1,length.out=nrow(d10))
    lines(x=x1, y=y1, col='grey',lwd=4,lty=2)
    # add example points - maybe remove 
    points(x=.68,y=.37,pch=15,cex=3,col=v1[2])
    points(x=.75,y=.29,pch=18,cex=4,col=v1[3])
    points(x=.87,y=.23,pch=17,cex=3,col=v1[4])
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
    plot(y=ss10$gini,x=ss10$len/1e6,pch=pch1,cex=cex1,xlab='',
         ylab='',ylim=c(.35,.85),xlim=c(0,70),axes=F,col=col1)
    axis(2,at=c(0.4,.6,.8),labels = c(0.4,.6,.8),cex.axis=axis.size,las=2)
    abline(h=0.68,lty=2)
    mtext('Gini Coefficient',2,2.6,cex = label.size)
    axis(1,at=c(0,25,50,75),labels=c(0,25,50,75),cex.axis=axis.size)
    mtext('Length (Mb)',1,2.6,cex = label.size)
    mtext('D',3,cex=2,at=0)
    plot(y=ss10$gini,x=ss10$pi,pch=pch1,cex=cex1,axes=F,xlab='',
         ylab='',ylim=c(.35,.85),xlim=c(0.0012,0.0027),col=col1)
    axis(2,at=c(0.4,.6,.8),labels = c(0.4,.6,.8),cex.axis=axis.size,las=2)
    abline(h=0.68,lty=2)
    axis(1,at=c(0.0012,0.0017,0.0022,0.0027),labels=c(1.2,1.7,2.2,2.7),cex.axis=axis.size)
    mtext('pi (x1000)',1,2.6,cex = label.size)
    mtext('E',3,at=1.2e-3,cex=2)
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
  #general 
  plotit = T
  if(plotit){
    pdf('plots/v1_figure3.pdf',width=12,height=12)
  }
  
  m1 <- matrix(c(1,7,7,
                 8,3,4,
                 2,5,6),
               byrow = T,nrow=3,ncol=3)
  layout(m1,
         heights = c(1,1,1),widths=c(1,1,1))
  axis.size = 1.2
  mtext.size = 1 
  viocol = viridis::mako(4)[c(4,3,2)]
  
  # GLOBAL HOTSPOTS
  hist(log(d10$CH_cM/mean(d10$CH_cM,na.rm=T)),101,xlab='',ylab='',xlim=c(-6,4),
       axes=F,main='',col=viocol[1],border = NA) -> h
  hist(log(d10$CH_cM[hot10]/mean(d10$CH_cM,na.rm=T)),h$breaks,add=T,col=viocol[3],
       border=NA)
  abline(v=log(5),col=viocol[3],lwd=2,lty=2)
  hist(log(d10$CH_cM[cold10]/mean(d10$CH_cM,na.rm=T)),h$breaks,add=T,col=viocol[2],
       border=NA)
  abline(v=log(.04),col=viocol[2],lwd=2,lty=2)
  axis(1,at=seq(-6,4,2),labels = seq(-6,4,2),cex.axis=axis.size)
  mtext('RRR (log)',1,2.5,cex=mtext.size)
  axis(2,at=seq(0,3e3,1e3),labels=c(0,1,2,3),cex.axis=axis.size,las=2)
  mtext('Counts (x1000)',2,2.5,cex=mtext.size)
  mtext('A',3,at=-6,cex=2)
  
  # TSS / CPGIS
  elcol <- viridis::mako(4)[c(2,3)]
  plot(colMeans(tsmat,na.rm=T),type='l',axes=F,xlab='',ylab='',lwd=4,
       xlim=c(1,81),ylim=c(.9,1.8),main='',
       col=elcol[1],lty=1)
  lines(colMeans(cpmat,na.rm=T),col=elcol[2],lwd=4,lty=2)
  axis(1,at=c(1,21,41,61,81),labels=c(-40,-20,0,20,40),cex.axis=axis.size)
  axis(2,at=c(.9,1.2,1.6),labels=c(.9,1.2,1.6),cex.axis=axis.size,las=2)
  mtext('Distance from element (kb)',1,2.5,cex=mtext.size)
  mtext('RRR80',side=2,line=2.5,cex=mtext.size)
  legend('topleft',lty=c(1,3),col=elcol,legend=c('TSS','CGIs'),lwd=4,bty='n',
         cex=1)
  abline(v=41,lwd=1,lty=1)
  mtext('G',3,at=1,cex=2)
  
  # VIOLIN PLOTS
  ## Distance to end 
  def <- c(5.1,4.1,4.1,2.1)
  par(mar=c(2.1,4.1,2.1,2.1))
  vioplot::vioplot(d10$dist[cold10],d10$dist[random10],d10$dist[hot10],col=viocol[c(2,1,3)],
                   ylim=c(0,.55), axes=FALSE, xaxt='n',yaxt='n',
                   main='')
  # axis(1, at = c(1,2,3), labels=c('deserts','random','jungles'),cex.axis=axis.size)
  axis(2,at = seq(0,0.5,0.1),labels=seq(0,50,10),cex.axis=axis.size,las=2)
  #mtext('***',line = -.8,at=2,cex=mtext.size)
  mtext('Distance from LG end (%)',2,2.5,cex=1.1)
  mtext('C',3,at=.4,cex=2)
  ### GC 
  vioplot::vioplot(100*d10$gc[cold10],d10$gc[random10]*100,100*d10$gc[hot10],col=viocol[c(2,1,3)],main='',
                   ylim=c(20,88), axes=FALSE, xaxt='n',yaxt='n')
  # axis(1, at = c(1,2,3), labels=c('deserts','random','jungles'),cex.axis=axis.size)
  axis(2,at = seq(20,80,20),labels=seq(20,80,20),cex.axis=axis.size,las=2)
  mtext('GC %',2,2.5,cex=mtext.size)
  mtext('D',3,at=.4,cex=2)
  par(mar=c(4.1,4.1,2.1,2.1))
  ### PI 
  vioplot::vioplot(d10$CH_pi[cold10],d10$CH_pi[random10],d10$CH_pi[hot10],col=viocol[c(2,1,3)],main='',
                   ylim=c(0,0.005), axes=FALSE, xaxt='n',yaxt='n')
  axis(1, at = c(1,2,3), labels=c('coldspots','random','hotspots'),cex.axis=axis.size)
  axis(2,at = seq(0,0.005,0.001),labels=seq(0,5,1),cex.axis=axis.size,las=2)
  mtext('PI x(1000)',2,2.5,cex=mtext.size)
  mtext('windows',1,2.5,cex=mtext.size)
  mtext('E',3,at=.4,cex=2)
  ### genes 
  vioplot::vioplot(d10$ngenes[cold10],d10$ngenes[random10],d10$ngenes[hot10],col=viocol[c(2,1,3)],main='',
                   ylim=c(0,8), axes=FALSE, xaxt='n',yaxt='n')
  axis(1, at = c(1,2,3), labels=c('coldspots','random','hotspots'),cex.axis=axis.size)
  axis(2,at = seq(0,8,2),labels=seq(0,8,2),cex.axis=axis.size,las=2)
  mtext('Gene count',2,2.5,cex=mtext.size)
  mtext('windows',1,2.5,cex=mtext.size)
  mtext('F',3,at=.4,cex=2)
  par(mar=def)
  
  ### FINE SCALE ### 
  
  ### example chr
  ss='Super-Scaffold_39'
  tmp1 <- d10[hot10,]
  tmp2 <- d10[cold10,]
  plot(y=d10$CH_cM[d10$ss==ss]/mean(d10$CH_cM,na.rm=T),
       x= d10$start[d10$ss==ss],pch=20, 
       col=add.alpha('grey70',.5),axes=F,xlab='',ylab='',
       ylim=c(-2,18),cex=1.5,xlim=c(0,15e6),main='')
  axis(1,at=seq(0,15e6,5e6),labels=seq(0,15,5),cex.axis=axis.size)
  mtext('Position (Mb)',1,2.5,cex=mtext.size)
  mtext('LG 32',3,0,cex=mtext.size)
  axis(2,at=c(0,5,10,15),labels=c(0,5,10,15),las=2,cex.axis=axis.size)
  mtext('RRR',2,2.5,cex=mtext.size)
  points(y=tmp1$CH_cM[tmp1$ss==ss]/mean(d10$CH_cM,na.rm=T),
         x=tmp1$start[tmp1$ss==ss], col=viocol[3],pch=20,
         cex=3)
  points(y=tmp2$CH_cM[tmp2$ss==ss]/mean(d10$CH_cM,na.rm=T), pch=20,
         x=tmp2$start[tmp2$ss==ss], col=add.alpha(viocol[2],.6),
         cex=3)
  mtext('B',3,at=0,cex=2)
  segments(x0=lm$pos[lm$ss==ss & lm$femaleCOs >0 ],
           x1=lm$pos[lm$ss==ss & lm$femaleCOs >0 ],
           y0=-2,y1=-1.1,lwd=2,col=add.alpha(femalecolor,.3))
  segments(x0=lm$pos[lm$ss==ss & lm$maleCOs >0 ],
           x1=lm$pos[lm$ss==ss & lm$maleCOs >0 ],
           y0=-1,y1=-.2,lwd=2,col=add.alpha(malecolor,.3))
  # abline(h=5,col=viocol[3],lwd=2)
  # mtext('D',3,at=0,cex=2)
  # l1 <- loess(d10$ngenes[d10$ss==ss]~d10$start[d10$ss==ss],span=.05)
  # lines(l1)
  # axis(4,at=seq(0,3,1))
  
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

####################
### Fig 4 part 2 ###
####################

# make pdf and merge with map in illustrator
{
  # global parameters 
  mtext.size=1.5
  axis.size=1.2
  popcol = viridis::mako(6,direction=1)[c(4,2,3,1,5)]
  plotit = F
  if(plotit){
    pdf('plots/v1_figure4.pdf',12,12)
  }
  layout(matrix(c(1,3,2,4),nrow=2,ncol=2),heights = c(.5,.5),widths = c(.5,.5))
  
  # 1 empty to fill with map 
  plot(0,pch='',xlab='',axes=F,ylab='')
  mtext('A',3,cex=2,adj=-.1)
  
  ## 2 Example LG32 
  {
    par(mar=c(4,3,2,1)+0.1)
    ss = 'Super-Scaffold_14'
    ## GB
    plot(d1$GB_cM[d1$ss==ss& d1$end < 5e6]~
           d1$end[d1$ss==ss& d1$end < 5e6],
         type='l',lwd=2,ylim=c(0,200),col=popcol[1],
         xlab='',ylab='',axes=F,pch=20,xlim=c(0,5.4e6),main=paste('LG',lms$lg[lms$ss==ss]))
    segments(x=0,x1=5e6,y0=10*mean(d1$CH13_cM),y1=10*mean(d1$CH13_cM),lty=2)
    tmp3 <- d1[gbhot,]
    points(y=rep(30,length(tmp3$end[tmp3$ss==ss& tmp3$end < 5e6])),
           x=tmp3$end[tmp3$ss==ss& tmp3$end < 5e6], col=popcol[1],
           pch=16,cex=1)
    ## PT
    lines(40+d1$PT_cM[d1$ss==ss& d1$end < 5e6]~
            d1$end[d1$ss==ss& d1$end < 5e6],
          type='l',lwd=2,col=popcol[2])
    segments(x=0,x1=5e6,y0=40+10*mean(d1$CH13_cM),y1=40+10*mean(d1$CH13_cM),lty=2)
    tmp2 <- d1[pthot,]
    points(y=rep(70,length(tmp2$end[tmp2$ss==ss& tmp2$end < 5e6])),
           x=tmp2$end[tmp2$ss==ss & tmp2$end < 5e6], col=popcol[2],
           pch=16, cex=1)
    # CH
    lines(80+d1$CH_cM[d1$ss==ss & d1$end < 5e6 ]~
            d1$end[d1$ss==ss& d1$end < 5e6],
          type='l',lwd=2,col=popcol[3])
    segments(x=0,x1=5e6,y0=80+10*mean(d1$CH13_cM),y1=80+10*mean(d1$CH13_cM),lty=2)
    tmp1 <- d1[hot1,]
    points(y=rep(110,length(tmp1$end[tmp1$ss==ss & tmp1$end < 5e6])),
           x=tmp1$end[tmp1$ss==ss& tmp1$end < 5e6], col=popcol[3],
           pch=16, cex=1)
    ## CH13 
    lines(120+d1$CH13_cM[d1$ss==ss & d1$end < 5e6]~
            d1$end[d1$ss==ss & d1$end < 5e6],
          type='l',lwd=2,col=popcol[4],lty=1)
    segments(x=0,x1=5e6,y0=120+10*mean(d1$CH13_cM),y1=120+10*mean(d1$CH13_cM),lty=2)
    tmp3 <- d1[ch13hot,]
    points(y=rep(150,length(tmp3$end[tmp3$ss==ss & tmp3$end < 5e6])),
           x=tmp3$end[tmp3$ss==ss & tmp3$end < 5e6], col=popcol[4],
           pch=16,cex=1)
    ## CH13_2
    lines(160+d1$CH13_2_cM[d1$ss==ss & d1$end < 5e6]~
            d1$end[d1$ss==ss & d1$end < 5e6],
          type='l',lwd=2,col=popcol[5],lty=1)
    segments(x=0,x1=5e6,y0=160+10*mean(d1$CH13_cM),y1=160+10*mean(d1$CH13_cM),lty=2)
    tmp3 <- d1[ch13_2hot,]
    points(y=rep(190,length(tmp3$end[tmp3$ss==ss& tmp3$end < 5e6])),
           x=tmp3$end[tmp3$ss==ss& tmp3$end < 5e6], col=popcol[5],
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
    call <- cor(d1[,c(13,11,5,7,9)],use='complete.obs')
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
  # {
  #   venn <- VennDiagram::venn.diagram(list(hot1,ch13hot,ch13_2hot),
  #                                     category.names=c('CH','PT','GB'),
  #                                     filename=NULL,
  #                                     disable.logging = T,
  #                                     label.col = 'white',
  #                                     cex=1.5,
  #                                     main='Local Hotspot overlap',
  #                                     fill=popcol[c(3,2,1)],
  #                                     col=rep('black',3),
  #                                     lwd=rep(2,3))
  #   vp <- viewport(x=.52,y=.4,width=.4,height=.4,just=c('left','top'))
  #   pushViewport(vp)
  #   grid.draw(venn)
  #   grid.draw(venn)
  #   popViewport()
  #   plot(0,pch='',xlab='',axes=F,ylab='')
  #   mtext('d',3,cex=2,adj=0.05)
  #   par(mar=c(5.1, 4.1, 4.1, 2.1))
  #   
  # }
  
  {
    par(mar=c(4,4,4,4))
    
    
    plot(colMeans(hotmat,na.rm=T),axes=F,xlab='',ylab='',xlim=c(21,61),
         ylim=c(0,10),col=popcol[3],lwd=3,cex=2,type='l',main='Hotspot sharing')
    lines(colMeans(gbhotmat,na.rm=T),col=popcol[1],lwd=3)
    lines(colMeans(ch13hotmat,na.rm=T),col='grey')
    lines(colMeans(ch13_2hotmat,na.rm=T),col='grey')
    lines(colMeans(ch13_3hotmat,na.rm=T),col='grey')
    lines(colMeans(ch13_4hotmat,na.rm=T),col='grey')
    lines(colMeans(ch13_5hotmat,na.rm=T),col='grey')
    lines(colMeans(pthotmat,na.rm=T),col=popcol[2],lwd=3)
    legend('topright',lwd=c(3,3,3,3),col=c(popcol[c(1,2,3)],'grey'),
           legend=c('GB','PT','CH','CH subsets'),bty='n',cex=1)
    abline(v=41,lwd=0.2,lty=2,col=add.alpha('grey',.5))
    axis(1,at=seq(21,61,10),labels=seq(-20,20,10))
    axis(4,at=seq(0,10,5))
    mtext('RRR80',4,2)
    mtext('Distance from CH hotspot (kb)',1,2)
    mtext('D',3,cex=2,adj=-.1)
    
  }
  if(plotit){
    dev.off()
  }
}


##### ENO OF MANUSCRIPT FIGURES ######


