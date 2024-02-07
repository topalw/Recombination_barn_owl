
{
  # only works in Rstudio : 
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
  library(tools)
  library(sf)
  colp=c('#FF0000','orange','violet',
         '#006600','#012169')
  mtext.size = 1.1 
  axis.size=1.1
}
# ACQUIRE DATA AND CURATE #
{ 
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
  lm <- read.table('../1.data/1.LepMAP/orders/fitted_values_gams.table',h=T)
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

#######################################
#######################################
###### Supplementary table XX #########
#######################################
#######################################
{# metadata table 
  meta <- read.table('../2.metadata/metadata_Bioproject_v4.tsv',h=T,sep='\t')
  meta <- meta[meta$REF_Panel==1,c(2,3,4,6,7,8,9,11,14,15)]
  meta$Bioproject[is.na(meta$Bioproject)] <- '666'
  names(meta) <- c('RingID', 'Sample_name', 'Population','Sex', 'Read_count',
                   'Coverage', 'Bioproject', 'dataset', 'Longitude', 'Latitude')
  lm3 <- read.table('../2.metadata/8_LM_250_RP502.txt')[2,-c(1,2)]
  lm.inds <- unique(unname(unlist(lm3)))
  meta$dataset <- ''
  meta$dataset[meta$RingID %in% lm.inds] <- 'LM3'
  p.s <- read.csv('../1.data/2.pyrho/smcpp/1.smcfiles/CH_samples.csv',h=F)
  p.s <- gsub('X', '', unname(unlist(p.s)))
  meta$dataset[meta$RingID %in% p.s] <- ifelse(meta$dataset[meta$RingID %in% p.s] == '',
                                               'CH', 
                                               paste(meta$dataset[meta$RingID %in% p.s],
                                                     'CH',sep=';'))
  p.s <- read.csv('../1.data/2.pyrho/smcpp/1.smcfiles/CH13_samples.csv',h=F)
  p.s <- gsub('X', '', unname(unlist(p.s)))
  meta$dataset[meta$RingID %in% p.s] <- ifelse(meta$dataset[meta$RingID %in% p.s] == '',
                                               'CH13', 
                                               paste(meta$dataset[meta$RingID %in% p.s],
                                                     'CH13',sep=';'))
  p.s <- read.csv('../1.data/2.pyrho/smcpp/1.smcfiles/CH13_2_samples.csv',h=F)
  p.s <- gsub('X', '', unname(unlist(p.s)))
  meta$dataset[meta$RingID %in% p.s] <- ifelse(meta$dataset[meta$RingID %in% p.s] == '',
                                               'CH13_2', 
                                               paste(meta$dataset[meta$RingID %in% p.s],
                                                     'CH13_2',sep=';'))
  p.s <- read.csv('../1.data/2.pyrho/smcpp/1.smcfiles/GB_samples.csv',h=F)
  p.s <- gsub('X', '', unname(unlist(p.s)))
  meta$dataset[meta$RingID %in% p.s] <- ifelse(meta$dataset[meta$RingID %in% p.s] == '',
                                               'GB', 
                                               paste(meta$dataset[meta$RingID %in% p.s],
                                                     'GB',sep=';'))
  p.s <- read.csv('../1.data/2.pyrho/smcpp/1.smcfiles/PT_samples.csv',h=F)
  p.s <- gsub('X', '', unname(unlist(p.s)))
  meta$dataset[meta$RingID %in% p.s] <- ifelse(meta$dataset[meta$RingID %in% p.s] == '',
                                               'PT', 
                                               paste(meta$dataset[meta$RingID %in% p.s],
                                                     'PT',sep=';'))
  write.csv(meta,'../2.metadata/Supp_table_samples.csv',quote=F,row.names = F)
}

#######################################
#######################################
###### Supplementary table XX #########
#######################################
#######################################
{ # Linkage group info supplementary table 
lms$min <- NA # make min coordinate 
lms$max <- NA # make max coordinate 
for(i in 1:nrow(lms)){
  lg <- lms$lg[i]
  lms$min[i] <- min(lm$pos[lm$lg == lg])
  lms$max[i] <- max(lm$pos[lm$lg == lg])
}
lms2 <- lms[,c(1,2,16,17,11,3,4,5,7,8,9)] # choose columns 
# tidy 
lms2$m <- round(lms2$m,2) # round 
lms2$f <- round(lms2$f,2) # round 
lms2$av <- round(lms2$av,2) # round 
lms2$lenN <- round(lms2$lenN,2) # round 
lms2$ss[lms2$ss == "ss2.2"] <- 'ss2' # remake ss2.2 
lms2$ss <- gsub('ss','',lms2$ss) # remove ss prefix 
names(lms2) <- c('Linkage group', 'Super-Scaffold', 'Starting position', 
                 'Ending position', 'Length', 'Male', 'Female', 'Sex-averaged',
                 'SNPs','COs in male', 'COs in female') # rename columns 
write.csv(lms2, '../2.metadata/supp_table_lgs.csv',row.names = F) 
}


###################################
###################################
###################################
#### SUPP FIGURE X ################
###################################
###################################
###################################
{ # Marey maps 
  pdf('plots/SUPP_mareyMaps_all.pdf',18,18)
  par(mfrow=c(7,6),mar=c(2,2,2,2))
  i = 1
  for(lg in unique(lm$lg)[order(unique(lm$lg))]){
    ymax = max(max(lm$m[lm$lg == lg],na.rm=T),max(lm$f[lm$lg == lg],na.rm=T))
    plot(y=lm$m[lm$lg == lg], x=lm$pos[lm$lg == lg], ylim=c(0,ymax),
        pch=20, xlab='bp', ylab='cM', main=paste('LG', lg),
        col=viridis::viridis(3,alpha = .5)[1], cex=1.5)
    points(y=lm$f[lm$lg == lg], x=lm$pos[lm$lg == lg], 
         pch=20, col=viridis::viridis(3,alpha = .5)[2], cex=1.5)
    if(i %in% c(39)){
      legend('bottomright',pch=20,col=viridis::viridis(3)[c(1,2)],legend=c('Male','Female'),cex=1.5,
             bty='n')
    }
      i = i + 1 
  }
  dev.off()
}

###################################
###################################
###################################
#### SUPP FIGURE X ################
###################################
###################################
###################################
{ # PLOT ALL RECOMBINATION LANDSCAPES IN CH IN 100KB 
  pdf('plots/SUPP_pyrho_CH_all100kb.pdf',18,18)
  par(mfrow=c(7,6),mar=c(2,2,2,2))
  # make plot of alternating scaffolds 
  for(ss in unique(gsub('ss','Super-Scaffold_',lms$ss[order(lms$len,decreasing = T)]))){
    plot(d100$CH_cM[d100$ss == ss]~d100$start[d100$ss == ss], type='l', col='black',lwd=2,
         axes=T,xlab='',ylab='',main=paste0(ss), ylim=c(0,25))
  }
  dev.off()
}

###################################
###################################
###################################
#### SUPP FIGURE X ################
###################################
###################################
###################################
{ # COMPARISON OF ALL POPS WITH LM3
  pdf('plots/SUPP_corr_pyrho_LM3_allpops.pdf',8,8)
  d1m <- read.table('../1.data/3.windows/total_1Mb.table',h=T)
  # remove inf 
  d1m <- d1m[-which(d1m$ped == -Inf),]
  d1m$ss <- rename.ss(d1m$ss) # match ss names 
  # calculate correlations 
  lms$corCH <- NA
  lms$corCH13 <- NA
  lms$corCH13_2 <- NA
  lms$corPT <- NA
  lms$corGB <- NA
  # only for lgs with sufficient length (here > 5Mb)
  for( ss in lms$ss[lms$len > 5]){
    lms$corCH[lms$ss==ss] <- cor(d1m$CH_cM[d1m$ss==ss],d1m$ped[d1m$ss==ss],use='complete.obs')
    lms$corCH13[lms$ss==ss] <- cor(d1m$CH13_cM[d1m$ss==ss],d1m$ped[d1m$ss==ss],use='complete.obs')
    lms$corCH13_2[lms$ss==ss] <- cor(d1m$CH13_2_cM[d1m$ss==ss],d1m$ped[d1m$ss==ss],use='complete.obs')
    lms$corPT[lms$ss==ss] <- cor(d1m$PT_cM[d1m$ss==ss],d1m$ped[d1m$ss==ss],use='complete.obs')
    lms$corGB[lms$ss==ss] <- cor(d1m$GB_cM[d1m$ss==ss],d1m$ped[d1m$ss==ss],use='complete.obs')
    }
# plot
boxplot(lms$corCH,lms$corCH13,lms$corCH13_2,lms$corPT,lms$corGB,xaxt='none',ylab='correlation',
        main='per scaffold correlation between pyrho-LM3',xlab='population',pch=' ',
        col=colp,axes=F, ylim=c(0.4,1))
# add points to boxplots 
points(lms$corCH~jitter(rep(1,nrow(lms)),amount=.2),pch=20)
points(lms$corCH13~jitter(rep(2,nrow(lms)),amount=.2),pch=20)
points(lms$corCH13_2~jitter(rep(3,nrow(lms)),amount=.2),pch=20)
points(lms$corPT~jitter(rep(4,nrow(lms)),amount=.2),pch=20)
points(lms$corGB~jitter(rep(5,nrow(lms)),amount=.2),pch=20)
# add means 
points(colMeans(lms[18:ncol(lms)]),pch=18,col='grey90')
# axes 
axis(2, at = seq(0.4,1,0.2), labels=seq(0.4,1,0.2))
axis(1, at =seq(1,5), labels=c('CH','CH13','CH13_2','PT','GB'))
dev.off()
}


###################################
###################################
###################################
#### SUPP FIGURE X ################
###################################
###################################
###################################
{ # CORRELATION OF POPULATIONS PER SUPER-SCAFFOLD
  lens <- d100 %>% group_by(ss) %>% summarize(ss=unique(ss),len=max(end)-min(start))
  pdf('plots/SUPP_compare_pops_100kb.pdf',10,10)
  plot(0,pch='',xlim=c(0,40), ylim=c(0,1),axes=F,xlab='',ylab='')
  for(ss in gsub('ss','Super-Scaffold_',lms$ss[order(lms$len,decreasing = T)])){
    pt.test <- cor.test(d100$CH_cM[d100$ss==ss],d100$PT_cM[d100$ss==ss])
    x = which(lens$ss[order(lens$len,decreasing = T)] == ss)
    points(y=pt.test$estimate,x=x,col=colp[4],pch=16,cex=1.5)
    segments(x0=x,y0=pt.test$conf.int[[1]],y1=pt.test$conf.int[[2]],col=colp[4])
    gb.test <- cor.test(d100$CH_cM[d100$ss==ss],d100$GB_cM[d100$ss==ss])
    points(y=gb.test$estimate,x=x,col=colp[5],pch=16,cex=1.5)
    segments(x0=x,y0=gb.test$conf.int[[1]],y1=gb.test$conf.int[[2]],col=colp[5])
    ch13.test  <- cor.test(d100$CH_cM[d100$ss==ss],d100$CH13_cM[d100$ss==ss])
    points(y=ch13.test$estimate,x=x,col=colp[2],pch=16,cex=1.5)
    segments(x0=x,y0=ch13.test$conf.int[[1]],y1=ch13.test$conf.int[[2]],col=colp[2])
    ch13_2.test  <- cor.test(d100$CH_cM[d100$ss==ss],d100$CH13_2_cM[d100$ss==ss])
    points(y=ch13_2.test$estimate,x=x,col=colp[3],pch=16,cex=1.5)
    segments(x0=x,y0=ch13_2.test$conf.int[[1]],y1=ch13_2.test$conf.int[[2]],col=colp[3])
  }
  for(i in seq(0.5,39.5)){
    abline(v=i,lwd=.2,lty=1,col='grey80')
  }
  axis(1,at=c(0,10,20,30,39),labels=c(0,10,20,30,39),cex.axis=1.2)
  axis(2,at=seq(0,1,.25),labels=seq(0,1,.25),cex.axis=1.2)
  mtext('Linkage group (ordered by physical length)',1,2.2,cex=1.2)
  mtext('Correlation with CH (100kb)',2,2.2,cex=1.25)
  legend('bottomleft',pch=16,col=colp[-1],legend=c('CH-CH13','CH-CH13_2',
                                                   'CH-PT','CH-GB'),
         bty='n')
  dev.off()
}

###################################
###################################
###################################
#### SUPP FIGURE X ################
###################################
###################################
###################################
{# HOTSPOT OVERLAP 
  pdf('plots/SUPP_hotspot_overlap.pdf',8,8)
  # true overlap
  hs <- as.matrix(read.table('../1.data/2.pyrho/scaled_datasets/hotspot_sharing.table'))
  par(mar=c(1,1,1,1))
  # plot corrplot of hotspots 
  # local hotspots above diagonal
  # global hotspots below diagonal
  corrplot(hs,method='color',is.corr=F,
           col=COL1('Blues'),
           type='full', number.digits=1,
           addCoef.col = 'black',cl.pos = 'n',
           diag=F, tl.cex = 1.7, number.cex = 1.5,
           tl.col=colp, tl.srt=45,na.label = ' ', tl.pos = 'd')
  dev.off()
}

###################################
###################################
###################################
#### SUPP FIGURE X ################
###################################
###################################
###################################
{ # pruned points and GAM 
  pdf('plots/SUPP_GAMS.pdf',18,18)
  par(mfrow=c(7,6),mar=c(2,2,2,2))
  for(i in 1:39){
  or1 <- read.table(paste0('../1.data/1.LepMAP/orders/1st_ordering/Call2_tf1_maf05_miss05_1kb_pruned_f.call_order_',i,'_nogp.txt.perf'))
  or1$av <- (or1$V3 + or1$V4)/2
  lm1 <- lm(or1$av~or1$V2)
  if(lm1$coefficients[2] < 0){
    or1$av <- abs( or1$av - max(or1$av, na.rm=T) )
  }
  plot(or1$av~or1$V2,col=add.alpha('grey',.8),pch=20, xlab='bp',ylab='cM',
       main=paste('LG',i))
  #points(all1$av[all1$lg==1]~all1$pos[all1$lg==1],col=add.alpha('blue',.5),pch=20)
  points(lm$av[lm$lg==i]~lm$pos[lm$lg==i],col=add.alpha('red',.4),pch=20)
  }
  dev.off()
}

###################################
###################################
###################################
#### SUPP FIGURE X ################
###################################
###################################
###################################
{ # SMC histories 
path1 <- '../1.data/2.pyrho/smcpp/3.plots/'
smc.f <- list.files(path1, '*.csv')
smc <- data.frame(label=c(),x=c(),y=c())
for(f in smc.f){
  tmp <- read.csv(paste0(path1,f))[,c(1,2,3)]
  smc <- rbind(smc,tmp)
}
names(smc) <- c('Population','Generations','Ne')
smc$Population <- factor(smc$Population,levels=c('CH','CH13','CH13_2','PT','GB'))
p1 <- ggplot(data=smc,mapping=aes(x=Generations,y=Ne,col=Population)) + 
  geom_line(size=1.5) + 
  scale_color_manual(values=colp) +
  scale_x_continuous(trans='log10') + 
  theme_minimal()
pdf('plots/SUPP_smc_all.pdf',8,8)
print(p1)
dev.off()
}

###################################
###################################
###################################
#### SUPP FIGURE X ################
###################################
###################################
###################################
{ # Random hotspot sharing 
# establish random expectation of overlap
s1 <- seq(1,nrow(d1))
overlappt <- rep(0,10000)
overlapgb <- rep(0,10000)
for(i in 1:10000){
  overlappt[i] <- sum(sample(s1, length(hot1)) %in% sample(s1, length(pthot)))/length(hot1)
  overlapgb[i] <- sum(sample(s1, length(hot1)) %in% sample(s1, length(gbhot)))/length(hot1)
  print(i)
}
overlaps <- cbind.data.frame(overlappt, overlapgb)
# plot 
pdf('plots/SUPP_random_hotspot_sharing.pdf',6,8)
par(mfrow=c(1,2))
hist(overlaps$overlappt * 100, main='random sharing between PT and CH (10000 samples)',
     xlab='Overlap percentage') 
hist(overlaps$overlapgb, main='random sharing between GB and CH (10000 samples)',
     xlab='Overlap percentage')
dev.off()
}

###################################
###################################
###################################
#### SUPP FIGURE X ################
###################################
###################################
###################################
{# pedigree verification of Gini and pi
gal <- data.frame(ss = unique(d10$ss),
                  ginid10 = rep(0,41),
                  giniped = rep(0,41),
                  ginigen = rep(0,41),
                  pi = rep(0,41)
                  )
for(i in 1:nrow(gal)){
  gal$ginid10[i] <- Gini(d10$CH_cM[d10$ss==gal$ss[i]])
  gal$giniped[i] <- Gini(d1m$ped[d1m$ss==gal$ss[i]])
  gal$ginigen[i] <- Gini(d10$tss[d10$ss==gal$ss[i]])
  gal$pi[i] <- mean(d10$CH_pi[d10$ss==gal$ss[i]],na.rm=T)
}
pdf('plots/SUPP_lm3_gini.pdf',10,8)
par(mfrow=c(1,2))
plot(gal$giniped~gal$ginid10,pch=16,xlab='Gini coefficient (LD - 10kb)',
     ylab='Gini coefficient (linkage map - 1Mb)',cex=2,xlim=c(0.2,.9),ylim=c(.2,.9))
mtext('A',3,2)
plot(gal$giniped~gal$pi,pch=16,xlab='pi (10kb)',
     ylab='Gini coefficient (linkage map - 1Mb)',cex=2,ylim=c(0,1))
mtext('B',3,2)
dev.off()
}


dfd$gene_gini <- rep(0,nrow(dfd))
for(i in 1:nrow(dfd)){
 dfd$gene_gini[i] <- Gini(d10$tss[d10$ss==dfd$ss[i]])
  
}
plot(dfd$gene_gini~dfd$gini)
cor(dfd$gene_gini,dfd$gini)
