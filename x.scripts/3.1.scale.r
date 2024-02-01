### TOOLS
library(tidyverse)
source('functions.r')
sizes <- c('10kb','100kb','1Mb','5Mb') 
size.translator <- c(10000,100000,1000000,5000000)
options(scipen=10000)

for(size in sizes){
  
tmp <- size.translator[sizes==size]
### DATA
# windows of 5mb - used for scaffold length
scaffs <- read_delim('../1.data/3.windows/genome_5mb_windows.bed',col_names = c('ss','start','end'))
# lepmap results 
lm <- read_delim('../1.data/1.LepMAP/orders/fitted_values_gams.table')
# pyrho ch results - raw 
ch <- read_delim(paste0('../1.data/2.pyrho/estimates/CH_m200_total.theta.intersect.bed_',tmp,'.thetaxN'))
pt <- read_delim(paste0('../1.data/2.pyrho/estimates/PT_26_total.theta.intersect.bed_',tmp,'.thetaxN'))
gb <- read_delim(paste0('../1.data/2.pyrho/estimates/GB_26_total.theta.intersect.bed_',tmp,'.thetaxN'))
ch13 <- read_delim(paste0('../1.data/2.pyrho/estimates/CH13_26_total.theta.intersect.bed_',tmp,'.thetaxN'))
ch13_2 <- read_delim(paste0('../1.data/2.pyrho/estimates/CH13_2_26_total.theta.intersect.bed_',tmp,'.thetaxN'))

### TRANSFORMATIONS
# summarize length of scaffolds
scaffs$ss <- rename.ss(scaffs$ss)
ss.s <- scaffs %>% group_by(ss) %>% summarize(ss = unique(ss), len = max(end) )
# summarize length of LGs 
lm <- lm[! is.na(lm$av),]
lm.s <- lm %>% group_by(lg) %>% summarize(lg = unique(lg),ss=unique(ss),  lm.cM = max(av))
lm.s$len <- ss.s$len[match(lm.s$ss,ss.s$ss)] # add length
# transform pyrho rho to cM (Haldane's mapping function)
ch$hald <- hald(ch)
pt$hald <- hald(pt)
gb$hald <- hald(gb)
ch13$hald <- hald(ch13)
ch13_2$hald <- hald(ch13_2)
# fix ss2_2
ch$ss[ch$ss == 'Super-Scaffold_2' & ch$end <= 28692481] <- 'Super-Scaffold_2.2'
pt$ss[pt$ss == 'Super-Scaffold_2' & pt$end <= 28692481] <- 'Super-Scaffold_2.2'
gb$ss[gb$ss == 'Super-Scaffold_2' & gb$end <= 28692481] <- 'Super-Scaffold_2.2'
ch13$ss[ch13$ss == 'Super-Scaffold_2' & ch13$end <= 28692481] <- 'Super-Scaffold_2.2'
ch13_2$ss[ch13_2$ss == 'Super-Scaffold_2' & ch13_2$end <= 28692481] <- 'Super-Scaffold_2.2'
# summarize pyrho and merge with LM
ch.s <- ch %>% group_by(ss) %>% summarize(ss = unique(ss), ch.cM = sum(hald,na.rm=T))
ch.s$ss <- rename.ss(ch.s$ss)
pt.s <- pt %>% group_by(ss) %>% summarize(ss = unique(ss), pt.cM = sum(hald,na.rm=T))
pt.s$ss <- rename.ss(pt.s$ss)
gb.s <- gb %>% group_by(ss) %>% summarize(ss = unique(ss), gb.cM = sum(hald,na.rm=T))
gb.s$ss <- rename.ss(gb.s$ss)
ch13.s <- ch13 %>% group_by(ss) %>% summarize(ss = unique(ss), ch13.cM = sum(hald,na.rm=T))
ch13.s$ss <- rename.ss(ch13.s$ss)
ch13_2.s <- ch13_2 %>% group_by(ss) %>% summarize(ss = unique(ss), ch13_2.cM = sum(hald,na.rm=T))
ch13_2.s$ss <- rename.ss(ch13_2.s$ss)
py.all <- plyr::join_all(list(ch.s,pt.s,gb.s,ch13.s,ch13_2.s), by='ss',type='left')
d <- left_join(lm.s,py.all, by='ss')

# scale pyrho
scales <- data.frame('pop'=c('ch','pt','gb','ch13_1','ch13_2'),
                     'scaling_factor'=rep(0,5))
# for each pop get scaling factor but make sure you dont use ss3 and 49 (cause its 2 ss in LM)
for(i in 5:9){
  scales$scaling_factor[i-4] <- sum(d$lm.cM) / sum(d[,i])
}
ch$cM.scaled <- ch$hald * scales$scaling_factor[1]
pt$cM.scaled <- pt$hald * scales$scaling_factor[2]
gb$cM.scaled <- gb$hald * scales$scaling_factor[3]
ch13$cM.scaled <- ch13$hald * scales$scaling_factor[4]
ch13_2$cM.scaled <- ch13_2$hald * scales$scaling_factor[5]

### PLOT
pdf(paste0('plots/',size,'scaling_factor.jpeg'),height=8,width=8)
boxplot(cbind(d$lm.cM/d$ch.cM,d$lm.cM/d$pt.cM,d$lm.cM/d$gb.cM,d$lm.cM/d$ch13.cM,d$lm.cM/d$ch13_2.cM),
        ylim=c(0,5), ylab='LepMAP3 cM / pyrho cM', main=paste0(size,' scaling factor per Super Scaffold'),xaxt='n')
points(y=scales$scaling_factor,x=1:5,pch=18,col='red')
axis(1,at=1:5,labels=c('CH','PT','GB','CH13','CH13_2'))
dev.off()

### WRITE 
write_delim(ch,paste0('../1.data/2.pyrho/scaled_datasets/CH_',size,'_scaled.txt'),quote = 'none',eol='\n')
write_delim(pt,paste0('../1.data/2.pyrho/scaled_datasets/PT_',size,'_scaled.txt'),quote = 'none',eol='\n')
write_delim(gb,paste0('../1.data/2.pyrho/scaled_datasets/GB_',size,'_scaled.txt'),quote = 'none',eol='\n')
write_delim(ch13,paste0('../1.data/2.pyrho/scaled_datasets/CH13_',size,'_scaled.txt'),quote = 'none',eol='\n')
write_delim(ch13_2,paste0('../1.data/2.pyrho/scaled_datasets/CH13_2_',size,'_scaled.txt'),quote = 'none',eol='\n')

write.table(scales,paste0('../1.data/2.pyrho/scaled_datasets/',size,'_scaling_factors.txt'),quote=F,row.names=F)

if(size=='10kb'){write.table(d,'../1.data/PerLG_lm_py_unscaled.txt',quote=F,row.names=F)}

}
