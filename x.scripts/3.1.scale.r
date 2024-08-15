### TOOLS
library(tidyverse)
source('functions.r')
sizes <- c('10kb','100kb','1Mb','5Mb') 
size.translator <- c(10000,100000,1000000,5000000)
options(scipen=10000)


# lepMap results 
lengths <- read_delim('../1.data/lengths.txt')
# lepmap results 
lm <- read_delim('../1.data/1.LepMAP/R.AllSNPs_physical/linkage_map_perLG.txt')
lm$ss[lm$lg==24] <- 'Super-Scaffold_2.2'
for(size in sizes){
  
tmp <- size.translator[sizes==size]
### DATA
# pyrho ch results - raw 
ch <- read_delim(paste0('../1.data/2.pyrho/estimates/CH_total.theta.intersect.bed_',tmp,'.thetaxN'))
pt <- read_delim(paste0('../1.data/2.pyrho/estimates/PT_total.theta.intersect.bed_',tmp,'.thetaxN'))
gb <- read_delim(paste0('../1.data/2.pyrho/estimates/GB_total.theta.intersect.bed_',tmp,'.thetaxN'))
ch13 <- read_delim(paste0('../1.data/2.pyrho/estimates/CH13_total.theta.intersect.bed_',tmp,'.thetaxN'))
ch13_2 <- read_delim(paste0('../1.data/2.pyrho/estimates/CH13_2_total.theta.intersect.bed_',tmp,'.thetaxN'))
ch13_3 <- read_delim(paste0('../1.data/2.pyrho/estimates/CH13_3_total.theta.intersect.bed_',tmp,'.thetaxN'))
ch13_4 <- read_delim(paste0('../1.data/2.pyrho/estimates/CH13_4_total.theta.intersect.bed_',tmp,'.thetaxN'))
ch13_5 <- read_delim(paste0('../1.data/2.pyrho/estimates/CH13_5_total.theta.intersect.bed_',tmp,'.thetaxN'))

### TRANSFORMATIONS
# dirty scale Z to make it in haldane's function 
ch$thetaxN[ch$ss %in% c('Super-Scaffold_13','Super-Scaffold_49')] <- ch$thetaxN[ch$ss %in% c('Super-Scaffold_13','Super-Scaffold_49')] * .1 
pt$thetaxN[pt$ss %in% c('Super-Scaffold_13','Super-Scaffold_49')] <- pt$thetaxN[pt$ss %in% c('Super-Scaffold_13','Super-Scaffold_49')] * .1
gb$thetaxN[pt$ss %in% c('Super-Scaffold_13','Super-Scaffold_49')] <- gb$thetaxN[gb$ss %in% c('Super-Scaffold_13','Super-Scaffold_49')] * .1 
# transform pyrho rho to cM (Haldane's mapping function)
{
ch$hald <- hald(ch)
pt$hald <- hald(pt)
gb$hald <- hald(gb)
ch13$hald <- hald(ch13)
ch13_2$hald <- hald(ch13_2)
ch13_3$hald <- hald(ch13_3)
ch13_4$hald <- hald(ch13_4)
ch13_5$hald <- hald(ch13_5)
}
{
# fix ss2_2
ch$ss[ch$ss == 'Super-Scaffold_2' & ch$end <= len2.2] <- 'Super-Scaffold_2.2'
pt$ss[pt$ss == 'Super-Scaffold_2' & pt$end <= len2.2] <- 'Super-Scaffold_2.2'
gb$ss[gb$ss == 'Super-Scaffold_2' & gb$end <= len2.2] <- 'Super-Scaffold_2.2'
ch13$ss[ch13$ss == 'Super-Scaffold_2' & ch13$end <= len2.2] <- 'Super-Scaffold_2.2'
ch13_2$ss[ch13_2$ss == 'Super-Scaffold_2' & ch13_2$end <= len2.2] <- 'Super-Scaffold_2.2'
ch13_3$ss[ch13_3$ss == 'Super-Scaffold_2' & ch13_3$end <= len2.2] <- 'Super-Scaffold_2.2'
ch13_4$ss[ch13_4$ss == 'Super-Scaffold_2' & ch13_4$end <= len2.2] <- 'Super-Scaffold_2.2'
ch13_5$ss[ch13_5$ss == 'Super-Scaffold_2' & ch13_5$end <= len2.2] <- 'Super-Scaffold_2.2'
}
# keep only lm scaffolds  
{
  ch <- ch[ch$ss %in% lm$ss,]
  pt <- pt[pt$ss %in% lm$ss,]
  gb <- gb[gb$ss %in% lm$ss,]
  ch13 <- ch13[ch13$ss %in% lm$ss,]
  ch13_2 <- ch13_2[ch13_2$ss %in% lm$ss,]
  ch13_3 <- ch13_3[ch13_3$ss %in% lm$ss,]
  ch13_4 <- ch13_4[ch13_4$ss %in% lm$ss,]
  ch13_5 <- ch13_5[ch13_5$ss %in% lm$ss,]
}
{
# summarize pyrho and merge with LM
ch.s <- ch %>% group_by(ss) %>% summarize(ss = unique(ss), ch.cM = sum(hald,na.rm=T))
#ch.s$ss <- rename.ss(ch.s$ss)
pt.s <- pt %>% group_by(ss) %>% summarize(ss = unique(ss), pt.cM = sum(hald,na.rm=T))
#pt.s$ss <- rename.ss(pt.s$ss)
gb.s <- gb %>% group_by(ss) %>% summarize(ss = unique(ss), gb.cM = sum(hald,na.rm=T))
#gb.s$ss <- rename.ss(gb.s$ss)
ch13.s <- ch13 %>% group_by(ss) %>% summarize(ss = unique(ss), ch13.cM = sum(hald,na.rm=T))
#ch13.s$ss <- rename.ss(ch13.s$ss)
ch13_2.s <- ch13_2 %>% group_by(ss) %>% summarize(ss = unique(ss), ch13_2.cM = sum(hald,na.rm=T))
#ch13_2.s$ss <- rename.ss(ch13_2.s$ss)
ch13_3.s <- ch13_3 %>% group_by(ss) %>% summarize(ss = unique(ss), ch13_3.cM = sum(hald,na.rm=T))
ch13_4.s <- ch13_4 %>% group_by(ss) %>% summarize(ss = unique(ss), ch13_4.cM = sum(hald,na.rm=T))
ch13_5.s <- ch13_5 %>% group_by(ss) %>% summarize(ss = unique(ss), ch13_5.cM = sum(hald,na.rm=T))

py.all <- plyr::join_all(list(ch.s,pt.s,gb.s,ch13.s,ch13_2.s,ch13_3.s,ch13_4.s,ch13_5.s), by='ss',type='left')

d <- left_join(lm,py.all, by='ss')
}
# scale pyrho
scales <- data.frame('pop'=c(),
                     'ss'=c(),
                     'scaling_factor'=c())
pops <- c('ch','pt','gb','ch13','ch13_2','ch13_3','ch13_4','ch13_5')
# work through ss and scale each with a unique scaling factor
for(i in 1:nrow(d)){
  lepmap.length <- d$av[i]
  ss <- d$ss[i]
  # get scaling factor for each ss per pop
  for(j in 11:18){
    sf <- as.numeric(lepmap.length / d[i,j])
    scales <- rbind.data.frame(scales,
                               data.frame('pop'=pops[j-10],
                                          'ss'=ss,
                                          'scaling_factor'=sf))
    
  }
}

# scale per ss
for(ss in unique(scales$ss)){
  tmp <- scales[scales$ss==ss,]
  ch$hald[ch$ss==ss] <- ch$hald[ch$ss==ss] * tmp$scaling_factor[1]
  pt$hald[pt$ss==ss] <- pt$hald[pt$ss==ss] * tmp$scaling_factor[2]
  gb$hald[gb$ss==ss] <- gb$hald[gb$ss==ss] * tmp$scaling_factor[3]
  ch13$hald[ch13$ss==ss] <- ch13$hald[ch13$ss==ss] * tmp$scaling_factor[4]
  ch13_2$hald[ch13_2$ss==ss] <- ch13_2$hald[ch13_2$ss==ss] * tmp$scaling_factor[5]
  ch13_3$hald[ch13_3$ss==ss] <- ch13_3$hald[ch13_3$ss==ss] * tmp$scaling_factor[6]
  ch13_4$hald[ch13_4$ss==ss] <- ch13_4$hald[ch13_4$ss==ss] * tmp$scaling_factor[7]
  ch13_5$hald[ch13_5$ss==ss] <- ch13_5$hald[ch13_5$ss==ss] * tmp$scaling_factor[8]
  
}

### PLOT
tiff(paste0('plots/',size,'scaling_factor.jpeg'),920,920,res=100)
par(mfrow=c(1,2))
boxplot(scales$scaling_factor~scales$pop,
        ylim=c(0,8), ylab='LepMAP3 cM / pyrho cM', main=paste0(size,' scaling factor per Super Scaffold'),xaxt='n')
#points(y=scales$scaling_factor,x=1:8,pch=18,col='red')
axis(1,at=1:8,labels=c('CH','PT','GB','CH13','CH13_2','CH13_3','CH13_4','CH13_5'))
abline(h=1,lty=2)
boxplot(scales$scaling_factor~scales$ss,
        ylim=c(0,10), ylab='LepMAP3 cM / pyrho cM', main=paste0(size,' scaling factor per Super Scaffold'),xaxt='n',
        xlab='Super-Scaffold')
axis(1,at=seq(1:length(unique(scales$ss))), labels=gsub('Super-Scaffold_','',unique(scales$ss)),las=2,cex.axis=.7)
abline(h=1,lty=2)
dev.off()

# subset only relevant
cols=c(1,2,3,5)

### WRITE 
write_delim(ch[,cols],paste0('../1.data/2.pyrho/scaled_datasets/CH_',size,'_scaled.txt'),quote = 'none',eol='\n')
write_delim(pt[,cols],paste0('../1.data/2.pyrho/scaled_datasets/PT_',size,'_scaled.txt'),quote = 'none',eol='\n')
write_delim(gb[,cols],paste0('../1.data/2.pyrho/scaled_datasets/GB_',size,'_scaled.txt'),quote = 'none',eol='\n')
write_delim(ch13[,cols],paste0('../1.data/2.pyrho/scaled_datasets/CH13_',size,'_scaled.txt'),quote = 'none',eol='\n')
write_delim(ch13_2[,cols],paste0('../1.data/2.pyrho/scaled_datasets/CH13_2_',size,'_scaled.txt'),quote = 'none',eol='\n')
write_delim(ch13_3[,cols],paste0('../1.data/2.pyrho/scaled_datasets/CH13_3_',size,'_scaled.txt'),quote = 'none',eol='\n')
write_delim(ch13_4[,cols],paste0('../1.data/2.pyrho/scaled_datasets/CH13_4_',size,'_scaled.txt'),quote = 'none',eol='\n')
write_delim(ch13_5[,cols],paste0('../1.data/2.pyrho/scaled_datasets/CH13_5_',size,'_scaled.txt'),quote = 'none',eol='\n')

write.table(scales,paste0('../1.data/2.pyrho/scaled_datasets/',size,'_scaling_factors.txt'),quote=F,row.names=F)

#if(size=='10kb'){write.table(d,'../1.data/PerLG_lm_py_unscaled.txt',quote=F,row.names=F)}

}
