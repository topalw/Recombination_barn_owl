setwd("C:/Users/topalw/Desktop/PhD/Analyses/Chapter_A/x.scripts/")
subset <- read.table('../1.data/1.LepMAP/R.LM3_physical/physical_orders_all.txt',h=T)
renamer <- read.table('../1.data/1.LepMAP/R.LM3_physical/linkage_map_perLG.txt',h=T)
source('functions.r')


#### PRELUDE #### 

# add ss3 & ss49 
# proof of end3_start3_start49_end49 ordering 
jpeg('plots/SUPP_1ss3_49.jpg',8,8,units = 'in',res=200,quality=200)
tmp <- read.table('../1.data/1.LepMAP/orders/1st_ordering/Call2_tf1_maf05_miss05_1kb_pruned_f.call_order_20_nogp.txt.perf')
plot(tmp$V3~tmp$V2,col=c('red','black')[as.numeric(factor(tmp$V1,levels=c('Super-Scaffold_3','Super-Scaffold_49')))],
     xlab='position',ylab='cM',pch=20)
legend('topright',legend=c('Super-Scaffold_3','Super-Scaffold_49'),col=c('red','black'),pch=20)
dev.off()

#################

#### FIX UNIQ ORDERS #####

#pdf('plots/R.comparing_orders.pdf',12,12)
jpeg('plots/SUPP_ALL_LEPMAP3.jpg',14,18,units = 'in',res=600,quality=800)

par(mfrow=c(7,6))
for(lg in c(1:19,21:39)){
  ss <- renamer$ss[renamer$lg==lg]
  d <- read.table(paste0('../1.data/1.LepMAP/R.AllSNPs_physical/perf_orders/',ss,'_f.call_lod5_map.phys_order.txt.perf'))
  names(d) <- c('ss','pos','m','f')
  d$pos <- gsub('\\*','',d$pos)
  d2 <- d 
  # prune ends 
  d2$m <- prune.ends(d2$m,100,2)
  d2$f <- prune.ends(d2$f,100,2)
  # make sex-averaged
  d2$av <- (d2$m + d2$f) /2 
  # add LG
  d2$lg <- lg 
  if(lg==1){
    dall <- d2[,c(6,1,2,3,4,5)]
  } else {
    dall <- rbind.data.frame(dall,d2[,c(6,1,2,3,4,5)])
  }
  # plot new order 
  plot(d2$av~d2$pos,pch=16,cex=.8,main=paste0('LG',lg),
       xlab='physical position (bp)',ylab='Genetic position (cM)',sub=ss)
  points(subset$av[subset$lg==lg]~subset$pos[subset$lg==lg],
         col=add.alpha('red',.2),cex=.8,pch=20)
  abline(v=min(subset$pos[subset$lg==lg]),lty=2)
  abline(v=max(subset$pos[subset$lg==lg]),lty=2)
  legend('topleft',pch=c(20,20),col=c('black','red'),
         legend=c('LOD5 SNPs - physical order ', 'LOD15 SNPs - physical order'),
         bty='n')
}


###### UNIQUE CASES ########


########### SS3 and SS49 = LG20 MERGE ##########

# fix 
d <- read.table('../1.data/1.LepMAP/R.AllSNPs_physical/perf_orders/merge_3_49_f.call_lod5_map.phys_order.txt.perf')
names(d)<- c('ss','pos','m','f')
#plot(d$av~d$pos,col=c('red','black')[as.numeric(factor(d$ss,levels=c('Super-Scaffold_3','Super-Scaffold_49')))])
# step 1 set end of 49 to be 0 cM 
d$f[d$ss=='Super-Scaffold_49'] <- abs(d$f[d$ss=='Super-Scaffold_49']  - max(d$f[d$ss=='Super-Scaffold_49'] ))
d$m[d$ss=='Super-Scaffold_49'] <- abs(d$m[d$ss=='Super-Scaffold_49']  - max(d$m[d$ss=='Super-Scaffold_49'] ))
# step 2 cM of 3 = 3 + max 49 
d$f[d$ss=='Super-Scaffold_3'] <- d$f[d$ss=='Super-Scaffold_3'] + max(d$f[d$ss=='Super-Scaffold_49'])
d$m[d$ss=='Super-Scaffold_3'] <- d$m[d$ss=='Super-Scaffold_3'] + max(d$m[d$ss=='Super-Scaffold_49'])
# step 3 positions of 49 come before 3 so 3 = 3 + 49 and 49 is inverted 
d$newpos <- d$pos
d$newpos[d$ss=='Super-Scaffold_3'] <- d$newpos[d$ss=='Super-Scaffold_3'] + 4323735
d$newpos[d$ss=='Super-Scaffold_49'] <- abs(d$newpos[d$ss=='Super-Scaffold_49'] - 4323735)
d <- d[order(d$newpos,decreasing = F),]
d$m <- prune.ends(d$m,100,2)
d$f <- prune.ends(d$f,100,2)
d$av <- (d$m + d$f) / 2
plot(d$av~d$newpos,col=c('purple','black')[as.numeric(factor(d$ss,levels=c('Super-Scaffold_3','Super-Scaffold_49')))],
     pch=20,xlab='Merged physical position', ylab='Genetic position (cM)', main='LG20')
legend('topleft',pch=20,col=c('black','purple'),legend=c('Super-Scaffold_49','Super-Scaffold_3'),
       bty='n')

dev.off()

dall$newpos <- dall$pos 
d$lg <- rep(20,nrow(d))

dall <- rbind.data.frame(dall,d[,c(7,1,2,3,4,6,5)])


########### SS13 and SS42 = Z MERGE ##########
len42 <- 41426968
len13 <- 48838493
PAR <- data.frame('start'=c(1),'end'=c(4415000))
# Z 
dm <- read.table('../1.data/1.LepMAP/R.AllSNPs_physical/perf_orders/Z2_f.call_lod5_map.phys_order_Male.txt.perf')
df <- read.table('../1.data/1.LepMAP/R.AllSNPs_physical/perf_orders/Z2_f.call_lod5_map.phys_order_Female.txt.perf')
names(dm)<- c('ss','pos','m','f')
names(df)<- c('ss','pos','m','f')
d <- dplyr::full_join(dm[,c(1,2,3)],df[,c(1,2,4)],by=c('ss','pos'))
d$pos <- gsub('\\*','',d$pos)
d$pos <- as.numeric(d$pos)
d$newpos <- d$pos
# 42 inverted 
d$newpos[d$ss=='Super-Scaffold_42'] <- abs(d$newpos[d$ss=='Super-Scaffold_42'] - len42 )
# 42 cM inverted 
d$m[d$ss=='Super-Scaffold_42'] <- abs(d$m[d$ss=='Super-Scaffold_42'] - max(d$m[d$ss=='Super-Scaffold_42'],na.rm=T)) + max(d$m[d$ss=='Super-Scaffold_13'],na.rm=T)
d$f[d$ss=='Super-Scaffold_42'] <- abs(d$f[d$ss=='Super-Scaffold_42'] - max(d$f[d$ss=='Super-Scaffold_42'],na.rm=T)) + max(d$f[d$ss=='Super-Scaffold_13'],na.rm=T)
# 42 at end 
d$newpos[d$ss=='Super-Scaffold_42'] <- d$newpos[d$ss=='Super-Scaffold_42'] + len13
d <- d[order(d$newpos,decreasing = F),]
d$m <- prune.ends(d$m,100,2)
# translate PAR to newpos 
newposPAR <- d$newpos[which(d$ss=='Super-Scaffold_42' & d$pos < PAR$end[1])]
# no female recombination anywhere except PAR 
d$f[! d$newpos %in% newposPAR] <- 0.0
d$f[ d$newpos %in% newposPAR] <- d$f[ d$newpos %in% newposPAR]  - min(d$f[ d$newpos %in% newposPAR] ,na.rm=T)
#d$f <- prune.ends(d$f,100,2)
jpeg('plots/SUPP_Z.jpg',8,8,units = 'in',res=200,quality=400)
plot(d$m~d$newpos,col=c('red','black')[as.numeric(factor(d$ss,levels=c('Super-Scaffold_13','Super-Scaffold_42')))],
     pch=16,ylim=c(0,max(c(d$m,d$f),na.rm=T)))
points(d$f~d$newpos,col=c('red','black')[as.numeric(factor(d$ss,levels=c('Super-Scaffold_13','Super-Scaffold_42')))],
       pch=17)
rect(xleft=min(newposPAR),ybottom=0,xright=max(newposPAR),ytop=100,col=add.alpha('blue',.2))

legend('topleft',pch=c(16,16),legend=c('Male Super-Scaffold_13','Male Super-Scaffold_42'),
       col=c('red','black'))
dev.off()

d$av <- (d$m + d$f)/ 2 
d$lg <- rep(40,nrow(d))
dall <- rbind.data.frame(dall,d[,c(7,1,2,3,4,6,5)])

# # write output
#options(scipen=9999)

# dall$ss[dall$lg==24] <- gsub('_2','_2.2',dall$ss[dall$lg==24])
#write.table(dall,'../1.data/1.LepMAP/R.AllSNPs_physical/physical_orders_full.txt',quote=F,row.names=F,col.names=T)
#tr <- dall[!is.na(dall$av),c(1,2,3,6)]
#tr$ss <- gsub('_2.2','_2',tr$ss)
#write.table(tr,'C:/Users/topalw/Desktop/linkage_map_snps.txt',quote=F,row.names=F)
