### TOOLS
library(tidyverse)
source('functions.r')

### DATA
# windows of 5mb - used for scaffold length
scaffs <- read_delim('../1.data/3.windows/genome_5mb_windows.bed',col_names = c('ss','start','end'))
# lepmap results 
lm <- read_delim('../1.data/1.LepMAP/orders/best_likelihood_pruned_orders.table')
# pyrho ch results - raw 
ch <- read_delim('../1.data/2.pyrho/estimates/CH_m200_total.theta',col_names = c('ss','a','b','rho'))
pt <- read_delim('../1.data/2.pyrho/estimates/PT_26_total.theta',col_names = c('ss','a','b','rho'))
gb <- read_delim('../1.data/2.pyrho/estimates/gb_26_total.theta',col_names = c('ss','a','b','rho'))
ch13_1 <- read_delim('../1.data/2.pyrho/estimates/CH13_26_total.theta',col_names = c('ss','a','b','rho'))
ch13_2 <- read_delim('../1.data/2.pyrho/estimates/CH13_2_26_total.theta',col_names = c('ss','a','b','rho'))

### TRANSFORMATIONS
# summarize length of scaffolds
scaffs$ss <- rename.ss(scaffs$ss)
ss.s <- scaffs %>% group_by(ss) %>% summarize(ss = unique(ss), len = max(end) )
# summarize length of LGs 
lm <- lm[! is.na(lm$av),]
lm.s <- lm %>% group_by(lg) %>% reframe(lg = unique(lg),ss=unique(ss),  lm.cM = max(av))
lm.s$len <- ss.s$len[match(lm.s$ss,ss.s$ss)] # add length
# Super-Scaffold_2 correction (first LG is all but start then only start - start= 1-28161559)
lm.s$len[lm.s$lg==1] <- lm.s$len[lm.s$lg==1] - min(lm$pos[lm$lg==1])
lm.s$len[lm.s$lg!=1 & lm.s$ss=='ss2'] <-  min(lm$pos[lm$lg==1])
lm.s$ss[lm.s$lg!=1 & lm.s$ss=='ss2'] <- 'ss2.2'
# transform pyrho rho to cM (Haldane's mapping function)
ch$hald <- hald(ch)
pt$hald <- hald(pt)
gb$hald <- hald(gb)
ch13_1$hald <- hald(ch13_1)
ch13_2$hald <- hald(ch13_2)
# fix ss2_2
ch$ss[ch$ss == 'Super-Scaffold_2' & ch$b <= 28692481] <- 'Super-Scaffold_2.2'
pt$ss[pt$ss == 'Super-Scaffold_2' & pt$b <= 28692481] <- 'Super-Scaffold_2.2'
gb$ss[gb$ss == 'Super-Scaffold_2' & gb$b <= 28692481] <- 'Super-Scaffold_2.2'
ch13_1$ss[ch13_1$ss == 'Super-Scaffold_2' & ch13_1$b <= 28692481] <- 'Super-Scaffold_2.2'
ch13_2$ss[ch13_2$ss == 'Super-Scaffold_2' & ch13_2$b <= 28692481] <- 'Super-Scaffold_2.2'
# summarize pyrho and merge with LM
ch.s <- ch %>% group_by(ss) %>% summarize(ss = unique(ss), ch.cM = sum(hald))
ch.s$ss <- rename.ss(ch.s$ss)
pt.s <- pt %>% group_by(ss) %>% summarize(ss = unique(ss), pt.cM = sum(hald))
pt.s$ss <- rename.ss(pt.s$ss)
gb.s <- gb %>% group_by(ss) %>% summarize(ss = unique(ss), gb.cM = sum(hald))
gb.s$ss <- rename.ss(gb.s$ss)
ch13_1.s <- ch13_1 %>% group_by(ss) %>% summarize(ss = unique(ss), ch13_1.cM = sum(hald))
ch13_1.s$ss <- rename.ss(ch13_1.s$ss)
ch13_2.s <- ch13_2 %>% group_by(ss) %>% summarize(ss = unique(ss), ch13_2.cM = sum(hald))
ch13_2.s$ss <- rename.ss(ch13_2.s$ss)
py.all <- plyr::join_all(list(ch.s,pt.s,gb.s,ch13_1.s,ch13_2.s), by='ss',type='left')
d <- left_join(lm.s,py.all, by='ss')

# scale pyrho
scales <- data.frame('pop'=c('ch','pt','gb','ch13_1','ch13_2'),
                     'scaling_factor'=rep(0,5))
for(i in 5:9){
  scales$scaling_factor[i-4] <- sum(d$lm.cM) / sum(d[,i])
}
ch$cM.scaled <- ch$hald * scales$scaling_factor[1]
pt$cM.scaled <- pt$hald * scales$scaling_factor[2]
gb$cM.scaled <- gb$hald * scales$scaling_factor[3]
ch13_1$cM.scaled <- ch13_1$hald * scales$scaling_factor[4]
ch13_2$cM.scaled <- ch13_2$hald * scales$scaling_factor[5]

### WRITE 
write_delim(ch,'../1.data/2.pyrho/scaled_datasets/ch_scaled_all.txt',quote = 'none',eol='\n')
write_delim(pt,'../1.data/2.pyrho/scaled_datasets/pt_scaled_all.txt',quote = 'none',eol='\n')
write_delim(gb,'../1.data/2.pyrho/scaled_datasets/gb_scaled_all.txt',quote = 'none',eol='\n')
write_delim(ch13_1,'../1.data/2.pyrho/scaled_datasets/ch13_1_scaled_all.txt',quote = 'none',eol='\n')
write_delim(ch13_2,'../1.data/2.pyrho/scaled_datasets/ch13_2_scaled_all.txt',quote = 'none',eol='\n')

write.table(scales,'../1.data/2.pyrho/scaled_datasets/scaling_factors.txt',quote=F,row.names=F)
write.table(d,'../1.data/PerLG_lm_py_unscaled.txt',quote=F,row.names=F)
