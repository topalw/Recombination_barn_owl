##################### PREAMPLE - THINGS WE USE IN HERE ########################
library(readr)
library(dplyr)
source('functions.r')
shit <- c('ss24','ss25','ss100000100078','ss200000178', 'ss100000100064', 'ss4')
p <- c('orange','blue')

run1 <- read_delim('../1.data/1.LepMAP/orders/likelihoods/run1.likelihoods',col_names = c('lg','L'))
run2 <- read_delim('../1.data/1.LepMAP/orders/likelihoods/run2.likelihoods',col_names = c('lg','L'))
run3 <- read_delim('../1.data/1.LepMAP/orders/likelihoods/run3.likelihoods',col_names = c('lg','L'))

ssmap <- read_delim('lg_ss_map.txt')
ssmap$ss <- rename.ss(ssmap$ss)
shit.lg <- ssmap[match(shit,ssmap$ss),1]

run1 <- run1[! run1$lg %in% shit.lg$lg,]
run2 <- run2[! run2$lg %in% shit.lg$lg,]
run3 <- run3[! run3$lg %in% shit.lg$lg,]

full <- full_join(run1,run2,by='lg')
full <- full_join(full,run3,by='lg')
names(full) <- c('lg','l1','l2','l3')
full$best <- 0
for(i in 1:nrow(full)){
full$best[i] <-  which(full[i,c(2,3,4)] == max(full[i,c(2,3,4)]))
}

files <- c('../1.data/1.LepMAP/orders/1st_ordering/full_orders_fixed349_nosmall_ascendingAll_mor.txt',
           '../1.data/1.LepMAP/orders/1st_ordering/full_orders_fixed349_nosmall_ascendingAll_mor.txt',
           '../1.data/1.LepMAP/orders/1st_ordering/full_orders_fixed349_nosmall_ascendingAll_mor.txt')
d1 <- read_delim(files[1])
d1$run <- 1
d2 <- read_delim(files[2])
d2$run <- 2
d3 <- read_delim(files[3])
d3$run <- 3
dall <- rbind(d1,d2,d3)

tmp <- left_join(d1,d2,by=c('ss','pos'))
tmp <- left_join(tmp,d3,by=c('ss','pos'))
tmp <- tmp[! is.na(tmp$av.x) | ! is.na(tmp$av.y) | ! is.na(tmp$av),]
# plot correlations 
pdf('plots/Marey_comparing_runs.pdf',width=12,height=12)
par(mfrow=c(3,3))
plot(0,pch='',axes=F,xlab='',ylab='')
plot(tmp$m.x~tmp$m.y,pch=20,col=add.alpha('black',.2),xlab='',ylim=c(0,80),ylab='')
plot(tmp$m.x~tmp$m,pch=20,col=add.alpha('black',.2),xlab='',ylim=c(0,80),ylab='')
plot(tmp$f.y~tmp$f.x,pch=20,col=add.alpha('black',.2),xlab='',ylim=c(0,80),ylab='')
plot(0,pch='',axes=F,xlab='',ylab='')
plot(tmp$m.y~tmp$m,pch=20,col=add.alpha('black',.2),xlab='',ylim=c(0,80),ylab='')
plot(tmp$f~tmp$f.x,pch=20,col=add.alpha('black',.2),xlab='',ylim=c(0,80),ylab='')
plot(tmp$f~tmp$f.y,pch=20,col=add.alpha('black',.2),xlab='',ylim=c(0,80),ylab='')
plot(0,pch='',axes=F,xlab='',ylab='')
dev.off()

best.orders <- data.frame('lg'=factor(),'ss'=character(),'pos'=integer(),
                          'm'=double(),'f'=double(),'av'=double())
for(i in 1:nrow(full)){
  lg <- full$lg[i]
  best <- full$best[i]
  best.orders <- rbind(best.orders,dall[dall$run==best & dall$lg==lg,c(5,1,2,3,4,6)])
}
unique(dall$lg)

write.table(best.orders,'../1.data/1.LepMAP/orders/best_likelihood_pruned_orders.table',quote=F,row.names = F)

