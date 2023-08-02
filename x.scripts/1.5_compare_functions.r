library(dplyr)
source('functions.r')

# read and format 
names1 <- c('ss','pos','m','f')
kos <- read.table('../1.data/1.LepMAP/orders/other_functions/full_orders_fixed349_nosmall_ascendingAll_kos.txt',h=T)[,c(1,2,6)]
mor <- read.table('../1.data/1.LepMAP/orders/other_functions/full_orders_fixed349_nosmall_ascendingAll_mor.txt',h=T)[,c(1,2,6)]
hal <- read.table('../1.data/1.LepMAP/orders/best_likelihood_pruned_orders.table',h=T)[,c(2,3,6)]

all <- left_join(hal,kos,by=c('ss','pos'))
names(all)[c(3,4)] <- c('hal','kos')
all <- left_join(all, mor,by=c('ss','pos'))
names(all)[5] <- c('mor')

pdf('plots/Comparing_functions.pdf',width=10,height=6)
par(mfrow=c(1,3))
plot(all$hal ~ all$kos, pch=20, col=add.alpha('black',.2), xlab="Kosambi's", ylab="Haldane's")
abline(0,1,lty=2, col='red')
text(x=60,y=5, labels=paste0('r=',round(cor(all$hal,all$kos,use='complete.obs'),2)))
plot(all$hal ~ all$mor, pch=20, col=add.alpha('black',.2), xlab="Morgan's", ylab="Haldane's")
abline(0,1,lty=2, col='red')
text(x=60,y=5, labels=paste0('r=',round(cor(all$hal,all$mor,use='complete.obs'),2)))
plot(all$kos ~ all$mor, pch=20, col=add.alpha('black',.2), xlab="Morgan's", ylab="Kosambi's")
abline(0,1,lty=2, col='red')
text(x=60,y=5, labels=paste0('r=',round(cor(all$kos,all$mor,use='complete.obs'),2)))
mtext('Comparing mapping functions',outer=T,line=-2.5)
dev.off()


