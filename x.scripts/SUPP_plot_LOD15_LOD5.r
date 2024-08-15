setwd('C:/Users/topalw/Desktop/PhD/Analyses/Chapter_A/x.scripts/')
library(tidyverse)
source('functions.r')

old <- read_delim('../1.data/1.LepMAP/orders/best_likelihood_pruned_orders.table')
phys <- read_delim('../1.data/1.LepMAP/R.AllSNPs_physical/physical_orders_full_COs.txt')

jpeg('plots/SUPP_LOD15_LOD5.jpg',14,18,units = 'in',res=600,quality=800)
par(mfrow=c(7,6))
for(lg in 1:39){
plot(phys$av[phys$lg==lg]~phys$newpos[phys$lg==lg],
     type='l',lwd=4,xlab='Physical position',main=paste('LG',lg),
     ylim=c(0,max(phys$av[phys$lg==lg],old$av[old$lg==lg],na.rm=T)),
     ylab='cM')
points(old$av[old$lg==lg]~old$pos[old$lg==lg],pch=20,cex=1,col='purple')
legend('topleft',lty=c(1,1),lwd=c(4,4),col=c('black','purple'),
       legend=c('LOD5 physical','LOD15 de-novo'),bty='n')
}
dev.off()
