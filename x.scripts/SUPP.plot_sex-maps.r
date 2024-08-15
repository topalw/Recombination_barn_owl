lm <- readr::read_delim('../1.data/1.LepMAP/R.AllSNPs_physical/physical_orders_full.txt')
options(scipen=0)
source('functions.r')
lm$newpos[lm$lg==1] <- lm$newpos[lm$lg==1]  - start2
jpeg('plots/SUPP_sex_orders.jpg',14,18,units = 'in',res=600,quality=800)
#tiff('plots/SUPP_sex_specific.tiff',1440,1920,res=100)
par(mfrow=c(7,6))
for(lg in unique(lm$lg)){
  tmp <- lm[lm$lg==lg,]
  ss = unique(lm$ss[lm$lg==lg])
  plot(tmp$m~tmp$newpos,type='l',col='purple',lwd=3,
       xlab='Physical position (bp)', ylab='Genetic position (cM)',
       main=paste0('LG',lg),ylim=c(0,max(c(tmp$m,tmp$f),na.rm=T)),sub=ifelse(length(ss) >  1, paste(ss[1],ss[2]),ss))
  lines(tmp$f ~ tmp$newpos,lwd=3,col='orange')
  legend('topleft',lwd=2,col=c('purple','orange'),legend=c('Male','Female'),bty='n')
  if(lg==40){
    rect(xleft=85850542,xright=90242836,ybot=0,ytop=100,col=add.alpha('blue',alpha=.3))
  }
}
dev.off()
