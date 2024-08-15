palette = add.alpha(c('purple','orange'),.6)

pdf('plots/SUPP_compare100kb.pdf',12,16)
par(mfrow=c(5,4))
lgs <- lms$lg[! lms$lg %in% c(20,40)]
lgs <- lgs[order(lgs,decreasing = F)]
for(lg in lgs){
  ss = lms$ss[lms$lg==lg]
  tmp2 <- d100[d100$ss==ss,]
  tmp <- d1m[d1m$ss==ss,]
  plot(tmp2$CH_cM~tmp2$start,type='l',col=palette[1],lwd=4,lty=2,
       ylim=c(0,max(c(tmp2$CH_cM,tmp$av,tmp$f,tmp$m),na.rm=T)),
       ylab='cM/Mb',xlab='Physical position',main=paste0('LG',lg))
  lines(tmp$av~tmp$start,col=palette[2],lwd=4,lty=1) 
  legend('topright',lwd=3,col=palette,lty=c(2,1),legend=c('pyrho','LepMAP3'),
         bty='n')
}
dev.off()