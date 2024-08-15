jpeg('plots/SUPP_10kb_pyrho.jpg',14,18,units = 'in',res=600,quality=800)
par(mfrow=c(7,6))
for(lg in unique(lm$lg)){
  if(lg==40){
    next
  } else if(lg !=20){
    ss=unique(lm$ss[lm$lg==lg])
    plot(y=d10$CH_cM[d10$ss==ss],
         x=d10$start[d10$ss==ss]/1e6,
         ylim=c(-2,50),
         col='black',
         type='l',
         lwd=1,
         xlab='position (Mb)',
         ylab='recombination rate cM/Mb',
         main=paste('LG',lg))
    segments(x0=lm$pos[lm$lg ==lg & lm$maleCOs > 0]/1e6,
             x1=lm$pos[lm$lg ==lg & lm$maleCOs > 0]/1e6,
             y0=-1,
             y1=-.2,
             col=malecolor,
             lwd=1)
    segments(x0=lm$pos[lm$lg ==lg & lm$femaleCOs > 0]/1e6,
             x1=lm$pos[lm$lg ==lg & lm$femaleCOs > 0]/1e6,
             y0=-2,
             y1=-1.1,
             col=femalecolor,
             lwd=1)
  } else {
    ss1='Super-Scaffold_49'
    ss2='Super-Scaffold_3'
    tmp <- d10[d10$ss %in% c(ss1,ss2),]
    tmp$newpos <- tmp$start
    tmp$newpos[tmp$ss==ss1] <- abs(tmp$newpos[tmp$ss==ss1] - lengths$len[lengths$ss==ss1])
    tmp$newpos[tmp$ss==ss2] <- tmp$newpos[tmp$ss==ss2] + lengths$len[lengths$ss==ss1]
    plot(y=tmp$CH_cM[order(tmp$newpos)],
         x=tmp$newpos[order(tmp$newpos)]/1e6,
         ylim=c(-2,50),
         col='black',
         type='l',
         lwd=1,
         xlab='position',
         ylab='recombination rate cM/Mb',
         main=paste('LG',lg))
    segments(x0=lm$newpos[lm$lg ==lg & lm$maleCOs > 0]/1e6,
             x1=lm$newpos[lm$lg ==lg & lm$maleCOs > 0]/1e6,
             y0=-1,
             y1=-.2,
             col=malecolor,
             lwd=1)
    segments(x0=lm$newpos[lm$lg ==lg & lm$femaleCOs > 0]/1e6,
             x1=lm$newpos[lm$lg ==lg & lm$femaleCOs > 0]/1e6,
             y0=-2,
             y1=-1.1,
             col=femalecolor,
             lwd=1)
  }
}
dev.off()