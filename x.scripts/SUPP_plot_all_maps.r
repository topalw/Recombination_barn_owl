d10 <- d10[!is.na(d10$ss),]
jpeg('plots/SUPP_all_maps.jpg',14,18,units = 'in',res=600,quality=800)
par(mfrow=c(7,6))
for(ss in unique(d10$ss)){
  l <- max(d10$end[d10$ss==ss],na.rm=T)
  
  plot(d10$GB_cM[d10$ss==ss]~d10$end[d10$ss==ss],
       type='l',lwd=2,ylim=c(0,200),col=popcol[1],
       xlab='',ylab='',axes=F,pch=20,main=paste('LG',lms$lg[lms$ss==ss]),
       sub=ss)
  
  ## PT
  lines(40+d10$PT_cM[d10$ss==ss]~
          d10$end[d10$ss==ss],
        type='l',lwd=2,col=popcol[2])
  # CH
  lines(80+d10$CH_cM[d10$ss==ss ]~
          d10$end[d10$ss==ss],
        type='l',lwd=2,col=popcol[3])
  ## CH13 
  lines(120+d10$CH13_cM[d10$ss==ss]~
          d10$end[d10$ss==ss] ,
        type='l',lwd=2,col=popcol[4],lty=1)
  ## CH13_2
  lines(160+d10$CH13_2_cM[d10$ss==ss] ~
          d10$end[d10$ss==ss] ,
        type='l',lwd=2,col=popcol[5],lty=1)
  axis(1,at=seq(0,l,5e6),labels=seq(0,l,5e6),cex.axis=axis.size)
  mtext('Mb',1,3,cex=mtext.size)
  mtext('cM/Mb',2,3,cex=mtext.size)
  axis(2,at=c(0,20),labels=c(0,20),cex.axis=axis.size,las=2)
  axis(2,at=c(40,60),labels=c(0,20),cex.axis=axis.size,las=2)
  axis(2,at=c(80,100),labels=c(0,20),cex.axis=axis.size,las=2)
  axis(2,at=c(120,140),labels=c(0,20),cex.axis=axis.size,las=2)
  axis(2,at=c(160,180),labels=c(0,20),cex.axis=axis.size,las=2)
  text('GB',x=l,y=10)
  text('PT',x=l,y=50)
  text('CH',x=l,y=90)
  text('CH13',x=l,y=130,cex=.8)
  text('CH13_2',x=l,y=170,cex=.8)
}
dev.off()