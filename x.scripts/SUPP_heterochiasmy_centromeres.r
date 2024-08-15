  
  jpeg('plots/SUPP_heterochiasmy_centros.jpg',12,16,units = 'in',res=400,quality=800)
  x <- c()
  rm <- c()
  rf <- c()
  layout(matrix(c(1,seq(1,19)),byrow = T,nrow=5,ncol=4))
  for(lg in lms2$lg[!is.na(lms2$centro)]){
    ss = lms$ss[lms$lg==lg]
    if(ss=='Super-Scaffold_2'){
      tmpx <- abs(d1m$start[d1m$ss==ss] - lms2$centro[lms2$lg==lg] - start2) / (lms2$len[lms2$lg==lg]*1e6)
    }else {
      tmpx <- abs(d1m$start[d1m$ss==ss] - lms2$centro[lms2$lg==lg] ) / (lms2$len[lms2$lg==lg]*1e6)
    }
    x <- c(x,tmpx)
    tmprm <- d1m$m[d1m$ss==ss] / mean(d1m$m[d1m$ss==ss],na.rm=T)
    tmprf <- d1m$f[d1m$ss==ss] / mean(d1m$f[d1m$ss==ss],na.rm=T)
    rm <- c(rm,tmprm)
    rf <- c(rf,tmprf)
  }
  cdf <- cbind.data.frame(x,rm,rf)
  lm <- loess(cdf$rm~cdf$x,span=0.1)
  lf <- loess(cdf$rf~cdf$x,span=0.1)
  plot(0,pch='',xlim=c(0,1),ylim=c(0,10),ylab='Recombination / mean',
       xlab='Relative position ',main='All')
  legend('topright',lwd=3,col=c(malecolor,femalecolor),legend=c('Male','Female'),bty='n',cex=.7)
  lines(lm,col=malecolor)
  lines(lf,col=femalecolor)
  for(lg in lms2$lg[!is.na(lms2$centro)]){
    ss = lms$ss[lms$lg==lg]
    if(ss=='Super-Scaffold_2'){
      tmpx <- abs(d1m$start[d1m$ss==ss] - lms2$centro[lms2$lg==lg] - start2) / (lms2$len[lms2$lg==lg]*1e6)
    }else {
      tmpx <- abs(d1m$start[d1m$ss==ss] - lms2$centro[lms2$lg==lg] ) / (lms2$len[lms2$lg==lg]*1e6)
    }
    tmprm <- d1m$m[d1m$ss==ss] / mean(d1m$m[d1m$ss==ss],na.rm=T)
    tmprf <- d1m$f[d1m$ss==ss] / mean(d1m$f[d1m$ss==ss],na.rm=T)
    plot(0,pch='',xlim=c(0,1),ylim=c(0,10),ylab='Recombination / mean ',
         xlab='Relative position ',main=paste0('LG',lg),
         sub=round(lms2$len[lms2$lg==lg],2))
    lm <- loess(tmprm~tmpx)
    lf <- loess(tmprf~tmpx)
    lines(lm,col=malecolor)
    lines(lf,col=femalecolor)
    legend('topright',lwd=3,col=c(malecolor,femalecolor),legend=c('Male','Female'),bty='n',cex=.7)
  }
  dev.off()
