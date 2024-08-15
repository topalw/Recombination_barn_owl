jpeg('plots/SUPP_CO_MF_allLG.jpg',10,10,units = 'in',res=600,quality=800)
par(mar=c(5.1,4.1,4.1,4.1))
lengths2 <- lengths %>% group_by(lg) %>% summarize(len=sum(len))
lgs <- lengths2$lg[order(lengths2$len,decreasing = F)]
lens <- round(lengths2$len[order(lengths2$len,decreasing = F)] / 1e6,2)
plot(0,pch='',ylim=c(0,140),xlim=c(0,1),
     axes=F,xlab='Proportion of total length',ylab='LG')
axis(1, at = seq(0,1,.25),labels = seq(0,100,25))
axis(2, at = seq(1,120,3),lgs,las=2,cex.axis=.8)
axis(4,at= seq(1,120,3),lens,las=2,cex.axis=.8)
counter = 0
lm$prop <- 0
for(lg in lgs){
  lm$prop[lm$lg==lg] <- lm$newpos[lm$lg==lg] / max(lm$newpos[lm$lg==lg])
  if(lg==1){
    lm$prop[lm$lg==lg] <- (lm$newpos[lm$lg==lg] - len2.2)/ (max(lm$newpos[lm$lg==lg]) - len2.2)
  }
  ss=lms$ss[lms$lg==lg]
  y=counter
  fco <- lm$prop[lm$lg==lg & lm$femaleCOs >0 ]
  mco <- lm$prop[lm$lg==lg & lm$maleCOs >0 ]
  points(x=fco,y=rep(y,length(fco)),
                     pch=20,cex=1,col=femalecolor)
  points(x=mco,y=rep(y+1,length(mco)),
         pch=20,cex=1,col=malecolor)
  counter=counter+3
  abline(h=counter-1,col=add.alpha('grey',.4))
}
legend('topright',pch=20,col=c(malecolor,femalecolor),legend=c('Male CO','Female CO'),
       bty='n')
dev.off()
