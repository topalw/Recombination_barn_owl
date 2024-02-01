{# 
i = 1 
hot2 <- hot1 
while(i <= length(hot2)){
  indx <- ifelse(hot2 %in% c(hot2[i]-1,hot2[i] + 1),F,T)
  hot2 <- hot2[indx]
  i <- i +1 
}
chhot2 <- hot2
#
i = 1 
hot2 <- pthot 
while(i <= length(hot2)){
  indx <- ifelse(hot2 %in% c(hot2[i]-1,hot2[i] + 1),F,T)
  hot2 <- hot2[indx]
  i <- i +1 
}
pthot2 <- hot2 
#
i = 1 
hot2 <- gbhot 
while(i <= length(hot2)){
  indx <- ifelse(hot2 %in% c(hot2[i]-1,hot2[i] + 1),F,T)
  hot2 <- hot2[indx]
  i <- i +1 
}
gbhot2 <- hot2 
# 
i = 1 
hot2 <- ch13hot
while(i <= length(hot2)){
  indx <- ifelse(hot2 %in% c(hot2[i]-1,hot2[i] + 1),F,T)
  hot2 <- hot2[indx]
  i <- i +1 
}
ch13hot2 <- hot2 
# 
i = 1 
hot2 <- ch13_2hot
while(i <= length(hot2)){
  indx <- ifelse(hot2 %in% c(hot2[i]-1,hot2[i] + 1),F,T)
  hot2 <- hot2[indx]
  i <- i +1 
}
ch13_2hot2 <- hot2 
}



### PLOTS 
{
jpeg('plots/venn1.jpg')
venn1 <- VennDiagram::venn.diagram(list(chhot2,pthot2),
                                  category.names=c('CH','PT'),
                                  filename=NULL,
                                  disable.logging = T,
                                  label.col = 'white',
                                  cex=1.5,
                                  main='Local Hotspot overlap',
                                  fill=c('#FF0000','#006600'),
                                  col=rep('white',2),
                                  lwd=rep(4,2))
grid.draw(venn1)
grid.draw(venn1)
dev.off()
jpeg('plots/venn2.jpg')
venn2 <- VennDiagram::venn.diagram(list(chhot2,gbhot2),
                                   category.names=c('CH','GB'),
                                   filename=NULL,
                                   disable.logging = T,
                                   label.col = 'white',
                                   cex=1.5,
                                   main='Local Hotspot overlap',
                                   fill=c('#FF0000','#012169'),
                                   col=rep('white',2),
                                   lwd=rep(4,2))
grid.draw(venn2)
grid.draw(venn2)
dev.off()
jpeg('plots/venn3.jpg')
venn3 <- VennDiagram::venn.diagram(list(ch13hot2,gbhot2),
                                   category.names=c('CH13','GB'),
                                   filename=NULL,
                                   disable.logging = T,
                                   label.col = 'white',
                                   cex=1.5,
                                   main='Local Hotspot overlap',
                                   fill=c('#FF2000','#012169'),
                                   col=rep('white',2),
                                   lwd=rep(4,2))
grid.draw(venn3)
grid.draw(venn3)
dev.off()
jpeg('plots/venn4.jpg')
venn4 <- VennDiagram::venn.diagram(list(ch13_2hot2,gbhot2),
                                   category.names=c('CH13_2','GB'),
                                   filename=NULL,
                                   disable.logging = T,
                                   label.col = 'white',
                                   cex=1.5,
                                   main='Local Hotspot overlap',
                                   fill=c('#FF2000','#012169'),
                                   col=rep('white',2),
                                   lwd=rep(4,2))
grid.draw(venn4)
grid.draw(venn4)
dev.off()
jpeg('plots/venn5.jpg')
venn5 <- VennDiagram::venn.diagram(list(ch13hot2,pthot2),
                                   category.names=c('CH13','GB'),
                                   filename=NULL,
                                   disable.logging = T,
                                   label.col = 'white',
                                   cex=1.5,
                                   main='Local Hotspot overlap',
                                   fill=c('#FF2000','#006600'),
                                   col=rep('white',2),
                                   lwd=rep(4,2))
grid.draw(venn5)
grid.draw(venn5)
dev.off()
jpeg('plots/venn6.jpg')
venn6 <- VennDiagram::venn.diagram(list(ch13_2hot2,pthot2),
                                   category.names=c('CH13_2','GB'),
                                   filename=NULL,
                                   disable.logging = T,
                                   label.col = 'white',
                                   cex=1.5,
                                   main='Local Hotspot overlap',
                                   fill=c('#FF2000','#006600'),
                                   col=rep('white',2),
                                   lwd=rep(4,2))
grid.draw(venn6)
grid.draw(venn6)
dev.off()
jpeg('plots/venn7.jpg')
venn7 <- VennDiagram::venn.diagram(list(chhot2,ch13hot2),
                                   category.names=c('CH','CH13_2'),
                                   filename=NULL,
                                   disable.logging = T,
                                   label.col = 'white',
                                   cex=1.5,
                                   main='Local Hotspot overlap',
                                   fill=c('#FF0000','#FF0000'),
                                   col=rep('white',2),
                                   lwd=rep(4,2))
grid.draw(venn7)
grid.draw(venn7)
dev.off()
jpeg('plots/venn8.jpg')
venn8 <- VennDiagram::venn.diagram(list(chhot2,ch13_2hot2),
                                   category.names=c('CH','CH13_2'),
                                   filename=NULL,
                                   disable.logging = T,
                                   label.col = 'white',
                                   cex=1.5,
                                   main='Local Hotspot overlap',
                                   fill=c('#FF0000','#FF0000'),
                                   col=rep('white',2),
                                   lwd=rep(4,2))
grid.draw(venn8)
grid.draw(venn8)
dev.off()
jpeg('plots/venn9.jpg')
venn9 <- VennDiagram::venn.diagram(list(ch13hot2,ch13_2hot2),
                                   category.names=c('CH13','CH13_2'),
                                   filename=NULL,
                                   disable.logging = T,
                                   label.col = 'white',
                                   cex=1.5,
                                   main='Local Hotspot overlap',
                                   fill=c('#FF0000','#FF0000'),
                                   col=rep('white',2),
                                   lwd=rep(4,2))
grid.draw(venn9)
grid.draw(venn9)
dev.off()
}
# hotspot overlaps
# make matrix 
tmp <- matrix(NA,5,5)
rownames(tmp) <- colnames(tmp) <- c('CH13_2','CH13','CH','PT','GB')

# fill upper part with local hotspot overlap 
# ch 
tmp[3,4] <- 100*sum(c(chhot2-1,chhot2,chhot2+1)%in% pthot2)/length(chhot2)
tmp[3,5] <- 100*sum(c(chhot2-1,chhot2,chhot2+1)%in% gbhot2)/length(chhot2)
tmp[2,3] <- 100*sum(c(chhot2-1,chhot2,chhot2+1)%in% ch13hot2)/length(chhot2)
tmp[1,3] <- 100*sum(c(chhot2-1,chhot2,chhot2+1)%in% ch13_2hot2)/length(chhot2)
# ch13
tmp[2,4] <- 100*sum(c(ch13hot2-1,ch13hot2,ch13hot2+1) %in% pthot2)/length(ch13hot2)
tmp[2,5] <- 100*sum(c(ch13hot2-1,ch13hot2,ch13hot2+1) %in% gbhot2)/length(ch13hot2)
tmp[1,2] <- 100*sum(c(ch13hot2-1,ch13hot2,ch13hot2+1) %in% ch13_2hot2)/length(ch13hot2)
# ch13_2 
tmp[1,4] <- 100*sum(c(ch13_2hot2-1,ch13_2hot2,ch13_2hot2+1) %in% pthot2)/length(ch13_2hot2)
tmp[1,5] <- 100*sum(c(ch13_2hot2-1,ch13_2hot2,ch13_2hot2+1) %in% gbhot2)/length(ch13_2hot2)
# pt - gb
tmp[4,5] <- 100*sum(c(pthot2 - 1,pthot2,pthot2+1) %in% gbhot2)/length(pthot2)
# fill lower part with global hotspot overlaps 
tmp[4,3] <- 100*sum(hot10 %in% pthot10)/length(hot10)
tmp[5,3] <- 100*sum(hot10 %in% gbhot10)/length(hot10)
tmp[3,2] <- 100*sum(hot10 %in% ch13hot10)/length(hot10)
tmp[3,1] <- 100*sum(hot10 %in% ch13_2hot10)/length(hot10)
# ch13
tmp[4,2] <- 100*sum(ch13hot10 %in% pthot10)/length(ch13hot10)
tmp[5,2] <- 100*sum(ch13hot10 %in% gbhot10)/length(ch13hot10)
tmp[2,1] <- 100*sum(ch13hot10 %in% ch13_2hot10)/length(ch13hot10)
# ch13_2 
tmp[4,1] <- 100*sum(ch13_2hot10 %in% pthot10)/length(ch13_2hot10)
tmp[5,1] <- 100*sum(ch13_2hot10 %in% gbhot10)/length(ch13_2hot10)
# pt - gb
tmp[5,4] <- 100*sum(pthot10 %in% gbhot10)/length(pthot10)
write.table(tmp,'../1.data/2.pyrho/scaled_datasets/hotspot_sharing.table')

tmp1 <- matrix(NA, ncol=length(hot2), nrow=4)


