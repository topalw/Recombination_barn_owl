setwd('C:/Users/topalw/Desktop/PhD/Analyses/Chapter_A/R.hotspots/quality_control/')

pdf('../../x.scripts/plots/compare_hotspots.pdf')
par(mfrow=c(2,2))
# look into HWE
d <- read.table('hotspots/CH.hotspots.012')
d <- d[,-1]
miss <- apply(d,2,FUN=function(c) return(sum(c==-1)) )
het <- apply(d,2,FUN=function(c) return(sum(c==1)/sum(c!=-1)) )
af <- apply(d,2,FUN=function(c) return(sum(c[c!=-1])/sum(c!=-1)/2) )

dn <- read.table('nonhotspots/CH.nonHotspots.012')
nmiss <- apply(dn,2,FUN=function(c) return(sum(c==-1)) )
hetn <- apply(dn,2,FUN=function(c) return(sum(c==1)/sum(c!=-1)) )
afn <- apply(dn,2,FUN=function(c) return(sum(c[c!=-1])/sum(c!=-1)/2) )


plot(het~af,pch=20,col='black',main='HWE')
points(y=hetn,x=afn,pch=20,col=add.alpha('red',.4),cex=.8)
lines(x=seq(0,1,0.01),y=seq(0,1,0.01)*(1-seq(0,1,0.01))*2,col='blue',lty=2,lwd=4)
legend('topright',pch=c(20,20,20),col=c('black','red','blue'),
       legend=c('hotspots','nonHotspots','2pq'))

# QUAL 
q <- read.table('hotspots/CH.hotspots.lqual',h=T)
nq <- read.table('nonhotspots/CH.nonHotspots.lqual',h=T)
boxplot(q$QUAL,nq$QUAL,names=c('hotspots','nonhotspots'),main='QUALITY',xlab='',ylab='QUAL')

# MISS 
boxplot(miss,nmiss,names=c('hotspots','nonhotspots'),main='MISSINGNESS',xlab='',ylab='F_MISS')

# DP
dp <- read.table('hotspots/CH.hotspots.ldepth',h=T)
ndp <- read.table('nonhotspots/CH.nonHotspots.ldepth',h=T)
boxplot(dp$SUM_DEPTH,ndp$SUM_DEPTH,names=c('hotspots','nonhotspots'),main='DEPTH',xlab='',ylab='DP')

dev.off()