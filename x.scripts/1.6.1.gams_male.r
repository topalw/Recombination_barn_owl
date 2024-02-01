# tools
source('functions.r')
library(tidyverse)
library(splines)
library(mgcv)
library(scam)
# data 
lm <- read.table('../1.data/1.LepMAP/orders/best_likelihood_pruned_orders.table',h=T)
size <- lm %>% group_by(lg) %>% summarize(lg = unique(lg),
                                          size = max(pos) - min(pos)) 
lm$lg <- factor(lm$lg, levels=size$lg[order(size$size, decreasing = T)])
# repair ss3-ss49
lm <- lm[lm$ss != 'ss49',]
lm$pos[lm$ss=='ss3'] <- lm$pos[lm$ss=='ss3'] - 4308849
# rm NA
lm <- lm[! is.na(lm$m), ]
# initiate plot
pdf('plots/linkage_maps_gams_male.pdf',12,12)
# loop through each lg
for(lg in levels(lm$lg)){
tmp <- lm[lm$lg == lg, ]
# reprune ends cause weird ! 
tmp$m <- prune.ends(tmp$m)
tmp <- tmp[! is.na(tmp$m),]
# loess 
lo1 <- loess(m~pos,data=tmp, span=0.1)
pos.grid <- data.frame('pos'=seq(min(tmp$pos),max(tmp$pos),1000))
lo1.p <- predict(lo1, newdata=pos.grid,se=T)
# monotonic gam spline 
ss.tmp <- smooth.spline(tmp$pos, tmp$m, cv=T)# find k from df 
gam2 <- scam(m ~ s(pos,k = round(ss.tmp$df,0), bs='mpi'), data = tmp) # long step 
saveRDS(gam2, file=paste0('rds/mono_gam_male_lg',lg,'_',round(ss.tmp$df,0),'df.rds'))
gam2.p <- predict(gam2, newdata=pos.grid)
# plot 
plot(m~pos, data=tmp, ylab='m cM', xlab='physical pos',
     main=paste0('LG=',unique(lg)))
lines(lo1.p$fit ~ pos.grid$pos, col='navy', lwd=2,lty=1)
lines(gam2.p ~ pos.grid$pos, col='magenta', lwd=2,lty=1)
legend('bottomright',lty=c(1,1), col=c('navy','magenta'), 
       legend=c('loess span=0.1', 'monotonic splines GAM'))
}
dev.off()
