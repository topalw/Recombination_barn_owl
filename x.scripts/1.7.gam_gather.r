source('functions.r')
library(tidyverse)

lm <- read.table('../1.data/1.LepMAP/orders/best_likelihood_pruned_orders.table',h=T)
lm <- lm[lm$ss != 'ss49',]
models <- list.files('rds/','*')
gam.lm <- data.frame()
## merge and plot 
pdf('plots/LM3_gams_fitted_all.pdf',width = 16, height = 16)
par(mfrow=c(2,2))
for(lg in unique(lm$lg)){
  # sex averaged 
  av <- readRDS(
    paste0('rds/',models[grep(paste0('mono_gam_lg',lg,'_'),models)])
  )
  plot(av$model$av ~ av$model$pos, col = add.alpha('black',.4),pch=16, xlab='Position',
       ylab='cM', main = paste0('sex-averaged map - LG ', lg),
       sub = paste0(unique(lm$ss[lm$lg == lg])))
  points(av$fitted.values ~ av$model$pos, col=add.alpha('red',0.4),pch=20)
  # males
  m <- readRDS(
    paste0('rds/',models[grep(paste0('mono_gam_male_lg',lg,'_'),models)])
  )
  plot(m$model$m ~ m$model$pos, col = add.alpha('black',.4),pch=16, xlab='Position',
       ylab='cM', main = paste0('male map - LG ', lg),
       sub = paste0(unique(lm$ss[lm$lg == lg])))
  points(m$fitted.values ~ m$model$pos, col=add.alpha('red',0.4),pch=20)
  # merge 
  tmp <- full_join(data.frame(m=m$fitted.values,pos=m$model[,2]),
                   data.frame(av=av$fitted.values,pos=av$model[,2]), by = 'pos')
  # females
  f <- readRDS(
    paste0('rds/',models[grep(paste0('mono_gam_female_lg',lg,'_'),models)])
  )
  plot(f$model$f ~ f$model$pos, col = add.alpha('black',.4),pch=16, xlab='Position',
       ylab='cM', main = paste0('female map - LG ', lg),
       sub = paste0(unique(lm$ss[lm$lg == lg])))
  points(f$fitted.values ~ f$model$pos, col=add.alpha('red',0.4),pch=20)
  # merge 2 
  tmp <- full_join(tmp,
                   data.frame(f=f$fitted.values,pos=f$model[,2]), by = 'pos')
  tmp$lg <- rep(lg,nrow(tmp))
  if(lg == 24) { 
    tmp$ss <- rep('ss2.2', nrow(tmp))
  } else { 
      tmp$ss <- rep(unique(lm$ss[lm$lg == lg]), nrow(tmp))
      }
  # make them start at 0 
  tmp$m <- tmp$m - min(tmp$m,na.rm=T)
  tmp$f <- tmp$f - min(tmp$f,na.rm=T)
  tmp$av <- tmp$av - min(tmp$av,na.rm=T)
  o1 <- order(tmp$pos)
  plot(tmp$av[o1]~tmp$pos[o1], col='black',lwd=2,lty=2, type='l',
       ylim=c(0,max(max(tmp$av,na.rm=T),max(tmp$m,na.rm=T),max(tmp$f,na.rm=T))), ylab='cM', xlab='position')
  lines(tmp$m[o1]~tmp$pos[o1], col='blue',lwd=2)
  lines(tmp$f[o1]~tmp$pos[o1], col='red',lwd=2)
  gam.lm <- rbind(tmp,gam.lm)
}
dev.off()

write.table(gam.lm, file='../1.data/1.LepMAP/orders/fitted_values_gams.table',
            quote=F, row.names = F)



