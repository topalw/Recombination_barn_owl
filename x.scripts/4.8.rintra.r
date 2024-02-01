library(mgcv)
library(dplyr)
library(scam)

source('functions.r')
lm <- read.table('../1.data/1.LepMAP/orders/fitted_values_gams.table',
                 h=T)
lm$ss[lm$lg == 24] <- 'ss2.2'
gams <- list.files('rds/','mono_gam_lg*')
sizes <- lm %>% group_by(lg) %>% summarize('size'=max(pos)-min(pos))
sizes$l2 <- sizes$size / sum(sizes$size)
sizes$rintra <- NA
for(lg in sizes$lg){
 l2 <- sizes$l2[sizes$lg == lg]
 ss <- unique(gsub('ss','Super-Scaffold_',lm$ss[lm$lg==lg]))
 gam1 <- readRDS(paste0('rds/',gams[grep(paste0('lg',lg,'_'),gams)]))
 pos <- seq(min(lm$pos[lm$lg==lg]),max(lm$pos[lm$lg==lg]),10000)
 b <- predict(gam1,newdata = as.data.frame(pos))
 bdist <- matrix(data=NA,nrow=length(b),ncol=length(b))
 for(i in seq(1,length(b)-1)){
   for(j in seq(i,length(b))){
      bdist[i,j] <- rev.haldane(b[j] - b[i])
    }
 }
 sizes$rintra[sizes$lg==lg] = mean(bdist[upper.tri(bdist)])*l2^2
}



gams <- list.files('rds/','mono_gam_female_lg*')
sizes$rintra_F <- NA
for(lg in sizes$lg){
  l2 <- sizes$l2[sizes$lg == lg]
  ss <- unique(gsub('ss','Super-Scaffold_',lm$ss[lm$lg==lg]))
  gam1 <- readRDS(paste0('rds/',gams[grep(paste0('lg',lg,'_'),gams)]))
  pos <- seq(min(lm$pos[lm$lg==lg]),max(lm$pos[lm$lg==lg]),10000)
  b <- predict(gam1,newdata = as.data.frame(pos))
  bdist <- matrix(data=NA,nrow=length(b),ncol=length(b))
  for(i in seq(1,length(b)-1)){
    for(j in seq(i,length(b))){
      bdist[i,j] <- rev.haldane(b[j] - b[i])
    }
  }
  sizes$rintra_F[sizes$lg==lg] = mean(bdist[upper.tri(bdist)])*l2^2
}
 

gams <- list.files('rds/','mono_gam_male_lg*')
sizes$rintra_M <- NA
for(lg in sizes$lg){
  l2 <- sizes$l2[sizes$lg == lg]
  ss <- unique(gsub('ss','Super-Scaffold_',lm$ss[lm$lg==lg]))
  gam1 <- readRDS(paste0('rds/',gams[grep(paste0('lg',lg,'_'),gams)]))
  pos <- seq(min(lm$pos[lm$lg==lg]),max(lm$pos[lm$lg==lg]),10000)
  b <- predict(gam1,newdata = as.data.frame(pos))
  bdist <- matrix(data=NA,nrow=length(b),ncol=length(b))
  for(i in seq(1,length(b)-1)){
    for(j in seq(i,length(b))){
      bdist[i,j] <- rev.haldane(b[j] - b[i])
    }
  }
  sizes$rintra_M[sizes$lg==lg] = mean(bdist[upper.tri(bdist)])*l2^2
}

write.table(sizes,'../1.data/1.LepMAP/orders/rintra.table',
            quote=F,row.names=F)
