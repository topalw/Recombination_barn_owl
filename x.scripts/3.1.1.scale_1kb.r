setwd('C:/Users/topalw/Desktop/PhD/Analyses/Chapter_A/x.scripts/')
source('functions.r')

pathe = '../1.data/2.pyrho/estimates/'
files <- list.files(pathe, '*1000.thetaxN')

lm <- read.table('../1.data/1.LepMAP/R.AllSNPs_physical/linkage_map_perLG.txt',h=T)
lm$ss[lm$lg==24] <- 'Super-Scaffold_2.2'

for(file in files){
  d <- read.table(paste0(pathe,file),h=T)
  d$hald <- hald(d)
  d$ss[d$ss == 'Super-Scaffold_2' & d$end <= len2.2] <- 'Super-Scaffold_2.2'
  d <- d[d$ss %in% lm$ss,]
  for(ss in lm$ss){
    sf <- as.numeric(lm$av[lm$ss==ss] / sum(d$hald[d$ss==ss],na.rm=T))
    d$hald[d$ss==ss] <-d$hald[d$ss==ss] * sf
  }
#  scale.factor <- lmlen / sum(d$hald,na.rm=T)
#  d$cM.scaled <- d$hald * scale.factor
  d <- d[,c(1,2,3,5,6)]
  pop <- gsub('_total.theta.intersect.bed_1000.thetaxN','',file)
  names(d)[c(4,5)] <- c(paste0(pop,'_hi'),paste0(pop,'_cM'))
  file2=paste0(file,'.scaled')
  write.table(d,paste0('../1.data/2.pyrho/scaled_datasets/',file2), quote=F,row.names=F)
}