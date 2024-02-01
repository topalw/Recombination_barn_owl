library(tidyverse)
source('functions.r')

snps <- read.table('../1.data/1.LepMAP/maps/Call2_tf1_maf05_miss05_1kb_pruned_f.call.snps.txt',h=T)
snps$name <- seq(1,nrow(snps))# if out of range use the non pruned  
snps$CHR[snps$CHR =='Super-Scaffold_2' & snps$POS <= 28692481] <- 'Super-Scaffold_2.2'
# this will be the guide for which run to use 
lik <- read.table('../1.data/1.LepMAP/orders/likelihoods/best_for_eachLG.txt',h=T)

# this will be the guide for which markers to keep 
final <- read.table('../1.data/1.LepMAP/orders/fitted_values_gams.table',h=T)

######################

name.v <- c('1st_ordering','2nd_ordering','3rd_ordering')
end.v <- c('_nogp.txt','_nogp_run2.txt','_nogp_run3.txt')

for(lg in unique(lik$lg)){
  keep <- snps$name[snps$CHR == unique(gsub('ss','Super-Scaffold_',final$ss[final$lg==lg])) 
                    & snps$POS %in% final$pos[final$lg==lg]]
  ordering <- as.numeric(lik$best[lik$lg == lg])
  d <- read.table(paste0('../1.data/1.LepMAP/orders/',name.v[ordering],
                         '/Call2_tf1_maf05_miss05_1kb_pruned_f.call_order_',lg,
                         end.v[ordering]))
  s1 <- seq(8,90,1)
  s2 <- seq(7,88,3)
  relevant <- c(1, s1[! s1 %in% s2])
  d <- d[d$V1 %in% keep,relevant]
  m <- seq(2,57,2)
  f <- seq(3,57,2)
  if(lg == unique(lik$lg)[1]){
  fin <- data.frame(
  'pos' = snps$POS[match(d$V1,snps$name)],
  'lg' = rep(lg,nrow(d)),
  'maleCO' = rowSums(d[,m]),
  'femaleCO' = rowSums(d[,f]))
  } else{
    fin <- rbind(fin,
                 data.frame(
                   'pos' = snps$POS[match(d$V1,snps$name)],
                   'lg' = rep(lg,nrow(d)),
                   'maleCO' = rowSums(d[,m]),
                   'femaleCO' = rowSums(d[,f])))
  }
}
sum(as.numeric(table(fin$maleCO))[-1]) # 2787 COs
sum(as.numeric(table(fin$femaleCO))[-1]) # 2887 COs 

write.table(fin,'../1.data/1.LepMAP/orders/COs_male_female.table',quote=F,row.names=F)
