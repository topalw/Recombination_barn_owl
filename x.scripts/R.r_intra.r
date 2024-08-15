setwd('C:/Users/topalw/Desktop/PhD/Analyses/Chapter_A/x.scripts/')
library(tidyverse)
source('functions.r')
####### prepare ####### 
lm <- read_delim('../1.data/1.LepMAP/R.AllSNPs_physical/physical_orders_full.txt')
lengths <- read.table('../2.metadata/Tyto_reference_Jan2020.bed')[,c(1,3)]

# get nonN length 
names(lengths) <- c('ss','len')
# modify LG24 and SS2 
lengths <- rbind.data.frame(lengths,
                            data.frame('ss'='Super-Scaffold_2.2',
                                       'len'=28222627))
# there are interspersed 467'039 Ns so remove 28222627 + 467039 from length of second part 
lengths$len[lengths$ss=='Super-Scaffold_2'] <- lengths$len[lengths$ss=='Super-Scaffold_2']  - start2
# keep only in LGs 
lengths <- lengths[lengths$ss %in% lm$ss,]
# annotate LGs 
lengths$lg <-lm$lg[match(lengths$ss,lm$ss)]
# group
lms <- lm %>% group_by(ss) %>% summarize(lg=unique(lg),
                                         m=max(m,na.rm=T),
                                         f=max(f,na.rm=T),
                                         av=max(av,na.rm=T),
                                         len=NA,
                                         loci = n())
lms$len <- lengths$len[match(lms$ss,lengths$ss)]

# l2 as proportion of length in each LG - same as interspersed markers 
lms$l2 <- (lms$len / sum(lms$len))^2
# correct lg 1 length 
lm$newpos[lm$lg==1] <- lm$newpos[lm$lg==1] - start2 
####### calculate rintra ####### 

# the calculation entails the following steps 
# STEP 1 make B a set of equally spaced pseudopoints
# STEP 2 for each b in B find proximate 'a' markers in linkage map 
# before b= a(b) & after b= a(b) + 1 calculate their distance z = pos(a(b)+1) - pos(a(b))
# and calculate their distance to the b marker x= pos(a(b)+1) - pos(b) & 
# y = pos(b) - pos(a(b)) and estimate the Db = x/z * D(a(b)) + y/z D(a(b)+1)
# where D(x) is the cM position of marker x 
# then calculate all pairwise distances in B set and get average * l2 which r_intra

# For pseudomarkers before the first a set 0 
# For pseudomarkers after last set max 
# assume cM = shuffling (Morgan's function)
# after estimating D(b) for each b then get mean of it and weigh by l2 (L^2)
# do it for both male and female 

# a function to find proximal a markers and calculate x,y,z,Db
# a(b) is named left and a(b)+1 is named right as relative positions to b
proximity_mine <- function(b,A,Da){
  # if before first or after last 
  if(b < A[1]){
    Db = 0
  } else if(b > A[length(A)]){
    Db = max(Da)
  } else { 
    # now we are talking 
    dif <- b - A 
    # maximum negative number 
    left <- which(dif == max(dif[dif<0]))
    # minimum positive number
    right <- which(dif == min(dif[dif>0]))
    x = A[right] - b
    y = b - A[left]
    z = A[right] - A[left]
    Db = (x/z) * Da[left] + (y/z)*Da[right]
  }
  return(Db)
}

step = 10000 # 10kb step between b pseudomarkers 

lms$rintra_m <- 0
lms$rintra_f <- 0 
# repeat for each LG
for(lg in 1:40){
  print(lg)
  # handle merged 
  if(lg %in% c(20,40)){
    # add them since length and l2 is based on scaffolds
    lg.l2 <- sum(lms$l2[lms$lg==lg])
    last = sum(lms$len[lms$lg==lg])
  } else {
    lg.l2 <- lms$l2[lms$lg==lg]
    last <- lms$len[lms$lg==lg]
  } 
  # STEP 1 - make pseudo markers b 
  first = 1 
  
  B <- seq(first,last,by=step)
  Db_m <- rep(0,length(B)) # male positions 
  Db_f <- rep(0,length(B)) # female positions 
  # STEP 2 - find proximate markers a - use newpos for merged lgs (20/40)
  A <- lm$newpos[lm$lg==lg]
  Da_m <- lm$m[lm$lg==lg] # male positions 
  Da_f <- lm$f[lm$lg==lg] # female positions 
  for(i in 1:length(B)){
    b <- B[i]
    Db_m[i] <- proximity_mine(b,A,Da_m)
    Db_f[i] <- proximity_mine(b,A,Da_f)
  } 
  # distance matrices
  dm <- as.matrix(dist(Db_m))
  df <- as.matrix(dist(Db_f))
  # mean * l2 
  lms$rintra_m[lms$lg==lg] <- mean(dm[upper.tri(dm)],na.rm=T) * lg.l2
  lms$rintra_f[lms$lg==lg] <- mean(df[upper.tri(df)],na.rm=T) * lg.l2 
}

# since we modified them lets write them:
write.table(lengths,'../1.data/lengths.txt',quote=F,row.names=F)
write.table(lms,'../1.data/1.LepMAP/R.AllSNPs_physical/linkage_map_perLG.txt',quote=F,row.names=F)
