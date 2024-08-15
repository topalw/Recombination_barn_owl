setwd("C:/Users/topalw/Desktop/PhD/Analyses/Chapter_A/x.scripts/")
source('functions.r')

########################################################
########################################################
########################################################
########################################################

# this will be the guide for which markers to keep 
final <- readr::read_delim('../1.data/1.LepMAP/R.AllSNPs_physical/physical_orders_full.txt')
final <- final[! final$lg %in% c(20,40),] # handle separately 
# how many kids are in each family 
kids.fams <- c(3,5,6,5,4,5,5,5,2,4,4,4,4,4,5,5,6,5,5,5,4,4,3,6,6,8,8,8)
final$ss[final$lg == 24] <- 'Super-Scaffold_2.2'

final$maleCOs <- NA
final$femaleCOs <- NA

for(male in c(TRUE, FALSE)){
  # a matrix with snps in rows and kids in columns
  for( lg in unique(final$lg) ){
    ss <- unique(final$ss[final$lg==lg])
    # read order to calculate COs 
    fs <- paste0('../1.data/1.LepMAP/R.AllSNPs_physical/orders/',ss,'_f.call_lod5_map.phys_order.txt')
    d <- read.table(fs, colClasses = 'character') 
    # these columns in the ordered output are:
    # ex 111000 -> paternal haplotypes 111 | maternal haplotyes 000 of 3 kids
    # and then 2 columns with 0/1 if CO happened in Male or Female 
    s1 <- seq(7,90,1) # all such columns
    s2 <- seq(7,88,3) # family info column index ex 111000 above 
    fam.co <- s1[! s1 %in% s2] # the parental crossover columns
    # make two dataframe with position and CO in each parent in each family (cos)
    # and one with the family phase in each kid in each family
    crossovers <- d[,fam.co] # crossovers info and positions 
    haplotypes <- d[,s2] # family phase info and positions
    # make numeric 
    for(i in 1:ncol(crossovers)){ 
      crossovers[,i] <- as.numeric(crossovers[,i])
    }
    # either keep only male or female crossovers 
    if(male){
      crossovers <- crossovers[,seq(1,ncol(crossovers),2)] # keep only male COs 
      final$maleCOs[final$lg==lg] <- as.numeric(rowSums(crossovers))
    }else { 
      crossovers <- crossovers[,seq(2,ncol(crossovers),2)] # keep only female COs 
      final$femaleCOs[final$lg==lg] <- as.numeric(rowSums(crossovers ))
    }
  }
}

# handle LG20 and LG40 separately cause they have unique ordering


##########################################
##########################################
############## handle LG 20 ##############
##########################################
##########################################
final2 <- readr::read_delim('../1.data/1.LepMAP/R.AllSNPs_physical/physical_orders_full.txt')
lm20 <- final2[final2$lg == 20,]
lm20$maleCOs <- NA
lm20$femaleCOs <- NA
for(male in c(TRUE, FALSE)){
  d <- read.table('../1.data/1.LepMAP/R.AllSNPs_physical/orders/merge_3_49_f.call_lod5_map.phys_order.txt')
  # need to order d based on newpos and then everything is smooth 
  tmp <- read.table('../1.data/1.LepMAP/R.AllSNPs_physical/perf_orders/merge_3_49_f.call_lod5_map.phys_order.txt.perf')
  indx <-  match(paste0(tmp$V1,tmp$V2),paste0(lm20$ss,lm20$pos))
  d <- d[indx,]
  # these columns in the ordered output are:
  # ex 111000 -> paternal haplotypes 111 | maternal haplotyes 000 of 3 kids
  # and then 2 columns with 0/1 if CO happened in Male or Female 
  s1 <- seq(7,90,1) # all such columns
  s2 <- seq(7,88,3) # family info column index ex 111000 above 
  fam.co <- s1[! s1 %in% s2] # the parental crossover columns
  # make two dataframe with position and CO in each parent in each family (cos)
  # and one with the family phase in each kid in each family
  crossovers <- d[,fam.co] # crossovers info and positions 
  haplotypes <- d[,s2] # family phase info and positions
  # make numeric 
  for(i in 1:ncol(crossovers)){ 
    crossovers[,i] <- as.numeric(crossovers[,i])
  }
  # either keep only male or female crossovers 
  if(male){
    crossovers <- crossovers[,seq(1,ncol(crossovers),2)] # keep only male COs 
    lm20$maleCOs <- as.numeric(rowSums(crossovers))
  }else { 
    crossovers <- crossovers[,seq(2,ncol(crossovers),2)] # keep only female COs 
    lm20$femaleCOs <- as.numeric(rowSums(crossovers))
  }
}

final <- rbind.data.frame(final, lm20)


##########################################
##########################################
############## handle LG 40 ##############
##########################################
##########################################
# Z wont work as above  so just take positions of COs in male and females 

final2 <- readr::read_delim('../1.data/1.LepMAP/R.AllSNPs_physical/physical_orders_full.txt')
lm20 <- final2[final2$lg == 40,]
lm20$maleCOs <- NA
lm20$femaleCOs <- NA
for(male in c(TRUE, FALSE)){
  if(male){
  d <- read.table('../1.data/1.LepMAP/R.AllSNPs_physical/orders/Z2_f.call_lod5_map.phys_order_Male.txt')
  tmp <- read.table('../1.data/1.LepMAP/R.AllSNPs_physical/perf_orders/Z2_f.call_lod5_map.phys_order_Male.txt.perf')
  } else{
  d <- read.table('../1.data/1.LepMAP/R.AllSNPs_physical/orders/Z2_f.call_lod5_map.phys_order_Female.txt')
  tmp <- read.table('../1.data/1.LepMAP/R.AllSNPs_physical/perf_orders/Z2_f.call_lod5_map.phys_order_Female.txt.perf')
  }
  tmp$V2 <- gsub('\\*','',tmp$V2)
  d[,1] <- tmp[,1]
  d[,2] <- tmp[,2]
  # make output table 
  tlg <- as.data.frame(matrix(data=0,ncol=3,nrow=nrow(d)))
  tlg[,1] <- as.numeric(rep(40,nrow(tlg)))
  tlg[,2] <- paste0(d[,1],d[,2]) # instead of newpos give old 
  # these columns in the ordered output are:
  # ex 111000 -> paternal haplotypes 111 | maternal haplotyes 000 of 3 kids
  # and then 2 columns with 0/1 if CO happened in Male or Female 
  s1 <- seq(7,90,1) # all such columns
  s2 <- seq(7,88,3) # family info column index ex 111000 above 
  fam.co <- s1[! s1 %in% s2] # the parental crossover columns
  # make two dataframe with position and CO in each parent in each family (cos)
  # and one with the family phase in each kid in each family
  crossovers <- d[,fam.co] # crossovers info and positions 
  haplotypes <- d[,s2] # family phase info and positions
  
  # make numeric 
  for(i in 1:ncol(crossovers)){ 
    crossovers[,i] <- as.numeric(crossovers[,i])
  }
  # either keep only male or female crossovers 
  if(male){
    crossovers <- crossovers[,seq(1,ncol(crossovers),2)] # keep only male COs 
    tlg[,3] <- as.numeric(rowSums(crossovers))
    lm20$maleCOs <- tlg[match(paste0(lm20$ss,lm20$pos),tlg[,2]),3]
  }else { 
    crossovers <- crossovers[,seq(2,ncol(crossovers),2)] # keep only female COs 
    tlg[,3] <- as.numeric(rowSums(crossovers))
    lm20$femaleCOs <- tlg[match(paste0(lm20$ss,lm20$pos),tlg[,2]),3]
    
  }
}

# set non PAR female COs to 0 
lm20$femaleCOs<- ifelse(lm20$f==0,0,lm20$femaleCOs)
# merge 
final <- rbind.data.frame(final,lm20)
# set all crossovers to 0 in pruned loci
final$maleCOs <- ifelse(is.na(final$m),NA,final$maleCOs)
final$femaleCOs <- ifelse(is.na(final$f),NA,final$femaleCOs)

############### WRITE ################

write.table(final,'../1.data/1.LepMAP/R.AllSNPs_physical/physical_orders_full_COs.txt',quote=F,row.names = F)


