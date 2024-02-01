### FUNCTION 1 - modify color transparency
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

### FUNCTION 2 - rename Super-Scaffold_## to ss##
rename.ss <- function(column){
  return(gsub('Super-Scaffold_','ss',column))
}

### FUNCTION 3 - make cM column from rho using Haldane's function
hald <- function(data){
  hald <- -0.5 * log(1 - (2*(data$thetaxN )), base=exp(1)) * 100
  return(hald)
}

rev.haldane <- function(cm){
  theta <- (1 - exp(-2*cm/100))/2
  return(theta)
}

### FUNCTION 4 - make windowed results from lepmap results  
window.ped <- function(windows, ped, sex.based = FALSE){
  if(sex.based){
    res <- data.frame('m'=rep(NA,nrow(windows)),
                      'f'=rep(NA,nrow(windows)))
    
  } else { res <- rep(NA,nrow(windows)) }
  for(i in 1:nrow(windows)){
    start <- windows$start[i]
    end <- windows$end[i]
    ss <- windows$ss[i]
    indx <- which(ped$pos < end & 
                    ped$pos >= start & 
                    ped$ss == ss)
    if(length(indx) > 1 ){
    if(sex.based){
      res$m[i] <- max(ped$m[indx], na.rm=T) - 
        min(ped$m[indx], na.rm=T)
      res$f[i] <- max(ped$f[indx], na.rm=T) - 
        min(ped$f[indx], na.rm=T)
    } else {
      res[i] <- max(ped$av[indx], na.rm=T) - 
        min(ped$av[indx], na.rm=T)
    }
    } else {next} 
  }
  return(res)
}


### FUNCTION 5 -  remove outlier markers from edges of LepMap3 orders 
prune.ends <- function(cMcol, nloc=20, cm.cutoff=2){
  nloc <- ifelse(nloc > length(cMcol),1 ,nloc) 
  # if the number of loci is bigger than the length then just take the first 

  # Indexes at start: 
  a1 <- max(which(diff(cMcol[1:nloc]) > cm.cutoff)) +1 
  # of those markers that show a jump bigger than cutoff, which one has the max pos
  # +1 because diff prints from 2
  
  # Indexes at end: 
  end.indx <- (length(cMcol)-nloc+1):length(cMcol) 
  # take the last nloc markers (sequence -> n-nloc +1  : n)
  b1 <- end.indx[which(abs(diff(cMcol[end.indx])) > cm.cutoff) ]
  # take the smallest of b1 (since it shows the first jump)
  
  # If -Inf give NA
  a <- ifelse(! is.na(a1 %% 1),a1,NA)[1] 
  b <- ifelse(! is.na(b1 %% 1), b1, NA)[1]  
  # index 1 is to remove similar values popping up

  # If NA give start or end
  a <- ifelse(! is.na(a), a, 1) 
  b <- ifelse(! is.na(b), b, length(cMcol) ) 

  # rescale to start at 0
  cMcol[a : b] <- cMcol[a : b]  - min(cMcol[a : b])
  
  # set all else to NA : 
  true.indx <- 1:length(cMcol)
  outside <- true.indx[! true.indx %in% seq(a,b)]
  cMcol[outside] <- NA
  
  # return pruned column
  return(cMcol)
}

cmmb <- function(row){
  return((row[3]/(row[2]-row[1])) * 1e6 )
}

# Missingness 
miss <- function(snp){
  return(sum(is.na(snp)))
}
maf <- function(snp){
  return(mean(snp,na.rm=T)/2)
}
het <- function(snp){
  return(sum(snp==1,na.rm=T) / length(!is.na(snp)))
}

# distance to end
distance <- function(dat){
  dat$dist <- NA
  for(ss in unique(dat$ss)){
    indx <- which(dat$ss == ss) # get scaffold 
    if(length(indx > 1)){ 
      dat$dist[indx] <- seq(0,length(indx)-1,1)/(length(indx)-1)
      #{ 0, 1, .... n-1 }/ n-1 = 0, 0.0..1, 1 
    }
  }
  # make 0 -> 1, 0 -> 0.5 -> 0
  dat$dist[which(dat$dist > 0.5)] <- abs(dat$dist[which(dat$dist > 0.5)] - 1) 
  return(dat$dist)
}
