#!/usr/bin/env Rscript
### 0 - Read arguments and sanity check 
args = commandArgs(trailingOnly=TRUE)
if(length(args) != 2 ) { 
	stop('Please specifiy ### 1 - the path and name of the total.theta.intersect.bed file ### 2 - the window bed file')
}
input <- args[1]  
input2 <- args[2]
### 1 - Library and file reading 
library(dplyr)
library(readr)
options(scipen=99999)
# read intersection file 
inters <- read_delim(input,
                     col_names = c('py_ss','py_start','py_end', 'py_theta',
                                   'n_ss','n_start','n_end','overlap'))
# read window file to create the new data frame 
nf <- read_delim(input2,
                 col_names = c('ss','start','end'))

### 2 - Calculations 
# make unique column for all
nf$uni <- paste(nf$ss, nf$start, '_') 
inters$nf_uni <- paste(inters$n_ss, inters$n_start, '_')
inters$thetaxN <- inters$py_theta * inters$overlap
# make sum of theta x N for each interval and add to nf df
tmp <- aggregate( inters$thetaxN, by=list(inters$nf_uni), FUN=sum)
# paste into window dataframe
nf$thetaxN <- tmp$x[ match(nf$uni, tmp$Group.1) ]
nf <- nf[,-4]
### 3 - output
# get window size to separate outputs of different Window sizes
w <- nf$end[1]-nf$start[1]
write_delim(nf,paste0(input,'_',w,'.thetaxN'),
            quote = 'none')

### END ### 
