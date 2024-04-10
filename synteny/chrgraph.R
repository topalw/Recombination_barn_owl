# invoke as :
# Rscript --vanilla .\tiny\utils\chrgraph.R <ref> <uniprotspeciescode> ...
# if none is provided, use CHICK

ref <- "TYTAL"
args <- "CHICK" # which species ? - always for one - 
library(tidyverse)
library(ggplotlyExtra)
library(plotly)

order <- read.csv('new_order.csv',h=F)
names(order) <- c('lg','ss','chicken')
l = args[1] 
  chrinfo <- read_delim("chrinfo.2.txt", delim = "\t")
  chrinfo <- group_by(chrinfo, speciescode) %>% arrange(krank) %>% mutate(offset = lag(seqlength), offset = replace_na(offset, 0), offset = cumsum(offset), xoffset = offset + ((krank - 1) * 2000000))
  dnazoo <- read_delim("dnazoo.scfinfo.txt", delim = "\t")
  dnazoo$krank[dnazoo$krank >= 24] <- dnazoo$krank[dnazoo$krank >= 24] + 1 
  dnazoo <- dnazoo %>% add_row(seqid = 'Super-Scaffold_2.2',
                               seqlength = 63464669, 
                               chrtype = 'macro' , 
                               krank = 24, 
                               speciescode = 'TYTAL', 
                               chrname = 's2.2',
                               cname = 'Tyto alba')
  dnazoo$seqlength[dnazoo$seqid=='Super-Scaffold_2'] <- 28222627
  renamer <- read.csv('../2.metadata/supp_table_lgs.csv')
  renamer[1,2] <- '2.2'
  r2 <- cbind.data.frame('old'=gsub('^','s',renamer[,2]),
                   'new'=gsub('^','LG',renamer$Linkage.group))
  r2$old[r2$new == 'LG34'] <- 's2042'
  r2$old[r2$new == 'LG14'] <- 's106'
  dnazoo$chrname[match(r2$old,dnazoo$chrname)] <- r2$new
  # reorder 
  dnazoo$krank[match(order$lg,dnazoo$chrname)] <- c(seq(1,22),seq(25,nrow(order)+2))
  # select useful lgs 
  #good <- c(which(dnazoo$chrname %in% r2$new), which(dnazoo$chrname %in% c('s13','s42') & dnazoo$speciescode =='TYTAL'))
  #dnazoo <- dnazoo[good,]
  
  # give useless scaffolds a higher krank (not present in linkage groups) to avoid plotting gaps 
  useless <- which(! dnazoo$chrname %in% c(r2$old,r2$new,c('s42','s13')) & dnazoo$speciescode =='TYTAL')
  dnazoo$krank[useless] <- seq(nrow(order)+3,nrow(order)+2 + length(useless)) 

  
  # make krank and offset based on chr
  dnazoo <- group_by(dnazoo, speciescode) %>% 
    arrange(krank) %>% 
    mutate(offset = lag(seqlength), offset = replace_na(offset, 0), offset = cumsum(offset), xoffset = offset + ((krank - 1) * 2000000))
  chrinfo <- bind_rows(dnazoo, select(chrinfo, colnames(dnazoo)))
  # species info 
  speciesinfo.orig <- read_delim("species.txt", delim = "\t")
  speciesinfo <- tribble(~uniprotspeciescode,l,ref)
  speciesinfo <- left_join(speciesinfo, select(speciesinfo.orig, uniprotspeciescode, "Common name"), by = "uniprotspeciescode")
  speciesinfo <- speciesinfo %>% mutate(rphyloorder = row_number())
  chrinfo <- left_join(chrinfo, select(speciesinfo, uniprotspeciescode, rphyloorder), by = c("speciescode"="uniprotspeciescode"))
  chrinfo <- filter(chrinfo, ! is.na(rphyloorder))
  # read alignment chainpairs 
  minalnlen <- 100000 #100kb min length
  mergedaln <- list()
  for (i in 1:2){
    if (i == 1){
      qcode <- l
      tcode <- ref
      aln <- paste0(qcode,".query.chainpairs.tab")
    }else{
      qcode <- ref
      tcode <- l
      aln <- paste0(tcode,".target.chainpairs.tab")
      # paste third line and update rphyloorder
      chrkept <- chrinfo %>% filter(rphyloorder == 1)
      chrinfo <- chrinfo %>% mutate(rphyloorder = replace(rphyloorder, rphyloorder == 1, 3))
      speciesinfo <- speciesinfo %>% add_row(speciesinfo[1,])
      speciesinfo[3,3] = 3
    }
    print(paste(qcode, tcode, aln))
    b <- NULL
    b <- read_delim(aln, "\t",col_names = F)
    colnames(b) <- c("tacc", "tstart", "tend", "qacc", "qstart", "qend", "strand")
      if(i==2){
        # remove ss2 
        b <- b[b$qacc != 'Super-Scaffold_2',] 
        # add new 
        tmp <- read_delim('chick.v.ss2.query.chainpairs.tab',
                          col_names = c("tacc", "tstart", "tend", "qacc", "qstart", "qend", "strand"))
        b <- rbind(tmp,b)
      } else { 
        b <- b[b$tacc != 'Super-Scaffold_2',] 
        tmp <- read_delim('chick.v.ss2.target.chainpairs.tab',
                         col_names = c("tacc", "tstart", "tend", "qacc", "qstart", "qend", "strand"))
        b <- rbind(tmp,b)
      }
    b <- left_join(b, dplyr::select(chrinfo, tacc = seqid, tchr = chrname, trank = krank, toffset=xoffset, torder = rphyloorder, tchrtype = chrtype, tcode = speciescode), by = "tacc") 
    b <- left_join(b, dplyr::select(chrinfo, qacc = seqid, qchr = chrname, qrank = krank, qoffset=xoffset, qorder = rphyloorder, qchrtype = chrtype, qcode = speciescode), by = "qacc")
    mergedaln[[i]] <- b 
  }
  b <- bind_rows(mergedaln)
  # extract Z 
  z <- b[b$qchr %in% c('s13','s42') | b$tchr %in% c('s13','s42'),]
  # extract linkage groups
  b <- b[ifelse(b$tcode == 'TYTAL', grep('LG',b$tchr), grep('LG',b$qchr)),]
  b <- rbind(b,z)
  c <- mutate(b, tlen=tend-tstart, qlen=qend-qstart) %>% 
    filter((qlen >= minalnlen | tlen >= minalnlen) & ! is.na(qchr) & ! is.na(tchr)) %>% 
    mutate(alnid = paste("aln", row_number(), sep = "")) %>%
    mutate(alntype = paste(tchrtype, qchrtype, sep = ":"), x1 = tstart, x2 = tend, x1 = x1 + toffset, x2 = x2 + toffset, x3 = if_else(strand == '-', qstart + qoffset, qend + qoffset), x4 = if_else(strand == '-', qend + qoffset, qstart + qoffset), y1 = torder - 0.01, y2 = torder - 0.01, y3 = qorder + 0.16, y4 = qorder + 0.16) %>%
    select(alnid,x1,x2,x3,x4,y1,y2,y3,y4,alntype)
  # flip assembly 
  tmpor <- which(c$x3 < c$x4)
  tmp <- c$x3[tmpor] 
  c$x3[tmpor] <- c$x4[tmpor]
  c$x4[tmpor] <- tmp
  
  t1 <- select(c, alnid, alntype, x1:x4) %>% pivot_longer(x1:x4, names_to = "ctype", values_to = "xcoords") %>% select(alnid, alntype, xcoords)
  t2 <- select(c, alnid, alntype, y1:y4) %>% pivot_longer(y1:y4, names_to = "ctype", values_to = "ycoords") %>% select(alnid, alntype, ycoords)
  t <- bind_cols(t1, t2)
  colnames(t) <- c("alnid","alntype", "xcoords","tmp1","tmp2", "ycoords")
  t <- t %>% mutate( alntype = case_when(alntype == "macro:macro" ~ "macro", alntype == "micro:micro" ~ "micro", alntype == "macro:micro" ~ "macro X micro", alntype == "micro:macro" ~ "macro X micro"))
  # which chr to plot 
  tytal <- chrinfo[chrinfo$speciescode == 'TYTAL',]
  z <- chrinfo[chrinfo$speciescode == 'TYTAL' & chrinfo$seqid %in% c('Super-Scaffold_13','Super-Scaffold_42'),]
  tmp <- chrinfo[chrinfo$speciescode != 'TYTAL',]
  tytal <- tytal[grep('LG',tytal$chrname),]
  chrinfo <- rbind(tmp,tytal,z)
  chrinfo <- bind_rows(chrkept,chrinfo)
  # colour Z 
  chrinfo$chrtype[chrinfo$chrname == 'chrZ'] <- 'Z'
  chrinfo$chrtype[chrinfo$chrname %in% c('s13','s42')] <- 'Z'
  chrinfo <- chrinfo[chrinfo$rphyloorder!=3,] # remove top row 
  t <- t[t$ycoords < 2,] # keep only barn owl to chicken 
  # plot 
  p <- ggplot() + 
    geom_rect(data = chrinfo, aes(xmin=xoffset, xmax=xoffset + seqlength, ymin=rphyloorder,ymax= rphyloorder+0.15, fill=chrtype), alpha = 0.4) +
    geom_text(data = chrinfo, aes(x=xoffset, y = jitter(rphyloorder + 0.075,amount=0.01), label = gsub('chr','',gsub('LG','',chrinfo$chrname))), hjust = 0, size = 2.5) + 
    geom_polygon(data = t, aes(x=xcoords,y=ycoords, fill = alntype, group = alnid), alpha = 0.4) +
    scale_y_continuous(breaks = 1:nrow(speciesinfo[c(1,2),])+ 0.075, labels = arrange(speciesinfo[c(1,2),], rphyloorder[c(1,2)]+0.075) %>% pull("Common name")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), 
	  text = element_text(size = 16), 
	  axis.title = element_blank(),
	  axis.ticks.x = element_blank(),
	  axis.text.x = element_blank())

 print(p)
  
pdf(paste0(ref,'vs',l,'.pdf'),width=12,height=10)
print(p) 	
dev.off()
 
jpeg(paste0(ref,'vs',l,'.jpeg'),width=12,height=10,unit='in',quality = 1000,res=800, antialias = "cleartype" )
print(p)
dev.off()

 
 
 #  htmlwidgets::saveWidget(partial_bundle(ggplotly(p)), paste0(ref,"vs",l,".html"))
