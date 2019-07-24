#!/usr/bin/env Rscript

library("tidyverse")
args<-commandArgs()
args
dir <- str_split(args[grep("^in",args)],"=",simplify=T)[2]

paths=list.dirs(path=dir,full.names=T)
tracking_files <- list.files(path=paths,pattern="genes.fpkm_tracking",full.names=T)
fpkms <- lapply(tracking_files,function(x) read.delim(x,header=T,stringsAsFactors=F))
names(fpkms) <- str_split(tracking_files,"/",simplify=T)[,3]
fpkm_df <- map2_df(fpkms,names(fpkms),function(x,y) {
dat <- x[,c("FPKM")]
return(dat)}
)

## add gene names- they're the same in all samples
fpkm_df$GeneID <- fpkms[[1]]$tracking_id
fpkm_df <- dplyr::select(fpkm_df,GeneID,everything())

##remove duplicate ids
fpkm_df <- distinct(fpkm_df)


write.csv(fpkm_df ,file=file.path(dir,"cufflink_fpkms_all_samples.csv"),row.names=F)
  
