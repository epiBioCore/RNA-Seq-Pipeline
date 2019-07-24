#!/usr/bin/env Rscript

library("stringr")
args<-commandArgs()
count_file<-str_split(args[grep("counts",args)],"=",simplify=T)[2]
out <- str_split(args[grep("--out",args)],"=",simplify=T)[2]
counts<-read.delim(count_file,header=TRUE,stringsAsFactors=F,skip=1,check.names=F)

rownames(counts)<-counts$Geneid
annotation <- counts[,2:6]
counts<-counts[,-c(1:6)]

colnames(counts)<-gsub("_properly_paired_sorted.bam|_sorted.bam","",basename(colnames(counts)))
write.csv(counts, file = file.path(out,"featureCounts_for_DESeq2.csv"))
write.csv(annotation,file = file.path(out,"gene_lengths.csv"))
