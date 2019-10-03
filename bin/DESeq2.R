#!/usr/bin/env Rscript

################################################################################
#####################  Packages  ###############################################
################################################################################


require(DESeq2) || stop("The DESeq2 library is not available!")
require(gplots) || stop("The gplots library is not available!")
require(RColorBrewer) || stop("The RcolorBrewer library is not available!")
#require(reshape2) || stop("The reshape2 library is not available!")
require(GenomicFeatures) || stop("The GenomicFeatures library is not available!")
require(ggplot2) || stop("The ggplot2 library is not available!")
require(tidyverse) || stop ("The tidyverse library is not available")
require(scales) || stop("The scales library is not available")
################################################################################
####################  Functions  ###############################################
################################################################################

## This Function will allow us to make a MA plot from a list of results objects

maPlot.lists <- function(x,i,file) {
  pdf(file = file)
  plotMA(x, main=i,alpha=0.05, ylim=c(-6,6))
  abline(h=c(2,-2), col='red')
  dev.off()
}


##function for adding the means of depth norm counts ,or fpkm per group the DEtable 
addMeans<-function(x,means,comp) {

  detable<-as.data.frame(x)
  
  ##get the two group names
  groups<-unlist(strsplit(comp,"vs"))
  group_means<-means[,groups]
  
  ##match the gene names in the mean df with the results df
  group_means<-group_means[match(rownames(detable),rownames(group_means)),]
  detable<-cbind(detable,group_means)
  return(detable)}


## modified DESeq2 plotPCA function, to plot addtional pc dimensions
plotPCA2 <- function (object, intgroup = "condition", ntop = 500, dim1= 1,dim2 = 2,returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame( pca$x[, dim1],pca$x[, dim2], group = group, 
                  intgroup.df, name = colnames(object))
  dim_names <- paste0("PC",c(dim1,dim2))
  colnames(d)[1:2] <- dim_names
  if (returnData) {
    attr(d, "percentVar") <- percentVar[c(dim1,dim2)]
    return(d)
  }
  ggplot(data = d, aes_string(x = dim_names[1], y = dim_names[2], color = "group")) + 
    geom_point(size = 3) + xlab(paste0(dim_names[1],": ", round(percentVar[dim1] * 
                                                        100), "% variance")) + ylab(paste0(dim_names[2],": ", round(percentVar[dim2] * 
                                                                                                            100), "% variance")) + coord_fixed()
}



################################################################################
theme_set(theme_bw())
args<-commandArgs()

args
count_file<-str_split(args[grep("--counts",args)],"=",simplify=T)[2]
sample_file<-str_split(args[grep("--samples",args)],"=",simplify=T)[2]
comparisons_file<-str_split(args[grep("--comparisons",args)],"=",simplify=T)[2]
gene_lengths_file <-str_split(args[grep("gene_lengths",args)],"=",simplify=T)[2]
out<-str_split(args[grep("--out",args)],"=",simplify=T)[2]
species<-str_split(args[grep("--species",args)],"=",simplify=T)[2]
fpkms <- str_split(args[grep("--fpkms",args)],"=",simplify=T)[2]


counts<-read.csv(count_file,header=TRUE,stringsAsFactors=F,row.names=1,check.names = F)

sample_info<-read.delim(sample_file,header=F,stringsAsFactors=F)
gene_lengths <- read.csv(gene_lengths_file,header=T,stringsAsFactor=F,row.names = 1)

##adds column names for single end or paired end sample files
if(ncol(sample_info) == 4) {
colnames(sample_info) <- c("sampleName","coreNumber","Group","read1")
} else {
colnames(sample_info) <- c("sampleName","coreNumber","Group","read1","read2")
}
comparisons<-read.delim(comparisons_file,header=F,stringsAsFactors=F,col.names=c("treat","control"))


####################  Algorithm   ##############################################

## Create data frame containing experimental design.  

colData <- data.frame(row.names=sample_info$coreNumber,
                      condition=sample_info$Group,
                      replicate=sample_info$sampleName)


##match rownames of colData with the colnames of count matrix
##check if rownames and colnames match

counts <- counts[,match(rownames(colData),colnames(counts))]
## Create DESeqDataSet
## Consolidates the data into an object that can be utilized by the program for 
## calculations and plotting.  Requires you to specify the experimental design 
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=colData,
                              design= ~ condition)


##add annotation, make sure the rownames of counts and annotation match although they should
gene_lengths <- gene_lengths[match(rownames(counts),rownames(gene_lengths)),]
colnames(gene_lengths)[5] <- "basepairs"

mcols(dds) <- DataFrame(gene_lengths)
## Call algorithm on the assembled data set
dds <- DESeq(dds,betaPrior= T)

## Plot the dispersion of the experiment to verify that the Algorithm's 
## assumptions are valid for this dataset.  This will also show us if 
## the variance is too LOW in the samples indicating an error in replication

pdf(file=file.path(out,'dispModel.pdf'))

plotDispEsts(dds)
dev.off()

####################  Count and FPKM data   ####################################

## Retrieve count data and clean up the data frame 
count_data <- counts(dds, normalized=T)
colnames(count_data) <- colData$replicate 
# 
# ## Here we will add mean counts for each sample to the count dataframe.
## The is a purely convience operation to make the output tables more useful.
## Then write those tables to file.s

count_means<-t(apply(count_data,1,function(x) tapply(x,colData(dds)$condition,mean,na.rm=T)))



##count means and counts should be in the same order, but just in case....
count_means<-count_means[match(rownames(count_data),rownames(count_means)),]

count_data_with_means <- data.frame(cbind(count_data,count_means),check.names=T)

write.csv(count_data_with_means, file=file.path(out,'depthNormCount.csv'), quote=F)


###This code uses fpkms from DESeq2
####same with fpkm 
#fpkm_data <- fpkm(dds, robust=T)
#colnames(fpkm_data) <- colData$replicate
## 
## ## Here we will add mean counts for each sample to the count dataframe.
### The is a purely convience operation to make the output tables more useful.
### Then write those tables to file.s
#
#fpkm_means<-t(apply(fpkm_data,1,function(x) tapply(x,colData(dds)$condition,mean,na.rm=T)))
#
#
###fpkm means and counts should be in the same order, but just in case....
#fpkm_means<-fpkm_means[match(rownames(fpkm_data),rownames(fpkm_means)),]
#
#fpkm_data_with_means <- data.frame(cbind(fpkm_data,fpkm_means))
#
#write.csv(fpkm_data_with_means, file=file.path(out,'fpkm_values.csv'), quote=F)

#####Currently using fpkms from cufflinks

##import table from cufflinks
fpkms <- read.csv(fpkms,header=T,stringsAsFactors=F,check.names=F)
idx <- fpkms$GeneID
fpkms <- as.matrix(fpkms[,-1])
rownames(fpkms) <- idx

##first match the gene names in dds object to cufflinks fpkm table
fpkm_data<-fpkms[match(rownames(count_data),rownames(fpkms)),]

##check that fpkm table is in same order as colData
fpkm_data <- fpkm_data[,match(rownames(colData(dds)),colnames(fpkm_data))]
colnames(fpkm_data) <- colData$replicate
dim(fpkm_data)
fpkm_means<-t(apply(fpkm_data,1,function(x) tapply(x,colData(dds)$condition,mean,na.rm=T)))


#fpkm means and counts should be in the same order, but just in case....
fpkm_means<-fpkm_means[match(rownames(fpkm_data),rownames(fpkm_means)),]

fpkm_data_with_means <- data.frame(cbind(fpkm_data,fpkm_means),check.names =F)

write.csv(fpkm_data_with_means, file=file.path(out,'fpkm_values.csv'), quote=F)


####################  QC   ########################################

## Perform rLog transformation on the data 
rld <- rlog(dds)

##save the rld
save(rld,file=file.path(out,"rld.rda"))



hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)

##sample hm
dists <- dist(t(assay(rld)))
mat <- as.matrix(dists)
rownames(mat) <- colnames(mat) <- as.character(colData$replicate)
hc <- hclust(dists)

## Plot HM and save to file.
pdf(file=file.path(out,'distClustering.pdf'))
hm <- heatmap.2(mat, Rowv=as.dendrogram(hc), symm=T, trace='none',
          col=rev(hmcol), margin=c(13,13))
dev.off()

###save distance matrix for multiqc
write.csv(hm$carpet,file = file.path(out,"sample_distance_matrix_for_mq.csv"),quote = F)

## a PCA plot is another way to look at the data in a similar way.  Here we 
## plot one and save it to file.
pdf(file=file.path(out,'PCA_pc1vspc2.pdf'))

plotPCA(rld, intgroup=c('condition'))
dev.off()

pca <- plotPCA(rld, intgroup=c('condition'),returnData = T)
col<-hue_pal()(nlevels(as.factor(pca$condition)))


pcaForMQC<-function(pca) {
pca$color<-col[as.numeric(pca$condition)]
##make the data section for yaml config file
## 4  spaces before line
## and a space between colon and value
paste0('    ',pca$name,': {x: ',pca$PC1,', y: ',pca$PC2,', color: "',pca$color,'"}')
}

for_mqc<-pcaForMQC(pca)


write.table(for_mqc,file = file.path(out,"pca_for_mq.yaml"),quote =F,row.names = F, col.names = F)

##plot additional dimensions
pdf(file=file.path(out,'PCA_pc1vspc3.pdf'))
plotPCA2(rld, dim1=1,dim2=3,intgroup=c('condition'))
dev.off()


pdf(file=file.path(out,'PCA_pc2vspc3.pdf'))
plotPCA2(rld, dim1=2,dim2=3,intgroup=c('condition'))
dev.off()

####################   DE Tests    #############################################

## Set up the results comparisons. 
## first make a list of contrasts using the comparisons file
## contrasts NEED to be in a list even for one contrasts

contrasts<-lapply(seq(1:nrow(comparisons)),function(x) {
  c("condition",comparisons[x,"treat"],comparisons[x,"control"])
})
names(contrasts)<-paste(comparisons$treat,comparisons$control,sep="vs")


Res <- lapply(contrasts,function(x) results(dds,contrast=x,alpha = 0.05))
file_names <- file.path(out,paste0(names(Res),"_maPlot.pdf"))
Map(maPlot.lists,Res,names(Res),file_names)


## Now let's order the result object, subset by adjusted p-value, print a 
## summary and rearrange things a bit to save these results to file
Res <- lapply(Res,function(x) x[order(x$padj), ])

Res_with_means<- Map(function(x,means,comp) {
        dat <- addMeans(x,means,comp)
        colnames(dat)[7:8] <- paste(colnames(dat)[7:8],"count.mean",sep = ".")
        return(dat)}
,x=Res,comp=names(Res),MoreArgs=list(means = count_means))



Res_with_means_fpkm <- Map(function(x,means,comp) {
        dat <- addMeans(x,means,comp)
        colnames(dat)[9:10] <- paste(colnames(dat)[9:10],"fpkm.mean",sep = ".")
        return(dat)}
,x=Res_with_means,comp=names(Res),MoreArgs=list(means = fpkm_means))


## make the gene ids a column
Res_with_means_fpkm <- lapply(Res_with_means_fpkm,function(x) {
		x$GeneID <- rownames(x)
		dat <- select(x,GeneID,everything())
		return(dat)
		})

###add annotation
if(species == "mm10"){
annot_file <- "/lower_bay/local_storage/annotation_db/Mus_musculus/UCSC/mm10/Annotation/Genes/mm10_gene_descriptions.txt"
} else if (species == "hg19") {
annot_file <- "/lower_bay/local_storage/annotation_db/Homo_sapiens/UCSC/hg19/Annotation/Genes/hg19_gene_descriptions.txt"
}

gene_descriptions <- read.delim(annot_file,header=T,stringsAsFactors=F)
head(gene_descriptions)

Res_with_annotations <- lapply(Res_with_means_fpkm,function(x) {
			dplyr::left_join(x,gene_descriptions,by=c("GeneID"="external_gene_name"))
			})

lapply(Res_with_annotations,head)
###export all results
file_names<-file.path(out,paste0(names(Res_with_annotations),"_DEtable_ALL_genes.csv"))

Map(write.csv,Res_with_annotations,file=file_names,MoreArgs = list(row.names = F))


## get significant genes
sigRes <- lapply(Res_with_annotations,function(x) subset(x, padj<=0.05))

##export a summary of the DEGS
summary<-data.frame(Down=sapply(sigRes,function(x) nrow(subset(x,log2FoldChange < 0))),
                    Up=sapply(sigRes,function(x) nrow(subset(x,log2FoldChange > 0))),
                    Total=sapply(sigRes,nrow))    
                                 
summary$Comparison<-names(sigRes)
summary <- select(summary,Comparison,everything())
write.table(summary,file=file.path(out,"DE_summary.txt"),sep="\t",quote=F,row.names=F)


###export Sig results
##remove any data frames with no results
sigRes <- sigRes[sapply(sigRes,nrow)>0]


if (length(sigRes)>0) {
file_names <- file.path(out,paste0(names(sigRes),"_DEtable_SIG_genes.csv"))

Map(write.csv,sigRes,file=file_names,MoreArgs = list(row.names = F))
}

## save DESeqDataSet
save(dds, file=file.path(out,"DESeqDataSet.rda"))

#save sessionInfo()
sink(file=file.path(out,"SessionInfo.txt"))
sessionInfo()
sink()

