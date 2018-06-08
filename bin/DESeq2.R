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
################################################################################
####################  Functions  ###############################################
################################################################################

## This Function will allow us to make a MA plot from a list of results objects

maPlot.lists <- function(x,i) {
  pdf(paste(i, '_maPlot.pdf', sep=''))
  plotMA(x, main=paste(i, 'alpha=0.01', sep=' '), 
         alpha=0.01, ylim=c(-6,6))
  abline(h=c(2,-2), col='red')
  dev.off()
}


##function for adding the means of depth norm counts per group the DEtable 
addMeans<-function(x,comp) {
  detable<-as.data.frame(x)
  
  ##get the two group names
  groups<-unlist(strsplit(comp,"vs"))
  groups<-paste(groups,"mean",sep=".")
  means<-count_means[,groups]
  
  ##match the gene names in the mean df with the results df
  means<-means[match(rownames(detable),rownames(means)),]
  detable<-cbind(detable,means)
  return(detable)}




################################################################################
args<-commandArgs()

count_file<-str_split(args[grep("counts",args)],"=",simplify=T)[2]
sample_file<-str_split(args[grep("samples",args)],"=",simplify=T)[2]
comparisons_file<-str_split(args[grep("comparisons",args)],"=",simplify=T)[2]

if (grepl("featureCounts", count_file)) {
counts<-read.delim(count_file,header=TRUE,stringsAsFactors=F,skip=1,check.names=F) 
rownames(counts)<-counts$Geneid
counts<-counts[,-c(1:6)]
 colnames(counts)<-gsub("_filtered_sortedByCoord.out.bam|Aligned.*","",colnames(counts))
} else {
counts<-read.csv(count_file,header=TRUE,stringsAsFactors=F,row.names=1)
}
head(counts)
#counts<-read.csv(count_file,header=TRUE,stringsAsFactors=F,row.names=1)
sample_info<-read.delim(sample_file,header=T,stringsAsFactors=F,col.names=c("filePrefix","coreNumber","sampleName","Treat"))
comparisons<-read.delim(comparisons_file,header=F,stringsAsFactors=F,col.names=c("treat","control"))
head(counts)

####################  Algorithm   ##############################################

## Create data frame containing experimental design.  
## This one will be set up with both condition and replicate variables.
## These variable names should be edited to match your experiment.
## Order matters!  Check your file_list variable!  
colData <- data.frame(row.names=sample_info$filePrefix,
                      condition=sample_info$Treat,
                      replicate=sample_info$sampleName)

head(colData)
##match rownames of colData with the colnames of count matrix
##check if rownames and colnames match

colData<-colData[match(colnames(counts),rownames(colData)),]
## Create DESeqDataSet
## Consolidates the data into an object that can be utilized by the program for 
## calculations and plotting.  Requires you to specify the experimental design 
dds <- DESeqDataSetFromMatrix(countData=counts,
                              colData=colData,
                              design= ~ condition)


## Ensure that reference level is set to the control before calling the 
## algorithm on the dds object.  

## Call algorithm on the assembled data set
dds <- DESeq(dds)

## Plot the dispersion of the experiment to verify that the Algorithm's 
## assumptions are valid for this dataset.  This will also show us if 
## the variance is too LOW in the samples indicating an error in replication
pdf('dispModel.pdf')
plotDispEsts(dds)
dev.off()

####################  Count and FPKM data   ####################################

## Retrieve count data and clean up the data frame 
count_data <- counts(dds, normalized=T)
colnames(count_data) <- colData$replicate 

## We need the gene lengths to calculte fpkm.  
## the following code will retrieve and calculate the gene lengths from the gtf
## file then add them to the dds object.  This may take some time...

# gr<-import('Assembly/merged.annotated.gtf')
# txdb<-makeTxDbFromGRanges(gr)
# txdb <- makeTxDbFromGFF('Assembly/merged.gtf')
# exons.list.per.gene <- exonsBy(txdb, by='gene')
# exonic.gene.sizes <- lapply(exons.list.per.gene,
#                             function(x){sum(width(reduce(x)))})
# gene_size <- data.frame(unlist(exonic.gene.sizes))
# gene_size <- data.frame(cbind(gene_size, 
#                               kb=apply(gene_size, 1, 
#                                        function(x){x/1000})))
# colnames(gene_size)[1] <- 'length'
# ## add to dds object
# mcols(dds)$basepairs <- gene_size[,'length']
# ## Get fpkm data
# fpkm_data <- fpkm(dds)
# colnames(fpkm_data) <- colnames(count_data)
# 
# ## Here we will add mean counts for each sample to the count and fpkm data frames.
## The is a purely convience operation to make the output tables more useful.
## Then write those tables to file.
## You will need to match these column names to your condition variable.
## Additional lines will need to be added if you have more than two conditions

count_means<-t(apply(count_data,1,function(x) tapply(x,colData(dds)$condition,mean,na.rm=T)))
head(count_means)
##count means and counts should be in the same order, but just in case....
count_means<-count_means[match(rownames(count_data),rownames(count_means)),]
colnames(count_means) <- paste(colnames(count_means),"mean",sep=".")

count_data_with_means <- data.frame(cbind(count_data,count_means))
# fpkm_data <- data.frame(cbind(fpkm_data, 
#                               CtlMean=apply(fpkm_data[,1:2],1,mean),
#                               TreatedMean=apply(fpkm_data[,3:4],1,mean)))
# 
write.csv(count_data, file='depthNormCount.csv', quote=F)
# write.csv(fpkm_data, file='fpkmData_151101.csv', quote=F)

####################  Conditional QC   ########################################

## Perform rLog transformation on the data 
rld <- rlog(dds)

##save the rld
rld_assay<-as.data.frame(assay(rld))
rld_assay$gene_id<-rownames(rld_assay)
rld_assay<-select(rld_assay,gene_id,everything())
write.csv(rld_assay,file="rld.csv",row.names=F,quote=F)


hmcol <- colorRampPalette(brewer.pal(9, 'GnBu'))(100)

##sample hm
dists <- dist(t(assay(rld)))
mat <- as.matrix(dists)
rownames(mat) <- colnames(mat) <- as.character(colData$replicate)
hc <- hclust(dists)

## Plot HM and save to file.  Run without the pdf() and dev.off() fxns to plot
## in the R window before saving to file.  There are several other file type 
## options besides pdf  
pdf('distClustering.pdf')
hm <- heatmap.2(mat, Rowv=as.dendrogram(hc), symm=T, trace='none',
          col=rev(hmcol), margin=c(13,13))
dev.off()

###save distance matrix for multiqc
write.csv(hm$carpet,file = "sample_distance_matrix_for_mq.csv",quote = F)

## a PCA plot is another way to look at the data in a similar way.  Here we 
## plot one and save it to file.
pdf('PCA.pdf')
pca <- plotPCA(rld, intgroup=c('condition'),returnData = T)
dev.off()

col<-brewer.pal(n=3,"Set1")

pcaForMQC<-function(pca) {
pca$color<-col[as.numeric(pca$condition)]
##make the data section for yaml config file
## 4  spaces before line
## and a space between colon and value
paste0('    ',pca$name,': {x: ',pca$PC1,', y: ',pca$PC2,', color: "',pca$color,'"}')
}

for_mqc<-pcaForMQC(pca)

write.table(for_mqc,file = "pca_for_mq.yaml",quote =F,row.names = F, col.names = F)
####################   DE Tests    #############################################

## Set up the results comparisons. 
## first make a list of contrasts using the comparisons file
## contrasts NEED to be in a list even for one contrasts

contrasts<-lapply(seq(1:nrow(comparisons)),function(x) {
  c("condition",comparisons[x,"treat"],comparisons[x,"control"])
})
names(contrasts)<-paste(comparisons$treat,comparisons$control,sep="vs")

resultsNames(dds)
Res <- lapply(contrasts,function(x) results(dds,contrast=x,alpha = 0.05))
Map(maPlot.lists,Res,names(Res))


## Now let's order the result object, subset by adjusted p-value, print a 
## summary and rearrange things a bit to save these results to file
Res <- lapply(Res,function(x) x[order(x$padj), ])




Res_with_means<-Map(addMeans,Res,names(Res))

###export all results
file_names<-paste0(names(Res_with_means),"_DEtable_ALL_genes.csv")

Map(write.csv,Res_with_means,file=file_names)


## get significant genes
sigRes <- lapply(Res_with_means,function(x) subset(x, padj<=0.05))

###export Sig results
file_names<-paste0(names(sigRes),"_DEtable_SIG_genes.csv")

Map(write.csv,sigRes,file=file_names)


##export a summary of the DEGS
summary<-data.frame(Down=sapply(sigRes,function(x) nrow(subset(x,log2FoldChange < 0))),
                    Up=sapply(sigRes,function(x) nrow(subset(x,log2FoldChange > 0))),
                    Total=sapply(sigRes,nrow))    
                                 
summary$Comparison<-names(sigRes)
summary<-select(summary,Comparison,everything())
write.table(summary,file="DE_summary.txt",sep="\t",quote=F,row.names=F)

##may meed to format differently for multiqc
save(dds, file="DESeqDataSet.rda")

#save sessionInfo()
sink(file="SessionInfo.txt")
sessionInfo()
sink()

