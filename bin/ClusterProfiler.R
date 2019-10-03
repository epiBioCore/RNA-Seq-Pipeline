#!/usr/bin/env Rscript

library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(rWikiPathways)
library(tidyverse)
library(SPIA)


theme_set(theme_bw())
args<-commandArgs()
args
DE_dir<-str_split(args[grep("DE",args)],"=",simplify=T)[2]
org <- str_split(args[grep("org",args)],"=",simplify=T)[2]
out <- str_split(args[grep("out",args)],"=",simplify=T)[2]

# obtaining data
lx=list.files(pattern ="DEtable_ALL_genes.csv",path=DE_dir,full.names=T)
DE=lapply(lx,read.csv)
names(DE)=gsub("_DEtable_ALL_genes.csv","",basename(lx))

# orgDb

if(org=="mm10") {
  orgDb <- "org.Mm.eg.db"
  kegg_ord <- "Mmu"
  organism <- "Mus Musculus"
}else if (org == "hg19") {
  orgDb <-"org.Hs.eg.db"
  kegg_ord <- "Hsa"
  organism <- "Homo Sapiens"
} else {
  stop("Stop! Either no organism provided or organism is not supported.")
}




#Get gene lists of DE genes with entrezid added

genelists <- map2(DE,names(DE),function(i,x){
  ids_all= bitr(i$GeneID, fromType="SYMBOL", toType=c("ENTREZID", "ENSEMBL"), OrgDb=orgDb)
  
  merge(i,ids_all,by.x="GeneID",by.y="SYMBOL")
})

sig_genes <-map(genelists,~subset(.x,padj <= 0.05))




###############GO term enrichment  for all three onologies, BP, CC, MF

ontologies <- c("CC","BP","MF")
GO_enrich=map2(genelists,sig_genes, function(genelist,sig_genes){
  map(ontologies,function(ont) {
    
    enrichGO(gene =sig_genes$ENTREZID,
             universe =genelists$ENTREZID,
             OrgDb=orgDb, ont=ont,
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.01,
             qvalueCutoff  = 0.05, readable =TRUE)
  })})

GO_enrich <- map(GO_enrich,function(x) {
  names(x) <- ontologies
  return(x)
})
GO_enrich_ul <- unlist(GO_enrich,recursive = F)

names(GO_enrich_ul) <- sub("(.{2}$)",paste0("GO_","\\1"),names(GO_enrich_ul))

map2(GO_enrich_ul,names(GO_enrich_ul),function(i,x){write.table(i,file=file.path(out,paste0(x,".txt")), sep="\t", quote=F,row.names=F)} )


#########KEGG pathway and module enrichment
###KEGG pathways


KEGG <-map(sig_genes,function(x){
  
  enrichKEGG(gene= x$ENTREZID,organism=kegg_org,pvalueCutoff = 0.05)
})

names(KEGG) <- paste0(names(KEGG),"_KEGG")
map2(KEGG,names(KEGG),function(i,x){write.table(i,file=file.path(out,paste0(x,".txt")), sep="\t", quote=F, row.names=F)} )


#KEGG module analysis

KEGG_mod=map(sig_genes,function(x){
  
  enrichMKEGG(gene = x$ENTREZID,organism = kegg_org)
  
})

names(KEGG_mod) <- paste0(names(KEGG_mod),"_KEGG_Module")
map2(KEGG_mod, names(KEGG_mod),function(i,x){write.table(i,file=file.path(out,paste0(x,".txt")), sep="\t", quote=F, row.names=F)} )

################ Gene set enrichment (GSEA)
GSEA=map(genelists,function(x){
  
  
  top<-x[!is.na(x$ENTREZID),]
  
  top<-top[!duplicated(top$ENTREZID),]
  
  GSEA_list=top$log2FoldChange
  
  names(GSEA_list)<-as.vector(top$ENTREZID)
  
  GSEA_list=sort(GSEA_list, decreasing=T)
  
  gseGO(geneList     = GSEA_list,
        
        OrgDb        = orgDb,
        
        ont          = "ALL",
        
        nPerm        = 1000,
        
        minGSSize    = 100,
        
        maxGSSize    = 500,
        
        pvalueCutoff = 0.05,
        
        verbose      = FALSE)
  
})

names(GSEA) <- paste0(names(GSEA),"_GSEA")
map2(GSEA,names(GSEA),function(i,x){write.table(i,file=file.path(out,paste0(x,".txt")), sep="\t", quote=F,row.names=F)} )




########### wikipathways
wp.Mm.gmt <- rWikiPathways::downloadPathwayArchive(organism="Mus musculus", format = "gmt")
wp2gene <- clusterProfiler::read.gmt(wp.Mm.gmt)
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid,gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid,name) #TERM2NAME

ewp <- map(sig_genes,~enricher(.x$ENTREZID, TERM2GENE = wpid2gene, TERM2NAME = wpid2name))

safe_setReadable <- safely(~setReadable(.,org.Mm.eg.db, keyType = "ENTREZID"))
ewp <- map(ewp,~safe_setReadable(.x))

ewp_to_keep <- map_lgl(ewp,function(x) is.null(x$error))

ewp<- ewp[ewp_to_keep]

ewp <- map(ewp,"result")
names(ewp) <- paste0(names(ewp),"_wikipathways")
map2(ewp,names(ewp),function(i,x){write.table(i,file=file.path(out,paste0(x,".txt")), sep="\t", quote=F,row.names=F)} )

############spia

safe_spia <- safely(~spia(de=.x,all=.y,organism=kegg_org,plots = F,beta = NULL,combine = "fisher",verbose = "FALSE"))
##to handle errors if there isn't enough DE genes i

top <- map(genelists, function(x) x[!is.na(x$ENTREZID),])
top<-map(top,~.x[!duplicated(.x$ENTREZID),])
tgl <- map(top,~subset(.x,padj<=0.05))
tgl <- tgl[sapply(tgl,nrow)>1]
DE <- map(tgl,"log2FoldChange")
DE <- map2(DE,tgl,function(d,t) {
  names(d) <- t$ENTREZID
  return(d)
})

ALL <- top[names(top) %in% names(tgl)]
ALL <- map(ALL,"ENTREZID")

res <- map2(DE,ALL,~safe_spia(de=.x,all=.y))
res_to_keep <- map_lgl(res,~is.null(.x$error))

res <- res[res_to_keep]
map2(res,names(res),function(r,name) {
  write.table(r$result,file = file.path(out,paste0(name,"_SPIA_enrichment.txt")),sep = "\t",row.names = F,quote = F)
})
###############Plotting
###first create "safe" functions to handle erros
safe_barplot <- safely(~barplot(.,showCategory=12))
safe_barplot2 <- safely(~barplot(.,showCategory=12,drop=T))
safe_emapplot <- safely(emapplot)
safe_cnetplot <- safely(cnetplot)
safe_dotplot <- safely(dotplot)

### then create a list the plotting functions 
safeplot_lists <- list(barplot=safe_barplot,bartplot2=safe_barplot2, dotplot=safe_dotplot,emapplot=safe_emapplot,cnetplot=safe_cnetplot)



### put all all objects to be plotted together.
plot_objects <- list(GO_enrich_ul,KEGG,KEGG_mod,GSEA,ewp)
plot_objects <- unlist(plot_objects,recursive = F)



all_the_plots <- map(safeplot_lists,function(f) {
  map(plot_objects,function(enrich) {
    exec(f,enrich)
  })})


all_the_plots <- transpose(all_the_plots)
all_the_plots <- unlist(all_the_plots,recursive = F)


plots_to_keep <- map_lgl(all_the_plots,function(x) is.null(x$error))

plots_filtered <- all_the_plots[plots_to_keep]

map2(plots_filtered,names(plots_filtered),function(x,y) {
  ggsave(x$result,file=file.path(out,paste0(y,".pdf")))
})


##################save all of the enrichment objects
map2(plot_objects,names(plot_objects), function(x,y) save(x,file=file.path(out,paste0(y,"_enrichment_object.rda"))))

sink(file=file.path(out,"SessionInfo.txt"))
sessionInfo()
sink()


