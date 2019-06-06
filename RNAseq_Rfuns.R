#*******************************************************************
#
#       Functions for RNASeq analysis with DESeq2 and other tools
#
#                       Vendula Svendova (2019)
#
#*******************************************************************
options(java.parameters = "- Xmx1024m") # for xlsx to work
library(readxl)
library(DESeq2)
library(WriteXLS)
library(ggplot2)
library(ggrepel)
library(plotly)
library(reshape)
library(biomaRt)
library(dplyr)
library(stringr)
library(VennDiagram)
library(pheatmap)
library(EnhancedVolcano)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(genefilter)
library(Glimma)
library(limma)

load.counts.tabular <- function(from){
  # from - directory path, e.g. "/home/vendula/Documents/Projects/FacialPalsy/"
  setwd(from)
  list.filenames = list.files(pattern=".tabular$")
  df.CountMatrix = read.table(list.filenames[1], header=TRUE)
  for (i in 2:length(list.filenames))
  {
    sample = read.table(list.filenames[i], header=TRUE)
    df.CountMatrix = merge(df.CountMatrix, sample, by="Geneid")
  }
  rownames(df.CountMatrix) = df.CountMatrix$Geneid
  return(df.CountMatrix)
}


volcano.plot.single <- function(result, FCthres, pvalthres, where.to, name){
  # FCthres - Fold Chance threshold
  # pvalthres - adjusted pvalue threshold
  descr = mcols(result)
  comparison = substr(descr$description[3], start=17, stop=str_length(descr$description[3]))
  res = as.data.frame(result)
  for (i in 1:nrow(res)){
    if (!is.na(res$padj[i]) & (abs(res$log2FoldChange[i]) > FCthres & res$padj[i] < pvalthres)) { 
      res$color[i] = paste0("padj<",pvalthres," & |LogFC|>", FCthres)
    } else if (!is.na(res$padj[i]) & abs(res$log2FoldChange[i]) < FCthres & res$padj[i] > pvalthres) {
      res$color[i] = "n.s."
    } else if (!is.na(res$padj[i]) & abs(res$log2FoldChange[i]) < FCthres & res$padj[i] < pvalthres) {
      res$color[i] = paste0("padj<",pvalthres)
    } else if (!is.na(res$padj[i]) & abs(res$log2FoldChange[i]) > FCthres & res$padj[i] > pvalthres) {
      res$color[i] = paste0("|LogFC|>", FCthres)
    }
  }
  res$color=as.factor(res$color)
  
  ggplot(res, aes(log2FoldChange, -log10(res$padj))) + 
    geom_point(aes(color=color), shape=16) + 
    scale_color_manual(values=c("forestgreen", "grey48", "blue2", "red")) +
    geom_text(aes(label=ifelse(padj<pvalthres & abs(log2FoldChange)>FCthres, as.character(GeneID),'')),hjust="inward", vjust="inward") +
    #xlim(xlim) + ylim(ylim) +
    ggtitle(paste(comparison)) +
    xlab(bquote(~Log[2]~ 'fold change')) + 
    ylab(bquote(~-Log[10]~adjusted~italic(p))) + 
    geom_hline(yintercept = -log10(pvalthres), color="black", linetype="dashed") +
    geom_vline(xintercept = c(-FCthres, FCthres), color="black", linetype="dashed") +
    theme_bw() + theme(legend.title = element_blank())
  

  ggsave(paste0(where.to,name,".png"), width = 10, height = 10)
}


volcano.plot.EnhancedVolcano <- function(DESeqDataSet, annotated, res1, res2, res3, where.to, ...){
  groups = levels(DESeqDataSet$condition)
  group.name=as.character(design(DESeqDataSet))[2]
  
    p1 = EnhancedVolcano(res1,
                             lab = res1$GeneID,
                             title = paste0(groups[1]," vs. ",groups[2]),
                             x = 'log2FoldChange',
                             y = 'padj',
                             xlim=c(-6,6),
                             ylim=c(0,150),
                             xlab = bquote(~Log[2]~ 'fold change'),
                             ylab = bquote(~-Log[10]~adjusted~italic(p)),
                             pCutoff = 10e-6,
                             FCcutoff = 1.0,
                             transcriptLabSize = 4.0,
                             colAlpha = 1)
  
    p2 = EnhancedVolcano(res2,
                         lab = res2$GeneID,
                         title = paste0(groups[1]," vs. ",groups[3]),
                         x = 'log2FoldChange',
                         y = 'padj',
                         xlim=c(-6,6),
                         ylim=c(0,150),
                         xlab = bquote(~Log[2]~ 'fold change'),
                         ylab = bquote(~-Log[10]~adjusted~italic(p)),
                         pCutoff = 0.0001,
                         FCcutoff = 1.0,
                         transcriptLabSize = 4.0,
                         colAlpha = 1)
    
    p3 = EnhancedVolcano(res3,
                         lab = res3$GeneID,
                         title = paste0(groups[2]," vs. ",groups[3]),
                         x = 'log2FoldChange',
                         y = 'padj',
                         xlim=c(-6,6),
                         ylim=c(0,150),
                         xlab = bquote(~Log[2]~ 'fold change'),
                         ylab = bquote(~-Log[10]~adjusted~italic(p)),
                         pCutoff = 0.0001,
                         FCcutoff = 1.0,
                         transcriptLabSize = 4.0,
                         colAlpha = 1)
    
  

  png(paste0(where.to, "volcanos_geneID.png"), width = 1800, height = 1000)
  grid.arrange(p1,p2,p3,
               ncol=3,
               top = textGrob(...,
                              just = c('center'),
                              gp = gpar(fontsize = 25)))
  grid.rect(gp=gpar(fill=NA))
  dev.off()
  
  
}


venn.diagram.2or3groups <- function(CountMatrix, samples, where.to){
    # replace all 0s with NAs (so these genes are considered as not present = not expressed)
    df = CountMatrix
    df[] <- replace(row.names(df)[row(df)], !df, 0)
    df[df=="0"] <- NA
    
    # group lanes into three groups
    grouping = as.data.frame(samples) 
    groups = unique(grouping$dg)
    
    for ( i in 1:length(groups))
    {
    assign(paste0(groups[i],'.index'),which(colnames(df) %in% subset(grouping, dg==groups[i])$laneID))
    assign(groups[i], unique(data.frame(genes = c(t(df[,eval(as.symbol(paste0(groups[i],'.index')))])))))
    }    
    if (length(groups)==2){
      overlap.12 = sum(!is.na(intersect(eval(as.symbol(groups[1]))$genes, eval(as.symbol(groups[2]))$genes)))
      
      venn.plot = draw.pairwise.venn(area1=nrow(eval(as.symbol(groups[1]))), area2=nrow(eval(as.symbol(groups[2]))), 
                                     cross.area=overlap.12,
                                     category = c(groups[1], groups[2]),
                                     fill = c("blue", "red"),
                                     lty = "blank", col="white",
                                     cex = 2,cat.cex = 2,cat.col = c("blue", "red"))
      png(paste0(where.to,"venn.diagram.2groups.png"))
      grid.draw(venn.plot)
      dev.off()
    }
    if (length(groups)==3){
      overlap.12 = sum(!is.na(intersect(eval(as.symbol(groups[1]))$genes, eval(as.symbol(groups[2]))$genes)))
      overlap.13 = sum(!is.na(intersect(eval(as.symbol(groups[1]))$genes, eval(as.symbol(groups[3]))$genes)))
      overlap.23 = sum(!is.na(intersect(eval(as.symbol(groups[2]))$genes, eval(as.symbol(groups[3]))$genes)))
      overlap.all = sum(!is.na(intersect(eval(as.symbol(groups[1]))$genes, intersect(eval(as.symbol(groups[2]))$genes, eval(as.symbol(groups[3]))$genes))))
  
      venn.plot = draw.triple.venn(area1=nrow(eval(as.symbol(groups[1]))), area2=nrow(eval(as.symbol(groups[2]))), area3=nrow(eval(as.symbol(groups[3]))), 
                               n12=overlap.12, n23=overlap.23, n13=overlap.13, n123=overlap.all, 
                               category = c(groups[1], groups[2], groups[3]),
                               fill = c("blue", "red", "green"),
                               lty = "blank", col="white",
                               cex = 2,cat.cex = 2,cat.col = c("blue", "red", "green"))
    
      png(paste0(where.to,"venn.diagram.3groups.png"))
      grid.draw(venn.plot)
      dev.off()
    }
}




PCA.plot.ZenLab <- function(normalized, samples, name, where.to){
    PCs = plotPCA(normalized, intgroup=c("dg"), returnData=TRUE)
    colnames(df.samples)[1]="laneID"
    colnames(PCs)[ncol(PCs)]="laneID"
    PCs = merge(PCs, samples, by="laneID")
    PCs = PCs[,c("laneID", "PC1", "PC2", "group", "patientID")]
    WriteXLS("PCs", paste0(where.to,"/PCs.xls"), row.names = FALSE)

    ggplot(PCs, aes(x=PC1, y=PC2, col=group)) + geom_point() + theme_bw()
    ggsave(paste0(where.to,name,".png"), width = 10, height = 10)

    ggplot(PCs, aes(x=PC1, y=PC2, col=group)) + geom_point() + geom_text_repel(aes(label=paste(patientID)),hjust="inward", vjust="inward") + theme_bw()
    ggsave(paste0(where.to, name,"_annotated.png"), width = 10, height = 10)

}


heatmap.gene.counts <- function(DESeqDataSet, normalized, how.many, annotated, name, where.to){
  # DESeqDataSet: DESeqDataSet dataset of raw counts
  # normalized: DESeqTransform object of normalized counts
  # annotated: GeneID annotated to Ensembl names
  # filename: path where to save the plot
  select <- order(rowMeans(counts(DESeqDataSet,normalized=TRUE)),
                  decreasing=TRUE)[1:how.many]
  val.mat = assay(normalized)[select,]
  rownames(val.mat) = annotated$external_gene_name[match(rownames(val.mat), annotated$ensembl_gene_id)]
  
  patID = data.frame(laneID=as.character(colData(DESeqDataSet)[,"laneID"]), patientID=as.character(colData(DESeqDataSet)[,"patientID"]), dg = as.factor(colData(DESeqDataSet)[,"dg"]), stringsAsFactors = FALSE)
  patID$patientID = ave(patID$patientID, patID$patientID, FUN = function(i) paste0(i, '.', seq_along(i)))
  colnames(val.mat) <- patID$patientID
  identical(colnames(val.mat), patID$patientID) # test
  
  #annotation
  df = data.frame(dg=patID$dg)
  rownames(df) = patID$patientID # patient IDs
  
  # save plot
  png(paste0(where.to, name,".clustered.", how.many, "genes.png"), width = 2000, height = 1000)
  pheatmap(val.mat, show_rownames=TRUE, show_colnames = TRUE,
           cluster_rows=FALSE, cluster_cols=TRUE,
           annotation_col=df, annotation_names_col = FALSE, main="Heatmap of log2(n+1) transformed counts for 20 most expressed genes")
  dev.off()
  
  
        # The same as above, but reorder so the dg groups are beside each other
         patID.new = patID[order(patID$dg),]
         val.mat.new = val.mat[,order(patID$dg)]
         identical(colnames(val.mat.new), patID.new$patientID) # test
  
        #annotation
         df = data.frame(dg=patID.new$dg)
         rownames(df) = patID.new$patientID # patient IDs
  
  
         png(paste0(where.to, name,".reordered.",how.many,"genes.png"), width = 2000, height = 1000)
         pheatmap(val.mat.new, show_rownames=TRUE, show_colnames = TRUE,
                  cluster_rows=FALSE, cluster_cols=FALSE,
                   annotation_col=df, annotation_names_col = FALSE, main="Heatmap of log2(n+1) transformed counts for 20 most expressed genes")
         dev.off()
  
}


heatmap.sample2sample.dist <- function(normalized, name, where.to){  
  sampleDists <- dist(t(assay(normalized)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- normalized$dg
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  png(paste0(where.to,name, ".png"), width = 2000, height = 2000)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors,
           main="Heatmap of sample-to-sample distances")
  dev.off()
}

heatmap.gene.deviations <- function(DESeqDataSet, normalized, how.many, annotated, name, where.to){  # unfinished
  topVarGenes <- head(order(-rowVars(assay(normalized))), how.many) # top genes with the highest variance (in counts) across samples
  mat <- assay(normalized)[ topVarGenes, ]
  mat <- mat - rowMeans(mat)
  rownames(mat) = annotated$external_gene_name[match(rownames(mat), annotated$ensembl_gene_id)]
  patID = data.frame(laneID=normalized$laneID, patientID=normalized$patientID, stringsAsFactors = FALSE)
  patID$patientID = ave(patID$patientID, patID$patientID, FUN = function(i) paste0(i, '.', seq_along(i)))
  
  df <- as.data.frame(colData(normalized)[,c("dg")])
  colnames(mat) = patID$patientID
  rownames(df) = patID$patientID
  png(paste0(where.to, name, ".png"), width = 2000, height = 1500)
  pheatmap(mat, annotation_col=df, main = paste0("Heatmap of top ",how.many," genes that, in a specific sample, deviate most from the gene's average across all samples"))
  dev.off()
  
}
