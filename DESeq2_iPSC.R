##Aim: To analyze the expression of hiPSC samples 
##data: the output file from htseq-count command
##installation
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
#Clean variable
rm(list = ls()) 
#load library
library("DESeq2")
library("dplyr")
library("stringr")
library("ggplot2")

#Set the working directory 
directory <- "~/iPSC_htseq_all"
file.set <- directory %>% list.files(full.names = T)
file.names <- file.set %>% str_extract("(?<=all/).*(?=_)")

#output directory
setwd("~/htseq_DESeq2/")

#Set the prefix for output file name ##
outputPrefix <- "_All_DESeq2_"

#create sampleTable for DESeq2
sampleFiles <- file.set[1:38] 
sampleFiles #check

sampleNames <- c("hiPSCs(1103)-1","hiPSCs(1103)-2","hiPSCs(1103)-3","hiPSC(1103)-CM-lac-1","hiPSC(1103)-CM-lac-2"," hiPSC(1103)-CM-lac-3",
                 "hiPSC(1103)-CM-TDI-1","hiPSC(1103)-CM-TDI-2"," hiPSC(1103)-CM-TDI-3","hiPSC(1103)-CM-1","hiPSC(1103)-CM-2"," hiPSC(1103)-CM-3",
                 "hAFSC-iPSC-CM-lac-1","hAFSC-iPSC-CM-lac-2","hAFSC-iPSC-CM-lac-3","hAFSC-iPSC-CM-TDI-1","hAFSC-iPSC-CM-TDI-2","hAFSC-iPSC-CM-TDI-3",
                 "hAFSC-iPSC-CM-1","hAFSC-iPSC-CM-2","hAFSC-iPSC-CM-3","hAFSC-iPSC-1","hAFSC-iPSC-2","hAFSC-iPSC-3","hAFSC-iPSC-4",
                 "hESC-CM(RUES2-CM)-lac-1","hESC-CM(RUES2-CM)-lac-2","hESC-CM(RUES2-CM)-lac-3","hESC-CM(RUES2-CM)-TDI-1","hESC-CM(RUES2-CM)-TDI-2","hESC-CM(RUES2-CM)-TDI-3",
                 "hESC-CM(RUES2-CM)-1","hESC-CM(RUES2-CM)-2","hESC-CM(RUES2-CM)-3","hESC(RUES2)-1","hESC(RUES2)-2","hESC(RUES2)-3","hESC(RUES2)-4") #with official name

sampleCondition <- c("SC","SC","SC","CM_lac","CM_lac","CM_lac","CM_TDI","CM_TDI","CM_TDI","CM","CM","CM",
                     "CM_lac","CM_lac","CM_lac","CM_TDI","CM_TDI","CM_TDI","CM","CM","CM","SC","SC","SC","SC",
                     "CM_lac","CM_lac","CM_lac","CM_TDI","CM_TDI","CM_TDI","CM","CM","CM","SC","SC","SC","SC")

sampleCell <- c("1103","1103","1103","1103","1103","1103","1103","1103","1103","1103","1103","1103",
                "AFSC","AFSC","AFSC","AFSC","AFSC","AFSC","AFSC","AFSC","AFSC","AFSC","AFSC","AFSC","AFSC",
                "RUES2","RUES2","RUES2","RUES2","RUES2","RUES2","RUES2","RUES2","RUES2","RUES2","RUES2","RUES2","RUES2")

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition, cell = sampleCell)

treatments = c("SC","CM","CM_lac","CM_TDI")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = "",
                                       design = ~cell+condition+cell:condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,
                                      levels = treatments)

#caculation
dds <- DESeq(ddsHTSeq)
res <- results(dds)
vsd <- varianceStabilizingTransformation(dds, blind=T)

ddsdata <- as.data.frame(counts(dds,normalized = TRUE))
ddsdata$gene_id <- rownames(ddsdata)

#add gene symbol for normalized counts dataframe to save the file
genesymbol.df <- read.table("~/ENSG_to_Genesymbol.csv", sep = ',', header = T)
ddsdata <- left_join(ddsdata, genesymbol.df, by = c("gene_id"="gene"))
ddsdata <- ddsdata %>% relocate("gene_id","Genesymbol", .before ="hiPSCs(1103)-1" )

#order the data with  c("SC","CM","CM_lac","CM_TDI") order before saving the normalized count file
neworder <- c("gene_id","Genesymbol",
              "hiPSCs(1103)-1","hiPSCs(1103)-2","hiPSCs(1103)-3","hiPSC(1103)-CM-1","hiPSC(1103)-CM-2"," hiPSC(1103)-CM-3",
              "hiPSC(1103)-CM-lac-1","hiPSC(1103)-CM-lac-2"," hiPSC(1103)-CM-lac-3","hiPSC(1103)-CM-TDI-1","hiPSC(1103)-CM-TDI-2"," hiPSC(1103)-CM-TDI-3",
              "hAFSC-iPSC-1","hAFSC-iPSC-2","hAFSC-iPSC-3","hAFSC-iPSC-4","hAFSC-iPSC-CM-1","hAFSC-iPSC-CM-2","hAFSC-iPSC-CM-3",
              "hAFSC-iPSC-CM-lac-1","hAFSC-iPSC-CM-lac-2","hAFSC-iPSC-CM-lac-3","hAFSC-iPSC-CM-TDI-1","hAFSC-iPSC-CM-TDI-2","hAFSC-iPSC-CM-TDI-3",
              "hESC(RUES2)-1","hESC(RUES2)-2","hESC(RUES2)-3","hESC(RUES2)-4", "hESC-CM(RUES2-CM)-1","hESC-CM(RUES2-CM)-2","hESC-CM(RUES2-CM)-3",
              "hESC-CM(RUES2-CM)-lac-1","hESC-CM(RUES2-CM)-lac-2","hESC-CM(RUES2-CM)-lac-3","hESC-CM(RUES2-CM)-TDI-1","hESC-CM(RUES2-CM)-TDI-2","hESC-CM(RUES2-CM)-TDI-3")

out_newOrder.df <- ddsdata[,neworder]
write.table(out_newOrder.df, paste0(Sys.Date(),outputPrefix,"result-normCounts.tsv"), sep = '\t', col.names = T, row.names = F, quote = F)

##----- create plots for samples -------
pdf(paste0(Sys.Date(),"_all_sample_tree.pdf"), height = 5, width = 8)
dists <- dist(t(assay(vsd)))
plot(hclust(dists))
graphics.off()

#plot distance by PCA
pdf(paste0(Sys.Date(),"_all_PCAplot.pdf"), height = 5, width = 6)
p1 <- plotPCA(vsd, intgroup=c("condition", "cell")) 
p2 <- p1 + theme_bw()
print(p2)

#plot distance by PCA with different presenting 
pcaData <- plotPCA(vsd, intgroup=c("condition", "cell"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=cell)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_bw()
print(p)
graphics.off()

#plot distance heatmap 
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library("pheatmap")

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$cell,vsd$condition, sep="_")
colnames(sampleDistMatrix) <- paste(vsd$cell,vsd$condition, sep="_")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pdf(paste0(Sys.Date(),"_all_sampledistance_rename.pdf"), height = 8, width = 8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, legend = F)
graphics.off()

##------ TDI verse Lac or CM verse SC to generate DEGs per group and valcano plot -----
## 1103CM-TDI verse 1103CM-lac for example 
sampleFiles <- file.set[c(4:9)]    #file.set[c(13:18)] -for ACM-TDI vs ACM-lac ;#file.set[c(26:31)] -for RCM-TDI vs RCM-lac
sampleNames <- file.names[c(4:9)]  #file.names[c(13:18)]  -for ACM-TDI vs ACM-lac ; #file.names[c(26:31)] -for RCM-TDI vs RCM-lac
sampleCondition <- c("CM_lac","CM_lac","CM_lac","CM_TDI","CM_TDI","CM_TDI") #or c("SC","SC","SC","CM","CM","CM") 

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = "",
                                       design = ~ condition)
colData(ddsHTSeq)$condition <- factor(colData(ddsHTSeq)$condition,)

#calculation
dds <- DESeq(ddsHTSeq)
res <- results(dds)

resdata <- merge(as.data.frame(res), as.data.frame(counts(dds,normalized =TRUE)), by = 'row.names', sort = FALSE)
resdata$color <- ifelse(resdata$pvalue<0.05 & abs(resdata$log2FoldChange)>= 1,ifelse(resdata$log2FoldChange > 1,"UP","DOWN"),"NO")

#add gene symbol 
resdata <- left_join(resdata, genesymbol.df, by = c('Row.names'= 'gene'))

#plot volcano plot
mycolor <- c("blue","gray","red")
pdf(paste0(Sys.Date(),"_volcano_plot_1103CM_lacvsTDI.pdf"))
p1 <- ggplot(resdata, aes(x=log2FoldChange, y=-log10(pvalue), col = color)) +
  geom_point(size = 0.2) + xlim(-20, 20) +
  theme_bw() +
  scale_colour_manual(values = mycolor) +
  labs(x="log2(TDI/lac)",y="-log10(p-value)") +
  geom_hline(yintercept = -log10(0.05), lty=8,col="black",lwd=0.8) +
  geom_vline(xintercept = c(-1, 1), lty=8,col="black",lwd=0.8) +
  theme(legend.position = "none",
        panel.grid=element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        text = element_text(size = 15)) + ggtitle("1103CM_TDI vs 1103CM_lac")
print(p1)
graphics.off()

#order results by padj value 
subset_res = subset(resdata, padj<0.05)
subset_res <- subset_res[order(subset_res$padj),]

#save the significant DEGs with normalized reads to csv
names(subset_res)[1] <- 'gene'
write.csv(subset_res, file = paste0(Sys.Date(),outputPrefix, "_1103CM_TDI_vs_lac-results-normCounts-padj5.csv"), row.names = F)


