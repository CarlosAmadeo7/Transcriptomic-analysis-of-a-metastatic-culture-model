### aCSF in vitrio experiment pipeline.

### Loading libraries
### suppressPackagesStart is suppressing all warnings from the packages 
## Installing packages 
## Simple installation 
####install.packages()
## from Bioconductor 
####install.packages("tidyplots")
####if (!require("BiocManager", quietly = TRUE))
####install.packages("BiocManager")

### Loading Libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(data.table)
  library(RColorBrewer)
  library(tidyverse)
  library(preprocessCore)
  library(future.apply)
  library(DESeq2)
  library(pheatmap)
  library(sva)
  library(viridis)
  library(limma)
  library(emmeans)
  library(broom)
  library(janitor)
  library(tidyplots)
  library(dplyr)
  library(writexl)
  library(scToppR)
  library(fgsea)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
})

## Establising a directory 
setwd("")
## Creating a dge directory
dir.create("deg")

### Loading feature counts
## genes that include protein_coding_genes, pseudogenes, etc.
exp <- read_csv("gene_count.csv")
dim(exp)
head(exp)

## Filtering for protein coding genes, 19934 protein coding genes
p_exp <- exp %>%
  filter(gene_biotype=="protein_coding") %>%
  dplyr::select(-gene_id,-gene_chr,-gene_start,-gene_end,-gene_strand,-gene_length,-gene_biotype,-gene_description,-tf_family) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), sum)) %>%   ## Adding up possible repetitive genes
  arrange(gene_name) %>%
  column_to_rownames("gene_name")


### Creating metadata from scratch, setting up factors and levels
metadata<- data.frame(
  Subject_ID = c("MB03_0_1","MB03_0_2","MB03_0_3",   
                 "MB03_aCSF_1","MB03_aCSF_2","MB03_aCSF_3"),
  Condition = c("Vehicle", "Vehicle", "Vehicle",
                 "aSCF_media", "aSCF_media", "aSCF_media"),
  Sex = c("Female", "Female", "Female", "Female"
          ,"Female", "Female")
)
rownames(metadata) <- metadata$Subject_ID
metadata$Condition<-factor(trimws(metadata$Condition))
metadata$Condition<-factor(metadata$Condition, levels = c("Vehicle","aSCF_media"))

#### Saving metadata
write.table(metadata,file = "dge/metadata.txt", sep = "\t")

## Normalizing counts for PCA and visualization purposes 
## Converting to cpms: Counts per million 
cpm <-future_apply(p_exp, 2, function(x) x/sum(as.numeric(x)) * 10^6) ### number 2 is the FUN value , there is no default 

## Filter low cpm counts for both conditions 
filter<-apply(cpm, 1, function(x) all(x[1:3]>=0.5) | all(x[4:6]>=0.5))
counts_filt <- p_exp[filter,]  ##12449 genes
cpm_filt <- cpm[filter,]   ##12449 genes

#### Saving cpm and filtered counts 
cpm_filt_xlsx<-cpm_filt %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Gene")

counts_filt_xlsx<-counts_filt %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Gene")

write_xlsx(cpm_filt_xlsx, "dge/gene_cpm.xlsx")
write_xlsx(counts_filt_xlsx, "dge/gene_filtered_counts.xlsx")


## Applying log_transformation
logCPM<- log2(cpm_filt +1)
boxplot(logCPM)
## Normalizing quantiles from cpm counts
p<- normalize.quantiles(as.matrix(logCPM))
rownames(p)<-rownames(logCPM)
colnames(p) <-colnames(logCPM)
boxplot(p)

## Performing PCA using the quantile normalization from the cpms
pca.Samples<-prcomp(t(p))
PCi<-data.frame(pca.Samples$x, Condition= metadata$Condition, Sex= metadata$Sex, ID= metadata$Subject_ID )
eig <- (pca.Samples$sdev)^2 ### Calculating the eigenvalue
variance <- eig*100/sum(eig) ## Calculating variance
pdf("dge/PCA_cpms.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(PCi,
          x = "PC1",
          y = "PC2",
          color = "Condition", palette=c("red","black"),
          shape = "Sex", size = 4,label = "ID") + 
  xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
  ylab(paste("PC2 (",round(variance[2],1),"% )"))+
  theme_classic()
dev.off()

#### Applying surrogate variables to clean dataset 
pd<- metadata 
pd$Sex<- NULL
pd$Subject_ID<- NULL

levels(pd$Condition)


mod <- model.matrix(~., pd)
mod0 <- model.matrix(~ 1, pd)

svaobj <- sva(as.matrix(p),mod,mod0,n.sv=NULL,B=100)


svaobj$sv <- data.frame(svaobj$sv)
colnames(svaobj$sv) = c(paste0('SV',seq(svaobj$n.sv)))
pdSv <- cbind(pd,svaobj$sv)

# Regression Expression
pd_sva <- pdSv %>%
  dplyr::select(-Condition) %>% #Removing variables that may cause cofounding
  droplevels()

betas<-future_lapply(1:nrow(p), function(x)
{
  lm(unlist(p[x,])~., data = pd_sva)
})

residuals<-future_lapply(betas, function(x)residuals(summary(x)))
residuals<-do.call(rbind, residuals)
p_regressed <- residuals+matrix(future_apply(p, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
rownames(p_regressed)<-rownames(p)
write.table(p_regressed,"dge/expression_regressed.txt",sep="\t",quote=F)

pdf("dge/PCA_AdjustedForConfound.pdf",width=8,height=8,useDingbats=FALSE)
pca.Sample<-prcomp(t(p_regressed))
PCi<-data.frame(pca.Sample$x,Condition=pd$Condition, ID = rownames(pd))
eig <- (pca.Sample$sdev)^2
variance <- eig*100/sum(eig)
ggscatter(PCi, x = "PC1", y = "PC2",
          color = "Condition",palette=c("red","black"), 
          size = 3,label = "ID")+
  xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
  ylab(paste("PC2 (",round(variance[2],1),"% )"))+
  theme_classic()
dev.off()


#######################################
# Emmeans post-hoc

# Modeling
model1 <- 'geneExpr ~ Condition * Sex' # Null model



## Vehicle 1:
### aCSF: 3
##########################################################################################################
## Surrogate variables (SV) application for correcting PCAs, I am going to do this later 
##########################################################################################################


## Differential expression analysis
### Making sure the columns and rows match before the Deseq2 object 
all(rownames(metadata)==colnames(counts_filt)) # cool
############## DESEQ2 ###############
#####################################
## Creating DESEq object 
dds<- DESeqDataSetFromMatrix(countData = counts_filt,
                             colData = metadata,
                             design = ~ Condition)
counts(dds)
## Estimating factors
dds<-estimateSizeFactors(dds)
## Performing DESEq
dds<-DESeq(dds, 
           full = design(dds),
           betaPrior = FALSE)

## Saving results to a variable, 12449 genes
res<-results(dds, contrast = c("Condition","aSCF_media","Vehicle"))
head(res[order(res$padj),]) ## order to FDR
deg<-subset(res, padj < 0.05 & abs(log2FoldChange) > 0.3 )

#### saving this as excel file 
CTX_FullTab <-as.data.frame(results(dds, contrast = c("Condition","aSCF_media","Vehicle"), cooksCutoff = F, independentFiltering = F))%>%
  rownames_to_column("Gene")
CTX_DGE <- CTX_FullTab %>%
  mutate(Abs = abs(log2FoldChange)) %>%
  filter(padj < 0.05 & Abs > 0.3) %>%   ### FDR<0.05 and absolute value of 0.3 ( 23% up or down)
  arrange(desc(Abs))

save(CTX_FullTab, CTX_DGE,metadata, file = "dge/CTX_Dge_Data.RData")

### Shrinkage just in case , 290 genes
res_shrink<-lfcShrink(dds, coef = "Condition_aSCF_media_vs_Vehicle", type = "ashr" )
degs_shrink<-subset(res_shrink, padj < 0.05 & abs(log2FoldChange) > 0.3 )


head(degs_shrink[order(degs_shrink$log2FoldChange, decreasing = FALSE),])
head(deg[order(deg$log2FoldChange, decreasing = F),])
head(CTX_FullTab[order(CTX_FullTab$log2FoldChange, decreasing = FALSE),])

CTX_DGE_shrink<-as.data.frame(degs_shrink) %>%
  rownames_to_column("Gene")

### Saving them as excel files
write_xlsx(CTX_DGE_shrink,"dge/CTX_Dge_shrink.xlsx" )
write_xlsx(CTX_FullTab, "dge/CTX_FullTab.xlsx")
write_xlsx(CTX_DGE, "dge/CTX_DGE.xlsx")



## Opening dataset
openxlsx::write.xlsx(CTX_FullTab,
                     file = "dge/CTX_FullTab.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(CTX_DGE,
                     file = "dge/CTX_DGE.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE,
                     borders = "columns",
                     sheetName="Stats")

##### edgeR ##################
##############################
metadata$Condition <- factor(metadata$Condition, levels = c("Vehicle", "aSCF_media"))
levels(metadata$Condition)

library(edgeR)
# Create DGEList
dge <- DGEList(counts = counts_filt, group = metadata$Condition)

# Filter: keep genes with CPM > 1 in at least 2 samples
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
## Normalization and inspection 
dge <- calcNormFactors(dge)
plotMDS(dge, labels = metadata$Condition)  # optional visualization
### Design matrix and dispersion estimation 
design <- model.matrix(~ Condition, data = metadata)
dge <- estimateDisp(dge, design)
## Fit model and testing 
fit <- glmFit(dge, design)
lrt <- glmLRT(fit, coef = 2)  # coef=2 tests "aSCF_media vs Vehicle"
## Obtaing the results 
res <- topTags(lrt, n = Inf)$table  # all results

# View top genes
head(res)
res[rownames(res) == "EPAS1",]
#Visualization plots
#### Add Log with the -log10 and the absolite value of logfc
#### Add threshold columns with padj and logfc > 0.3
#### Add direction
#### Full table 
df <- CTX_FullTab %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>%
  mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg")) ### why 0.3
dim(df)

### Top genes based on Direction and padj
### These genes are the most statistically significance
top_labelled <- df %>%
  group_by(Direction) %>%
  na.omit() %>%
  arrange(padj) %>%
  top_n(n = 7, wt = LOG)

head(top_labelled)

### Volcano plots 
pdf("dge/Volcano_Plot_CTX.pdf",width=6,height=6,useDingbats=FALSE)
ggscatter(df,
          x="log2FoldChange",
          y="LOG",
          color = "Threshold",
          palette=c("grey","red"),
          size = 1,
          alpha=0.3,
          shape=19)+
  xlab("log2(Fold Change)")+
  ylab("-log10(FDR)")+
  geom_vline(xintercept = 0, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = 0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_vline(xintercept = -0.3, colour = "black",linetype="dotted",size=1,alpha=0.5) +
  geom_hline(yintercept = 1.3, colour = "grey",linetype="dotted",size=1,alpha=0.5) +
  geom_text_repel(data = top_labelled,
                  mapping = aes(label = Gene),
                  size = 5,
                  box.padding = unit(0.4, "lines"),
                  point.padding = unit(0.4, "lines"))+
  theme(legend.position="none")+
  ylim(0,30) + xlim(-5,+5)
dev.off()

### Heatmaps
## Filtering the genes from p that are differentially expressed
mat <- p[rownames(p)%in% CTX_DGE$Gene,] 
anno<-metadata
Genotype<- c("red", "black")
names(Genotype) <- c("Vehicles", "aSCF media")
anno_colors <- list(Genotype = Genotype)
pdf("dge/Heatmap_CTX_2.pdf",width=6,height=6,useDingbats=FALSE)
pheatmap(mat,
         scale="row",
         show_rownames = F,
         annotation=anno,
         annotation_colors = anno_colors)
dev.off()







