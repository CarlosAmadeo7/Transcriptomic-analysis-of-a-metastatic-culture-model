### aCSF in vitrio experiment pipeline.
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
pdf("dge/PCA_cpms.pdf",width=7,height=6,useDingbats=FALSE)
ggscatter(PCi,
          x = "PC1",
          y = "PC2",
          color = "Condition", palette=c("red","black"),
          shape = "Sex", size = 4,label = "ID") + 
  xlab(paste("PC1 (",round(variance[1],1),"% )"))+ 
  ylab(paste("PC2 (",round(variance[2],1),"% )"))+
  theme_classic()
dev.off()


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

CTX_FullTable_shrink<-as.data.frame(res_shrink) %>%
  rownames_to_column("Gene")
CTX_DGE_shrink<-as.data.frame(degs_shrink) %>%
  rownames_to_column("Gene")

### Saving them as excel files
write_xlsx(CTX_FullTab, "dge/CTX_FullTab.xlsx")
write_xlsx(CTX_DGE, "dge/CTX_DGE.xlsx")
write_xlsx(CTX_FullTable_shrink,"dge/CTX_FullTab_shrink.xlsx" )
write_xlsx(CTX_DGE_shrink,"dge/CTX_Dge_shrink.xlsx" )


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

# Visualization plots
#### Add Log with the -log10 and the absolite value of logfc
#### Add threshold columns with padj and logfc > 0.3
#### Add direction
#### Full table 
df <- CTX_FullTable_shrink %>%
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
  top_n(n = 5, wt = LOG)

head(top_labelled)

### Volcano plots 
pdf("dge/Volcano_Plot_CTX.pdf",width=8,height=6,useDingbats=FALSE)
ggplot(df, aes(x = log2FoldChange, y = LOG)) +
  theme_classic(base_size = 13) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    axis.ticks   = element_line(linewidth = 0.8),
    axis.ticks.length = unit(0.25, "cm")
  ) +
  geom_point(aes(color = Direction), size = 1, alpha = 0.8, show.legend = FALSE) +
  scale_color_manual(values = c(DownReg = "#5CACDB", NS = "grey", UpReg = "#EA7FA3")) +
  ggnewscale::new_scale_color() +
  geom_text_repel(
    data = top_labelled,
    aes(label = Gene, color = Direction),
    max.overlaps = Inf,
    size = 4,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(DownReg = "#24658C", UpReg = "#A93D5D")) +
  geom_vline(xintercept = 0, linetype = 2, color = "red") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red") +
  xlab("avg_log2FC") +
  ylab("-log10(adj p-value)") +
  scale_x_continuous(limits = c(-5, 5))
dev.off()

### Heatmaps
## Filtering the genes from p that are differentially expressed
mat <- p[rownames(p)%in% CTX_DGE_shrink$Gene,] 
anno<-metadata
Genotype<- c("red", "black")
names(Genotype) <- c("Vehicles", "aSCF media")
anno_colors <- list(Genotype = Genotype)
pdf("dge/Heatmap_CTX_2.pdf",width=6,height=6,useDingbats=FALSE)
pheatmap(mat,
         scale="row",
         show_rownames = F,
         annotation=anno,
         annotation_colors = anno_colors,
         color = viridis(100))
dev.off()

## Boxplots 
#  boxplots
mat <- p[rownames(p)%in% top_labelled$Gene,] %>%
  t() %>%
  as.data.frame() %>%
  mutate(Genotype = metadata$Condition) %>%
  pivot_longer(!Genotype, names_to = "Gene", values_to="Exp")

ggboxplot(mat, "Genotype", "Exp", color = "Genotype",
          palette = c("red", "black")) +
  xlab("")+ 
  ylab("log2(Expression Adjusted)")+
  theme_classic() + 
  facet_wrap(.~Gene,scales="free",ncol=4,nrow=3) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) 

### Gene ontology, KEEG, etc
tmp1 <- CTX_FullTable_shrink %>%
  mutate(LOG = -log10(padj), ABS = abs(log2FoldChange)) %>%
  mutate(Threshold = if_else(padj < 0.05 & ABS > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg")) %>%
  filter(Threshold == TRUE)

markers_list <- split(tmp1$Gene, tmp1$Direction)

for(i in 1:length(markers_list)){
  markers_list[[i]] <- bitr(markers_list[[i]], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
}
genelist <- list()
for(i in 1:length(markers_list)){
  genelist[[i]] <- markers_list[[i]]$ENTREZID
}
names(genelist) <- names(markers_list)

# GO_BP
markers_GOBP <- compareCluster(geneClusters = genelist, fun = "enrichGO", ont="BP", OrgDb="org.Hs.eg.db")
# KEGG
markers_kegg <- compareCluster(geneClusters = genelist, fun="enrichKEGG", organism="hsa", keyType="kegg")
# WikiPathway
markers_WP <- compareCluster(geneClusters = genelist, fun="enrichWP", organism="Homo sapiens")
# Reactomeでclusterprofiler
markers_Reactome <- compareCluster(genelist, fun = "enrichPathway")

## For visualization...........
## go_bp
col <- c("#EA7FA3", "#5CACDB")
names(col) <- unique(tmp1$Direction)
p<- enrichplot::cnetplot(markers_GOBP,
                         showCategory = 10,
                         size_category = 2,
                         node_label= 'category', size_edge = 0.01)
pdf("dge/ORA/DEGs_GOBP.pdf",height = 8, width = 12)
p + scale_fill_manual(values = col) + theme(
  text = element_text(size = 10),        
  legend.text = element_text(size = 10)
)
dev.off()

## go_bp
col <- c("#EA7FA3", "#5CACDB")
names(col) <- unique(tmp1$Direction)
p<- enrichplot::cnetplot(markers_kegg,
                         showCategory = 10,
                         size_category = 2,
                         node_label= 'category', size_edge = 0.01)
pdf("dge/ORA/DEGs_KEEG.pdf",height = 8, width = 12)
p + scale_fill_manual(values = col) + theme(
  text = element_text(size = 10),        
  legend.text = element_text(size = 10)
)
dev.off()
## go_bp
col <- c("#EA7FA3", "#5CACDB")
names(col) <- unique(tmp1$Direction)
p<- enrichplot::cnetplot(markers_Reactome,
                         showCategory = 10,
                         size_category = 2,
                         node_label= 'category', size_edge = 0.01)
pdf("dge/ORA/DEGs_REACTOME.pdf",height = 8, width = 12)
p + scale_fill_manual(values = col) + theme(
  text = element_text(size = 10),        
  legend.text = element_text(size = 10)
)
dev.off()
## go_bp
col <- c("#EA7FA3", "#5CACDB")
names(col) <- unique(tmp1$Direction)
p<- enrichplot::cnetplot(markers_WP,
                         showCategory = 10,
                         size_category = 2,
                         node_label= 'category', size_edge = 0.01)
pdf("dge/ORA/DEGs_WP.pdf",height = 8, width = 12)
p + scale_fill_manual(values = col) + theme(
  text = element_text(size = 10),        
  legend.text = element_text(size = 10)
)
dev.off()


save(markers_GOBP, markers_kegg, markers_WP, markers_Reactome, file = "dge/ORA/ORA_outputs.RData")
#TEST<- markers_GOBP@compareClusterResult

####### GSEA and GO ontology 
library(aPEAR)
conv <- bitr(CTX_FullTable_shrink$Gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)
df_rank <- CTX_FullTable_shrink %>%
  inner_join(conv, by = c("Gene" = "SYMBOL")) %>%
  filter(!is.na(log2FoldChange)) %>%
  distinct(ENTREZID, .keep_all = TRUE) %>%   # avoid duplicate ENTREZ IDs
  arrange(desc(log2FoldChange))
geneList <- df_rank$log2FoldChange
names(geneList) <- df_rank$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

## BP
gse_list_BP <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType  = "ENTREZID",
                     ont = "BP", minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, verbose = FALSE)

TEST<- gse_list_BP@result |> filter (qvalue < 0.1)
enrichmentNetwork(TEST, drawEllipses = TRUE, fontSize = 2.5, repelLabels = T)

## REACTOME
library(ReactomePA)
gse_list_REACTOME <- gsePathway(geneList = geneList, organism = "human", minGSSize = 10,
                                maxGSSize = 500, pvalueCutoff  = 1, verbose = FALSE)

TEST<- gse_list_REACTOME@result |> filter (qvalue < 0.1)
enrichmentNetwork(TEST, drawEllipses = TRUE, fontSize = 2.5, repelLabels = T)

## Barplot 
library(forcats)
up <- TEST %>%dplyr::filter(NES > 0) %>%dplyr::arrange(desc(NES)) %>%dplyr::slice_head(n = 5)
down <- TEST %>%dplyr::filter(NES < 0) %>%dplyr::arrange(NES) %>% dplyr::slice_head(n = 5)
top10 <- dplyr::bind_rows(up, down)
top10 <- top10 %>%mutate(Description = fct_reorder(Description, NES))
ggplot(top10, aes(x = NES, y = Description, fill = NES)) +
  geom_col() +
  scale_fill_viridis_c() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "NES", y = "GO BP Pathways")


#### Running Transcription factor activity 
# Using decoupleR
library(decoupleR)
p
metadata
deg_table<- CTX_DGE%>%
  dplyr::select(Gene, log2FoldChange, stat, padj) %>% 
  dplyr::filter(!is.na(stat)) %>% 
  tibble::column_to_rownames(var = "Gene") %>%
  as.matrix()

net <- decoupleR::get_collectri(organism = 'human', split_complexes = FALSE)

contrast_acts <- decoupleR::run_ulm(mat = deg_table[, 'stat', drop = FALSE], 
                                         net = net, 
                                         .source = 'source', 
                                         .target = 'target',
                                         .mor='mor', 
                                         minsize = 5) 

# Filter top TFs in both signs
n_tfs = 25
f_contrast_acts <- contrast_acts %>%
  dplyr::mutate(rnk = NA)
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))
tfs <- f_contrast_acts %>%
  dplyr::arrange(rnk) %>%
  head(n_tfs) %>%
  dplyr::pull(source)
f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)

# viz
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
p <- ggplot2::ggplot(data = f_contrast_acts, 
                     mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                            y = score)) + 
  ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                    color = "black",
                    stat = "identity") +
  ggplot2::scale_fill_gradient2(low = colors[1], 
                                mid = "whitesmoke", 
                                high = colors[2], 
                                midpoint = 0) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                 axis.text.x = ggplot2::element_text(angle = 45, 
                                                     hjust = 1, 
                                                     size = 10, 
                                                     face = "bold"),
                 axis.text.y = ggplot2::element_text(size = 10, 
                                                     face = "bold"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) +
  ggplot2::xlab("TFs")

p

### interpretation 
# aCSF looks more like a stress-adaptation state with reduce proliferation, redox, metabolic rewriring.
# NQO1 is high : NQO1 is a canonical antioxidant / quinone detox gene and a common readout of NRF2 pathway activation
# GGT5 is also strongly up induced and fits glutathione handling / extracellular redox adaptation.

# Cell cylce is broadly lower, suggesting a more quiescent-like adapattion? , most likely a stress-adapation
# Metabolic rewriring to adaptation: HK2, SLC2A1/GLUT1, ME1, FLT1, BHLHE41 are up.
# BDH1 and PLIN5 are up, which hints at altered lipid / ketone handling.
# TXNIP is strongly down, TXNIP is a negative regulator of glucose uptake and often opposes glycolytic adaptation, its suppression is consistent with increased glucose-scavenging 

# This looks like cells trying to cope with a nutrient-poor or atypical fluid environment.
# The pattern suggests a shift toward stress-resistant metabolism, including more glucose uptake and possibly alternate fuel handling.

# Neuronal differentaition also down
# NPTX1 strongly down., SHISA9, CACNG2, CNIH3, EPHA7, SEMA3A, SEMA3G, TNR, PTPRS, GAD1, NEFH are lower.

## ECM/Adhesion is messy:
#LAMA1 down strongly , ITGA5 down, ELN, FBLN1, CDH2, CDH6, SDC4, FSCN1, SPOCK1, 
# COL5A1 and COLGALT2 are up.

# STEM features 
# PROM1/CD133: Up log2FC 0.364, padj 0.0349.

# KLF rewiring 
#KLF10, KLF11 down, KLF15 up,
# KLF family shifts often reflect changes in growth control, stress response, differentiation, and metabolism.
# KLF15 in particular fits metabolic adaptation logic better than a simple proliferation story.


sessionInfo()


