# This is an RNA-Seq analysis script released as a part of Keiser et al., 2021 Nature Medicine
# This script was used to perform analysis of Study1 RNA-Seq data which includes the following samples.
# Study1_DCN01
# Study1_DCN02
# Study1_DCN03
# Study1_DCN04
# Study1_DCN05
# Study1_DCN06
# Study1_DCN07
# Study1_DCN08
# Study1_DCN09
# Study1_DCN10
# These samples are archived at NCBI Gene Expression Omnibus (GEO), accession number: GSE182666.

setwd("~/Documents/Davidson_Lab/Megan/NatureMed_GithubRepo/Keiser_NatMed_2021/")

#Import dependency packages
library(DESeq2)
library(readr)

# Read in the merged counts matrix and metadata table provided with this script from the Keiser_NatMed_2021 github repository.
countdata <- readRDS(file = "Study1_merged_countsdata.rds")
coldata <- readRDS(file = "Study1_merged_metadata.rds")

# Create the DESeq2 set from the input countdata and coldata
dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  colData = coldata,
  design = ~ Treatment)

###############
### Compare miS1 treated "Mid" samples to control "Buffer" samples.  
# Run Deseq2 differential expression analysis
dds <- DESeq(dds)
res <- results(dds, name="Treatment_Mid_vs_Buffer")
res <- res[order(res$padj),]

# Log fold change shrinkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="Treatment_Mid_vs_Buffer", type="apeglm")
resLFC[order(resLFC$pvalue),]

# Visualize results of "res" and apeglm shrinkage estimator adjusted results "resLFC"
plotMA(res, ylim=c(-8,20))
plotMA(resLFC, ylim=c(-14,14))

# Get the summary of differentially expressed genes
summary(res)
summary(resLFC)

# Plot the normalized counts of features of interest
# First look at the expression constructs
plotCounts(dds, gene=c("GpKFBextmU6miS1newStfr"), intgroup="Treatment")

# This can also be used to look for genes of interest by Rhesus macaque Ensembl ID.
plotCounts(dds, gene="ENSMMUG00000011153", intgroup="Treatment")

###############
### Perform visualizations
# PCA plot
plotPCA(vsd, intgroup=c("Patient", "Tissue", "Treatment"))
plotPCA(vsd, intgroup=c("Treatment"))


######################
# GO analysis
######################
# Install goseq package
#BiocManager::install("goseq")
#BiocManager::install("org.Mmu.eg.db")
library(goseq)
setwd("~/Documents/Davidson_Lab/Megan/")
# Format deseq2 output for GO analysis
res_ordered_Mid_vs_Buffer <- data.frame(res)
rownames(res_ordered_Mid_vs_Buffer) <- res_ordered_Mid_vs_Buffer$X1
genes=as.integer(p.adjust(res_ordered_Mid_vs_Buffer$pvalue[res_ordered_Mid_vs_Buffer$padj!=0],
                          method="BH")<.05)
names(genes)=row.names(res_ordered_Mid_vs_Buffer[res_ordered_Mid_vs_Buffer$padj!=0,])
table(genes)

# Check supported organism
supportedOrganisms()
supportedOrganisms()[supportedOrganisms()$Genome=="rheMac10",]

# Get gene length for rheMac10
library(readr)
mMul_GeneLength <- read_delim("mMul_GeneLength.txt", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)
mMul_GeneLength <- data.frame(mMul_GeneLength)
mMul_GeneLength["length"] <- mMul_GeneLength[,3] - mMul_GeneLength[,2]
mMul_GeneLength <- mMul_GeneLength[,c(1,4)]
colnames(mMul_GeneLength) <- c("GeneID","Length")

genesDF <- data.frame(genes)
#genesDF["GeneID"] <- row.names(res)[1:1000]
genesDF["GeneID"] <- row.names(res) # This line was added because the line above gave an error
genesDFmerged <- merge(genesDF, mMul_GeneLength, by = "GeneID", all = TRUE)
genesDFmerged[is.na(genesDFmerged)] <- 0
genes <- as.integer(genesDFmerged$genes)
names(genes) <- genesDFmerged$GeneID
#genesDFmerged <- genesDFmerged[genesDFmerged$genes!=0,]

# Fitting the probability weighting function (PWF)
pwf=nullp(genes,"rheMac10","ensGene", bias.data = genesDFmerged$Length)
head(pwf)

GO.wall=goseq(pwf,"rheMac10","ensGene", use_genes_without_cat=TRUE)

# Get only the significant categories
sigCats <- GO.wall[which(GO.wall[,2] < 0.50),] # getting the sig categories

# Annotate the GO categories
terms <- stack(lapply(mget(cats, GO.wall), MF))




#####################################
# Pull genes out of a GO category 
#####################################
# Get the gene lists of "numDFinCat" in GO.wall report
getGeneLists <- function(pwf, goterms, genome, ids){
  gene2cat <- getgo(rownames(pwf), genome, ids)
  cat2gene <- split(rep(names(gene2cat), sapply(gene2cat, length)),
                    unlist(gene2cat, use.names = FALSE))
  out <- list()
  for(term in goterms){
    tmp <- pwf[cat2gene[[term]],]
    tmp <- rownames(tmp[tmp$DEgenes > 0, ])
    out[[term]] <- tmp
  }
  out
}

goterms <- GO.wall$category
goList <- getGeneLists(pwf, goterms, "rheMac10", "ensGene")

goList_ImmuneResponse1 <- data.frame(goList[1], stringsAsFactors = FALSE)
goList_ImmuneResponse2 <- data.frame(goList[2], stringsAsFactors = FALSE)
ImmuneResponse1_2 <- c(goList_ImmuneResponse1$GO.0002376)

# Set up data.frame of GO terms in ImmuneResponse1
ImmuneResponse_genes_golist <- data.frame(ImmuneResponse1_2)
colnames(ImmuneResponse_genes_golist) <- c("GeneNames")

# Set up data.frame of DE results from SO vs EC test
DE_Results_SOvsEC <- data.frame(resLFC)
DE_Results_SOvsEC["GeneNames"] <- rownames(DE_Results_SOvsEC)

ImmuneDE_GenesList <- merge(DE_Results_SOvsEC, ImmuneResponse_genes_golist, by = "GeneNames")

# Get human readble gene names
library(readr)
mart_export <- read_table2("mart_export.txt", 
                           skip = 0)


mart_export <- mart_export[,1:2]
colnames(mart_export) <- c("GeneNames","StableID")

ImmuneDE_GenesList_geneNames <- merge(ImmuneDE_GenesList, mart_export, by = "GeneNames")


##########################
# Generate Heatmaps
##########################
# Create a heatmap
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)

library("pheatmap")
df <- as.data.frame(colData(dds)[,c("Tissue","Treatment")])

# Generate a heatmap from the top 100 most 
pheatmap(assay(ntd)[head(rownames(res), n = 100),], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, scale = "row")
pheatmap(assay(vsd)[ImmuneResponse1_2,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, scale = "row")
