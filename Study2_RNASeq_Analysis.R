# This is an RNA-Seq analysis script released as a part of Keiser et al., 2021 Nature Medicine
# This script was used to perform analysis of Study2 RNA-Seq data which includes the following samples.
# Study2_EmptyCapsid_1_DCN
# Study2_EmptyCapsid_2_DCN
# Study2_EmptyCapsid_3_DCN
# Study2_intmiS1_1_DCN
# Study2_intmiS1_2_DCN
# Study2_intmiS1_3_DCN
# Study2_miSCA7Stuffer_1_DCN
# Study2_miSCA7Stuffer_2_DCN
# Study2_miSCA7Stuffer_3_DCN
# Study2_StufferOnly_1_DCN
# Study2_StufferOnly_2_DCN
# Study2_StufferOnly_3_DCN
# These samples are archived at NCBI Gene Expression Omnibus (GEO), accession number: GSE182666.

#Import dependency packages
library(DESeq2)
library(readr)

# Load counts matrix files
merged_countsdata <- readRDS(file = "~/Documents/Davidson_Lab/Megan/NatureMed_GithubRepo/Study2_merged_countsdata.rds")
merged_metadata <- readRDS(file = "~/Documents/Davidson_Lab/Megan/NatureMed_GithubRepo/Study2_merged_metadata.rds")

# Create DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = merged_countsdata,
  colData = merged_metadata,
  design = ~ Treatment)

dds <- DESeq(dds)


###############
### First lets compare treatment "Stuffer Only (Symptomatic samples)" to "Empty Capsid (non-symptomatic)".  
res <- results(dds, contrast=c("Treatment","SO","EC"), name="Treatment_EC_iMIS1")
resordered <- res[order(res$padj),]

#Log fold change shirnkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="Treatment_SO_vs_EC", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
#Explore results using visualization
plotMA(res, ylim=c(-8,8))
plotMA(resLFC, ylim=c(-14,14))

# Plot the normalized counts of features of interest
# First look at the expression constructs
plotCounts(dds, gene=c("GpKFBextmU6miS1newStfr"), intgroup="Treatment")
plotCounts(dds, gene=c("pKAAV_EF1a_hATXN1L_intron2_miS1"), intgroup="Treatment")
# This can also be used to look for genes of interest by Rhesus macaque Ensembl ID.
plotCounts(dds, gene="ENSMMUG00000011153", intgroup="Treatment")



###############
### Second lets compare treatment "Stuffer Only (Symptomatic samples)" to "Empty Capsid (non-symptomatic)".  
res <- results(dds, contrast=c("Treatment","IM","EC"), name="Treatment_EC_iMIS1")
resordered <- res[order(res$padj),]

#Log fold change shirnkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="Treatment_IM_vs_EC", type="apeglm")
resLFC <- resLFC[order(resLFC$padj),]
#Explore results using visualization
plotMA(res, ylim=c(-8,8))
plotMA(resLFC, ylim=c(-14,14))


#########################
# Create a heatmap
#########################
# Set the comparisson you are intersted in (as shown above)
res <- results(dds, name="Treatment_SO_vs_EC", cooksCutoff=FALSE)
res <- res[order(res$padj),]
head(res)

#Log fold change shirnkage for visualization and ranking
resLFC <- lfcShrink(dds, coef="Treatment_SO_vs_EC", type="apeglm", lfcThreshold = 1)
resLFCoredered <- resLFC[order(resLFC$svalue),]

# Create MA plots as shown above
plotMA(res, ylim = c(-15, 30), main="res")
plotMA(resLFC, ylim = c(-15, 15), main="LFC")

# Get summary stats
summary(res)
summary(resLFC)

# make a list of highly expresed miRNAs (This only works if you convert Ensembl gene IDs to gene names)
NeuronalMiRNAs <- SamplesCombined[grep("mml-mir", rownames(SamplesCombined)),]
NeuronalMiRNAs["rowMeans"] <- rowMeans(NeuronalMiRNAs)
NeuronalMiRNAs <- NeuronalMiRNAs[order(NeuronalMiRNAs$rowMeans, decreasing = TRUE),]
NeuronalMiRNAs <- rownames(NeuronalMiRNAs)[1:30]

# make a list of neuronal miRNAs
NeuronalMiRNAs <- c("mml-mir-128b", "mml-mir-145", "mml-mir-124a-2", "mml-mir-124a-1",
                    "mml-mir-125b-1", "mml-mir-125b-2", "mml-mir-99a", "mml-mir-30b",
                    "mml-mir-9-2", "mml-mir-9-3", "mml-let-7a-2")

NeuronalMiRNAs_Ensembl <- c("ENSMMUG00000027112","ENSMMUG00000027024","ENSMMUG00000026853",
                            "ENSMMUG00000026854","ENSMMUG00000026946","ENSMMUG00000028627",
                            "ENSMMUG00000026997","ENSMMUG00000027123","ENSMMUG00000027085",
                            "ENSMMUG00000027085","ENSMMUT00000048499","ENSMMUG00000028635")

# Create a heatmap
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
#rld <- rlog(dds, blind=FALSE)

library("pheatmap")
df <- as.data.frame(colData(dds)[,c("Tissue","Treatment")])
pheatmap(assay(ntd)[head(rownames(res), n = 100),], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, scale = "row")
pheatmap(assay(ntd)[ImmuneResponse1_2,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, scale = "none")



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
