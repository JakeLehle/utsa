# This script is concered with work I did in my postdoc with the MCCarrey Lab trying to identify the SSC signature and id epigentic marks from FOXC2 that could be establishing this cell state early in male development
# The first focus was to determine the SSC signature using single cell data collected by Lorena Roa from her work in the Hermann Lab during her doctoral studies.
# Here are the DEGs when comparing P3 dim to P6 Bright cells the degs are spcifically for the P3 Dim we know that is the case since the list indicates that Id4 was downregulated which we know has to go up in P6 bright vs P3 dim

P6_vs_P3 <- read.csv("/mnt/f/Golden_boy_project/Master_Dir/Single_cell/Files/excel/Jake/P6B_vs_P3D.csv")

P6_vs_P3_all <- as.data.frame(P6_vs_P3)

table(grepl("Id4", P6_vs_P3_all_UP$X))
class(P6_vs_P3_all$avg_log2FC)
length(P6_vs_P3_all$X)
P6_vs_P3_all_UP <- P6_vs_P3_all %>% filter(avg_log2FC < 0) 
length(P6_vs_P3_all_UP$X)



P6_vs_P3_LFC <- P6_vs_P3_all %>% filter(avg_log2FC <= -1.5 | avg_log2FC >= 1.5)

P6_vs_P3_LFC # 1 gene! Have to just use the p_val filtering



# We next look at the single cells P6 bright data vs P6 dim data

SC_DEG_P6B_vs_P6D <- read.csv("/mnt/e/Golden_boy_project/Master_Dir/Single_cell/Files/excel/Jake/P6B_vs_P6D.csv")
SC_DEG_P6B_vs_P6D <- as.data.frame(SC_DEG_P6B_vs_P6D)
length(SC_DEG_P6B_vs_P6D$X)
SC_DEG_P6B_vs_P6D_UP <- SC_DEG_P6B_vs_P6D %>% filter(avg_log2FC < 0) 

SC_DEG_P6B_vs_P6D_UP
length(SC_DEG_P6B_vs_P6D_UP$X)


SC_DEG_P6B_vs_P6D  %>% filter(avg_log2FC <= -1.5 | avg_log2FC >= 1.5)
# There are no genes that meet the 1.5 LFC cutoff



SSC_overlapping_SC <- intersect(P6_vs_P3_all$X, SC_DEG_P6B_vs_P6D$X)
SSC_overlapping_SC

SSC_overlapping_SC_UP <- intersect(SC_DEG_P6B_vs_P6D_UP$X, P6_vs_P3_all_UP$X)
SSC_overlapping_SC_UP

# There are limited numbers of DEGs within this comparison. John wanted me to look at Keren's bulk data for P6 samples to see if there would be more information about DEGs


#////////////////////////////////////////////////////////////////////////////////
# BULK RNA-seq work flow

# Get the reference genome within the terminal
#$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz
#$ wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/latest_assembly_versions/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.gtf.gz

library(DESeq2)
library(Rsubread)

setwd("/work/sdz852/goldenboy/bulk_rna_seq/fastq/")

# Step 1: prepare your own reference genome index.
buildindex(basename = "index", reference = "/work/sdz852/goldenboy/bulk_rna_seq/fastq/GCF_000001635.27_GRCm39_genomic.fna")


# Step 2: read mapping

# Functional Programming approach to aligning samples
# Here I create a custom lazy_align function that I will then use to lapply() a list of samples across in much more compact code

lazy_align <- function(sample) {
        align(index = "index",
                readfile1 = paste0("/work/sdz852/goldenboy/bulk_rna_seq/fastq/", sample ,"_1.fastq.gz"),
                readfile2 = paste0("/work/sdz852/goldenboy/bulk_rna_seq/fastq/", sample ,"_2.fastq.gz"),
                type = "rna",
                input_format = "FASTQ",
                output_format = "BAM",
                sortReadsByCoordinates = TRUE,
                useAnnotation = TRUE,
                annot.ext = "/work/sdz852/goldenboy/bulk_rna_seq/fastq/GCF_000001635.27_GRCm39_genomic.gtf",
                isGTF = TRUE,
                GTF.featureType = "exon",
                GTF.attrType = "gene_id",
                chrAliases = NULL,
                nthreads = 64)
}

filesfastq_2023 <- list.files(pattern = "\\_1.fastq.gz$")
filesfastq_2023
filesfastq_2023_cut <- substr(filesfastq_2023, 1,nchar(filesfastq_2023)-11)
filesfastq_2023_cut

lapply(filesfastq_2023_cut, lazy_align)


bam_files <- list.files(pattern = "\\.subread.BAM$")


compact_featureCounts <- function(sample) {
			      featureCounts(sample,
                              
                              # annotation
                              annot.ext="/work/sdz852/goldenboy/bulk_rna_seq/fastq/GCF_000001635.27_GRCm39_genomic.gtf",
                              isGTFAnnotationFile=TRUE,
                              GTF.featureType="exon",
                              GTF.attrType="gene_id",
                              chrAliases=NULL,
                              
                              # level of summarization
                              useMetaFeatures=TRUE,
                              
                              # overlap between reads and features
                              allowMultiOverlap=FALSE,
                              minOverlap=1,
                              largestOverlap=FALSE,
                              readExtension5=0,
                              readExtension3=0,
                              read2pos=NULL,
                              
                              # multi-mapping reads
                              countMultiMappingReads=FALSE,
                              fraction=FALSE,
                              
                              # read filtering
                              minMQS=0,
                              splitOnly=FALSE,
                              nonSplitOnly=FALSE,
                              primaryOnly=FALSE,
                              ignoreDup=FALSE,
                              
                              # strandness
                              strandSpecific=0,
                              
                              # exon-exon junctions
                              juncCounts=FALSE,
                              genome=NULL,
                              
                              # parameters specific to paired end reads
                              isPairedEnd=TRUE,
                              requireBothEndsMapped=TRUE,
                              checkFragLength=FALSE,
                              minFragLength=50,
                              maxFragLength=600,
                              countChimericFragments=TRUE,
                              autosort=TRUE,

			      # Threads
			      nthreads=32
)}


fc_SSC <- compact_featureCounts(bam_files)



setwd("/work/sdz852/goldenboy/bulk_rna_seq/")

save.image("Master_R_Env_1-5-24.RData")

q()

getwd()
#///////////////////////////////////////////////////////////////////////////////////////////////////

#Start the RNA work flow 

setwd("/mnt/f/Golden_boy_project/Master_Dir/Single_cell/Scripts/R/Jake")
load("Master_R_Env_1-5-24.RData")

bam_files


#Install Packages
BiocManager::install("edgeR")
install.packages("rlang")
BiocManager::install("tidyverse")
BiocManager::install("DESeq2")
install.packages("devtools")
BiocManager::install("Rsubread")
install.packages("stats")
install.packages("statmod")
install.packages("tibble")
install.packages("dplyr")

#Call libraries

library(edgeR)
library(rlang)
library(tidyverse)
library(DESeq2)
library(devtools)
library(Rsubread)
library(stats)
library(statmod)
library(tibble)
library(dplyr)
library(ggplot2)

 
# You can check that your fc object was set correctly using head() and specificng to return the information in the counts 
head(fc_SSC$counts)

ncol(fc_SSC$counts)

fc_SSC$targets
#edgeR analysis
group_all <- factor(c(rep("ID4_Bright", 3), rep("ID4_Dim", 3)))

y_all <- DGEList(counts=fc_SSC$counts,group=group_all)


#Filter out genes that have no expression
keepy_all <- filterByExpr(y_all)


y_all <- y_all[keepy_all,,keep.lib.sizes=FALSE]


#Normalize the library size using TMM normalization
y_all <- calcNormFactors(y_all, method = "TMM")


y_all$samples


##The design matrix: to determine your control and treatment
design_all <- model.matrix(~ 0 + group_all)

colnames(design_all) <- levels(group_all)

## Estimating the dispersion
y_all <- estimateDisp(y_all, design_all, robust=TRUE)

y_all$common.dispersion
# Common dispersion is low this is good for getting a high number of DEGs

## QL dispersions
fit_all <- glmQLFit(y_all, design_all, robust=TRUE)

head(fit_all$coefficients)


## Differential expression analysis
# the left is treatment group, right one is control.

EvsSSC_all <- makeContrasts(ID4_Bright - ID4_Dim, levels=design_all)


DEG_P6B_vs_P6D_All <- glmQLFTest(fit_all, contrast=EvsSSC_all)

summary(decideTests(DEG_P6B_vs_P6D_All))

DEG_P6B_vs_P6D_All <- as.data.frame(decideTests(DEG_P6B_vs_P6D_All))

colnames(DEG_P6B_vs_P6D_All) <- "DEGs"
DEG_P6B_vs_P6D_All
DEG_P6B_vs_P6D_UP <- DEG_P6B_vs_P6D_All %>% filter(DEGs == 1)
table(grepl("Id4", rownames(DEG_P6B_vs_P6D_UP)))

DEG_P6B_vs_P6D_All <- DEG_P6B_vs_P6D_All %>% filter(DEGs != 0) 
 

# to filter out only genes with 1.5 fold change or more in DEG
DEG_P6B_vs_P6D_LFC <- glmTreat(fit_all, contrast=EvsSSC_all, lfc = log2(1.5))

tmp <- DEG_P6B_vs_P6D_LFC$table 
tmp <- tmp %>% filter(logFC >= 1.5 | logFC <= -1.5)
tmp <- tmp %>% filter(PValue <= 0.01)
length(tmp$logFC)
summary(decideTests(DEG_P6B_vs_P6D_LFC))

DEG_P6B_vs_P6D_LFC <- as.data.frame(decideTests(DEG_P6B_vs_P6D_LFC))

colnames(DEG_P6B_vs_P6D_LFC) <- "DEGs"
DEG_P6B_vs_P6D_LFC
DEG_P6B_vs_P6D_LFC <- DEG_P6B_vs_P6D_LFC %>% filter(DEGs != 0) 

DEG_P6B_vs_P6D_LFC

#Compare this DEG list to what was published in the 2020 paper
DEGs_SSC_Bulk_Keren_2020 <- read.csv("/mnt/e/Golden_boy_project/Master_Dir/Single_cell/Files/excel/Jake/DEGs_SSC_Bulk_Keren_2020.csv")
DEGs_SSC_Bulk_Keren_2020 <- DEGs_SSC_Bulk_Keren_2020 %>% filter(logFC >= 1.5 | logFC <= -1.5)
DEGs_SSC_Bulk_Keren_2020 <- DEGs_SSC_Bulk_Keren_2020 %>% filter(PValue <= 0.01)
length(DEGs_SSC_Bulk_Keren_2020$Symbol)
Keren_sanity_check <- intersect(DEGs_SSC_Bulk_Keren_2020$Symbol, rownames(tmp)) 
length(Keren_sanity_check)
Keren_sanity_check

###############################################################################
# Okay let's compare this to the data we have collected for the P6 to P3 single comparison to the bulk data
###############################################################################


################################################################################
# Here I do a quick check to see to what extent that there is good overlap between the P6 bulk and cingle cell DEGs
DEG_P6B_vs_P6D_All



SSC_overlapping_SC_BULK <- intersect(rownames(DEG_P6B_vs_P6D_All), SC_DEG_P6B_vs_P6D$X)
length(SSC_overlapping_SC_BULK)
# Looks like there is good overlap between both of the datasets. ~ 80% of the genes from the cingle cell P6 DEGs are in the bulk data for P6

# Now lets do the comparison to find the overlap for the P6BvP3D SC and P6BvP6D Bulk which we will call SSC_overlapping

SSC_overlapping <- intersect(P6_vs_P3_all$X, rownames(DEG_P6B_vs_P6D_All))
SSC_overlapping

#John was interested in seeing the genes that are uniquely upregulated in SSCs

SSC_overlapping_UP <- intersect(P6_vs_P3_all_UP$X, rownames(DEG_P6B_vs_P6D_UP))
SSC_overlapping_UP # Very interesting 184 genes that contain genes tat have overlaps with FOXC2 ChIP binding

# Just another check lets compare the overlapping list to see how similar they are

Overlapping_2 <- intersect(SSC_overlapping, SSC_overlapping_SC)

Overlapping_2
# Not bad looks like 1/3 of the genes are overlapping and it contains Id4 and Egr1 which we know later on are important genes for the FOXC2 story 

#//////////////////////////////////////////////////////////////////////////////
# Let's see how many cell death or DNA repair genes are in these SSC lists we have generated

Cell_Death <- as.data.frame(read.csv("/mnt/f/Golden_boy_project/Master_Dir/Single_cell/Files/excel/Jake/Cell_Death.csv"))
Cell_Death
DNA_Repair <- as.data.frame(read.csv("/mnt/f/Golden_boy_project/Master_Dir/Single_cell/Files/excel/Jake/DNA_Repair.csv"))
DNA_Repair

# Here we do the SC comparison
SSC_overlapping_SC_CD <- intersect(SSC_overlapping_SC, Cell_Death$Final.list.cell.death)
length(SSC_overlapping_SC_CD)


SSC_overlapping_SC_DR <- intersect(SSC_overlapping_SC, DNA_Repair$Final.list.cell.repair)
length(SSC_overlapping_SC_DR)


SSC_overlapping_SC_CD_UP <- intersect(SSC_overlapping_SC_UP, Cell_Death$Final.list.cell.death)
length(SSC_overlapping_SC_CD_UP)


SSC_overlapping_SC_DR_UP <- intersect(SSC_overlapping_SC_UP, DNA_Repair$Final.list.cell.repair)
length(SSC_overlapping_SC_DR_UP)

# Now we do the SC nand bulk comparison
SSC_overlapping_CD <- intersect(SSC_overlapping, Cell_Death$Final.list.cell.death)
length(SSC_overlapping_CD)
SSC_overlapping_CD

SSC_overlapping_DR <- intersect(SSC_overlapping, DNA_Repair$Final.list.cell.repair)
length(SSC_overlapping_DR)
SSC_overlapping_DR


# Now we do the SC nand bulk comparison
SSC_overlapping_UP_CD <- intersect(SSC_overlapping_UP, Cell_Death$Final.list.cell.death)
length(SSC_overlapping_UP_CD)
SSC_overlapping_UP_CD


SSC_overlapping_UP_DR <- intersect(SSC_overlapping_UP, DNA_Repair$Final.list.cell.repair)
length(SSC_overlapping_UP_DR)
SSC_overlapping_UP_DR


# GO

install.packages("enrichR")
library(enrichR)


websiteLive <- getOption("enrichR.live")
websiteLive
if (websiteLive) {
  listEnrichrSites()
  setEnrichrSite("Enrichr") # Human or mouse genes   
}

if (websiteLive) dbs <- listEnrichrDbs()
dbs
dbs %>% filter(grepl("Panth", libraryName))
dbs <- c('GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023', 'KEGG_2019_Mouse', 'WikiPathways_2019_Mouse', "Panther_2016")


enriched <- enrichr(unique(SSC_overlapping), dbs)

enriched[['GO_Biological_Process_2023']]
enriched[['GO_Cellular_Component_2023']]
enriched[['GO_Molecular_Function_2023']]
enriched[['KEGG_2019_Mouse']]
enriched[['WikiPathways_2019_Mouse']]
enriched[['Panther_2016']]


tmp <- enriched[['GO_Biological_Process_2023']]


tmp$Term
tmp %>% filter(grepl("Death", Term))

x <- head(enriched[["GO_Biological_Process_2023"]])
x[,1] <- gsub("GO:", "GO_", x[,1])
kable(x)

plotEnrich(enriched[[1]], showTerms = 15, numChar = 90, y = "Count", orderBy = "P.value") +
  ggtitle("P6BvsP3 SC Data vs P6BvsP6D Bulk Data") +
  theme(
    axis.text = element_text(size=12, face="bold"),
    plot.title = element_text(size=35, face="bold", hjust = 1),
    axis.text.y = element_text(size=18, face="bold"),
    axis.title.x = element_text(size=12, face="bold"),
    axis.title.y = element_text(size=12, face="bold"),
    legend.title = element_text(size=12, face="bold")
  )

# If you wanted to peak through the raw pathway output not in a graphical format
enriched[[3]]


#///////////////////////////////////////////////////////////////////////////////

#Okay I have 320 genes that are the "SSC Barcode"
#Lets check these against any of the peaks from the FOXC2 P6 Bright ChIP-seq data and see if we can find overlaps
# First things first we need to get the FOXC2 data.

library(plyranges)
library(data.table)
library(GenomicRanges)
BiocManager::install("Mus.musculus")
library(Mus.musculus)

FOXC2_Peaks_df <- as.data.frame(read.csv("/mnt/e/Golden_boy_project/Master_Dir/Single_cell/Files/excel/Jake/FOXC2_P6B_Peaks.csv"))
head(FOXC2_Peaks_df)
FOXC2_Peaks_df[,2]
gr_FOXC2 <- GRanges(
  seqnames = FOXC2_Peaks_df[,2],
  ranges = IRanges(FOXC2_Peaks_df[,3], end = FOXC2_Peaks_df[,4]),
  strand = Rle(strand(c( "*")), c(length(FOXC2_Peaks_df[,2]))))



Mus.musculus

Mm_genes <- transcriptsBy(Mus.musculus, by="gene", columns=c("SYMBOL", "ENTREZID", "TXCHROM", "TXSTRAND"))
Mm_genes <- unlist(Mm_genes)
Mm_genes

Mm_genes_promoters <- promoters(trim(Mm_genes), upstream = 900, downstream = 400)


overlap_gene_FOXC2 <- gr_FOXC2 %>% 
  plyranges::join_overlap_intersect(Mm_genes) %>% unique()


unique(unlist(overlap_gene_FOXC2$SYMBOL))


#Here I do a comparison to see if any of the SSC barcode genes come up in the FOXC2 annotated peaks
intersect(unique(unlist(overlap_gene_FOXC2$SYMBOL)), SSC_overlapping)


overlap_promoters_FOXC2 <- gr_FOXC2 %>% 
  plyranges::join_overlap_intersect(Mm_genes_promoters)

overlap_promoters_FOXC2 <- as.data.frame(overlap_promoters_FOXC2$SYMBOL)


#Same comparison as above but with gene promoter regions
intersect(unique(overlap_promoters_FOXC2$value), SSC_overlapping)


#///////////////////////////////////////////////////////////////////////////////
# John had some questions about why I didn't see a peak in the Id4 promoter so here I put together some code to dine out that I was a little too conservative with my promoter regions and the peak is only 189 downstream of the Id4 promoter region I establish so nearly there

#Here is where we see overlaps with the Id4 gene and the FOXC2 ChIP-seq peaks

overlap_gene_FOXC2_Id4 <- overlap_gene_FOXC2 %>% filter(unlist(overlap_gene_FOXC2$SYMBOL) == 'Id4')
overlap_gene_FOXC2_Id4

# This is Id4 gene.

table(grepl('Id4', unlist(Mm_genes$SYMBOL)))
Mm_genes_Id4 <- Mm_genes %>% filter(unlist(Mm_genes$SYMBOL) == 'Id4') 
Mm_genes_Id4

# we can see the peak is ~600 bp downstream of the TSS and I cut of my "promoter" regions at 400 bp downstream of the promoter region 

#///////////////////////////////////////////////////////////////////////////////

# Okay John was very happy with all of this but was surprised that there were not as many genes that had a FOXC2 peak 
# Our hypothesis is that FOXC2 is helping epeigeneticaly mark the genes that will later come on marking the fully functional SSCs capable of recolonizing a testes if transplanted which we know ID4 bright are.
# Thus by Identifying the genes that are unique to Id4 bright SSCs we predict that we would see a larger number of these genes having FOXC2 bound in their promoter region or at least more than 4.
# I personally disagree with this assessment as that is too simple I think it is likely a combinatorial TF response that establishes cellular id and FOXC2 is likely a final piece which only has to bind to a few ver crucial genes one of which happens to be Id4
# To get an idea of what other TF might be playing a role I am going to look relative expression of all of the TFs in P6 bright cells from out bulk data and also do some motif discovery at the promoter regions of the genes we identified as part of the ssc barcode

# Let's start off with the heatmap because motif discovery takes forever
# I need a way to find all TF gene names. I found a database called TF link that gives out information collected about TF and their validated binding sites I'll start there and scrape toegther a good testing list to use for a heatmap comparison to see which are being expressed

TF_Link_DB <- as.data.frame(read.csv("/mnt/f/Golden_boy_project/Master_Dir/FOXC2/TF_Link_DB/TFLink_Mus_musculus_bindingSites_All_annotation_v1.0.csv"))

unique(TF_Link_DB$Name.TF)

# Okay cool that gave us 140 TF to start off with this definetly isn't all of the TF out there though.

install.packages("rvest")
library(rvest)

# retrieving the target web page 
riken_TF <- read_html("http://gerg.gsc.riken.jp/TFdb/tf_list.html")

riken_TF_tbl <- as.data.frame(html_table(riken_TF, fill = TRUE))
riken_TF_gene <- riken_TF_tbl$X3[2:length(riken_TF_tbl$X3)]
length(riken_TF_gene)

#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
######
###### Heatmap Section
######
# Rational for this section. In the enrichment of the dCpGs for the PGCLCs we saw that there was a significant enrichment of transcription factor binding sites.
# We first wanna see if that TF is even expressed in the cell type we are interested in and then from there take the TFs that are expressed and look to see the number of EREs that are located proximal to those regions. 

#Lets take this object and add it to our DGEListgroup_1 object as a data frame and then add on the ENTREZID as a ne column
BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)
require(org.Mm.eg.db)
GeneSymbol <- mapIds(org.Mm.eg.db, keys = rownames(y_all), keytype="SYMBOL", column="ENTREZID")
GeneSymbol

y_all$gene <- data.frame(ENTREZID = GeneSymbol)
head(y_all$gene)
SYMBOL <- rownames(y_all$gene)
ENTREZID <- y_all$gene[,1] 
y_all$gene <- data.frame(Symbol = SYMBOL, ENTREZID = ENTREZID)

head(y_all$counts)

Mm_genes <- transcriptsBy(Mus.musculus, by="gene", columns=c("SYMBOL", "ENTREZID", "TXCHROM", "TXSTRAND"))
BiocManager::install("TxDb.Mmusculus.UCSC.mm39.refGene")
library(TxDb.Mmusculus.UCSC.mm39.refGene)
gene.length <- transcriptsBy(TxDb.Mmusculus.UCSC.mm39.refGene, by = "gene")

gene.length_2 <- unlist(gene.length)
gene.length_2$ENTREZID <- names(gene.length_2)

names(gene.length_2) <- gene.length_2$tx_name
gene.length <- relist(gene.length_2, gene.length)
gene.length.df <- as.data.frame(gene.length)
gene.length.df
gene.length.df <- gene.length.df[ -c(1:2) ]
install.packages("dplyr")
library(dplyr)
gene.length.df.2 = gene.length.df %>% group_by(ENTREZID) %>% top_n(n = 1, wt = width) %>% distinct( ENTREZID, .keep_all = TRUE)
gene.length.df.2
gene.length.df.2$length = gene.length.df.2$width

y_all$gene = y_all$gene %>% left_join(dplyr::select(gene.length.df.2, c("length", "ENTREZID")), by = "ENTREZID")

length(unique(y_all$gene$Symbol))

head(y$gene)
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#*************************************************************************************
install.packages("pheatmap")
library(pheatmap)
install.packages("RColorBrewer")
library(RColorBrewer)
# to plot heatmap you have to scale the data using the edgeR cpm() function so the difference in reads will be indicated by log(cpm).
DGEList.cpm <- cpm(y_all, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
DGEList.cpm
# genes to plot : 
#genes.x = unique(TF_Link_DB$Name.TF)
genes.x = unique(riken_TF_gene)

df.cpm = as.data.frame(DGEList.cpm)
df.cpm$Symbol <- y_all$gene$Symbol


#####


df.cpm.subset = filter(df.cpm, Symbol %in% genes.x)
y_all$samples

mex.cpm.subset = as.matrix(df.cpm.subset[, 1:6])
rownames(mex.cpm.subset) = df.cpm.subset$Symbol
colnames(mex.cpm.subset) = group_all
colnames(mex.cpm.subset)

# 
pdf("Heatmap_TM4_mm39.pdf", width = 10, height = 10)
pheatmap(mex.cpm.subset, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
         border_color = NA,
         cellwidth = NA, 
         cellheight = NA, 
         scale = "none", 
         display_numbers = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE, 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", 
         clustering_method = "complete",
         show_rownames = T, 
         show_colnames = T, 
         main = NA)
dev.off()
#///////////////////////////////////////////////////////////////////////////////

length(riken_TF_gene)

TF_SSC <- intersect(unique(TF_Link_DB$Name.TF), SSC_overlapping)

TF_SSC <- intersect(unique(riken_TF_gene), SSC_overlapping)
TF_SSC

TF_SSC_UP <- intersect(unique(riken_TF_gene), SSC_overlapping_UP)
TF_SSC_UP
#intersect(unique(TF_Link_DB$Name.TF), unique(unlist(overlap_gene_FOXC2$SYMBOL)))

intersect(TF_SSC, unique(unlist(overlap_gene_FOXC2$SYMBOL)))
# I can compare this to the JASPAR database of all TF binding sites to see the overlap

#///////////////////////////////////////////////////////////////////////////////
Mm_genes <- unlist(Mm_genes)

unlist(Mm_genes$SYMBOL)
gr_SSC_barcode_genes <- Mm_genes[unlist(Mm_genes$SYMBOL) %in% SSC_overlapping]

gr_SSC_barcode_genes <- Mm_genes[unlist(Mm_genes$SYMBOL) %in% SSC_overlapping_UP]

gr_SSC_barcode_genes
gr_SSC_pormoters <- promoters(gr_SSC_barcode_genes, upstream = 1000, downstream = 1000)
gr_SSC_pormoters
mm.genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
mm.genome

library(memes)

gr_SSC_pormoter_seq <- gr_SSC_pormoters %>%
  get_sequence(mm.genome)

gr_SSC_pormoter_seq
seq <- create_sequences(rng.seed = 100)
# Shuffle sequences preserving dinucleotide frequency
shuffle <- shuffle_sequences(seq, k = 2, rng.seed = 100)



dreme_out <- runDreme(gr_SSC_pormoter_seq, control = "shuffle", dna = TRUE, sec = 3600)
dreme_out <- runDreme(gr_SSC_pormoter_seq, control = "shuffle", dna = TRUE, nmotifs = 20)
dreme_out
dreme_out %>% 
  to_list() %>% 
  view_motifs()





JASPAR <- read_meme("/mnt/c/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/JASPAR.meme")
options(meme_db = read_meme("/mnt/c/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/JASPAR.meme"))



tomtomout <- runTomTom(dreme_out, database = JASPAR)
names(tomtomout)
tomtomout$best_match_altname
tomtomout$best_match_pval
tomtomout$best_match_qval
view_motifs(tomtomout$best_match_motif)
tomtomout$best_match_altname
tomtomout[36,]

################################################################################
# Okay based on the Single cell ananlysis of E18.5 We found there was a population of cells that had an elevated expression of a subset of DNA repair genes so I'm going to put that together and check them for motifs and pathways
DNA_Repair

DNA_repair_sub_E18.5 <- c('Parp2', 'Bard1', 'Brca1', 'Brca2', 'Chek1', 'Gtf2h4', 'Msh2', 'Nthl1', 'Pola1', 'Rad1', 'Rad50', 'Rfc1', 'Usp10', 'Parg', 'Chaf1a', 'Sfpq', 'Ercc8', 'Dtl', 'Fancm', 'Fanci', 'Usp1', 'Fancl', 'Atm', 'Huwe1', 'Topbp1', 'Smc1a', 'Atrx')

gr_SSC_barcode_genes <- Mm_genes[unlist(Mm_genes$SYMBOL) %in% DNA_repair_sub_E18.5]

gr_SSC_pormoters <- promoters(gr_SSC_barcode_genes, upstream = 1000, downstream = 1000)


gr_SSC_pormoter_seq <- gr_SSC_pormoters %>%
  get_sequence(mm.genome)

gr_SSC_pormoter_seq


dreme_out <- runDreme(gr_SSC_pormoter_seq, control = "shuffle", dna = TRUE, sec = 3600)
dreme_out <- runDreme(gr_SSC_pormoter_seq, control = "shuffle", dna = TRUE, nmotifs = 20)
dreme_out <- runDreme(gr_SSC_pormoter_seq, control = "shuffle", dna = TRUE)
dreme_out
dreme_out %>% 
  to_list() %>% 
  view_motifs()





JASPAR <- read_meme("/mnt/c/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/JASPAR.meme")
options(meme_db = read_meme("/mnt/c/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/JASPAR.meme"))



tomtomout <- runTomTom(dreme_out, database = JASPAR)
names(tomtomout)
tomtomout$best_match_altname
tomtomout$best_match_pval
tomtomout$best_match_qval
tomtomout
view_motifs(tomtomout$best_match_motif)
tomtomout$best_match_altname
tomtomout[36,]

enriched <- enrichr(unique(DNA_Repair$Final.list.cell.repair), dbs)


tmp <- enriched[['GO_Biological_Process_2023']]
tmp$Genes
write.csv(as.data.frame(decideTests(tmp)), 
          file=".csv")

#Something to try tomorrow

library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#gets gene symbol, transcript_id and go_id for all genes annotated with GO:0007507
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_transcript_id', 'go_id'),
                   filters = 'go_id', values = 'GO:0007507', mart = ensembl)
