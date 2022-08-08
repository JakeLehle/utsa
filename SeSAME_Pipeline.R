getwd()
setwd("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/")
setwd("/mnt/c/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/")
# Clean up the environment and remove all the sample object files
rm(list = ls())
# make sure you are using the most recent version of R
install.packages("installr")
library(installr)
updateR()
#/////////////////////////////////////////////////////////////////////////

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sesame")
browseVignettes("sesame")

## ----message=FALSE, warning=FALSE, include=FALSE------------------------------
options(rmarkdown.html_vignette.check_title = FALSE)

## -----------------------------------------------------------------------------
sapply(c("sesame","sesameData","ExperimentHub"),
       function(x) as.character(packageVersion(x)))
with(R.Version(), paste0(major, ".", minor))
# This outputs the version of sesame, sesame data, and experiment hub as well as R 
## ----message=FALSE------------------------------------------------------------
library(sesame)
# You need to run the sesameCache() command once after installation to make sure the data is retrieved 
# This caches the sesameData package with uses the ExperimentHub infrastructure.
library(knitr)
library(SummarizedExperiment)
library(ggrepel)
library(pals)
library(wheatmap)
install.packages("tidyr")
library(tidyr)
library(dplyr)
library(RPMM)
library(ggplot2)
if (!require(devtools)) install.packages("devtools")
devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")
if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")
library("ggvenn")
install.packages('gprofiler2')  
library(gprofiler2)

# If you are picking up from running this pipeline previous session.
load("IVERR_Data.RData")
## -----------------------------------------------------------------------------
tools::R_user_dir("ExperimentHub", which="cache")
# This shows you location of the cached data on your computer. 

sesameDataCacheAll()

######## EXPERIMENTAL DATA FROM MOUSE

# First create an object with the path to the IDAT file
idat_dir = c("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/206011040054/IDAT")
idat_dir
IDATprefixes <- searchIDATprefixes(idat_dir)
basename(IDATprefixes)

######## QC 

# This section will show you how to look at certain aspects of the pre-processing steps and produce graphs that visualize the changes that you are making.
# These graphs are what you should use in a supplemental section to show the corrections that you are doing to the data.
# This is what the final processing step will look like lets go over what each of the letters in the "prep" step indicate.

sdf_preped = openSesame(IDATprefixes, BPPARAM = BiocParallel::SnowParam(8), prep="TQCDPB", func=NULL)
sdf_preped
# T: Infer strain, set mask.
# Q: Mask probes of poor design.
# C: Infer infinium channel intensity.
# D: Dye bias correction (non-linear)
# P: Detection of p-value masking using oob
# B: Equivalent to noob() background subtraction using oob

# Now that you know what each step of the final process will be lets go through steps and see the graphs that you can make from each 

# Open the file and produce a beta matrix and sdfs list

betas = openSesame(IDATprefixes, BPPARAM = BiocParallel::SnowParam(8))
sdfs = openSesame(IDATprefixes, BPPARAM = BiocParallel::SnowParam(8), func = NULL)

# Infer mouse strains from single sdf IDAT pairs. This is equivalent to the prep = "T" step.

sample_strains <- list()

for (strain in IDATprefixes) {
  sdf <- readIDATpair(strain, platform = "MM285")
  inferStrain(sdf, return.strain = TRUE)
  sample_strains <- append(sample_strains, inferStrain(sdf, return.strain = TRUE))
}
sample_strains

# If any of the strains look weird to you then you can manually check them and plot the probability of each strain in the sample.
# Lets pick two samples (1 control and 1 experimental) and check them individually.

sdf_1 <- readIDATpair("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/206011040054/IDAT/TM4_EtOH_1", platform = "MM285")
sdf_2 <- readIDATpair("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/206011040054/IDAT/TM4_BPS_200uM_1", platform = "MM285")

p = inferStrain(sdf_1, return.probability = TRUE)
df_sdf1 = data.frame(strain=names(p), probs=p)
ggplot(data = df_sdf1,  aes(x = strain, y = probs)) +
  geom_bar(stat = "identity", color="gray") +
  ggtitle("Strain Probabilities") +
  ylab("Probability") + xlab("") +
  scale_x_discrete(position = "top") +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=0),
        legend.position = "none")

# You can also use this comparison tool but this is much harder to tell what is going on here the BALB cj color lines up both on the column and row indicating goo alignment. However other strains alos have alignment so it's not the best way to display this information.
compareMouseStrainReference(getBetas(sdf_1))

#Looks like BALB_cj is the primary background for all of these samples which is good that it is consistent.
# Fascinating stuff to keep in mind for your own studies to try and maintain high similarities in the backgrounds of your mouse samples or at leats bring it up when discussing results.
# You can use this structure to also determine the ethnicity of Human samples which should be discussed in projects especially in the context of studying methylation patterns in underrepresented populations. 

# You can even take this information further and predict the average age of the samples in months from the beta values
# Haha these immortalized cells have markers that indicate they are pretty old which isn't a stretch by any means (15 months). 

predictMouseAgeInMonth(betas[,1])

# Sample preprocessing and quality control

sdf_preped = openSesame(IDATprefixes, BPPARAM = BiocParallel::SnowParam(8), prep="TQCDPB", func=NULL)
sdf_preped


# Using our sdf's from sample 1 and 2 lets look at the probe signal clarity and subtract background. 
# Fist let's infer the signal intensity from the sample to avoid bias coming from different iscan devices. This is equivalent to prep="C"
sdf_1.InfICorrected = inferInfiniumIChannel(sdf_1, verbose=TRUE)

# Let's look at the dye bias for this experiment

##### Dye Bias Correction 
# Residual dye bias can be corrected using nonlinear quantile interpolation with type-1 probes.
# Note non linear scaling does shift the beta values of type 1 probes keep that in mind if someone wants to see the raw un-scaled data.
par(mfrow=c(1,2))
sesameQC_plotRedGrnQQ(sdf_1.InfICorrected, main="Before")   
sesameQC_plotRedGrnQQ(dyeBiasL(sdf_1.InfICorrected), main="After")  # linear correction

# Go with nonlinear for better adjustment
par(mfrow=c(1,2))
sesameQC_plotRedGrnQQ(sdf_1, main="Before")
sesameQC_plotRedGrnQQ(dyeBiasNL(sdf_1.InfICorrected), main="After")  # nonlinear correction

sdf_1.InfICorrected <- dyeBiasNL(sdf_1.InfICorrected)

###### Background Subtraction
# Let's now look at background subtraction using oob. Equivalent to prep="B" 

par(mfrow=c(2,1), mar=c(3,3,2,1))
sesameQC_plotBetaByDesign(sdf_1.InfICorrected, main="Before", xlab="\beta")
sesameQC_plotBetaByDesign(noob(sdf_1.InfICorrected), main="After", xlab="\beta")

# Good the data has a very clean signal with only minor sites that have a mixture of both green and red probes.
# Let's look at some other quality control metrics that can give us an idea of how good our data is.

#QC of detection metric
qcs = openSesame(IDATprefixes, BPPARAM = BiocParallel::SnowParam(8), prep="", func=sesameQC_calcStats, funs="detection")
qcs_prep = openSesame(IDATprefixes, BPPARAM = BiocParallel::SnowParam(8), prep="TQCDPB", func=sesameQC_calcStats, funs="detection")

qcs
qcs_prep

sesameQC_plotBar(qcs_prep)

# The detection metrics are good as well around 89%. That's exactly what you want to see.

########### DML/DMR Analysis

# Let's start our DML(Differential Methylation Level) analysis
# To do this we need to create a summarized experiment object
# The object itself will be built on our beta values calculated earlier
# We will create a coldata dataframe object that will divide our object by IDAT file name, treatment group of the sample, and treatment concentration. 
# Finally I reassign the name of the assay so it will indicate it is recording beta values.
betas_prep = openSesame(IDATprefixes, BPPARAM = BiocParallel::SnowParam(8), prep="TQCDPB")
betas_prep
IDAT <- basename(IDATprefixes)
IDAT
treatment <- c("BPS_100uM", "BPS_100uM","BPS_100uM", 
               "BPS_1uM", "BPS_1uM", "BPS_1uM", 
               "BPS_200uM", "BPS_200uM", "BPS_200uM",
               "EtOH", "EtOH", "EtOH")
#condition <- c("treatment", "treatment", "treatment", "treatment", "treatment", "treatment", "treatment", "treatment", "treatment", "control", "control", "control")
df <- data.frame(IDAT, treatment)
se <- SummarizedExperiment(assays = betas_prep, colData = df)
names(assays(se)) <- c("betas")
cd = as.data.frame(colData(se)); rownames(cd) = NULL

# To check DML you have to exclude NA's use the checkLevels() function

se_ok = (checkLevels(assay(se), colData(se)$treatment))
sum(se_ok)                      # the number of CpGs that passes
se = se[se_ok,]

# Here they set the control treatment as that with only the vehicle EtOH
colData(se)$treatment <- relevel(factor(colData(se)$treatment), "EtOH")

#Now we model DNA methylaton variation treating treatment and condition as covariates. 
#The DML will fit the DNA methylation reading to a linear model and perform the corresponding slope test and goodness-of-fit test (F-test holding out each contrast variable)
#All of these results are returned in an object of class DMLSummary
smry = DML(se, ~treatment, mc.cores = BiocParallel::SnowParam(8))
save.image("IVERR_Data.RData")

#### TEST INERPRETATION
test_result = summaryExtractTest(smry)
# Rows are CpG/Loci and columns contain the slopes and p-values for each variable

colnames(test_result) # the column names, show four groups of statistics
#EST The slope estimate (aka the β coefficient, not to be confused with the DNA methylation β-value though) for continuous variable. DNA methylation difference of the current level with respect to the reference level for nominal contrast variables. Each suffix is concatenated from the contrast variable name (e.g., treatment) and the level name if the contrast variable is discrete (e.g, EtOH, BPS_1uM, BPS_100uM, BPS_200uM). For example, Est_BPS1uM should be interpreted as the estimated methylation level difference of cells treated with 1uM BPS compared to the reference group (which is EtOH, as set above). There is a special column named Est_`(Intercept)`. It corresponds to the base-level methylation of the reference (in this case a EtOH sample)
#PVAL The unadjusted p-values of t-testing the slope. This represents the statistical significance of the methylation difference. For example, Pval_treatmentBPS_1uM tests whether the cells treated with 1uM BPS is significantly different from the cells treated only with EtOH (the reference level) in DNA methylation. The Pval_`(Intercept)` tests whether the reference level is significantly different from zero.
#FPVAL The unadjusted p-value of the F-test contrasting the full model against a reduced model with the labeled contrast variable held out. Note that “Pval_” and “FPval_” are equivalent when the contrast variable is a 2-level factor, i.e., in the case of a pairwise comparison.
#EFF The effect size of each normial contrast variable. This is equivalent to the maximum slope subtracted by the minimum level including the reference level (0).
head(test_result)

### For your project PVAL will be the most important metric to see here.
###### GOODNESS Of FIT
# Another way to say this is.
# Is the CpG methylation treatment-specific? Rather than. Is the CpG more methylated in treated TM4 compared to normal TM4?
# Here we used 0.05 as the size threshold. This means a difference less than 5% is not considered biologically relevant and excluded. 
# We can go further and do a side by comparison of probes that are treatment specific and ones that are condition specific. 

# Here I use some dplyr functions to find out the number of probes that are have differences in DNA methylation with low p-values and high effect values. Because my group size is so small I expect this to only be a few probes.
test_result %>%
  mutate(treatment_specific =
           ifelse(FPval_treatment < 0.05 & Eff_treatment > 0.05, TRUE, FALSE)) %>%
  select(treatment_specific) %>% table
# Yep looks like only 700 probes in the whole data set pass this threshold. Beacuse of that I'll simply focus on probes that have a smaller than 0.05 p-value untill I can include more data into the experiment and increase the effect score.

# Now we can plot this difference with a volcano plot to show the differences in probes for treated vs control

ggplot(test_result)+ geom_point(aes(Est_treatmentBPS_1uM, -log10(Pval_treatmentBPS_1uM)))
ggplot(test_result)+ geom_point(aes(Est_treatmentBPS_100uM, -log10(Pval_treatmentBPS_100uM)))
ggplot(test_result)+ geom_point(aes(Est_treatmentBPS_200uM, -log10(Pval_treatmentBPS_200uM)))

# Let's make it pretty but focus on plotting the false discovery rate or FPval for stage specific changes as it relates to changes in beta values.
test_result
test_result_pretty <- test_result
test_result_pretty$Methylation
# Let's convert the value of the beta over into Hyper or Hypo and store it in this Methylation column and assign a significance to the change as well in the significance column.

test_result_pretty <- test_result_pretty %>%
  mutate(
    Significance = case_when(
      FPval_treatment <= 0.05 & FPval_treatment > 0.01 ~ "FDR 0.05",
      FPval_treatment <= 0.01 & FPval_treatment > 0.001 ~ "FDR 0.01",
      FPval_treatment <= 0.001 ~ "FDR 0.001",
      TRUE ~ "Unchanged"),
    Methylation = case_when(
      Est_treatmentBPS_100uM > 0 ~ "Hyper",
      Est_treatmentBPS_100uM < 0 ~ "Hypo",
      TRUE ~ "Unchanged"),
    #Zero = case_when(
    # Est_treatmentBPS_100uM < 0.01 & Est_treatmentBPS_100uM > 0.01 ~ "zero",
    #TRUE ~ "Changed")
  )

test_result_pretty
###### Side note. Quick way to find out how many of the probes are either hyper or hypo methylated in treated samples compared to controls.
test_result_pretty %>% 
  count(Methylation)
# Now with significance added for significant changes specific changes 
test_result_pretty %>%
  count(Methylation, Significance)

install.packages("rlang")
library("rlang")
devtools::install_github("cardiomoon/moonBook")
devtools::install_github("cardiomoon/webr")
require(ggplot2)
require(moonBook)
require(webr)
PieDonut(test_result_pretty, aes(pies=Methylation,donuts=Significance), start=3*pi/3, title="DNA METHYLATION CHANGE IN TM4 CELLS 
         FOLLOWING 100uM BPS", ratioByGroup=FALSE)

# 3D Exploded Pie Chart 
install.packages('plotrix')
library(plotrix)
slices <- c(57, 43)
lbls <- c("Hypo: 57%", "Hyper: 43%")
pie3D(slices,labels=lbls,explode=0.35,
      theta=pi/3,
      main="BPS 200uM",
      col=c('cornflowerblue', "brown2"),
      mar = c(0.5,1,3,0.5))

# Let's find out the higest and lowest changes in methylation in our groups

top_10 <- 10
top_probes <- bind_rows(
  test_result_pretty %>%
    filter(Methylation == 'Hyper') %>%
    filter(Est_treatmentBPS_100uM > 0.025) %>%
    arrange(Pval_treatmentBPS_100uM) %>%
    head(top_10),
  test_result_pretty %>%
    filter(Methylation == 'Hypo') %>%
    filter(Est_treatmentBPS_100uM < -0.05) %>%
    arrange(Pval_treatmentBPS_100uM) %>%
    head(top_10)
)

top_probes

ggplot(test_result_pretty, aes(Est_treatmentBPS_100uM, -log10(Pval_treatmentBPS_100uM))) +
  geom_point(aes(color = Significance), size = 2/5) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  xlim(-0.25, 0.2) +
  geom_label_repel(data = top_probes,
                   mapping = aes(Est_treatmentBPS_100uM, -log(Pval_treatmentBPS_100uM,10), label = Est_Probe_ID),
                   size = 2, max.overlaps = 40)

# Dr McCarrey wanted me to write in a section where I test the package and generate some false datasets that have no significant differences between the groups to see what that plot would look like so here is that.
rm(betas_prep_scramble)
betas_prep
betas_prep_scramble <- betas_prep[1:500,]
betas_prep_scramble[1:250,] <- runif(12, min=0.85, max=1)
betas_prep_scramble[251:500,] <- runif(12, min=0.25, max=0.31)
betas_prep_scramble
se_scramble <- SummarizedExperiment(assays = betas_prep_scramble, colData = df)
names(assays(se_scramble)) <- c("betas")
cd = as.data.frame(colData(se_scramble)); rownames(cd) = NULL
se_ok_scramble = (checkLevels(assay(se_scramble), colData(se_scramble)$treatment))
sum(se_ok_scramble)                      # the number of CpGs that passes
se_scramble = se_scramble[se_ok_scramble,]
colData(se_scramble)$treatment <- relevel(factor(colData(se_scramble)$treatment), "EtOH")
smry_scramble = DML(se_scramble, ~treatment, mc.cores = BiocParallel::SnowParam(8))
test_result_scramble = summaryExtractTest(smry_scramble)
test_result_scramble
test_result_scramble <- test_result_scramble %>%
  mutate(
    Significance = case_when(
      FPval_treatment <= 0.05 & FPval_treatment > 0.01 ~ "FDR 0.05",
      FPval_treatment <= 0.01 & FPval_treatment > 0.001 ~ "FDR 0.01",
      FPval_treatment <= 0.001 ~ "FDR 0.001",
      TRUE ~ "Unchanged"),
    Methylation = case_when(
      Est_treatmentBPS_100uM > 0 ~ "Hyper",
      Est_treatmentBPS_100uM < 0 ~ "Hypo",
      TRUE ~ "Unchanged")
  )
ggplot(test_result_scramble, aes(Est_treatmentBPS_100uM, -log10(Pval_treatmentBPS_100uM))) +
  geom_point(aes(color = Significance), size = 2/5) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  xlim(-0.25, 0.2) +
  ylim(0, 6)


#Okay now that you have the 10 most significant hyper and hypo methylated probes you can visualize them on the genome using this function

visualizeProbes(c("cg42515012_BC21"), betas_prep, platform='MM285')

# Oh wow that probe is located in the promoter of Bptf the bromodomain protein that has been linked to Alzhymers 
# Let's call that region out by gene now.

# Hyper
visualizeGene('Mmrn2', betas_prep, platform='MM285')
visualizeGene('Trim72', betas_prep, platform='MM285')
visualizeGene('Septin9', betas_prep, platform='MM285')
# Hypo
visualizeGene('Elovl6', betas_prep, platform='MM285')
visualizeGene('Mat2b', betas_prep, platform='MM285')
visualizeGene('Otud6a', betas_prep, platform='MM285')
visualizeGene('Fam219a', betas_prep, platform='MM285')
visualizeGene('Vmn1r72', betas_prep, platform='MM285')
visualizeGene('Dram2', betas_prep, platform='MM285')


# Another way to visulize this trend is as a dot plot here I show you how to do that so we can get something similar to the box and wisker plots I was making with the methylflash data.

df = data.frame(treatment = colData(se)$treatment,
                BetaValue = assay(se)[rownames(test_result)[nrow(test_result)],])


ggplot(df, aes(treatment, BetaValue)) + geom_smooth(method="lm") + geom_point()

# Based on all of these findings it looks like 100uM is the treatment group with the most significant changes and what we should use going forward for analysis.


###### Use this result to make a Venn diagram with ggvenn

test_result_probes_df_1uM <- test_result %>%
  mutate(treatment_specific =
           ifelse(Pval_treatmentBPS_1uM < 0.05, TRUE, FALSE)) %>%
  select(treatment_specific) %>% DataFrame()

test_result_probes_df_100uM <- test_result %>%
  mutate(treatment_specific =
           ifelse(Pval_treatmentBPS_100uM < 0.05, TRUE, FALSE)) %>%
  select(treatment_specific) %>% DataFrame()

test_result_probes_df_200uM <- test_result %>%
  mutate(treatment_specific =
           ifelse(Pval_treatmentBPS_200uM < 0.05, TRUE, FALSE)) %>%
  select(treatment_specific) %>% DataFrame()

# Remove NA's from data frames

treatment_specific_rn_1uM <- na.omit(rownames(test_result_probes_df_1uM)[test_result_probes_df_1uM[["treatment_specific"]] == TRUE])
treatment_specific_rn_100uM <- na.omit(rownames(test_result_probes_df_100uM)[test_result_probes_df_100uM[["treatment_specific"]] == TRUE])
treatment_specific_rn_200uM <- na.omit(rownames(test_result_probes_df_200uM)[test_result_probes_df_200uM[["treatment_specific"]] == TRUE])

# Put everything together and make the graph

specific_list <- list(
  "BPS 1uM" = treatment_specific_rn_1uM,
  "BPS 100uM" = treatment_specific_rn_100uM,
  "BPS 200uM" = treatment_specific_rn_200uM
) 

ggvenn(specific_list,
       fill_color = c("#0073C2FF", "#868686FF", "#DC0000B2"),
       set_name_size = 8,
       text_size = 6)

######## End of ggvenn section

#Use this list to find out if any class or probes are enriched in the population
treatment_specific_100uM <- test_result %>% dplyr::filter(Pval_treatmentBPS_100uM < 0.05) %>% select(Pval_treatmentBPS_100uM)

# Let's see if there an enrichment of probes in the stage_specific differential methylated loci regions within one of the prob groups
query <- rownames(treatment_specific_100uM)
query
test_result
results_treatment_specific <- testEnrichment(query, "chromHMM", platform = "MM285")
KYCG_plotEnrichAll(results_treatment_specific)
# Oh this is interesting looks like there are a number of groups from either the design group or the chromHMM group that we should look more into.
# Before that let's rank the enrichment between groups.
KYCG_plotDot(results_treatment_specific, n_max=20)
KYCG_plotBar(results_treatment_specific, n_max=20)

treatment_specific_100uM_beta <- test_result %>% dplyr::filter(Pval_treatmentBPS_100uM < 0.05) %>% select(Est_treatmentBPS_100uM)
treatment_specific_100uM$probe <- rownames(treatment_specific_100uM)
treatment_specific_100uM_beta$probe <- rownames(treatment_specific_100uM_beta)
treatment_specific_100uM <- merge(x = treatment_specific_100uM, y = treatment_specific_100uM_beta, by = "probe", all = TRUE)
treatment_specific_100uM_FPval <- test_result %>% dplyr::filter(Pval_treatmentBPS_100uM < 0.05) %>% select(FPval_treatment)
treatment_specific_100uM_FPval$probe <- rownames(treatment_specific_100uM_FPval)
treatment_specific_100uM <- merge(x = treatment_specific_100uM, y = treatment_specific_100uM_FPval, by = "probe", all = TRUE)
treatment_specific_100uM
treatment_specific_100uM <- treatment_specific_100uM %>%
  mutate(
    Significance = case_when(
      FPval_treatment <= 0.05 & FPval_treatment > 0.01 ~ "FDR 0.05",
      FPval_treatment <= 0.01 & FPval_treatment > 0.001 ~ "FDR 0.01",
      FPval_treatment <= 0.001 ~ "FDR 0.001",
      TRUE ~ "Unchanged"))

anno_design <- KYCG_annoProbes(query, "designGroup", silent = TRUE)
anno_design <- as.data.frame(anno_design)
colnames(anno_design) <- c("Annotation")
anno_design$probe <- rownames(anno_design)

treatment_specific_100uM <- merge(x = treatment_specific_100uM, y = anno_design, by = "probe", all = TRUE)
anno_chrom <- KYCG_annoProbes(query, "chromHMM", silent = TRUE)
anno_chrom <- as.data.frame(anno_chrom)
colnames(anno_chrom) <- c("Annotation")
anno_chrom$probe <- rownames(anno_chrom)
treatment_specific_100uM <- merge(x = treatment_specific_100uM, y = anno_chrom, by = "probe", all = TRUE)

treatment_specific_100uM

probe_anno <- read.csv(file="MouseMethylation-12v1-0_A2.csv")
probe_anno[7,]
colnames(probe_anno) <- c("probe", "Name", "AddressA_ID", "AlleleA_ProbeSeq", "AddressB_ID", "AlleleB_ProbeSeq", "Next_Base", "Color_Channel", "Col", "Probe_Type", "Strand",
                          "Strand_TB", "Strand_CO", "Infinium_Design_Type", "Rep_Num",  "CHR", "MAPINFO", "Species", "Genome_Build", "SourceSeq", "Forward_Sequence", "Top_Sequence",
                          "Genome_Build_NCBI", "N_Shelf", "N_Shore", "CpG_Island", "CpG_Island_chrom", "CpG_Island_chromStart", "CpG_Island_chromEnd", "CpG_Island_length",
                          "CpG_Island_cpgNum", "CpG_Island_gcNum", "CpG_Island_perCpg", "CpG_Island_perGc", "CpG_Island_obsExp", "S_Shore", "S_Shelf", "MFG_Change_Flagged")

treatment_specific_100uM <- merge(x = treatment_specific_100uM, y = probe_anno, by = "probe", all = TRUE)

probe_anno_gene <- read.csv(file="MouseMethylation-12v1-0_A1_Annotation_Mus_musculus.csv")
colnames(probe_anno_gene) <- c("Name", "Gene", "Transcript", "chrom", "chromStart", "chromEnd", "chromLength", "chromStrand", "Feature", "Source")

treatment_specific_100uM <- merge(x = treatment_specific_100uM, y = probe_anno_gene, by = "Name", all = TRUE)

treatment_specific_100uM_final <- treatment_specific_100uM %>%
  drop_na(Pval_treatmentBPS_100uM)

treatment_specific_100uM_final <- treatment_specific_100uM_final %>% distinct()
length(treatment_specific_100uM_final$Name)
# Okay now that we have built our annotation table for each of our treatment specific probes that were changed compared to the controls lets look at and filter oput each type of probe and start getting an idea of each of our groups that are enriched for each probe.
print(unique(treatment_specific_100uM_final$Annotation.x))
print(unique(treatment_specific_100uM_final$Annotation.y))

Enh_design <- treatment_specific_100uM_final %>%
  filter(Annotation.x == "Enhancer")

Enh <- treatment_specific_100uM_final %>%
  filter(Annotation.y == "Enh")
EnhPois <- treatment_specific_100uM_final %>%
  filter(Annotation.y == "EnhPois")
Enh_Joined <- rbind(Enh,EnhPois)
Enh_Joined
Gene_Enh_Joined <- print(unique(Enh_Joined$Gene))
write.table(Gene_Enh_Joined, file = "Gene_Enh_Joined.txt", sep = "\t",
            row.names = FALSE)
print(unique(treatment_specific_100uM_final$Annotation.x))

treatment_specific_100uM_final
FDR_0.05 <- treatment_specific_100uM_final %>%
  filter(Significance == "FDR 0.05")
FDR_0.01 <- treatment_specific_100uM_final %>%
  filter(Significance == "FDR 0.01")
FDR_0.001 <- treatment_specific_100uM_final %>%
  filter(Significance == "FDR 0.001")

# Okay now let's start building some chromosome maps so we can display all the sites where there are changes. 

install.packages("RIdeogram")
require(RIdeogram)
# Load in some sample data.

# Check the data
mouse_karyotype <- data.frame(
  Chr = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"),
  Start = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  End = c(195154279, 181755017, 159745316, 156860686, 151758149, 149588044, 144995196, 130127694, 124359700, 130530862, 121973369, 120092757, 120883175, 125139656, 104073951, 98008968, 95294699, 90720763, 61420004, 169476592, 91455967),
  CE_start = rep(110000, 21),
  CE_end = rep(300000, 21)
)


probes <- data.frame(
  Name = treatment_specific_100uM_final$Name,
  Chr = treatment_specific_100uM_final$CHR,
  Start = as.integer(treatment_specific_100uM_final$MAPINFO),
  End = as.integer(treatment_specific_100uM_final$MAPINFO) + 1,
  Value = treatment_specific_100uM_final$Est_treatmentBPS_100uM
)


probes_Enh_design <- data.frame(
  Name = Enh_design$Name,
  Chr = Enh_design$CHR,
  Start = as.integer(Enh_design$MAPINFO),
  End = as.integer(Enh_design$MAPINFO) + 1,
  Value = Enh_design$Est_treatmentBPS_100uM
)

probes_Enh <- data.frame(
  Name = Enh_Joined$Name,
  Chr = Enh_Joined$CHR,
  Start = as.integer(Enh_Joined$MAPINFO),
  End = as.integer(Enh_Joined$MAPINFO) + 1,
  Value = Enh_Joined$Est_treatmentBPS_100uM
)

probes <- probes %>% distinct()
probes
probes_no_name <- select(probes, -c("Name"))
probes_no_name <- probes_no_name %>% arrange(probes_no_name$End, by_group = TRUE)
print(unique(probes_no_name$Chr))
probes_no_name <- probes_no_name[probes_no_name$Chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"),]
probes_hypo <- probes_no_name %>% filter(Value < 0)
probes_hyper <- probes_no_name %>% filter(Value > 0)
length(probes_no_name$Chr)
probes_no_name$Value

# Hypo
ideogram(karyotype = mouse_karyotype, overlaid = probes_hypo, colorset1 = c("#0d4fd4"))
convertSVG("chromosome.svg", device = "png")
svg2pdf("chromosome.svg")
# Hyper
ideogram(karyotype = mouse_karyotype, overlaid = probes_hyper, colorset1 = c("#f21b1b"))
convertSVG("chromosome.svg", device = "png")
svg2pdf("chromosome.svg")


probes_Enh <- probes_Enh %>% distinct()
probes_Enh
probes_no_name_Enh <- select(probes_Enh, -c("Name"))
probes_no_name_Enh <- probes_no_name_Enh %>% arrange(probes_no_name_Enh$End, by_group = TRUE)
print(unique(probes_no_name_Enh$Chr))
probes_no_name_Enh <- probes_no_name_Enh[probes_no_name_Enh$Chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"),]
length(probes_no_name_Enh$Chr)

ideogram(karyotype = mouse_karyotype, overlaid = probes_no_name_Enh, colorset1 = c("#11f0e8"))
convertSVG("chromosome.svg", device = "png")
svg2pdf("chromosome.svg")

probes_Enh_design <- probes_Enh_design %>% distinct()
probes_no_name_Enh_design <- select(probes_Enh_design, -c("Name"))
probes_no_name_Enh_design <- probes_no_name_Enh_design %>% arrange(probes_no_name_Enh_design$End, by_group = TRUE)
print(unique(probes_no_name_Enh_design$Chr))
probes_no_name_Enh_design <- probes_no_name_Enh_design[probes_no_name_Enh_design$Chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"),]
length(probes_no_name_Enh_design$Chr)

ideogram(karyotype = mouse_karyotype, overlaid = probes_no_name_Enh_design, colorset1 = c("#f5ae0a"))
convertSVG("chromosome.svg", device = "png")
svg2pdf("chromosome.svg")



FDR_0.001 <- data.frame(
  Name = FDR_0.001$Name,
  Chr = FDR_0.001$CHR,
  Start = as.integer(FDR_0.001$MAPINFO),
  End = as.integer(FDR_0.001$MAPINFO) + 1,
  Value = FDR_0.001$Est_treatmentBPS_100uM
)
FDR_0.001 <- FDR_0.001 %>% distinct()
FDR_0.001 <- select(FDR_0.001, -c("Name"))
probes_no_name_FDR <- FDR_0.001 %>% arrange(FDR_0.001$End, by_group = TRUE)
print(unique(probes_no_name_FDR$Chr))
probes_no_name_FDR <- probes_no_name_FDR[probes_no_name_FDR$Chr %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "X", "Y"),]
length(probes_no_name_FDR$Chr)

ideogram(karyotype = mouse_karyotype, overlaid = probes_no_name_FDR, colorset1 = c("#25075c"))
convertSVG("chromosome.svg", device = "png")
svg2pdf("chromosome.svg")

# Very interesting this looks like there enrichment from enhancers, enhancer and enhancer poised chromatin, gene body, H3K36me3. 
# Let's try to pull out the designGroup and chromHMM data sets and try to annotate them to find out some more information about them.
results_treatment_specific_designGroup <- testEnrichment(query, "chromHMM", platform = "MM285")
results_treatment_specific_designGroup

cg_lists <- KYCG_getDBs("MM285.designGroup")
cg_lists
queries <- cg_lists[(sapply(cg_lists, length) > 40000)]
result_list <- lapply(cg_lists, results_treatment_specific)


# How cools is that this shows that there is a enrichment of methylation diffrences at Imprinted DMRs which has ben shown previously as well as a number of other sites like H19, Gene Bodies, Enhancers, Pesudogene TSS's, etc
results_treatment_specific
results_treatment_specific_filtered <- results_treatment_specific %>% dplyr::filter(overlap>10)
KYCG_plotEnrichAll(results_treatment_specific_filtered)
# Let's rank these by significance
KYCG_plotWaterfall(results_treatment_specific)
KYCG_plotWaterfal
# To dig more into each group we need to see what groups are available which can be done using the list database groups function.
KYCG_listDBGroups("MM285")
# Now that we can look at individual groups we can start pulling apart some of the data we are seeing from our enrich all data
# First let's do a dot plot of the design group to see which of the design groups are most differentially methylated
results_treatment_specific_designGroup <- testEnrichment(query, "designGroup", platform = "MM285")


# Next to find out more about if a region is hyper or hypo methylated lets to a volcano plot this time using the chromHMM group.
results_2tailed <- testEnrichment(query, alternative = "two.sided", platform = "MM285")
KYCG_plotVolcano(results_2tailed)

# You can see from this that there is quite a bit of data that you can already get and this might help you to determine what group to look into more. 
# However let's go back and look at all of our probes in the DML group and find out the ones that have the largest differences.




# Here we can model these functions as continuous. This shows CpGs that are positively associated with length of culture conditions.
smry2 = DML(se, ~condition + treatment, mc.cores = BiocParallel::SnowParam(8))
test_result2 = summaryExtractTest(smry2) %>% arrange(Est_stageblastocyst)

test_result2 %>% dplyr::select(Est_stageblastocyst, Pval_stageblastocyst) %>% tail # positive assoc.


###### DMR
# Okay so the way that this works is if there are significant differences in a group of probes that are closely located together they can be merged into a region that they call a DMR.
# You will then have to annotate these regions seperately to figure out if they are in promoter regions or gene regions use the wg-blimp script for this.
dmContrasts(smry)                       # pick a contrast from below
merged = DMR(se, smry, "treatmentBPS_100uM")       # merge CpGs to regions / segments
DMR <- merged %>% dplyr::filter(Seg_Pval_adj < 0.01) %>% arrange(Seg_Pval_adj)
DMR

# Here is where you can call out specific regions and get the probe names
visualizeRegion(
  'chr3',129540509,129621090, betas_prep, platform='MM285',
  show.probeNames = TRUE)
#Turns out this is a DMR for a zinc finger binding protein that is important for mRNA processing.

######## GO Analysis 

#Start with a probe list. For this example let's use the probes that are specific for blastocysts.

BPS_100um <- sesameData_getGenesByProbes(treatment_specific_rn_100uM)
BPS_100um

#  ## use gene name
gostres <- gost(BPS_100um$gene_name, organism = "mmusculus")
gostres$result[order(gostres$result$p_value),]
gostplot(gostres)
#  
## use Ensembl gene ID, note we need to remove the version suffix
gene_ids <- sapply(strsplit(names(blastocysts_genes),"\\."), function(x) x[1])
gostres <- gost(gene_ids, organism = "mmusculus")
gostres$result[order(gostres$result$p_value),]
gostplot(gostres)

# Both will give you the same thing.

# Pulling in the regionAnnotation and dmrAnnotation script from wg-blimp

library(data.table)
library(GenomicRanges)
BiocManager::install("rtracklayer")
library(rtracklayer)
library(tidyverse)

annotation.isExistingFileOrNone <- function(fileName) {
  
  if (is.null(fileName)) {
    return(FALSE)
  }
  
  if (length(fileName) == 1 && file.exists(fileName)) {
    return(TRUE)
  }
  
  stop(paste0("Annotation file not found: ", fileName))
  
}

annotation.annotateOverlap <- function (regionRanges, overlapRanges, geneTable) {
  
  intersection <- as.list(findOverlaps(regionRanges, overlapRanges))
  isOverlapping <- sapply(intersection, length) > 0
  overlappingGeneNames <- sapply(intersection, function (indices) paste(unique(geneTable$external_gene_name[indices]), collapse = ","))
  
  return(data.table(
    isOverlapping,
    overlappingGeneNames
  ))
}

annotation.annotateGenes <- function (regions, regionRanges, gtfRanges) {
  
  geneGtfRanges = gtfRanges[gtfRanges$type == "gene"]
  
  geneNames <- data.table(external_gene_name = geneGtfRanges$gene_name)
  
  regions[,c("gene_overlap", "gene_name")] <- annotation.annotateOverlap(regionRanges, geneGtfRanges, geneNames)
  
  return(regions)
}

annotation.annotateCGIslands <- function (regions, regionRanges, cgis) {
  
  cgiRanges  <- GRanges(seqnames = cgis$chrom, ranges = IRanges(cgis$chromStart, cgis$chromEnd))
  
  regions$cgi_overlap <- countOverlaps(regionRanges, cgiRanges) > 0
  
  return(regions)
}

annotation.annotatePromoters <- function (regions, regionRanges, gtfRanges, promoterTSSDistances) {
  
  transcriptRanges = gtfRanges[gtfRanges$type == "transcript"]
  geneNames = data.table(external_gene_name = transcriptRanges$gene_name)
  
  promoterStarts <- ifelse(
    strand(transcriptRanges) == "+",
    start(transcriptRanges) + promoterTSSDistances[1],
    end(transcriptRanges) - promoterTSSDistances[2]
  )
  
  promoterEnds <- ifelse(
    strand(transcriptRanges) == "+",
    start(transcriptRanges) + promoterTSSDistances[2],
    end(transcriptRanges) - promoterTSSDistances[1]
  )
  
  promoterRanges <- GRanges(seqnames = seqnames(transcriptRanges), ranges = IRanges(promoterStarts, promoterEnds))
  
  regions[,c("promoter_overlap", "promoter_name")] <- annotation.annotateOverlap(regionRanges, promoterRanges, geneNames)
  
  return(regions)
}

annotation.annotateRepeats <- function (regions, regionRanges, repeats, classes) {
  
  relevantRepeats <- repeats[repClass %in% classes]
  
  repeatRanges <- GRanges(relevantRepeats$genoName, IRanges(relevantRepeats$genoStart, relevantRepeats$genoEnd))
  
  regions$num_repeats <- countOverlaps(regionRanges, repeatRanges)
  
  return(regions)
}

annotation.annotateRegions <- function (regionTable, gzippedCgiFile, gzippedGTFFile, gzippedRepeatMaskerAnnotationFile, allowedBiotypes, promoterTSSDistances) {
  
  regionTable$length <- regionTable$end - regionTable$start
  
  # remove any trailing 'chr's for range intersection
  regionRanges  <- GRanges(seqnames = str_remove(regionTable$chr, "^chr"), ranges = IRanges(regionTable$start, regionTable$end))
  
  if (annotation.isExistingFileOrNone(gzippedCgiFile)) {
    
    cgis <- fread(cmd = paste("zcat", gzippedCgiFile))
    cgis$chrom <- str_remove(cgis$chrom, "^chr")
    regionTable <- annotation.annotateCGIslands(regionTable, regionRanges, cgis)
    
  }
  
  if (annotation.isExistingFileOrNone(gzippedGTFFile)) {
    
    gtfRanges <- import(gzippedGTFFile)
    gtfMetaData <- values(gtfRanges)
    
    noHgncNameIndices <- is.na(gtfRanges$gene_name)
    gtfRanges$gene_name[noHgncNameIndices] <- gtfRanges$gene_id[noHgncNameIndices]
    
    biotypeColumn <- str_which(colnames(gtfMetaData), "gene_(bio)?type")
    
    gtfRanges <- gtfRanges[gtfMetaData[,biotypeColumn] %in% allowedBiotypes]
    gtfRanges <- renameSeqlevels(gtfRanges, str_remove(seqlevels(gtfRanges), "^chr"))
    
    regionTable <- annotation.annotateGenes(regionTable, regionRanges, gtfRanges)
    regionTable <- annotation.annotatePromoters(regionTable, regionRanges, gtfRanges, promoterTSSDistances)
  }
  
  if (annotation.isExistingFileOrNone(gzippedRepeatMaskerAnnotationFile)) {
    
    repeats <- fread(cmd = paste("zcat", gzippedRepeatMaskerAnnotationFile))
    repeats$genoName <- str_remove(repeats$genoName, "^chr")
    regionTable <- annotation.annotateRepeats(regionTable, regionRanges, repeats, c("LINE", "SINE", "LTR", "DNA"))
    
  }
  
  return(regionTable)
}

###### NEXT SCRIPT

wgbs.annotateDMRs <- function (dmrFile, gzippedCgiFile, gzippedGeneFile, gzippedRepeatMaskerAnnotationFile, allowedBiotypes, promoterTSSDistances, annotatedDmrFile) {
  
  dmrs <- fread(dmrFile)
  
  dmrs <- annotation.annotateRegions(dmrs, gzippedCgiFile, gzippedGeneFile, gzippedRepeatMaskerAnnotationFile, allowedBiotypes, promoterTSSDistances)
  
  write.table(dmrs, file = annotatedDmrFile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = ";")
  
  return(dmrs)
}

/mnt/c/
  
  write.csv(New_DMR, "C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/DMR.csv", row.names = FALSE)
DMR_link <- c("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/DMR.csv.gz")
cgi_link <- c("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/cpi-GRCm38.csv.gz")
anno_link <- c("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/gencode.vM30.annotation.gtf.gz")
rptmsk_link <- c("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/rptmsk-GRCm38.csv.gz")
biotypes <- c('protein_coding',
              'lncRNA',
              'miRNA')
tss_range <- c(-1000, 500)
output <- c("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/annotated_DMR.csv")

wgbs.annotateDMRs(DMR_link, cgi_link, anno_link, rptmsk_link, biotypes, tss_range, output)
DMR_table 

DMR_table <- fread("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/DMR.table.csv")
write.csv(DMR_table, "C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/DMR.table.csv", row.names = FALSE)
promoter_genes <- DMR_table %>%
  filter(promoter_overlap == TRUE)
DMR_table

promoter_genes  <- promoter_genes$gene_name 

promoter_genes <- sapply(promoter_genes, as.character)
paste(promoter_genes, collapse = ", ")

install.packages("rlang")
library("rlang")
devtools::install_github("cardiomoon/moonBook")
devtools::install_github("cardiomoon/webr")
require(ggplot2)
require(moonBook)
require(webr)
PieDonut(DMR_table, aes(pies=Annotation,donuts=ret_name), start=3*pi/95, title="DMR Annotation")

Pie
#Repeat Masker selection
#Download the Repeat Masker database file from UCSC at https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1219084465_HRDjHL1J6LoO4NKjTqqdzM3L9RS0&clade=mammal&org=Mouse&db=mm39&hgta_group=allTracks&hgta_track=rmsk&hgta_table=0&hgta_regionType=genome&position=chr12%3A56%2C741%2C761-56%2C761%2C390&hgta_outputType=primaryTable&hgta_outFileName=rptmsk-GRCm38.csv
#Read the file into R
rptmsk <- read.csv(file="rptmsk-GRCm38.csv")
rptmsk

DMR
New_DMR <- DMR
New_DMR <- select(New_DMR, -c("Seg_Pval", "Probe_ID", "Estimate", "Std. Error", "t value", "Pr(>|t|)"))
New_DMR
library(plyr)
New_DMR <- ddply(New_DMR, .(Seg_ID), summarize,
                 chr = paste(unique(Seg_Chrm)),
                 start = paste(unique(Seg_Start)),
                 end = paste(unique(Seg_End)),
                 num_cpg = length(Seg_Chrm),
                 diff = paste(unique(Seg_Est)),
                 pValue = paste(unique(Seg_Pval_adj)),
                 Tool = c("Infinium"))
New_DMR <- select(New_DMR, -Seg_ID)


rptmsk_DMR <- DMR
rptmsk_DMR <- arrange(rptmsk_DMR, rptmsk_DMR$Seg_Chrm)
rptmsk_DMR
rptmsk_ordered <- rptmsk
rptmsk_ordered
rptmsk_ordered <- arrange(rptmsk_ordered, rptmsk_ordered$genoName)
rptmsk_DMR[1,3]
rptmsk_ordered[1,7]
#### REMOVE ROWS WITH NA
rptmsk_DMR <- rptmsk_DMR %>% drop_na(Seg_Start)
rptmsk_DMR[1,8]
rptmsk_ordered


ERV_result <- data.frame()
for (i in 1:nrow(rptmsk_DMR)) {
  for (j in 1:nrow(rptmsk_ordered)) {
    if (rptmsk_DMR[i,2] == rptmsk_ordered[j,6] & rptmsk_DMR[i,3] >= rptmsk_ordered[j,7] & rptmsk_DMR[i,4] <= rptmsk_ordered[j,8]) {
      ERV_rslt_hit <- data.frame(
        "chr" = rptmsk_DMR[i,2],
        "start" = rptmsk_DMR[i,3],
        "end" = rptmsk_DMR[i,4],
        "erv" = rptmsk_ordered[j,11])
      ERV_result <- rbind.data.frame(ERV_result, ERV_rslt_hit)
      print("Hooray")
    }
  }}

ERV_result

data(Random_RNAs_500, package="RIdeogram")
head(Random_RNAs_500)
DMR
library(stringr)
probes_dmr <- data.frame(
  Type = rep("DMR", 363),
  Shape = rep("box", 363),
  Chr = str_remove_all(DMR$Seg_Chrm, "[chr]"),
  Start = DMR$Seg_Start,
  End = DMR$Seg_End,
  color = rep("#0a0a0a", 363)
)

probes_dmr <- probes_dmr %>% distinct()

probes_dmr
print(unique(probes_dmr$Chr))

# DMR
ideogram(karyotype = mouse_karyotype, label = probes_dmr, label_type = "marker")
convertSVG("chromosome.svg", device = "png")
svg2pdf("chromosome.svg")

# Motif analysis of DMRs
# This has to be done on a LINUX vitrual machine 
BiocManager::install("memes")
library(memes)
library(universalmotif)
library(dplyr)
suppressPackageStartupMessages(library(GenomicRanges))
BiocManager::install("plyranges")
suppressPackageStartupMessages(library(plyranges))
DMR
DMR_table
DMR_table <- DMR_table %>%
  mutate(
    Annotation = case_when(
      promoter_overlap == TRUE & cgi_overlap == TRUE ~ "CpG_Island_Promoter",
      promoter_overlap == TRUE & cgi_overlap == FALSE ~ "Promoter",
      gene_overlap == TRUE ~ "Gene_Body",
      promoter_overlap == FALSE & cgi_overlap == FALSE & gene_overlap == FALSE ~ "Other"
    )
  )

gr <- GRanges(
  seqnames = DMR_table$chr,
  ranges = IRanges(as.numeric(DMR_table$start), end = as.numeric(New_DMR$end)),
  strand = Rle(strand(c( "*")), c(41)),
  annotation = DMR_table$Annotation)
gr

# NOTE: beware system-specific differences. As of MEME v5, dreme will compile using the default python installation on a system (either python2.7 or python3). The random number generator changed between python2.7 and python3, so results will not be reproducible between systems using python2.7 vs 3 even if setting the same random seed.
# One way to overcome this is to manually shuffle the sequences within R. This can be done easily using universalmotif::shuffle_sequences(). Set k = 2 to preserve dinucleotide frequency (similar to dreme’s built-in shuffle), and set rng.seed to any number to create a reproducible shuffle. The output of this function can be used directly as the control sequences.
# Create random sequences to use for this example
seq <- create_sequences(rng.seed = 100)
# Shuffle sequences preserving dinucleotide frequency
shuffle <- shuffle_sequences(seq, k = 2, rng.seed = 100)

BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")
mm.genome <- BSgenome.Mmusculus.UCSC.mm39::BSgenome.Mmusculus.UCSC.mm39

#by_anno <- gr %>% 
#  anchor_center() %>% 
#  mutate(width = 100) %>% 
#  split(mcols(.)$annotation)

#seq_by_annotation <- by_anno %>% 
#  get_sequence(mm.genome)

CRE_table <- fread("/mnt/c/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/enhancer.csv")
CRE_table
#length(CRE_table$annotation)
#CRE_table1 <- CRE_table %>% filter(annotation == "Gene_Body")
#CRE_table2 <- CRE_table %>% filter(annotation == "Other")
#CRE_table <- rbind(CRE_table1, CRE_table2)
gr <- GRanges(
  seqnames = CRE_table$chrom,
  ranges = IRanges(as.numeric(CRE_table$start), end = as.numeric(CRE_table$start + CRE_table$len)),
  strand = Rle(strand(c( "*")), c(512)),
  annotation = CRE_table$annotation)
gr

seq_by_annotation <- gr %>%
  get_sequence(mm.genome)

dreme_out <- runDreme(seq_by_annotation, control = "shuffle")

class(dreme_out)
dreme_out %>% 
  to_list() %>% 
  view_motifs()

# Let's compare our cis-regulatory element motifs to known transcripton factors found on the JASPAR database.  

JASPAR <- read_meme("/mnt/c/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/JASPAR.meme")
options(meme_db = read_meme("/mnt/c/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/Lehle-UnivTexas_MouseMethylation_20220627/Lehle-UnivTexas_MouseMethylation_20220627/JASPAR.meme"))
tomtomout <- runTomTom(dreme_out, database = JASPAR)
names(tomtomout)
tomtomout$best_match_qval
tomtomout
view_motifs(tomtomout$best_match_motif)
tomtomout$tomtom[[1]] %>% head(6)
