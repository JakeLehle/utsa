getwd()
setwd("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array")
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


load("Sesame_Mouse_Bartolomei_Data.RData")
## -----------------------------------------------------------------------------
tools::R_user_dir("ExperimentHub", which="cache")
# This shows you location of the cached data on your computer. 

######## SAMPLE DATA FROM MOUSE

sesameDataCacheAll()

# First create an object with the path to the IDAT file
idat_dir = c("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array")
idat_dir
IDATprefixes <- searchIDATprefixes(idat_dir)
IDATprefixes
basename(IDATprefixes)

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
# For instance all samples look like BUB_BnJ except sample 22?
# Lets plot first sample 1 then sample 22 and determine what's going on.

sdf_1 <- readIDATpair("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/GSM5778405_Female-blastocyst-1358-6", platform = "MM285")
sdf_22 <- readIDATpair("C:/Users/jakel/Documents/UTSA/Lab/IVERR/Infinium_Array/GSM5778431_Male-blastocyst-1620-2", platform = "MM285")

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

#Looks like BUB_BnJ only represents 0.2 of sample 22 the other 0.8 is C57BL_6
# Fascinating stuff to keep in mind for your own studies to try and maintain high similarities in the backgrounds of your mouse samples or at leats bring it up when discussing results.
# You can use this structure to also determine the ethnicity of Human samples which should be discussed in projects especially in the context of studying methylation patterns in underrepresented populations. 

# You can even take this infomation futher and predict the average age of the samples in months from the beta values
# Haha it's a negative number because these are samples are at day 2.5 to 4.5

predictMouseAgeInMonth(betas[,1])

# Sample preprocessing and quality control

sdf_preped = openSesame(IDATprefixes, BPPARAM = BiocParallel::SnowParam(8), prep="TQCDPB", func=NULL)
sdf_preped


# Using our sdf's from sample 1 and 22 lets look at the probe signal clarity and subtract background. 
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
sesameQC_plotRedGrnQQ(sdf_1.InfICorrected, main="Before")
sesameQC_plotRedGrnQQ(dyeBiasNL(sdf_1.InfICorrected), main="After")  # nonlinear correction

sdf_1.InfICorrected <- dyeBiasNL(sdf_1.InfICorrected)

###### Background Subtraction
# Let's now look at background subtraction using oob. Equivalent to prep="B" 

par(mfrow=c(2,1), mar=c(3,3,2,1))
sesameQC_plotBetaByDesign(sdf_1.InfICorrected, main="Before", xlab="\beta")
sesameQC_plotBetaByDesign(noob(sdf_1.InfICorrected), main="After", xlab="\beta")

# Wow that data looks rough even cleaned up there is a lot of ambiguity at many of these probe sites leading to intermediate probe values.
# However, this isn't my data, so not my problem. You now have an example of how to go through pre-processing and displayt these changes visually.
# Let's look at some other quality control metrics that can give us an idea of how good our data is (or in this case how bad haha).

#QC of detection metric
qcs = openSesame(IDATprefixes, BPPARAM = BiocParallel::SnowParam(8), prep="", func=sesameQC_calcStats, funs="detection")
qcs_prep = openSesame(IDATprefixes, BPPARAM = BiocParallel::SnowParam(8), prep="TQCDPB", func=sesameQC_calcStats, funs="detection")

qcs
qcs_prep

sesameQC_plotBar(qcs_prep)

# The detection metrics were also disappointing at or around 75%

########### DMR Analysis
# Let's start our DMR analysis
# To do this we need to create a summarized experiment object
# The object itself will be built on our beta values calculated earlier
# We will create a coldata dataframe object that will divide our object by IDAT file name, stage/age of the sample, and sex 
# Fianlly I throw in a little hack to reassign the name of the assay so it will indicate it is recording beta values.
betas_prep = openSesame(IDATprefixes, BPPARAM = BiocParallel::SnowParam(8), prep="TQCDPB")
IDAT <- basename(IDATprefixes)
stage <- c("blastocyst", "blastocyst", "blastocyst", "blastocyst", "blastocyst", "blastocyst", 
           "morula", "morula", "morula", "morula", "morula", "morula", 
           "natrual", "natrual", "natrual", "natrual", "natrual",
           "blastocyst", "blastocyst", "blastocyst", "blastocyst", "blastocyst", "blastocyst",
           "morula", "morula", "morula", "morula", "morula", "morula",
           "natrual", "natrual", "natrual", "natrual", "natrual")
sex <- c("Female", "Female", "Female", "Female", "Female", "Female", "Female", "Female", "Female", "Female", "Female", "Female", "Female", "Female", "Female", "Female", "Female", 
         "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Male", "Male")
df <- data.frame(IDAT, stage, sex)
df
se <- SummarizedExperiment(assays = betas_prep, colData = df)
names(assays(se)) <- c("betas") #hack hack hack
se
cd = as.data.frame(colData(se)); rownames(cd) = NULL
cd

# To check DMR you have to exclude NA's use the checkLevels() function

se_ok = (checkLevels(assay(se), colData(se)$sex) &
           checkLevels(assay(se), colData(se)$stage))
sum(se_ok)                      # the number of CpGs that passes
se = se[se_ok,]

# Here they set the control cell as the colon and the control sex as female
colData(se)$stage <- relevel(factor(colData(se)$stage), "natrual")
colData(se)$sex <- relevel(factor(colData(se)$sex), "Female")

#Now we model DNA methylaton variation treating stage and sex as covariates. 
#The DML will fit the DNA methylation reading to a linear model and perform the corresponding slope test and goodness-of-fit test (F-test holding out each contrast variable)
#All of these results are returned in an object of class DMLSummary
# This took ~4 hours on 6-15-22
smry = DML(se, ~stage + sex, mc.cores = BiocParallel::SnowParam(8))
smry
save.image("Sesame_Mouse_Bartolomei_Data.RData")

#### TEST INERPRETATION
test_result = summaryExtractTest(smry)
# Rows are CpG/Loci and columns contain the slopes and p-values for each variable
colnames(test_result) # the column names, show four groups of statistics
#EST The slope estimate (aka the β coefficient, not to be confused with the DNA methylation β-value though) for continuous variable. DNA methylation difference of the current level with respect to the reference level for nominal contrast variables. Each suffix is concatenated from the contrast variable name (e.g., tissue, sex) and the level name if the contrast variable is discrete (e.g, Cecum, Esophagus, Fat). For example, Est_tissueFat should be interpreted as the estimated methylation level difference of Fat compared to the reference tissue (which is Colon, as set above). There is a special column named Est_`(Intercept)`. It corresponds to the base-level methylation of the reference (in this case a Female Colon sample)
#PVAL The unadjusted p-values of t-testing the slope. This represents the statistical significance of the methylation difference. For example, Pval_stageblasticyst tests whether the blastocyst is significantly different from the normal embryo (the reference level) in DNA methylation. The Pval_`(Intercept)` tests whether the reference level is significantly different from zero.
#FPVAL The unadjusted p-value of the F-test contrasting the full model against a reduced model with the labeled contrast variable held out. Note that “Pval_” and “FPval_” are equivalent when the contrast variable is a 2-level factor, i.e., in the case of a pairwise comparison.
#EFF The effect size of each normial contrast variable. This is equivalent to the maximum slope subtracted by the minimum level including the reference level (0).
head(test_result)

### For your project PVAL will be the most important metric to see here.
###### GOODNESS Of FIT
# Another way to say this is.
# Is the CpG methylation tissue-specific? Rather than. Is the CpG more methylated in cultured blastocysts compared to normal embryos?
# Here we used 0.1 as the size threshold. This means a difference less than 10% is not considered biologically relevant and excluded. 
# We can go further and do a side by comparison of probes that are sex specific and ones that are cell type specific. 

test_result %>%
  mutate(sex_specific =
           ifelse(FPval_sex < 0.05 & Eff_sex > 0.1, TRUE, FALSE)) %>%
  mutate(stage_specific =
           ifelse(FPval_stage < 0.05 & Eff_stage > 0.1, TRUE, FALSE)) %>%
  select(sex_specific, stage_specific) %>% table

#Use this result to make a Venn diagram

test_result_probes_df <- test_result %>%
  mutate(sex_specific =
           ifelse(FPval_sex < 0.05 & Eff_sex > 0.1, TRUE, FALSE)) %>%
  mutate(stage_specific =
           ifelse(FPval_stage < 0.05 & Eff_stage > 0.1, TRUE, FALSE)) %>%
  select(sex_specific, stage_specific) %>% DataFrame()


test_result_probes_df$stage_specific
sex_specific_rn <- na.omit(rownames(test_result_probes_df)[test_result_probes_df[["sex_specific"]] == TRUE])
stage_specific_rn <- na.omit(rownames(test_result_probes_df)[test_result_probes_df[["stage_specific"]] == TRUE])

specific_list <- list(
  "Sex Specific" = sex_specific_rn,
  "Stage Specific" = stage_specific_rn
) 

ggvenn(specific_list,
       fill_color = c("#0073C2FF", "#868686FF"))


#Use this list to find out if any class or probes are enriched in the population
stage_specific <- test_result %>% dplyr::filter(FPval_stage < 0.05, Eff_stage > 0.1) %>% select(FPval_stage, Eff_stage)
stage_specific
# Let's see if there an enrichment of probes in the stage_specific differential methylated loci regions within one of the prob groups
query <- rownames(stage_specific)
query
results_stage_specific <- testEnrichment(query)
KYCG_plotEnrichAll(results_stage_specific)
# How cools is that this shows that there is a enrichment of methylation diffrences at Imprinted DMRs which has ben shown previously as well as a number of other sites like H19, Gene Bodies, Enhancers, Pesudogene TSS's, etc
results_stage_specific
results_stage_specific_filtered <- results_stage_specific %>% dplyr::filter(overlap>10)
KYCG_plotEnrichAll(results_stage_specific_filtered)
# Let's rank these by significance
KYCG_plotWaterfall(results_stage_specific_filtered)

# To dig more into each group we need to see what groups are available which can be done using the list database groups function.
KYCG_listDBGroups("MM285")
# Now that we can look at individual groups we can start pulling apart some of the data we are seeing from our enrich all data
# First let's do a dot plot of the design group to see which of the design groups are most differentially methylated
results_stage_specific_designGroup <- testEnrichment(query, "designGroup", platform = "MM285")
KYCG_plotDot(results_stage_specific_designGroup, n_max=20)
# Next to find out more about if a region is hyper or hypo methylated lets to a volcano plot this time using the chromHMM group.
results_2tailed <- testEnrichment(query, "chromHMM", alternative = "two.sided", platform = "MM285")
KYCG_plotVolcano(results_2tailed)

# You can see from this that there is quite a bit of data that you can already get and this might help you to determine what group to look into more. 
# However let's go back and look at all of our probes in the DML group and find out the ones that have the largest differences.

# Now we can plot this difference with a volcano plot to show the differences in probes for men compared to females
library(ggplot2)
ggplot(test_result) + geom_point(aes(Est_sexMale, -log10(Pval_sexMale)))

# Here they are showing the same thing but with blastocysts compared to normal embryos
ggplot(test_result)+ geom_point(aes(Est_stageblastocyst, -log10(Pval_stageblastocyst)))
# Let's make it pretty but focus on plotting the false discovery rate or FPval for stage specific changes as it relates to changes in beta values.


test_result_pretty <- test_result
test_result_pretty$Methylation

test_result_pretty <- test_result_pretty %>%
  mutate(
    Significance = case_when(
      FPval_stage <= 0.05 & FPval_stage > 0.01 ~ "FDR 0.05",
      FPval_stage <= 0.01 & FPval_stage > 0.001 ~ "FDR 0.01",
      FPval_stage <= 0.001 ~ "FDR 0.001",
      TRUE ~ "Unchanged"),
    Methylation = case_when(
      Est_stageblastocyst > 0 ~ "Hyper",
      Est_stageblastocyst < 0 ~ "Hypo",
      TRUE ~ "Unchanged")
  )

# Side note. Quick way to find out how many of the probes are either hyper or hypo methylated in blastocysts samples compared to controls.
test_result_pretty %>%
  count(Methylation)
# Now with significance added for stage specific changes 
test_result_pretty %>%
  count(Methylation, Significance)

# Let's find out the higest and lowest changes in methylation in our groups

top_10 <- 10
top_probes <- bind_rows(
  test_result_pretty %>%
    filter(Methylation == 'Hyper') %>%
    arrange(FPval_stage, desc(abs(Pval_stageblastocyst))) %>%
    head(top_10),
  test_result_pretty %>%
    filter(Methylation == 'Hypo') %>%
    arrange(FPval_stage, desc(abs(Pval_stageblastocyst))) %>%
    head(top_10)
)

test_result_pretty$Est_Probe_ID

ggplot(test_result_pretty, aes(Est_stageblastocyst, -log10(Pval_stageblastocyst))) +
  geom_point(aes(color = Significance), size = 2/5) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) +
  geom_label_repel(data = top_probes,
                   mapping = aes(Est_stageblastocyst, -log(Pval_stageblastocyst,10), label = Est_Probe_ID),
                   size = 2, max.overlaps = 20)

#Okay now that you have the 10 most significant hyper and hypo methylated probes you can visualize them on the genome using this function

visualizeProbes(c("cg30118034_TC11"), betas_prep, platform='MM285')

# Oh wow that probe is located in the promoter of Bptf the bromodomain protein that has been linked to Alzhymers 
# Let's call that region out by gene now.

visualizeGene('Bptf', betas_prep, platform='MM285')

# Here we can model these functions as continuous. This shows CpGs that are positively associated with length of culture conditions.
smry2 = DML(se, ~ stage + sex, mc.cores = BiocParallel::SnowParam(8))
test_result2 = summaryExtractTest(smry2) %>% arrange(Est_stageblastocyst)

test_result2 %>% dplyr::select(Est_stageblastocyst, Pval_stageblastocyst) %>% tail # positive assoc.
df = data.frame(Stage = colData(se)$stage,
                BetaValue = assay(se)[rownames(test_result2)[nrow(test_result2)],])
ggplot(df, aes(Stage, BetaValue)) + geom_smooth(method="lm") + geom_point()

###### DMR
# Okay so the way that this works is if there are significant differences in a group of probes that are closely located together they can be merged into a region that they call a DMR.
# You will then have to annotate these regions seperately to figure out if they are in promoter regions or gene regions use the wg-blimp script for this.
dmContrasts(smry)                       # pick a contrast from below
merged = DMR(se, smry, "stageblastocyst")       # merge CpGs to regions / segments
DMR <- merged %>% dplyr::filter(Seg_Pval_adj < 0.01) %>% arrange(Seg_Pval_adj)
DMR

# Here is where you can call out specific regions and get the probe names
visualizeRegion(
  'chr11',22971841,22972953, betas_prep, platform='MM285',
  show.probeNames = TRUE)
#Turns out this is a DMR for a zinc finger binding protein that is important for mRNA processing.

######## GO Analysis 

#Start with a probe list. For this example let's use the probes that are specific for blastocysts.

blastocysts_genes <- sesameData_getGenesByProbes(stage_specific_rn)
blastocysts_genes

#  ## use gene name
gostres <- gost(blastocysts_genes$gene_name, organism = "mmusculus")
gostres$result[order(gostres$result$p_value),]
gostplot(gostres)
#  
## use Ensembl gene ID, note we need to remove the version suffix
gene_ids <- sapply(strsplit(names(blastocysts_genes),"\\."), function(x) x[1])
gostres <- gost(gene_ids, organism = "mmusculus")
gostres$result[order(gostres$result$p_value),]
gostplot(gostres)

# Both will give you the same thing.
