#------------------------------------------------#
# Format example data to include in MSomicsDATA packages #
# Kelly Stratton #
# May 10, 2016 #
#------------------------------------------------#

library(MSomicsQC)

### Without a time component ###

# peptide #
mers_pep <- read.csv("O:/MSomicsQC_skeleton/example_data/WithoutTime/MCL003_peptides.csv", header=TRUE)
dim(mers_pep) # 17525    18

meta <- mers_pep[, 2]
mers_pep <- mers_pep[, c(1, 7:18)]
dim(mers_pep) # 17525    13

sampleinfo <- read.csv("O:/MSomicsQC_skeleton/example_data/WithoutTime/Sample_Info_Peptides.csv", header=TRUE)
dim(sampleinfo) # 12 3
names(sampleinfo)
sampleinfo_pep <- data.frame(SampleID = sampleinfo$SampleID, Condition = sampleinfo$Condition)

emeta <- read.csv("O:/MSomicsQC_skeleton/example_data/WithoutTime/MCL003_T_Mass_Tag_to_Protein_Map.csv", header=TRUE)
# need to map emeta to the peptide data #
emeta <- emeta[, 1:4]
inds <- match(mers_pep$Mass_Tag_ID, emeta$Mass_Tag_ID)
emeta2 <- emeta[inds, ]
dim(emeta2)
all(emeta2$Mass_Tag_ID == mers_pep$Mass_Tag_ID) #TRUE

# remove contaminant proteins #
inds <- unique(c(grep(pattern = "Contaminant_", emeta2$Protein), grep(pattern = "HORSE", emeta2$Protein), grep(pattern = "RABIT", emeta2$Protein), grep(pattern = "BOV", emeta2$Protein)))
length(inds)
emeta3 <- emeta2[-inds, ]
dim(emeta3)
mers_pep <- mers_pep[-inds, ]
dim(mers_pep)

pep_pepData <- as.pepData(e_data = mers_pep, f_data = sampleinfo_pep, e_meta = emeta3, edata_cname = "Mass_Tag_ID", samp_cname = "SampleID", emeta_cname = "Protein")
attr(pep_pepData, "data_info")$data_scale <- "abundance"

pep_edata <- mers_pep
pep_emeta <- emeta3
pep_fdata <- sampleinfo_pep

summary(pep_pepData)
# $num_samps
# [1] 12
#
# $num_edata
# [1] 17407
#
# $num_emeta
# [1] 2774
#
# $num_miss_obs
# [1] 82833
#
# $frac_missing
# [1] 0.3660464


saveRDS(pep_pepData, "O:/MSomicsQC_skeleton/example_data/WithoutTime/pep.RDS")
devtools::use_data(pep_pepData, pkg = "MSomicsDATA", overwrite = TRUE)
devtools::use_data(pep_edata, pkg = "MSomicsDATA", overwrite = TRUE)
devtools::use_data(pep_fdata, pkg = "MSomicsDATA", overwrite = TRUE)
devtools::use_data(pep_emeta, pkg = "MSomicsDATA", overwrite = TRUE)


# prot #
mers_prot <- read.csv("O:/MSomicsQC_skeleton/example_data/WithoutTime/MCL003_proteins_bpquant.csv", header=TRUE)
dim(mers_prot) # 2731   12
names(mers_prot)[1] <- "Reference"

sampleinfo <- read.csv("O:/MSomicsQC_skeleton/example_data/WithoutTime/Sample_Info_Proteins.csv", header=TRUE)
dim(sampleinfo) # 11 3
sampleinfo_pro <- data.frame(SampleID = sampleinfo$SampleID, Condition = sampleinfo$Condition)

pro_proData <- as.proData(e_data = mers_prot, f_data = sampleinfo_pro, edata_cname = "Reference", samp_cname = "SampleID")
attr(pro_proData, "data_info")$data_norm <- TRUE
attr(pro_proData, "data_info")$norm_method <- "median_center"
attr(pro_proData, "data_info")$data_scale <- "log2"
pro_edata <- pro_proData$e_data
pro_fdata <- pro_proData$f_data

summary(pro_proData)
# $num_samps
# [1] 11
#
# $num_edata
# [1] 2731
#
# $num_emeta
# NULL
#
# $num_miss_obs
# [1] 5242
#
# $frac_missing
# [1] 0.1599536

saveRDS(pro_proData, "O:/MSomicsQC_skeleton/example_data/WithoutTime/prot.RDS")
devtools::use_data(pro_proData, pkg = "MSomicsDATA", overwrite = TRUE)
devtools::use_data(pro_edata, pkg = "MSomicsDATA", overwrite = TRUE)
devtools::use_data(pro_fdata, pkg = "MSomicsDATA", overwrite = TRUE)



# metab #
mers_metab <- read.csv("O:/MSomicsQC_skeleton/example_data/WithoutTime/MCL003_Metabolites.csv", header=TRUE)
dim(mers_metab) # 80 13

sampleinfo <- read.csv("O:/MSomicsQC_skeleton/example_data/WithoutTime/Sample_Info_Metabolomics.csv", header=TRUE)
dim(sampleinfo) # 12 3
sampleinfo_metab <- data.frame(SampleID = sampleinfo$SampleID, Condition = sampleinfo$Condition)

metab_metabData <- as.metabData(e_data = mers_metab, f_data = sampleinfo_metab, edata_cname = "Metabolite", samp_cname = "SampleID")
attr(metab_metabData, "data_info")$data_scale <- "abundance"

# replace 0s with NAs
metab_metabData <- edata_replace(mintR_data = metab_metabData, x = 0, y = NA)
#  "148 instances of 0 have been replaced with NA"

metab_edata <- metab_metabData$e_data
metab_fdata <- metab_metabData$f_data

summary(metab_metabData)
# $num_samps
# [1] 12
#
# $num_edata
# [1] 80
#
# $num_emeta
# NULL
#
# $num_miss_obs
# [1] 148
#
# $frac_missing
# [1] 0.1423077

saveRDS(metab_metabData, "O:/MSomicsQC_skeleton/example_data/WithoutTime/metab.RDS")
devtools::use_data(metab_metabData, pkg = "MSomicsDATA", overwrite = TRUE)
devtools::use_data(metab_edata, pkg = "MSomicsDATA", overwrite = TRUE)
devtools::use_data(metab_fdata, pkg = "MSomicsDATA", overwrite = TRUE)


# lipid - pos #
mers_lipid_pos <- read.csv("O:/MSomicsQC_skeleton/example_data/WithoutTime/MCL003_lipids_pos.csv", header=TRUE)
dim(mers_lipid_pos) # 148  12

length(unique(mers_lipid_pos$LipidCommonName)) # 146

mers_lipid_pos <- mers_lipid_pos[-c(49, 84), ]
dim(mers_lipid_pos) # 148  12

sampleinfo <- read.csv("O:/MSomicsQC_skeleton/example_data/WithoutTime/Sample_Info_Lipids_Pos.csv", header=TRUE)
dim(sampleinfo) # 11 3

sampleinfo_lipid_pos <- data.frame(Sample_Name = sampleinfo$SampleID, Condition = sampleinfo$Condition)

lipid_lipidData <- as.lipidData(e_data = mers_lipid_pos, f_data = sampleinfo_lipid_pos, edata_cname = "LipidCommonName", samp_cname = "Sample_Name")
attr(lipid_lipidData, "data_info")$data_scale <- "abundance"
lipid_edata <- lipid_lipidData$e_data
lipid_fdata <- lipid_lipidData$f_data

summary(lipid_lipidData)
# $num_samps
# [1] 11
#
# $num_edata
# [1] 146
#
# $num_emeta
# NULL
#
# $num_miss_obs
# [1] 884
#
# $frac_missing
# [1] 0.5045662

saveRDS(lipid_lipidData, "O:/MSomicsQC_skeleton/example_data/WithoutTime/lipid_pos.RDS")
devtools::use_data(lipid_lipidData, pkg = "MSomicsDATA", overwrite = TRUE)
devtools::use_data(lipid_edata, pkg = "MSomicsDATA", overwrite = TRUE)
devtools::use_data(lipid_fdata, pkg = "MSomicsDATA", overwrite = TRUE)


# lipid - neg #
mers_lipid_neg <- read.csv("O:/MSomicsQC_skeleton/example_data/WithoutTime/MCL003_lipids_neg.csv", header=TRUE)
dim(mers_lipid_neg) # 144  12
length(unique(mers_lipid_neg$LipidCommonName)) # 143

sampleinfo <- read.csv("O:/MSomicsQC_skeleton/example_data/WithoutTime/Sample_Info_Lipids_Neg.csv", header=TRUE)
dim(sampleinfo) # 11 3

sampleinfo_lipid <- data.frame(SampleID = sampleinfo$SampleID, Condition = sampleinfo$Condition)

lipidN_lipidData <- as.lipidData(e_data = mers_lipid_neg, f_data = sampleinfo_lipid, edata_cname = "LipidCommonName", samp_cname = "SampleID", data_scale = "abundance")
lipid_edata <- lipidN_lipidData$e_data
lipid_fdata <- lipidN_lipidData$f_data

saveRDS(lipidN_lipidData, "O:/MSomicsQC_skeleton/example_data/WithoutTime/lipid_neg.RDS")
# devtools::use_data(lipidN_lipidData, pkg = "MSomicsDATA", overwrite = TRUE)
# devtools::use_data(lipidN_edata, pkg = "MSomicsDATA", overwrite = TRUE)
# devtools::use_data(lipidN_fdata, pkg = "MSomicsDATA", overwrite = TRUE)

#------------------------------------------------#

### With a time component ###

# metab #
icl_metab <- read.csv("O:/MSomicsQC_skeleton/example_data/WithTime/OMICS_ICL102_Metabolomics.csv", header=TRUE)
dim(icl_metab) # 80 61

sampleinfo_metab <- read.csv("O:/MSomicsQC_skeleton/example_data/WithTime/Sample_Info_Metabolomics.csv", header=TRUE)
dim(sampleinfo_metab) # 60 3

metab <- as.metabData(e_data = icl_metab, f_data = sampleinfo_metab, edata_cname = "Metabolite", samp_cname = "MetabSampleName")

# replace 0s with NAs
metab_time <- edata_replace(mintR_data = metab, x = 0, y = NA)
# "168 instances of 0 have been replaced with NA"
saveRDS(metab_time, "O:/MSomicsQC_skeleton/example_data/WithTime/metab.RDS")
devtools::use_data(metab_time, pkg = "MSomicsDATA")


# lipid - pos #
icl_lipid_pos <- read.csv("O:/MSomicsQC_skeleton/example_data/WithTime/OMICS_ICL102_LipidomicPos.csv", header=TRUE)
dim(icl_lipid_pos) # 245  61

sampleinfo_lipid_pos <- read.csv("O:/MSomicsQC_skeleton/example_data/WithTime/Sample_Info_Lipids_Pos.csv", header=TRUE)
dim(sampleinfo_lipid_pos) # 60 3

lipid_pos_time <- as.lipidData(e_data = icl_lipid_pos, f_data = sampleinfo_lipid_pos, edata_cname = "LipidCommonName", samp_cname = "LipidPosSampleName")

saveRDS(lipid_pos_time, "O:/MSomicsQC_skeleton/example_data/WithTime/lipid_pos.RDS")
devtools::use_data(lipid_pos_time, pkg = "MSomicsDATA")


# lipid - neg #
icl_lipid_neg <- read.csv("O:/MSomicsQC_skeleton/example_data/WithTime/OMICS_ICL102_LipidomicNeg.csv", header=TRUE)
dim(icl_lipid_neg) # 251  61

sampleinfo_lipid_neg <- read.csv("O:/MSomicsQC_skeleton/example_data/WithTime/Sample_Info_Lipids_Neg.csv", header=TRUE)
dim(sampleinfo_lipid_neg) # 60 3

lipid_neg_time <- as.lipidData(e_data = icl_lipid_neg, f_data = sampleinfo_lipid_neg, edata_cname = "LipidCommonName", samp_cname = "LipidNegSampleName")

saveRDS(lipid_neg_time, "O:/MSomicsQC_skeleton/example_data/WithTime/lipid_neg.RDS")
#devtools::use_data(lipid_neg_time, pkg = "MSomicsDATA")
