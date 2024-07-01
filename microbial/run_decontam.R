.libPaths(new = ".Rlibs")
# Add Cran mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Installation
# Bioconductor version 3.19 (BiocManager 1.30.23), R 4.4.0 (2024-04-24)

if (!requireNamespace("BiocManager", quietly = TRUE, force=FALSE))
install.packages("BiocManager")
BiocManager::install("decontam")
install.packages("ggplot2")
#install.packages("phyloseq")

# Call libraries
#library(ggplot2); packageVersion("ggplot2") # ‘3.5.1’
#library(phyloseq); packageVersion("phyloseq") # ‘1.48.0’
library(decontam); packageVersion("decontam") # ‘1.24.0’


## Create output mother_path
mother_path <-'/Users/lchen/00_projects/FBI/output/02_tables/03_intermediate/'
input_path <- '/Users/lchen/00_projects/FBI/output/02_tables/03_intermediate/'

# Perform decontamination isContaminant detection with non-batched information
# Read in taxonomy cfDNA relative fraction kraken counts.
data <- read.csv(paste(input_path, "contam_all_taxid_NO_batch.csv", sep=''))
df <-as.matrix(data)
# Read SRSLY input volume of each sample (contam type B)
concB <- read.csv(paste(input_path, "contamB_srsly_input_volume_NO_batch.csv", sep=''))
concB <- concB[,1]
contamdf.freq <- isContaminant(df, method='frequency', conc = concB)
write.csv(contamdf.freq, file=paste(mother_path, 'decontam_output_concB_no_batch_X_taxid.csv', sep = ''))

# Read inverse isolated input DNA volume * inverse srsly input volume
concA <- read.csv(paste(input_path, "contamA_input_volume_NO_batch.csv", sep=''))
concA <- concA[,1]
contamdf.freq2 <- isContaminant(df, method='frequency', conc = concA)
write.csv(contamdf.freq2, file=paste(mother_path, 'decontam_output_concA_no_batch_X_taxid.csv', sep = ''))

# Read inverse isolation input DNA volume -- this term is independent from SRSLY prep
contamA_alt_no_batch <- read.csv(paste(input_path, "ALT_contamA_input_volume_NO_batch.csv", sep = ''))
contamA_alt_no_batch <- contamA_alt_no_batch[,1]
contamdf.freq3 <- isContaminant(df, method='frequency', conc = contamA_alt_no_batch)
write.csv(contamdf.freq3, file=paste(mother_path, 'decontam_output_ALT_conA_no_batch_X_taxid.csv', sep=''))


# Here we execute the same 3 analysis but with batched-specific tests. The test report 1 p-value that is of the smallest p-value of all batches (minimum).
# concA
# Make data into right format
data_batched_concA <- read.csv(paste(input_path, "contamA_iso_all_taxid_BATCHED.csv",sep=''))
data_batched_concA <-as.matrix(data_batched_concA)
# Make data into right format
concA_batched <- read.csv(paste(input_path,'contamA_input_volume_BATCHED.csv',sep=''))
concA_batched <- concA_batched[,1]
# Make data into right format
batch_index_concA <- read.csv(paste(input_path, "contamA_iso_BATCHES.csv", sep = ''))
batch_index_concA <- factor(batch_index_concA[,1])
# Perform decontam, include batch index
contamdf.freq4 <- isContaminant(data_batched_concA, conc = concA_batched, batch=batch_index_concA)
write.csv(contamdf.freq4, file=paste(mother_path, 'decontam_output_contamA_BATCHED_X_taxid.csv', sep = ''))

# concB
# Make data into right format
data_batched_concB <- read.csv(paste(input_path, "contamB_prep_all_taxid_BATCHED.csv", sep = ''))
data_batched_concB <- as.matrix(data_batched_concB)

# Make data into right format
concB_batched <- read.csv(paste(input_path,'contamB_srsly_input_volume_BATCHED.csv',sep=''))
concB_batched <- concB_batched[,1]

# Make data into right format
batch_index_concB <- read.csv(paste(input_path,'contamB_prep_BATCHES.csv',sep=''))
batch_index_concB <- factor(batch_index_concB[,1])

# Perform decontam
contamdf.freq5 <- isContaminant(data_batched_concB, conc = concB_batched, batch=batch_index_concB)
write.csv(contamdf.freq5, file=paste(mother_path, 'decontam_output_contamB_BATCHED_X_taxid.csv', sep = ''))

## ALT contamA
ALT_concA_batched <- read.csv(paste(input_path,'ALT_contamA_input_volume_BATCHED.csv',sep=''))
ALT_concA_batched <- ALT_concA_batched[,1]
contamdf.freq6 <- isContaminant(data_batched_concA, conc = ALT_concA_batched, batch=batch_index_concA)
write.csv(contamdf.freq6, file=paste(mother_path, 'decontam_output_ALT_contamA_BATCHED_X_taxid.csv', sep = ''))
