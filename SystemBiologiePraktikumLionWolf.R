# Setup ########################################################################

# Install the required package from CRAN 
cran_packages <- c("readr", "dplyr", "bedr", "data.table", "tictoc", "progress")

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) { # !require still installs the packages
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# install from Bioconductor
bc_packages <- c("rtracklayer", "BiocGenerics", "S4Vectors", "GenomicRanges")
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
for (pkg in bc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE)
  library(pkg, character.only = TRUE)
}

rm(cran_packages, bc_packages, pkg)

data_path = "/home/lion/Nextcloud/Computational Biology/2. FS/Systembiologie/Praktikum/data"


# Data Aquisition ##############################################################
# Preselection on ENCODE database yields a list of 8746 experiments for cell lines K562 and GM12878

src <- read_tsv(
  "/home/lion/SysbioDataStorage/data-selection/selected-files-metadata.tsv", 
  col_names = T, 
  quote = "\"") 

experiments <- src %>%
  filter(`File assembly` %in% c("hg19")) %>%
  filter(`Output type` %in% c("optimal IDR thresholded peaks")) %>%
  filter(!is.na(`Biological replicate(s)`)) %>%
  filter(is.na(`Audit ERROR`))


# filter for cell line
# K562
experiments_K562 <- experiments %>% 
  filter(`Biosample term name` %in% c("K562"))


# GM12878
experiments_GM12878 <- experiments %>% 
  filter(`Biosample term name` %in% c("GM12878"))


# Data Transformation ##########################################################

# 1) defining function to cut bed regions to size around peak ##################

queensize <- function(folder_path, window_size, overwrite = FALSE) {
  tic("time for trimming: ")
  
  if (!overwrite) {
    # List to store the modified data frames
    modified_data_frames <- list()
  }
  # Check if folder_path is a directory or a single file
  if (file.info(folder_path)$isdir) {
    # If folder_path is a directory, get all the .bed files in it
    bed_files <- list.files(path = folder_path, pattern = "\\.bed$", full.names = TRUE)
  } else {
    # If folder_path is a single file, process it directly
    bed_files <- folder_path
  }
  
  # Create a progress bar
  pb <- progress_bar$new(total = length(bed_files), format = "[:bar] :percent :eta")
  
  # Loop through each BED file
  for (bed_file in bed_files) {
    # Read the BED file into a data frame
    bed_data <- read.table(bed_file, header = FALSE, sep = "\t")
    
    # Perform the operations on each row
    bed_data[, 2] <- floor((bed_data[, 2] + bed_data[, 10]) - (window_size/2))
    bed_data[, 3] <- bed_data[, 2] + window_size
    bed_data[, 10] <- round(window_size/2)
    
    if (overwrite) {
      # Overwrite the original BED file with the updated data
      write.table(bed_data, file = bed_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    } else {
      # Save the modified data frame to the list
      modified_data_frames[[basename(bed_file)]] <- bed_data
    }
    
    # Increment the progress bar
    pb$tick()
  }
  toc()
  # Return a message indicating the operation is complete
  return("Trimming complete. Original BED files have been updated.")
}


kingsize <- function(dataframe, windowsize) {
  tic("time for trimming: ")
  # Calculate the midpoint of the original start and end coordinates
  midpoint <- floor(dataframe$start + (dataframe$stop - dataframe$start) / 2)
  
  # Calculate the new start and end coordinates based on the midpoint and windowsize
  dataframe$start <- floor(midpoint - (windowsize / 2))
  dataframe$stop <- floor(midpoint + (windowsize / 2))
  
  toc()
  # Return the updated dataframe
  return(dataframe)
}


# 2) read in TFBS data #########################################################
TFBS_union <- fread("./data/TRN_TFBS_union.csv", sep = ",", header = TRUE)

# change the tf column from "V$name_*" to "name"
TFBS_union <- TFBS_union %>% 
  mutate(tf = sub("^V\\$", "", tf)) %>%
  mutate(tf = sub("^(.*?)_.*", "\\1", tf))
  

# 3) Retrieving common TF targets present in TFBS dataset and ChIPseq datasets #####################
ChIP_targets <- sort(unique(sub("-human", "", experiments$`Experiment target`)))
TFBS_targets <- sort(unique(TFBS_union[["tf"]]))
common_TF <- intersect(ChIP_targets, TFBS_targets)

# There are 76 common TF targets present in both datasets

# 4a) Filtering TFBS data for common TF targets ####################################
filtered_TFBS_union <- subset(TFBS_union, tf %in% common_TF)


# 4b) Filtering ChIPseq data for common TF targets #################################
# K562
filtered_experiments_K562 <- experiments_K562 %>% 
  mutate(`Experiment target` = sub("-human", "", `Experiment target`)) %>% 
  filter(`Experiment target` %in% common_TF)

write_tsv(filtered_experiments_K562, file.path(data_path, "K562-sharedTF-experiments-metadata.tsv"))
write_lines(filtered_experiments_K562$`File download URL`, file.path(data_path, "K562-sharedTF-experiments.txt"))

#sort(table(filtered_experiments_K562$`Experiment target`))
# all TF-targets have between 1 and 5 corresponding experiments


# GM12878
filtered_experiments_GM12878 <- experiments_GM12878 %>% 
  mutate(`Experiment target` = sub("-human", "", `Experiment target`)) %>% 
  filter(`Experiment target` %in% common_TF)

write_tsv(filtered_experiments_GM12878, file.path(data_path, "GM12878-sharedTF-experiments-metadata.tsv"))
write_lines(filtered_experiments_GM12878$`File download URL`, file.path(data_path, "GM12878-sharedTF-experiments.txt"))


# 5a) download appropriate files ################################################

# cd into K562 or GM12878 directory and download with: 

# xargs -L 1 curl -O -J -L < ../../K562-sharedTF-experiments.txt
# xargs -L 1 curl -O -J -L < ../../GM12878-sharedTF-experiments.txt


# 5b) execute bash script to sort and decompress all files #######################




# 6) create Results table #######################################################

# source bed2BWscores function
source(file.path(data_path, "bed2BWscore.R"))


# 6a) Creating the TFBS profiles for each common transcription factor ############

TFBS_profile_df <- data.frame(matrix(0, nrow = 501, ncol = 0))

i=1
# BaTcH pRoCeSsInG
target <- common_TF[i]
print(c("target is:", target))
TFBS_regions <- filtered_TFBS_union[filtered_TFBS_union[["tf"]] == target]

TFBS_regions <- as.data.frame(kingsize(TFBS_regions, windowsize = 500))
write.table(TFBS_regions, file = file.path(data_path, "TFBS_regions.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
rm(TFBS_regions)

TFBS_scores <- bed2BWscore(file = file.path(data_path, "TFBS_regions.bed"), header = FALSE)
file.remove(file.path(data_path, "TFBS_regions.bed"))

TFBS_scores <- as.matrix(TFBS_scores)
TFBS_profile <- colMeans(TFBS_scores)
rm(TFBS_scores)

TFBS_profile_df <- cbind(TFBS_profile_df, TFBS_profile)
colnames(TFBS_profile_df)[i] <- target

rm(TFBS_profile)
gc()
i = i + 1

write.csv(TFBS_profile_df, file = file.path(data_path, "TFBS_profiles.csv"))



# 6b) creating ChIPseq experiment profiles for each file accession -------------

# 6b-1) cell line K562 ---------------------------------------------------------

# Set the directory path
K562_path <- file.path(data_path, "download/K562")


# cut all regions in all ChIPseq files to appropriate length
queensize(K562_path, window_size = 501, overwrite = TRUE)


# Get a list of all files in the directory
K562_files <- list.files(K562_path)


# Initialize an empty dataframe with 501 columns for profile positions
K562_profile_df <- data.frame(matrix(ncol = 501, nrow = 0))
names(K562_profile_df) <- c(-250:250)
K562_profile_df$target <- character(0)
K562_profile_df$accession <- character(0)

for (file in K562_files){
  # store filepath
  filepath <- file.path(K562_path, file)
  
  # store file accession
  file_accession <- sub("\\.sorted\\.bed$", "", file)
  cat("current file:", file_accession, "\n")
  
  # store transcription factor target
  target <- filtered_experiments_K562$`Experiment target`[filtered_experiments_K562$`File accession` == file_accession]
  cat("target is:", target, "\n\n")
  
  # calculate profile for the file accessions transcription factor target
  ChIPseq_scores <- bed2BWscore(file = filepath, header=FALSE)
  ChIPseq_scores <- as.matrix(ChIPseq_scores)
  ChIPseq_profile <- colMeans(ChIPseq_scores)
  rm(ChIPseq_scores)
  gc()
  
  # Create a new row with profile positions as separate columns
  new_row <- data.frame(t(ChIPseq_profile))
  names(new_row) <- c(-250:250)
  new_row$target <- target
  new_row$accession <- file_accession
  
  # Append the new row to the dataframe
  K562_profile_df <- rbind(K562_profile_df, new_row)
}





# clean up dataframe K562_profile_df
K562_profile_df <- t(K562_profile_df)
colnames(K562_profile_df) <- K562_profile_df[nrow(K562_profile_df),]
K562_profile_df <- K562_profile_df[-(nrow(K562_profile_df) - 2), ]
rownames(K562_profile_df) <- c(seq(-250, 250), "target", "accession")


# clean up environment
rm(file, filepath, file_accession, K562_files, K562_path, new_row)
gc()


# 6b-2) cell line GM12878 ------------------------------------------------------
# Set the directory path
GM12878_path <- file.path(data_path, "download/GM12878")


# cut all regions in all ChIPseq files to appropriate length
queensize(GM12878_path, window_size = 501, overwrite = TRUE)


# Get a list of all files in the directory
GM12878_files <- list.files(GM12878_path)

# Initialize an empty dataframe with 501 columns for profile positions
GM12878_profile_df <- data.frame(matrix(ncol = 501, nrow = 0))
names(GM12878_profile_df) <- c(-250:250)
GM12878_profile_df$target <- character(0)
GM12878_profile_df$accession <- character(0)

for (file in GM12878_files){
  # store filepath
  filepath <- file.path(GM12878_path, file)
  
  # store file accession
  file_accession <- sub("\\.sorted\\.bed$", "", file)
  cat("current file:", file_accession, "\n")
  
  # store transcription factor target
  target <- filtered_experiments_GM12878$`Experiment target`[filtered_experiments_GM12878$`File accession` == file_accession]
  cat("target is:", target, "\n\n")
  
  # calculate profile for the file accessions transcription factor target
  ChIPseq_scores <- bed2BWscore(file = filepath, header=FALSE)
  ChIPseq_scores <- as.matrix(ChIPseq_scores)
  ChIPseq_profile <- colMeans(ChIPseq_scores)
  rm(ChIPseq_scores)
  gc()
  
  # Create a new row with profile positions as separate columns
  new_row <- data.frame(t(ChIPseq_profile))
  names(new_row) <- c(-250:250)
  new_row$target <- target
  new_row$accession <- file_accession
  
  # Append the new row to the dataframe
  GM12878_profile_df <- rbind(GM12878_profile_df, new_row)
}

# clean up dataframe GM12878_profile_df
GM12878_profile_df <- t(GM12878_profile_df)
colnames(GM12878_profile_df) <- GM12878_profile_df[nrow(GM12878_profile_df),]
GM12878_profile_df <- GM12878_profile_df[-(nrow(GM12878_profile_df) - 2), ]
rownames(GM12878_profile_df) <- c(seq(-250, 250), "target", "accession")

# clean up environment
rm(file, filepath, file_accession, GM12878_files, GM12878_path, new_row)
gc()



# 7) calculating the correlation coefficients for all file accessions and their respective TFBS profile -----------------------------------------------------------
# 7a) cell line K562 
# create three new rows for the correlation values, Mann-Whitney U and p-value
K562_profile_df <- rbind(K562_profile_df, rep(NA, 93))
K562_profile_df <- rbind(K562_profile_df, rep(NA, 93))
K562_profile_df <- rbind(K562_profile_df, rep(NA, 93))
K562_profile_df <- rbind(K562_profile_df, rep(NA, 93))
rownames(K562_profile_df)[504:507] <- c("correlation", "correlation p-value", "Mann-Whitney U", "p-value")


for (i in 1:dim(K562_profile_df)[2]) {
  # extract tf target
  target_row_index <- which(rownames(K562_profile_df) == "target")
  
  target <- K562_profile_df[target_row_index, i]
  
  # compute correlation coefficient
  correlation <- cor.test(as.numeric(K562_profile_df[1:501, i]), TFBS_profile_df[[target]], method = "pearson")
  correlation_row_index <- which(rownames(K562_profile_df) == "correlation")
  
  K562_profile_df[correlation_row_index, i] <- correlation$estimate # add current correlation coefficient
  K562_profile_df[correlation_row_index+1, i] <- correlation$p.value # add pvalue of correlation
  
  # compute mann-whitney U and p-value
  MWUtest <- wilcox.test(as.numeric(K562_profile_df[1:501, i]), TFBS_profile_df[[target]], alternative = "two.sided")
  MWU_row_index <- which(rownames(K562_profile_df) == "Mann-Whitney U")
  
  K562_profile_df[MWU_row_index, i] <- MWUtest$statistic
  
  pvalue_row_index <- which(rownames(K562_profile_df) == "p-value")
  K562_profile_df[pvalue_row_index, i] <- MWUtest$p.value
}

# clean up environment
rm(i, target_row_index, target, correlation, correlation_row_index, MWUtest, MWU_row_index, pvalue_row_index)


# 7b) cell line GM12878 --------------------------------------------------------
# create three new rows for the correlation values, Mann-Whitney U and p-value
GM12878_profile_df <- rbind(GM12878_profile_df, rep(NA, 54))
GM12878_profile_df <- rbind(GM12878_profile_df, rep(NA, 54))
GM12878_profile_df <- rbind(GM12878_profile_df, rep(NA, 54))
GM12878_profile_df <- rbind(GM12878_profile_df, rep(NA, 54))
rownames(GM12878_profile_df)[504:507] <- c("correlation", "correlation p-value", "Mann-Whitney U", "p-value")


for (i in 1:dim(GM12878_profile_df)[2]) {
  # extract tf target
  target_row_index <- which(rownames(GM12878_profile_df) == "target")
  
  target <- GM12878_profile_df[target_row_index, i]
  
  # compute correlation coefficient
  correlation <- cor.test(as.numeric(GM12878_profile_df[1:501, i]), TFBS_profile_df[[target]], method = "pearson")
  correlation_row_index <- which(rownames(GM12878_profile_df) == "correlation")
  
  GM12878_profile_df[correlation_row_index, i] <- correlation$estimate # add current correlation coefficient
  GM12878_profile_df[correlation_row_index+1, i] <- correlation$p.value # add pvalue of correlation
  
  # compute mann-whitney U and p-value
  MWUtest <- wilcox.test(as.numeric(GM12878_profile_df[1:501, i]), TFBS_profile_df[[target]], alternative = "two.sided")
  MWU_row_index <- which(rownames(GM12878_profile_df) == "Mann-Whitney U")
  
  GM12878_profile_df[MWU_row_index, i] <- MWUtest$statistic
  
  pvalue_row_index <- which(rownames(GM12878_profile_df) == "p-value")
  GM12878_profile_df[pvalue_row_index, i] <- MWUtest$p.value
}


# clean up environment
rm(i, target_row_index, target, correlation, correlation_row_index, MWUtest, MWU_row_index, pvalue_row_index)


# 8)    Bonferroni p-value correction ------------------------------------------
# multiply p-values with the number of experiments
K562_profile_df[505,] <- as.numeric(K562_profile_df[505,]) * 93 # pvalue correction for the correlation pvalue
K562_profile_df[507,] <- as.numeric(K562_profile_df[507,]) * 93 # pvalue correction for the wilcoxon pvalue

GM12878_profile_df[505,] <- as.numeric(GM12878_profile_df[505,]) * 54 # pvalue correction for the correlation pvalue
GM12878_profile_df[507,] <- as.numeric(GM12878_profile_df[507,]) * 54 # pvalue correction for the wilcoxon pvalue


# 9)    Extract 10 experiments with positive and negative correlation ----------
# 9a)   K562 ---------------------------------------------------------------------

# Extract the correlation coefficients from row 504 
K562_correlation_coefficients <- as.numeric(K562_profile_df[504, ])

# Get the indices for the 10 most negatively and positively correlated samples
neg_corr_indices <- order(K562_correlation_coefficients)[1:11] # delete not significant result in col 3
pos_corr_indices <- order(-K562_correlation_coefficients)[1:10]

# Extract the samples
K562_neg_corr_samples <- K562_profile_df[, neg_corr_indices]
K562_pos_corr_samples <- K562_profile_df[, pos_corr_indices]

rm(neg_corr_indices, pos_corr_indices, K562_correlation_coefficients)


# 9b)   GM12878 ------------------------------------------------------------------

# Extract the correlation coefficients from row 504
GM12878_correlation_coefficients <- as.numeric(GM12878_profile_df[504, ])

# Get the indices for the 10 most negatively and positively correlated samples
neg_corr_indices <- order(GM12878_correlation_coefficients)[1:10]
pos_corr_indices <- order(-GM12878_correlation_coefficients)[1:12] # non significant results in col 1 and 12

# Extract the samples
GM12878_neg_corr_samples <- GM12878_profile_df[, neg_corr_indices]
GM12878_pos_corr_samples <- GM12878_profile_df[, pos_corr_indices]

rm(neg_corr_indices, pos_corr_indices, GM12878_correlation_coefficients)




# 10)   global profiles of NF scores -----------------------------------------------
# 9a) TFBS
globalTFBS_profile <- rowMeans(TFBS_profile_df)
names(globalTFBS_profile) <- seq(-250, 250)

# 9b) K562
k562cut <- as.data.frame(K562_profile_df[1:501,])

# Get the column names
col_names <- names(k562cut)

# Convert each column to numeric
for (col in col_names) {
  k562cut[[col]] <- as.numeric(k562cut[[col]])
}

globalK562_profile <- rowMeans(k562cut)
rm(col_names, k562cut)


# 9c) GM12878
GM12878cut <- as.data.frame(GM12878_profile_df[1:501,])

# Get the column names
col_names <- names(GM12878cut)

# Convert each column to numeric
for (col in col_names) {
  GM12878cut[[col]] <- as.numeric(GM12878cut[[col]])
}

globalGM12878_profile <- rowMeans(GM12878cut)
rm(col_names, GM12878cut)

# 11)   plot the global trends ---------------------------------------------------
# Generate positions from -250 to 250
positions <- seq(-250, 250)

# Create a data frame with positions and values from the vectors
df <- data.frame(
  Position = rep(positions, 3),
  Value = c(globalTFBS_profile, globalK562_profile, globalGM12878_profile),
  Group = rep(c("NF-score on predicted binding sites", "NF-score over real binding sites in K562", "NF-score over real binding sites in GM12878"), each = length(positions))
)

# Load ggplot2 package
library(ggplot2)

# Create the plot
ggplot(df, aes(x = Position, y = Value, color = Group)) +
  geom_line() +
  labs(x = "Position", y = "NF-score", color = "Data source") +
  theme_minimal()

rm(df, positions)








# 12)   plot one positive correlation example ----------------------------------
# done in plot_positiveCorrelation.R



# 13)   plot 10 negative correlation examples ----------------------------------
# done in plot_negativeCorrelation.R


# 14)   Tables of correlations -------------------------------------------------
# K562
K562_positive_corr_table <- K562_pos_corr_samples[502:507,]
K562_negative_corr_table <- K562_neg_corr_samples[502:507,]

  
# GM12878
write.csv(t(K562_neg_corr_samples[502:507,]), file = file.path(data_path, "K562_negativeCorrelation.csv"))
write.csv(t(K562_pos_corr_samples[502:507,]), file = file.path(data_path, "K562_positiveCorrelation.csv"))
write.csv(t(GM12878_pos_corr_samples[502:507,]), file = file.path(data_path, "GM12878_positiveCorrelation.csv"))
write.csv(t(GM12878_neg_corr_samples[502:507,]), file = file.path(data_path, "GM12878_negativeCorrelation.csv"))


# 15)   plotting the correlation distribution  ---------------------------------------
# K562
correlation_dataK562 <- data.frame(Correlation = as.numeric(K562_profile_df[504,]))

K562corr <- ggplot(correlation_dataK562, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.1, fill = "white", color = "black") +
  labs(title = "Correlation Values Distribution for K562",
       x = "Correlation",
       y = "Frequency")

# GM12878
correlation_dataGM12878 <- data.frame(Correlation = as.numeric(GM12878_profile_df[504,]))

GM12878corr <- ggplot(correlation_dataGM12878, aes(x = Correlation)) +
  geom_histogram(binwidth = 0.1, fill = "white", color = "black") +
  labs(title = "Correlation Values Distribution for GM12878",
       x = "Correlation",
       y = "Frequency")

# plot together
library(gridExtra)
test <- grid.arrange(K562corr, GM12878corr, ncol=2)
ggsave("./correlationHistogram_all.jpg", plot=test, dpi=300)


 
# 16)   differences between cell lines -----------------------------------------------
K562_MAFK <- as.numeric(K562_profile_df[1:501,1])
K562_ESRRA <- as.numeric(K562_profile_df[1:501,17])

GM12878_MAFK <- as.numeric(GM12878_profile_df[1:501,41])
GM12878_ESRRA <- as.numeric(GM12878_profile_df[1:501,13])

# Generate positions from -250 to 250
positions <- seq(-250, 250)

# Create a data frame with positions and values from the vectors
df_cellline_comparison <- data.frame(
  Position = rep(positions, 4),
  Value = c(K562_MAFK, GM12878_MAFK, K562_ESRRA, GM12878_ESRRA),
  Group = rep(c("MAFK in K562", "MAFK in GM12878", "ESRRA in K562", "ESRRA in GM12878"), each = length(positions))
)

# Load ggplot2 package
library(ggplot2)

# Create the plot
ggplot(df_cellline_comparison, aes(x = Position, y = Value, color = Group)) +
  geom_line(size=0.9) +
  labs(x = "Position", y = "NF-score", color = "Data source") +
  ggtitle("NF-score profiles of MAFK and ESRRA associated regions in K562 and GM12878") +
  theme_minimal() +
  theme(legend.position = "bottom",  # Set the legend below the plot
        legend.box = "vertical",   # Display legend items in a vertical stack
        legend.direction = "vertical")+
  guides(color = guide_legend(ncol = 2))  # Set legend to 2 columns



