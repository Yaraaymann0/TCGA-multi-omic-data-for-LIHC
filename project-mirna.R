###############################################################################
# PIPELINE STEP 1: Load Required Packages
###############################################################################

library(BiocManager)
library(DESeq2)
library(NMF)
library(readxl)
library(ggplot2)
library(pheatmap)

###############################################################################
# PIPELINE STEP 2: Import miRNA Expression and Clinical Data
###############################################################################

# Load miRNA expression matrix
data = as.matrix(read.table("D:/C/Desktop/NU/Integrative/Project/mirna", header = TRUE, row.names = 1, sep = ""))
head(data)

# Loading the Clinical data file
install.packages("readxl")
library(readxl)

clinical_df <- read_excel(file.choose(), sheet = 1)

###############################################################################
# PIPELINE STEP 3: Clinical Data Cleaning and Sample Harmonization
###############################################################################

# Converting the excel file into a standard data frame 
clinical <- as.data.frame(clinical_df)

# Matching the barcodes: Converting dashes (-) in clinical sampleID to periods (.)
clinical$sampleID_dot = gsub("-", ".", clinical$sampleID)

# Identifying the common samples between both files
common_samples = intersect(colnames(data), clinical$sampleID_dot)

# Filtering and aligning both datasets so the samples are in the exact same order
mirna_clean = data[, common_samples]
clinical_clean = clinical[match(common_samples, clinical$sampleID_dot), ]

# Verify sample correspondence
all(colnames(mirna_clean) == clinical_clean$sampleID_dot)
dim(mirna_clean)
###############################################################################
# PIPELINE STEP 4: Data Quality Assessment and Normalization
###############################################################################

# Visualize raw expression distribution
png("miRNA_distribution_plots.png", width = 1000, height = 500, res = 120)

# Setting layout to 1 row, 2 columns
par(mfrow = c(1, 2))

# Drawing the first histogram (Raw data)
hist(mirna_clean, main = "Raw miRNA Expression", 
     xlab = "Expression Value", col = "skyblue", breaks = 50)

# Log2 transformation

mirna_log = log2(mirna_clean + 1)

# Drawing the second histogram into the file
hist(mirna_log, main = "Log2 Transformed miRNA", 
     xlab = "Log2(Expression + 1)", col = "salmon", breaks = 50)

# Closing and savi the file
dev.off()
# Visualize normalized distribution

# Normality assessment
# - Boxplots
boxplot(log2(mirna_clean[1:5,]+1))
# - Q-Q plots
qqnorm(mirna_clean[3,])
qqline(mirna_clean[3,])
# - Shapiro-Wilk test
shapiro.test(mirna_clean[100,])

###############################################################################
# PIPELINE STEP 5: Principal Component Analysis (PCA)
###############################################################################

# Calculating the variance for each microRNA row
row_variances = apply(mirna_log, 1, var)

# # Filter non-variable miRNAs (Keeping only miRNAs where variance is greater than 0)
variable_mirnas = mirna_log[row_variances > 0, ]

# Running PCA on the transposed filtered matrix
pca_result = prcomp(t(variable_mirnas), scale. = TRUE)

# Creating the PCA coordinates data frame
pca_df = data.frame(PC1 = pca_result$x[, 1], 
                    PC2 = pca_result$x[, 2], 
                    Vital_Status = clinical_clean$vital_status,
                    Gender = clinical_clean$gender)


# Calculate variance explained
var_explained = (pca_result$sdev^2) / sum(pca_result$sdev^2) * 100
# Visualize PCA by:
# Step A (Part 1) Plotting PC1 vs PC2 colored by Vital Status
ggplot(pca_df, aes(x = PC1, y = PC2, color = Vital_Status)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA of miRNA Expression Colored by Vital Status",
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)")) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

#Step B (Part 2): Re-color by Gender
ggplot(pca_df, aes(x = PC1, y = PC2, color = Gender)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "PCA of miRNA Expression Colored by Gender",
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)")) +
  theme_minimal() +
  scale_color_brewer(palette = "Set2")

###############################################################################
# PIPELINE STEP 6: Heatmap Visualization of Highly Variable miRNAs
###############################################################################

# Calculate the variance

row_variances <- apply(mirna_log, 1, var)

# Selecting the top 50 most variable miRNAs
top50 <- names(sort(row_variances, decreasing = TRUE))[1:50]

# Creating sample annotations
annotation <- data.frame(
  Vital_Status = clinical_clean$vital_status,
  Gender = clinical_clean$gender
)

rownames(annotation) <- clinical_clean$sampleID_dot

# Plotting the heatmap
pheatmap(
  mirna_log[top50, ],
  scale = "row",
  annotation_col = annotation,
  show_rownames = TRUE,
  show_colnames = FALSE
)

###############################################################################
# PIPELINE STEP 7: Differential Expression Analysis Using DESeq2
###############################################################################

# Creating a metadata table with vital status as the experimental condition
coldata <- data.frame(
  row.names = colnames(mirna_clean),
  vital_status = factor(clinical_clean$vital_status)
)

# Checking the groups
table(coldata$vital_status)

# Ensuring "LIVING" is the reference level
coldata$vital_status <- relevel(coldata$vital_status, ref = "LIVING")

# Creating the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(mirna_clean),
  colData = coldata,
  design = ~ vital_status
)

# Removing very lowly expressed miRNAs
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Running the differential expression analysis
dds <- DESeq(dds)

# Extracting results
res <- results(dds)

# Converting to a data frame
de_results <- as.data.frame(res)
# Adding miRNA names
de_results$miRNA <- rownames(de_results)

# Removing missing values
de_results <- na.omit(de_results)

# Renaming columns
colnames(de_results)[which(colnames(de_results)=="log2FoldChange")] <- "Log2FC"
colnames(de_results)[which(colnames(de_results)=="pvalue")] <- "PValue"
colnames(de_results)[which(colnames(de_results)=="padj")] <- "Padj"
###############################################################################
# PIPELINE STEP 8: Identification of Differentially Expressed miRNAs
###############################################################################

# Classify significant miRNAs
de_results$Significance <- "Not Significant"

de_results$Significance[
  de_results$Padj < 0.05 &
    de_results$Log2FC > 0.5
] <- "Up-regulated (Deceased)"

de_results$Significance[
  de_results$Padj < 0.05 &
    de_results$Log2FC < -0.5
] <- "Down-regulated (Deceased)"

# Generate volcano plot
library(ggplot2)

ggplot(de_results,
       aes(x = Log2FC,
           y = -log10(Padj),
           color = Significance)) +
  
  geom_point(size = 2, alpha = 0.8) +
  
  scale_color_manual(values = c(
    "Down-regulated (Deceased)" = "skyblue",
    "Not Significant" = "black",
    "Up-regulated (Deceased)" = "firebrick"
  )) +
  
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  
  geom_vline(xintercept = c(-0.5,0.5),
             linetype = "dashed") +
  
  labs(
    title = "Volcano Plot: Differentially Expressed miRNAs",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-value"
  ) +
  
  theme_minimal()
###############################################################################
# PIPELINE STEP 9: Export Significant Differentially Expressed miRNAs
###############################################################################
significant_mirnas <- subset(
  de_results,
  Padj < 0.05 &
    abs(Log2FC) > 0.5
)

significant_mirnas <- significant_mirnas[
  order(significant_mirnas$Padj),
]

head(significant_mirnas)

write.csv(
  significant_mirnas,
  "significant_differential_miRNAs_DESeq2.csv",
  row.names = FALSE
)

