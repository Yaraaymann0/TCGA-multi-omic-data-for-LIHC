# ============================================================
#      Transcriptome:  RNA-seq Using Limma
#          Dataset: Survival vs Decrease
# ============================================================
# WORKFLOW OVERVIEW:
#   1. Tools and R Packages used
#   2. Load count data & phenotype data
#   3. Harmonize sample IDs and Align the sample order
#   4. Explore & visualize the data distribution by histogram, Boxplot, heatmaps 
#   5. PCA plot 
#   6. Check normality by QQ plot & Shapiro test
#   7. Run limma for differential expression analysis
#   8. Filter significant DEGs and Export the table results
#   9. Volcano plot
#   10. Gene Enrichment & Pathway Analysis — gProfiler

# ============================================================


# ============================================================
# STEP 1: Tools and R Packages used 
# ============================================================

install.packages("ggrepel")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("gprofiler2")


#loading packages:
library(ggrepel) #for PCA
library(ggplot2) 
library(pheatmap)
library(RColorBrewer)
library(limma)  #for DE analysis 
library(gprofiler2)
library(ggplot2)

# ============================================================
# STEP 2: Load Count Data & Phenotype Data
# ============================================================

dataa <- read.csv2("~/Yara/Bioinformatics/Integrative Bioinformatics/project/exp.txt",
                   header=TRUE, 
                   row.names=1, 
                   sep="",
                   check.names = FALSE) 

pheno_data <- read.delim("~/Yara/Bioinformatics/Integrative Bioinformatics/project/liver.txt", 
                         row.names=1, 
                         header= TRUE,
                         sep= "\t",
                         check.names = FALSE)

str(dataa)
#first 5 genes and first 5 samples
print(dataa[1:5, 1:5])
# Values = FPKM/RPKM (continuous, pre-normalized — NOT raw counts)

table(pheno_data$vital_status)
table(pheno_data$gender)
table(pheno_data$sample_type)
#convert the data into numeric
dataa[] <- lapply(dataa, function(x) as.numeric(as.character(x)))
dataa <- as.matrix(dataa) # convert to matrix


# =============================================================================
# STEP 3: Harmonize sample IDs and Align the sample order 
# =============================================================================

colnames(dataa)
rownames(pheno_data)
# the samples name in the count table and phenotype are different in  barcode formatting; 
# the count table samples has barcodes written with dots, 
# while the phenotype data uses the same IDs with hyphens, 
# standardize them before matching samples by: 
colnames(dataa) <- gsub("\\.", "-", colnames(dataa))

#  --- Align the sample order --- 
phenoo= pheno_data[colnames(dataa), ]
all(colnames(dataa) == rownames(phenoo))



# =============================================================================
# STEP 4: Data Exploration & Visualization
# =============================================================================

# --- Data Dimensions ---
#data  → rows = genes= 20531,   columns = samples=423
#pheno → rows = samples= 438 , columns = metadata variables
dim(dataa)
dim(pheno_data)

# --- Check for missing values ---
sum(is.na(dataa))
sum(is.null(dataa))
sum(is.nan(dataa))
#no missing values 

# --- Histogram of LOG2-TRANSFORMED counts ---
hist(log2(dataa+1),col = "#1565C0", main = "Histogram — Log2 Counts")


# --- Data Distribution ---
#the data is negative binomial distribution

# --- Boxplot of the first 3 genes (log2-transformed) ---
# Each box represents one gene across all samples.
# Useful for spotting outlier genes or samples with unusual distributions.
boxplot(log2(dataa[1:3,] + 1),
        main = "Boxplot — First 3 Genes (log2 scale)",
        ylab = "log2(count + 1)",
        col  = "lightblue")


# --- Heatmap and Sample-to-Sample correlation  ---
#assess the correlation between expression of samples and cluster the data by similarity.

data_log <- log2(dataa + 1)
sampleDists<-as.dist(1 -cor(data_log))
pheatmap(sampleDists,
         main = "Sample-to-Sample Correlation",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         border_color = NA)



# =============================================================================
# STEP 5: PCA Plot to see sample structure 
# =============================================================================

# Remove zero-variance genes
zero_var <- apply(data_log, 1, var) == 0   # TRUE for genes with no variation
data_log_clean <- data_log[!zero_var, ]

pca_res <- prcomp(t(data_log_clean), scale. = TRUE, center = TRUE) # prcomp expects samples as ROWS, genes as COLUMNS → transpose
percentVar <- round(100 * summary(pca_res)$importance[2, 1:2])  # % variance for PC1, PC2

# --- PCA for Vital status --
pca_data <- data.frame(
  PC1          = pca_res$x[, 1],
  PC2          = pca_res$x[, 2],
  name         = rownames(pca_res$x),                         
  vital_status = phenoo[rownames(pca_res$x), "vital_status"]
)

ggplot(pca_data, aes(PC1, PC2, color = vital_status)) +
  geom_point(size = 3) +
  stat_ellipse() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()


# --- PCA for gender --

pca_data <- data.frame(
  PC1          = pca_res$x[, 1],
  PC2          = pca_res$x[, 2],
  name         = rownames(pca_res$x),                         
  gender = phenoo[rownames(pca_res$x), "gender"]
)

ggplot(pca_data, aes(PC1, PC2, color = gender)) +
  geom_point(size = 3) +
  stat_ellipse() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()

# --- PCA for sample type --

pca_data <- data.frame(
  PC1          = pca_res$x[, 1],
  PC2          = pca_res$x[, 2],
  name         = rownames(pca_res$x),                         
  sample_type = phenoo[rownames(pca_res$x), "sample_type"]
)

ggplot(pca_data, aes(PC1, PC2, color = sample_type)) +
  geom_point(size = 3) +
  stat_ellipse() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()


# =============================================================================
# STEP 6: Normality Check — QQ Plots & Shapiro-Wilk Test
# =============================================================================

# --- Q-Q Plot (visual normality check) --
#the 7rd gene in the matrix
qqnorm(data_log[7, ])
qqline(data_log[7, ], col = "red", lwd = 2)

# --- Shapiro Test ---
shapiro.test(as.numeric(data_log[100, ]))

#the p-value of first 100th gene is < 0.05, the data is not normally distributed  


# =============================================================================
# STEP 7: Run differential expression analysis using limma
# =============================================================================

# design matrix
phenoo$vital_status <- factor(phenoo$vital_status, levels = c("LIVING", "DECEASED"))
design <- model.matrix(~ vital_status, data = phenoo)

# fit model
fit <- lmFit(data_log_clean, design)
fit <- eBayes(fit)

# results
res <- topTable(fit, coef = "vital_statusDECEASED", number = Inf, adjust.method = "BH")


# =============================================================================
# STEP 8: Filter Significant DEGs and Export Results
# =============================================================================


#--- based on the significance padj < 0.05 and  significance FoldChange > 0.5 ---
res.filter = res[res$adj.P.Val < 0.05 & abs(res$logFC) > 0.5, ]
#log2FC > 0.5  means the gene is ~2.3x higher in decreased samples vs survival "Upregulated"
#log2FC < -0.5 means the gene is ~2.3x lower in decreased samples vs survival "Downregulated"

# --- Export the results of DEGs of decease/survive  ---
write.csv(res.filter, file = "DEG.csv", quote = F, row.names = T)

# =============================================================================
# STEP 9: Volcano Plot 
# =============================================================================


# --- Volcano plot  ---

#   X-axis → FoldChange (effect size: how much the gene changed)
#   Y-axis → -log10(padj)   (significance: higher = more significant)

#
# ONLY two genes in the top-right  → significantly UP in death/decrease   (red)
# RESET of Genes  → not significantly changed       (black)

# Reset plot layout to a single panel
par(mfrow = c(1, 1))

# Plot all genes as a background scatter (grey cloud)
with(res,
     plot(logFC, -log10(adj.P.Val),
          pch = 20,   
          main = "Decease vs Survive DEGs"))

# Highlight significantly UPREGULATED genes in red 
# (padj < 0.05 AND FC > +0.5 → higher in death/decrease)
with(subset(res, adj.P.Val < 0.05 & logFC > 0.5),
     points(logFC, -log10(adj.P.Val), pch = 20, col = "red"))

# Highlight significantly DOWNREGULATED genes in blue
# (padj < 0.05 AND log2FC < -0.5 → lower in death/decrease / higher in survive)
with(subset(res, adj.P.Val < 0.05 & logFC < -0.5),
     points(logFC, -log10(adj.P.Val), pch = 20, col = "blue"))

# Add a legend to explain the color coding
legend(x = -3.7, y = 9.5,
       c("upregulated", "downregulated"),
       cex = .8, bty = "n",
       col = c("red", "blue"),
       pch = 19)





# ============================================================
# STEP 10: Gene Enrichment & Pathway Analysis — gProfiler
# ============================================================

res$EntrezID <- sub("^X\\.+", "", rownames(res))
up_genes   <- res$EntrezID[res$adj.P.Val < 0.05 & res$logFC >  0.5]
down_genes <- res$EntrezID[res$adj.P.Val < 0.05 & res$logFC < -0.5]

# ── Run gProfiler on upregulated genes ──────────────────────
gp_up <- gost(
  query         = up_genes,
  organism      = "hsapiens",       # human
  ordered_query = FALSE,
  sources       = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
  significant   = TRUE,
  evcodes       = TRUE              # include evidence codes
)

# ── Run gProfiler on downregulated genes ────────────────────
gp_down <- gost(
  query         = down_genes,
  organism      = "hsapiens",
  ordered_query = FALSE,
  sources       = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
  significant   = TRUE,
  evcodes       = TRUE
)

# No result beacsue of the DEGs is only 2, and gProfiler needs a large gene list 
