# ═══════════════════════════════════════════════════════════════
##  TCGA-LIHC  |  Full 450K DNA Methylation Analysis
# ═══════════════════════════════════════════════════════════════


## Import Packages
required_pkgs <- c(
  "data.table", "limma", "minfi",
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "RColorBrewer", "missMethyl", "ggplot2",
  "ggrepel", "pheatmap", "stringr",
  "dplyr", "matrixStats"
)

library(matrixStats)
library(pheatmap)
library(patchwork)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg, ask = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

pal8 <- brewer.pal(8, "Dark2")

setwd("/Users/hossam/Downloads/integrative_project/methylation")

## Data Loading
METHY_FILE <- "/Users/hossam/Downloads/TCGA-LIHC-methylation450.tsv"

cat("Reading column names from Xena TSV")
header_line <- fread(METHY_FILE, nrows = 0)
all_cols    <- colnames(header_line)

# Probe-ID column is the non-TCGA column
probe_col <- all_cols[!grepl("^TCGA", all_cols)][1]
cat(sprintf("  Probe column name : '%s'\n", probe_col))

# Identify sample columns by tissue code before stripping the vial letter
# -01A / -01B = primary tumour; -11A/B = adjacent normal
tumour_cols <- all_cols[grepl("-01[AB]$", all_cols)]
normal_cols <- all_cols[grepl("-11[AB]$", all_cols)]

cat(sprintf("  Total samples     : %d\n", length(all_cols) - 1))
cat(sprintf("  Primary tumour    : %d\n", length(tumour_cols)))
cat(sprintf("  Adjacent normal   : %d\n", length(normal_cols)))

# Load only probe column + tumour columns
cols_to_load <- c(probe_col, tumour_cols)
cat("  Loading tumour columns …\n")
dt <- fread(METHY_FILE, select = cols_to_load, showProgress = FALSE)

probe_ids <- dt[[probe_col]]
betaMat   <- as.matrix(dt[, -1, with = FALSE])
rownames(betaMat) <- probe_ids
rm(dt, header_line); gc()

# Remove the ending letter A/B
colnames(betaMat) <- sub("[AB]$", "", colnames(betaMat))

# Remove duplicate barcodes (keep first vial if both -01A and -01B present)
dup_mask <- duplicated(colnames(betaMat))
if (any(dup_mask)) {
  cat(sprintf("  Removing %d duplicate barcodes (kept first vial)\n", sum(dup_mask)))
  betaMat <- betaMat[, !dup_mask]
}

cat(sprintf("Matrix loaded: %d probes × %d tumour samples\n",
            nrow(betaMat), ncol(betaMat)))
cat("Sample name preview:\n")
print(head(colnames(betaMat), 5))

## Clinical Matching
cat("Loading clinical & survival data")

clinical_raw <- read.table(
  "liver.txt", header = TRUE, sep = "\t",
  row.names = 1, check.names = FALSE,
  na.strings = c("", "NA")
)

surv_raw <- read.table(
  "survival.txt", header = TRUE, sep = "\t",
  row.names = 1, na.strings = c("", "NA")
)
surv_raw$Survival <- as.numeric(surv_raw$Survival)
surv_raw$Death    <- as.numeric(surv_raw$Death)

# Match matrix columns to clinical rows
shared_ids <- intersect(colnames(betaMat), rownames(clinical_raw))
betaMat    <- betaMat[, shared_ids]
clinical   <- clinical_raw[shared_ids, ]

cat(sprintf("  Tumour samples matched to clinical: %d\n", ncol(betaMat)))

# Key metadata
vital_status <- factor(clinical$vital_status,
                       levels = c("LIVING", "DECEASED"))
gender       <- factor(clinical$gender)
age          <- as.numeric(clinical$age_at_initial_pathologic_diagnosis)

stage_raw    <- clinical$pathologic_stage
stage_simple <- ifelse(
  grepl("Stage I[^IV]|Stage I$|Stage II[^IV]|Stage II$", stage_raw, perl = TRUE),
  "Early",
  ifelse(grepl("Stage III|Stage IV", stage_raw), "Late", NA)
)
stage_simple <- factor(stage_simple, levels = c("Early", "Late"))

cat("  vital_status distribution:\n"); print(table(vital_status, useNA = "ifany"))
cat("  gender distribution:\n");       print(table(gender,       useNA = "ifany"))
cat("  stage distribution:\n");        print(table(stage_simple, useNA = "ifany"))

## Probe Filtering
cat(sprintf("Probes before filtering: %d\n", nrow(betaMat)))

## 3a. Remove probes with > 20% missing values
na_frac <- rowMeans(is.na(betaMat))
betaMat <- betaMat[na_frac <= 0.20, ]
cat(sprintf("  After NA filter (>20%% missing): %d\n", nrow(betaMat)))

## 3b. Impute remaining NAs with row median
row_medians <- matrixStats::rowMedians(betaMat, na.rm = TRUE)
na_idx      <- which(is.na(betaMat), arr.ind = TRUE)
if (nrow(na_idx) > 0) {
  betaMat[na_idx] <- row_medians[na_idx[, 1]]
}
cat(sprintf("  NAs imputed with probe median: %d cells\n", nrow(na_idx)))

## 3c. Sex-chromosome probes
ann450k    <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
sex_probes <- rownames(ann450k)[ann450k$chr %in% c("chrX", "chrY")]
betaMat    <- betaMat[!rownames(betaMat) %in% sex_probes, ]
cat(sprintf("  After sex-chr removal: %d\n", nrow(betaMat)))

## 3d. SNP probes
snp_probes <- rownames(ann450k)[
  !is.na(ann450k$Probe_rs) | !is.na(ann450k$CpG_rs)
]
betaMat <- betaMat[!rownames(betaMat) %in% snp_probes, ]
cat(sprintf("  After SNP-probe removal: %d\n", nrow(betaMat)))

## Beta to M-Value Conversion
betaC <- pmax(pmin(betaMat, 0.999), 0.001)
mVals <- log2(betaC / (1 - betaC))
rm(betaMat); gc()

cat(sprintf("M-value matrix: %d probes × %d samples\n",
            nrow(mVals), ncol(mVals)))

# Colour vector aligned to vital_status (used in base-R plots below)
grp_col <- ifelse(vital_status == "DECEASED", pal8[1], pal8[2])

## Global Beta-Value Histogram
hist(betaC, breaks = 150,
     col    = "steelblue", border = "white",
     main   = "Global Beta-Value Distribution — TCGA-LIHC (all probes, tumour)",
     xlab   = "Beta Value",
     ylab   = "Frequency")

## Beta-Value Density: DECEASED vs LIVING
densityPlot(betaC, sampGroups = vital_status,
            main   = "Beta-Value Density: DECEASED vs LIVING",
            legend = FALSE, xlab = "Beta Value")
legend("topright", legend = levels(vital_status),
       text.col = pal8[1:2], bty = "n")

## M-Value Density: DECEASED vs LIVING
densityPlot(mVals, sampGroups = vital_status,
            main   = "M-Value Density: DECEASED vs LIVING",
            legend = FALSE, xlab = "M Value")
legend("topright", legend = levels(vital_status),
       text.col = pal8[1:2], bty = "n")

## Per-Sample QC Boxplot
set.seed(42)
# box_idx should be defined earlier in your Rmd; keep same logic if present
# Example (if it was there): box_idx <- sample(seq_len(ncol(betaC)), 80)
box_colors <- grp_col[box_idx]

dec_idx <- box_idx[vital_status[box_idx] == "DECEASED"]
liv_idx <- box_idx[vital_status[box_idx] == "LIVING"]
ordered_idx <- c(dec_idx, liv_idx)
ordered_col <- grp_col[ordered_idx]

par(mar = c(7, 4, 3, 1))
boxplot(betaC[, ordered_idx],
        col = ordered_col,
        las = 2, cex.axis = 0.3, outline = FALSE,
        main = "Per-Sample Beta Distribution QC (random 80-sample subset, grouped by status)",
        ylab = "Beta Value")
abline(v = sum(vital_status[ordered_idx] == "DECEASED") + 0.5,
       col = "black", lty = 2, lwd = 1.5)
legend("topright", legend = levels(vital_status),
       fill = pal8[1:2], bty = "n")
par(mar = c(5, 4, 4, 2))

## Heatmap: Top 1000 Most Variable CpGs
top1k <- order(matrixStats::rowVars(mVals), decreasing = TRUE)[1:1000]
heat_mat  <- t(scale(t(mVals[top1k, ])))

annot_col <- data.frame(
  VitalStatus = vital_status,
  Gender      = gender,
  Stage       = stage_simple,
  row.names   = colnames(heat_mat)
)
ann_colors <- list(
  VitalStatus = c(LIVING = pal8[2], DECEASED = pal8[1]),
  Gender      = c(MALE = "#4393C3", FEMALE = "#D6604D"),
  Stage       = c(Early = "#66C2A5", Late = "#FC8D62")
)

order_df <- annot_col
order_df$sample <- rownames(order_df)

order_df <- order_df[order(order_df$VitalStatus,
                           order_df$Gender,
                           order_df$Stage), ]

col_order <- order_df$sample

heat_mat   <- heat_mat[, col_order, drop = FALSE]
annot_col  <- annot_col[col_order, , drop = FALSE]

pheatmap(heat_mat,
         annotation_col    = annot_col,
         annotation_colors = ann_colors,
         show_colnames     = FALSE, show_rownames = FALSE,
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         main  = "Top 1000 Variable CpGs — Row-Scaled M-Values",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

## Heatmap: Top 50 Most Variable CpGs
top50   <- order(matrixStats::rowVars(mVals), decreasing = TRUE)[1:50]
heat50  <- t(scale(t(mVals[top50, ])))

annot_col$VitalStatus <- factor(annot_col$VitalStatus,
                                levels = c("DECEASED", "LIVING"))
annot_col$Gender      <- factor(annot_col$Gender,
                                levels = c("MALE", "FEMALE"))
annot_col$Stage       <- factor(annot_col$Stage,
                                levels = c("Early", "Late"))

order_df <- annot_col
order_df$sample <- rownames(order_df)

order_df <- order(order_df$VitalStatus,
                  order_df$Gender,
                  order_df$Stage)

col_order <- rownames(annot_col)[order_df]

heat50    <- heat50[, col_order, drop = FALSE]
annot_col <- annot_col[col_order, , drop = FALSE]

pheatmap(heat50,
         annotation_col    = annot_col,
         annotation_colors = ann_colors,
         show_colnames     = FALSE,
         show_rownames     = TRUE,
         fontsize_row      = 6,
         cluster_rows      = TRUE,
         cluster_cols      = FALSE,
         main  = "Top 50 Variable CpGs — Row-Scaled M-Values",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

## PCA / MDS — Sample Structure
top5k   <- order(rowVars(mVals), decreasing = TRUE)[1:5000]
pca_res <- prcomp(t(mVals[top5k, ]), scale. = FALSE)
pca_df  <- as.data.frame(pca_res$x[, 1:6])
pca_df$VitalStatus <- vital_status
pca_df$Gender      <- gender
pca_df$Stage       <- stage_simple
pca_df$Age         <- age

var_exp <- round(summary(pca_res)$importance[2, 1:6] * 100, 2)
cum_exp <- round(summary(pca_res)$importance[3, 1:6] * 100, 2)

pc_pairs <- list(
  c(1, 2), c(1, 3), c(1, 4),
  c(2, 3), c(2, 4), c(3, 4)
)

cat("Variance explained by PC1–PC5:\n")
print(var_exp)

## PCA — Vital Status
vs_cols <- c(LIVING = pal8[2], DECEASED = pal8[1])

p1 <- ggplot(pca_df, aes(PC1, PC2, colour = VitalStatus)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = vs_cols, na.value = "grey80") +
  labs(
    title  = "PC1 vs PC2",
    x      = paste0("PC1 (", var_exp[1], "%)"),
    y      = paste0("PC2 (", var_exp[2], "%)"),
    colour = "Vital Status"
  ) +
  theme_bw(base_size = 12)

p2 <- ggplot(pca_df, aes(PC1, PC3, colour = VitalStatus)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = vs_cols, na.value = "grey80") +
  labs(
    title  = "PC1 vs PC3",
    x      = paste0("PC1 (", var_exp[1], "%)"),
    y      = paste0("PC3 (", var_exp[3], "%)"),
    colour = "Vital Status"
  ) +
  theme_bw(base_size = 12)

p3 <- ggplot(pca_df, aes(PC1, PC4, colour = VitalStatus)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = vs_cols, na.value = "grey80") +
  labs(
    title  = "PC1 vs PC4",
    x      = paste0("PC1 (", var_exp[1], "%)"),
    y      = paste0("PC4 (", var_exp[4], "%)"),
    colour = "Vital Status"
  ) +
  theme_bw(base_size = 12)

p4 <- ggplot(pca_df, aes(PC2, PC3, colour = VitalStatus)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = vs_cols, na.value = "grey80") +
  labs(
    title  = "PC2 vs PC3",
    x      = paste0("PC2 (", var_exp[2], "%)"),
    y      = paste0("PC3 (", var_exp[3], "%)"),
    colour = "Vital Status"
  ) +
  theme_bw(base_size = 12)

p5 <- ggplot(pca_df, aes(PC2, PC4, colour = VitalStatus)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = vs_cols, na.value = "grey80") +
  labs(
    title  = "PC2 vs PC4",
    x      = paste0("PC2 (", var_exp[2], "%)"),
    y      = paste0("PC4 (", var_exp[4], "%)"),
    colour = "Vital Status"
  ) +
  theme_bw(base_size = 12)

p6 <- ggplot(pca_df, aes(PC3, PC4, colour = VitalStatus)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = vs_cols, na.value = "grey80") +
  labs(
    title  = "PC3 vs PC4",
    x      = paste0("PC3 (", var_exp[3], "%)"),
    y      = paste0("PC4 (", var_exp[4], "%)"),
    colour = "Vital Status"
  ) +
  theme_bw(base_size = 12)

p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(ncol = 3, guides = "collect")

## PCA — Gender
p1 <- ggplot(pca_df, aes(PC1, PC2, colour = Gender)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = c(MALE = "#4393C3", FEMALE = "#D6604D"),
                      na.value = "grey80") +
  labs(title  = "PC1 vs PC2",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "Gender") +
  theme_bw(base_size = 12)

p2 <- ggplot(pca_df, aes(PC1, PC3, colour = Gender)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = c(MALE = "#4393C3", FEMALE = "#D6604D"),
                      na.value = "grey80") +
  labs(title  = "PC1 vs PC3",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC3 (", var_exp[3], "%)"),
       colour = "Gender") +
  theme_bw(base_size = 12)

p3 <- ggplot(pca_df, aes(PC1, PC4, colour = Gender)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = c(MALE = "#4393C3", FEMALE = "#D6604D"),
                      na.value = "grey80") +
  labs(title  = "PC1 vs PC4",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC4 (", var_exp[4], "%)"),
       colour = "Gender") +
  theme_bw(base_size = 12)

p4 <- ggplot(pca_df, aes(PC2, PC3, colour = Gender)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = c(MALE = "#4393C3", FEMALE = "#D6604D"),
                      na.value = "grey80") +
  labs(title  = "PC2 vs PC3",
       x = paste0("PC2 (", var_exp[2], "%)"),
       y = paste0("PC3 (", var_exp[3], "%)"),
       colour = "Gender") +
  theme_bw(base_size = 12)

p5 <- ggplot(pca_df, aes(PC2, PC4, colour = Gender)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = c(MALE = "#4393C3", FEMALE = "#D6604D"),
                      na.value = "grey80") +
  labs(title  = "PC2 vs PC4",
       x = paste0("PC2 (", var_exp[2], "%)"),
       y = paste0("PC4 (", var_exp[4], "%)"),
       colour = "Gender") +
  theme_bw(base_size = 12)

p6 <- ggplot(pca_df, aes(PC3, PC4, colour = Gender)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = c(MALE = "#4393C3", FEMALE = "#D6604D"),
                      na.value = "grey80") +
  labs(title  = "PC3 vs PC4",
       x = paste0("PC3 (", var_exp[3], "%)"),
       y = paste0("PC4 (", var_exp[4], "%)"),
       colour = "Gender") +
  theme_bw(base_size = 12)

p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(ncol = 3, guides = "collect")

## PCA — Pathologic Stage
p1 <- ggplot(pca_df, aes(PC1, PC2, colour = Stage)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = c(Early = "#66C2A5", Late = "#FC8D62"),
                      na.value = "grey80") +
  labs(title  = "PC1 vs PC2",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "Stage") +
  theme_bw(base_size = 12)

p2 <- ggplot(pca_df, aes(PC1, PC3, colour = Stage)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = c(Early = "#66C2A5", Late = "#FC8D62"),
                      na.value = "grey80") +
  labs(title  = "PC1 vs PC3",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC3 (", var_exp[3], "%)"),
       colour = "Stage") +
  theme_bw(base_size = 12)

p3 <- ggplot(pca_df, aes(PC1, PC4, colour = Stage)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = c(Early = "#66C2A5", Late = "#FC8D62"),
                      na.value = "grey80") +
  labs(title  = "PC1 vs PC4",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC4 (", var_exp[4], "%)"),
       colour = "Stage") +
  theme_bw(base_size = 12)

p4 <- ggplot(pca_df, aes(PC2, PC3, colour = Stage)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = c(Early = "#66C2A5", Late = "#FC8D62"),
                      na.value = "grey80") +
  labs(title  = "PC2 vs PC3",
       x = paste0("PC2 (", var_exp[2], "%)"),
       y = paste0("PC3 (", var_exp[3], "%)"),
       colour = "Stage") +
  theme_bw(base_size = 12)

p5 <- ggplot(pca_df, aes(PC2, PC4, colour = Stage)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = c(Early = "#66C2A5", Late = "#FC8D62"),
                      na.value = "grey80") +
  labs(title  = "PC2 vs PC4",
       x = paste0("PC2 (", var_exp[2], "%)"),
       y = paste0("PC4 (", var_exp[4], "%)"),
       colour = "Stage") +
  theme_bw(base_size = 12)

p6 <- ggplot(pca_df, aes(PC3, PC4, colour = Stage)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_manual(values = c(Early = "#66C2A5", Late = "#FC8D62"),
                      na.value = "grey80") +
  labs(title  = "PC3 vs PC4",
       x = paste0("PC3 (", var_exp[3], "%)"),
       y = paste0("PC4 (", var_exp[4], "%)"),
       colour = "Stage") +
  theme_bw(base_size = 12)

p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(ncol = 3, guides = "collect")

## PCA — Age at Diagnosis
p1 <- ggplot(pca_df, aes(PC1, PC2, colour = Age)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_gradient(low = "#FEE0D2", high = "#A50F15", na.value = "grey80") +
  labs(title = "PC1 vs PC2",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)"),
       colour = "Age") +
  theme_bw(base_size = 12)

p2 <- ggplot(pca_df, aes(PC1, PC3, colour = Age)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_gradient(low = "#FEE0D2", high = "#A50F15", na.value = "grey80") +
  labs(title = "PC1 vs PC3",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC3 (", var_exp[3], "%)"),
       colour = "Age") +
  theme_bw(base_size = 12)

p3 <- ggplot(pca_df, aes(PC1, PC4, colour = Age)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_gradient(low = "#FEE0D2", high = "#A50F15", na.value = "grey80") +
  labs(title = "PC1 vs PC4",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC4 (", var_exp[4], "%)"),
       colour = "Age") +
  theme_bw(base_size = 12)

p4 <- ggplot(pca_df, aes(PC2, PC3, colour = Age)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_gradient(low = "#FEE0D2", high = "#A50F15", na.value = "grey80") +
  labs(title = "PC2 vs PC3",
       x = paste0("PC2 (", var_exp[2], "%)"),
       y = paste0("PC3 (", var_exp[3], "%)"),
       colour = "Age") +
  theme_bw(base_size = 12)

p5 <- ggplot(pca_df, aes(PC2, PC4, colour = Age)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_gradient(low = "#FEE0D2", high = "#A50F15", na.value = "grey80") +
  labs(title = "PC2 vs PC4",
       x = paste0("PC2 (", var_exp[2], "%)"),
       y = paste0("PC4 (", var_exp[4], "%)"),
       colour = "Age") +
  theme_bw(base_size = 12)

p6 <- ggplot(pca_df, aes(PC3, PC4, colour = Age)) +
  geom_point(size = 1.8, alpha = 0.75) +
  scale_colour_gradient(low = "#FEE0D2", high = "#A50F15", na.value = "grey80") +
  labs(title = "PC3 vs PC4",
       x = paste0("PC3 (", var_exp[3], "%)"),
       y = paste0("PC4 (", var_exp[4], "%)"),
       colour = "Age") +
  theme_bw(base_size = 12)

p1 + p2 + p3 + p4 + p5 + p6 +
  plot_layout(ncol = 3, guides = "collect")

## Scree Plot
scree_df <- data.frame(
  PC     = paste0("PC", 1:5),
  VarExp = var_exp[1:5]
)
scree_df$PC <- factor(scree_df$PC, levels = scree_df$PC)

ggplot(scree_df, aes(PC, VarExp)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = paste0(VarExp, "%")), vjust = -0.4, size = 4) +
  labs(
    title = "PCA Scree Plot",
    x     = NULL,
    y     = "% Variance Explained"
  ) +
  theme_bw(base_size = 12)

## MDS — Vital Status & Gender
par(mfrow = c(1, 2))
grp_col_full <- ifelse(vital_status == "DECEASED", pal8[1], pal8[2])
plotMDS(mVals[top5k, ], top = 1000, gene.selection = "common",
        col = grp_col_full, pch = 19, cex = 0.6,
        main = "MDS — Vital Status")
legend("topright", legend = levels(vital_status),
       col = pal8[1:2], pch = 19, bty = "n")

plotMDS(mVals[top5k, ], top = 1000, gene.selection = "common",
        col = ifelse(gender == "MALE", "#4393C3", "#D6604D"),
        pch = 19, cex = 0.6, main = "MDS — Gender")
legend("topright", legend = c("MALE", "FEMALE"),
       col = c("#4393C3", "#D6604D"), pch = 19, bty = "n")
par(mfrow = c(1, 1))

## Differential Methylation: DECEASED vs LIVING
valid_vs   <- !is.na(vital_status)
mVals_dea  <- mVals[, valid_vs]
betaC_dea  <- betaC[, valid_vs]
vital_dea  <- droplevels(vital_status[valid_vs])
gender_dea <- droplevels(gender[valid_vs])

cat(sprintf("Samples used: %d  (DECEASED=%d, LIVING=%d)\n",
            ncol(mVals_dea),
            sum(vital_dea == "DECEASED"),
            sum(vital_dea == "LIVING")))

design_dea <- model.matrix(
  ~ 0 + vital_dea + gender_dea,
  data = data.frame(vital_dea, gender_dea)
)
colnames(design_dea) <- gsub("vital_dea|gender_dea", "", colnames(design_dea))

fit1  <- lmFit(mVals_dea, design_dea)
conts <- makeContrasts(DECEASED_vs_LIVING = DECEASED - LIVING,
                       levels = design_dea)
fit2  <- contrasts.fit(fit1, conts)
fit2  <- eBayes(fit2)

cat("Significant CpGs (FDR < 0.05):\n")
print(summary(decideTests(fit2, p.value = 0.05)))

## Annotation & Effect Size
DMPs_raw      <- topTable(fit2, num = Inf,
                          coef = "DECEASED_vs_LIVING", sort.by = "p")
DMPs_raw$Name <- rownames(DMPs_raw)

ann_cols_wanted <- c("Name", "chr", "pos", "strand",
                     "UCSC_RefGene_Name", "Relation_to_Island")
ann_sub <- ann450k[ann450k$Name %in% DMPs_raw$Name, ann_cols_wanted]
ann_sub <- as.data.frame(ann_sub)

DMPs_full <- merge(DMPs_raw, ann_sub, by = "Name", all.x = TRUE, sort = FALSE)
DMPs_full <- DMPs_full[order(DMPs_full$P.Value), ]
rownames(DMPs_full) <- DMPs_full$Name

beta_DEC           <- rowMeans(betaC_dea[, vital_dea == "DECEASED"], na.rm = TRUE)
beta_LIV           <- rowMeans(betaC_dea[, vital_dea == "LIVING"],   na.rm = TRUE)
delta_beta         <- beta_DEC - beta_LIV
DMPs_full$delta_beta <- delta_beta[DMPs_full$Name]

sig_mask <- !is.na(DMPs_full$adj.P.Val) &
  DMPs_full$adj.P.Val < 0.05  &
  abs(DMPs_full$logFC) > 0.5
sig_DMPs <- DMPs_full[sig_mask, ]

cat(sprintf("\nSignificant DMPs (FDR<0.05, |logFC|>0.5): %d\n", nrow(sig_DMPs)))
cat(sprintf("  Hypermethylated in DECEASED (logFC > +0.5): %d\n",
            sum(sig_DMPs$logFC >  0.5, na.rm = TRUE)))
cat(sprintf("  Hypomethylated  in DECEASED (logFC < -0.5): %d\n",
            sum(sig_DMPs$logFC < -0.5, na.rm = TRUE)))
cat(sprintf("  Δβ range of sig. DMPs: %.3f to %.3f\n",
            min(sig_DMPs$delta_beta, na.rm = TRUE),
            max(sig_DMPs$delta_beta, na.rm = TRUE)))

## P-Value Histogram
hist(DMPs_full$P.Value, breaks = 50, col = "steelblue", border = "white",
     main = "P-Value Histogram",
     xlab = "Raw P-Value")
abline(v = 0.05, col = "red", lty = 2)

## Save Result Tables
write.csv(DMPs_full,
          file = "DMP_full_results_DECEASED_vs_LIVING.csv",
          row.names = FALSE)

sig_cols <- c("Name", "chr", "pos", "strand", "UCSC_RefGene_Name",
              "Relation_to_Island", "logFC", "AveExpr", "t",
              "P.Value", "adj.P.Val", "delta_beta")
write.csv(sig_DMPs[, sig_cols[sig_cols %in% colnames(sig_DMPs)]],
          file = "DMP_significant_DECEASED_vs_LIVING.csv",
          row.names = FALSE)

cat("Tables saved:\n")
cat("  DMP_full_results_DECEASED_vs_LIVING.csv\n")
cat("  DMP_significant_DECEASED_vs_LIVING.csv\n")

## Volcano Plot — Primary
gene_col <- if ("UCSC_RefGene_Name" %in% colnames(DMPs_full)) {
  sapply(str_split(DMPs_full$UCSC_RefGene_Name, ";"), `[`, 1)
} else rep(NA_character_, nrow(DMPs_full))

vol_df <- data.frame(
  CpG        = DMPs_full$Name,
  logFC      = DMPs_full$logFC,
  delta_beta = DMPs_full$delta_beta,
  neglog10p  = pmin(-log10(DMPs_full$P.Value), 50),
  FDR        = DMPs_full$adj.P.Val,
  Gene       = gene_col,
  stringsAsFactors = FALSE
)
vol_df$Group <- "Not Significant"
vol_df$Group[vol_df$FDR < 0.05 & vol_df$logFC >  0.5] <- "Hyper in DECEASED"
vol_df$Group[vol_df$FDR < 0.05 & vol_df$logFC < -0.5] <- "Hypo in DECEASED"

label_df <- vol_df %>%
  filter(Group != "Not Significant", !is.na(Gene), Gene != "") %>%
  slice_min(FDR, n = 20)

ggplot(vol_df, aes(logFC, neglog10p, colour = Group)) +
  geom_point(size = 0.7, alpha = 0.45) +
  scale_colour_manual(values = c("Hyper in DECEASED" = "firebrick",
                                 "Hypo in DECEASED"  = "royalblue",
                                 "Not Significant"   = "grey75")) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "black", linewidth = 0.4) +
  geom_text_repel(data = label_df, aes(label = Gene),
                  size = 2.5, max.overlaps = 25,
                  colour = "black", segment.color = "grey50",
                  box.padding = 0.35, min.segment.length = 0.2) +
  annotate("text", x =  Inf, y = Inf,
           label = paste0("Hyper: ", sum(vol_df$Group == "Hyper in DECEASED")),
           hjust = 1.1, vjust = 1.5, colour = "firebrick", size = 3.5) +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("Hypo: ", sum(vol_df$Group == "Hypo in DECEASED")),
           hjust = -0.1, vjust = 1.5, colour = "royalblue", size = 3.5) +
  labs(title    = "Volcano — DECEASED vs LIVING (TCGA-LIHC Full 450K)",
       subtitle = "FDR < 0.05  |  |logFC (M-value)| > 0.5",
       x        = "logFC [M-value; DECEASED − LIVING]",
       y        = "−log10(p-value)",
       colour   = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank())

## Top 6 Significant CpGs - Boxplots
avail_sig <- sig_DMPs$Name[sig_DMPs$Name %in% rownames(betaC_dea)]
top6_sig  <- head(avail_sig, 6)

if (length(top6_sig) > 0) {
  par(mfrow = c(2, 