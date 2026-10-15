###############################################################################
# PIPELINE STEP 1: Load Required Packages
###############################################################################

library(MOFA2)
library(readxl)
library(survival)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

###############################################################################
# PIPELINE STEP 2: Import Multi-Omics Datasets
###############################################################################

# mRNA
mRNA_data <- as.matrix(
  read.table("exp", header = TRUE, row.names = 1, check.names = FALSE)
)




# DNA methylation
methylation_data <- as.matrix(
  read.table("methy", header = TRUE, row.names = 1, check.names = FALSE)
)
# miRNA
miRNA_data <- as.matrix(
  read.table("mirna", header = TRUE, row.names = 1, check.names = FALSE)
  
)
###############################################################################
# PIPELINE STEP 3: Data Preprocessing and Normalization
###############################################################################

# Convert methylation β-values to M-values
methylation_data[methylation_data == 0] <- 1e-6
methylation_data[methylation_data == 1] <- 1 - 1e-6

methylation_data <- log2(methylation_data / (1 - methylation_data))

# Log2-transform mRNA
mRNA_data <- log2(mRNA_data + 1)

# Log2-transform miRNA
miRNA_data <- log2(miRNA_data + 1)

###############################################################################
# PIPELINE STEP 4: Harmonize Samples Across Omics Layers
###############################################################################

# Identify common samples
common_samples <- Reduce(
  intersect,
  list(
    colnames(mRNA_data),
    colnames(methylation_data),
    colnames(miRNA_data)
  )
)

cat("Common samples:", length(common_samples), "\n")

# Match sample order
mRNA_data <- mRNA_data[, common_samples]
methylation_data <- methylation_data[, common_samples]
miRNA_data <- miRNA_data[, common_samples]

# Verify identical sample order
identical(colnames(mRNA_data), colnames(methylation_data))
identical(colnames(mRNA_data), colnames(miRNA_data))

###############################################################################
# PIPELINE STEP 5: Feature Selection
###############################################################################

# Keep 5000 most variable genes
gene_var <- apply(mRNA_data, 1, var)
mRNA_data <- mRNA_data[
  order(gene_var, decreasing = TRUE)[1:min(5000, nrow(mRNA_data))],
]

# Keep 10000 most variable CpGs
meth_var <- apply(methylation_data, 1, var)
methylation_data <- methylation_data[
  order(meth_var, decreasing = TRUE)[1:min(10000, nrow(methylation_data))],
]

# Keep 500 most variable miRNAs
mirna_var <- apply(miRNA_data, 1, var)
miRNA_data <- miRNA_data[
  order(mirna_var, decreasing = TRUE)[1:min(500, nrow(miRNA_data))],
]

###############################################################################
# PIPELINE STEP 6: Construct the MOFA Input Object
###############################################################################

# Create data list
data_list <- list(
  mRNA = mRNA_data,
  Methylation = methylation_data,
  miRNA = miRNA_data
)

lapply(data_list, dim)

# Create MOFA object

MOFAobject <- create_mofa(data_list)

plot_data_overview(MOFAobject)

###############################################################################
# PIPELINE STEP 7: Configure MOFA Model
###############################################################################

# Data options

data_opts <- get_default_data_options(MOFAobject)

data_opts$scale_views <- TRUE

# Model options
model_opts <- get_default_model_options(MOFAobject)

model_opts$likelihoods <- c(
  "gaussian",
  "gaussian",
  "gaussian"
)

# Training options
train_opts <- get_default_training_options(MOFAobject)

train_opts$seed <- 42
train_opts$verbose <- TRUE
train_opts$convergence_mode <- "slow"
train_opts$maxiter <- 10000

###############################################################################
# PIPELINE STEP 8: Train the MOFA Model
###############################################################################

# Prepare model

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
# Run MOFA and Save trained model

outfile <- file.path(getwd(), "TCGA_MOFA_model.hdf5")

MOFAobject <- run_mofa(
  MOFAobject,
  outfile = outfile
)

###############################################################################
# PIPELINE STEP 9: Load Trained Model and Integrate Clinical Data
###############################################################################

# Load model

MOFAobject <- load_model("TCGA_MOFA_model.hdf5")
#Load clinical data
clinical <- read_excel("Clinical_data.xlsx")
clinical <- as.data.frame(clinical)

clinical$sample <- gsub("-", ".", clinical$sampleID)

rownames(clinical) <- clinical$sample
clinical$sampleID <- NULL

# Match samples
samples <- samples_names(MOFAobject)[[1]]
samples <- gsub("-", ".", samples)

common <- intersect(samples, rownames(clinical))

length(common)

MOFAobject <- subset_samples(MOFAobject,
                             samples=common)

clinical <- clinical[common, ]
# Add metadata
#Add metadata
clinical$sample <- rownames(clinical)

clinical <- clinical[
  ,
  c("sample",
    setdiff(colnames(clinical),"sample"))
]

samples_metadata(MOFAobject) <- clinical

#Verify

head(samples_metadata(MOFAobject))
###############################################################################
# PIPELINE STEP 10: Explore Latent Factors
###############################################################################

#Create survival variables
clinical$event <-
  ifelse(clinical$vital_status=="Deceaced",1,0)

clinical$time <-
  ifelse(
    !is.na(clinical$days_to_death),
    clinical$days_to_death,
    clinical$days_to_last_followup
  )

clinical$time <- as.numeric(clinical$time)
clinical$event <- as.numeric(clinical$event)

samples_metadata(MOFAobject) <- clinical
#Data overview
plot_data_overview(MOFAobject)
#Variance explained
r2 <- get_variance_explained(MOFAobject)

r2$r2_total

head(r2$r2_per_factor)

plot_variance_explained(
  MOFAobject,
  plot_total=TRUE
)


#Factor correlations
plot_factor_cor(MOFAobject)

# Factor visualization
plot_factor(
  MOFAobject,
  factors=1:10
)

plot_factor(
  MOFAobject,
  factors=1:10,
  color_by="pathologic_stage"
)

plot_factor(
  MOFAobject,
  factors=1:10,
  color_by="vital_status"
)

plot_factor(
  MOFAobject,
  factors=1:10,
  color_by="gender"
)
plot_factor(
  MOFAobject,
  factors = 1:10,
  color_by = "age_at_initial_pathologic_diagnosis"
)


#UMAP
set.seed(123)

MOFAobject <- run_umap(MOFAobject)

plot_dimred(
  
  MOFAobject,
  
  method="UMAP",
  
  color_by="pathologic_stage"
  
)

###############################################################################
# PIPELINE STEP 11: Clinical Association Analysis
###############################################################################

correlate_factors_with_covariates(
  
  MOFAobject,
  
  covariates=c(
    
    "age_at_initial_pathologic_diagnosis",
    
    "gender",
    
    "pathologic_stage",
    
    "vital_status"
    
  ),
  
  plot="r"
  
)

###############################################################################
# PIPELINE STEP 12: Identify Key Molecular Drivers
###############################################################################

# Top mRNA weights
plot_top_weights(
  
  MOFAobject,
  
  view="mRNA",
  
  factor=1,
  
  nfeatures=30
  
)
# Top methylation weights
plot_top_weights(
  
  MOFAobject,
  
  view="Methylation",
  
  factor=1,
  
  nfeatures=30
  
)
# Top miRNA weights
plot_top_weights(
  
  MOFAobject,
  
  view="miRNA",
  
  factor=1,
  
  nfeatures=30
  
)

###############################################################################
# PIPELINE STEP 13: Visualize Feature Weight Landscapes
###############################################################################

# Weight heatmaps
plot_weights_heatmap(
  MOFAobject,
  view="mRNA",
  factors=1:10
)

plot_weights_heatmap(
  MOFAobject,
  view="Methylation",
  factors=1:10
)

plot_weights_heatmap(
  MOFAobject,
  view="miRNA",
  factors=1:10
)
# Data heatmaps
plot_data_heatmap(
  
  MOFAobject,
  
  view="mRNA",
  
  factor=1,
  
  features=50,
  
  cluster_rows=TRUE,
  
  cluster_cols=FALSE
  
)
#Scatter plots
plot_data_scatter(
  
  MOFAobject,
  
  view="mRNA",
  
  factor=1,
  
  features=20,
  
  add_lm=TRUE
  
)

###############################################################################
# PIPELINE STEP 14: Export MOFA Results
###############################################################################
#Extract latent factors
Factors <- get_factors(
  
  MOFAobject,
  
  as.data.frame=TRUE
  
)

write.csv(
  
  Factors,
  
  "MOFA_Factors.csv",
  
  row.names=FALSE
  
)
#Extract weights
Weights <- get_weights(
  
  MOFAobject,
  
  as.data.frame=TRUE
  
)

write.csv(
  
  Weights,
  
  "MOFA_Weights.csv",
  
  row.names=FALSE
  
)
#Extract original data
Omics <- get_data(
  
  MOFAobject,
  
  as.data.frame=TRUE
  
)

write.csv(
  
  Omics,
  
  "MOFA_Data.csv",
  
  row.names=FALSE
  
)
#Export top genes of every factor
weights_mRNA <-
  
  get_weights(
    
    MOFAobject,
    
    view="mRNA",
    
    as.data.frame=TRUE
    
  )
for(i in 1:15){
  
  tmp <-
    
    subset(
      
      weights_mRNA,
      
      factor==paste0("Factor",i)
      
    )
  
  tmp <-
    
    tmp[
      
      order(abs(tmp$value),
            
            decreasing=TRUE),
      
    ]
  
  write.csv(
    
    head(tmp,100),
    
    paste0("Factor",i,"_Top100Genes.csv"),
    
    row.names=FALSE
    
  )
}
#Extract factors

fac <- get_factors(MOFAobject)[[1]]

#Loop through all factors

library(survival)

cox.results <- data.frame()

for(i in 1:ncol(fac)){
  
  clinical$Factor <- fac[,i]
  
  fit <- coxph(
    
    Surv(time,event)~Factor,
    
    data=clinical
    
  )
  
  s <- summary(fit)
  
  cox.results <-
    
    rbind(
      
      cox.results,
      
      data.frame(
        
        Factor=i,
        
        HR=s$coef[2],
        
        Pvalue=s$coef[5]
        
      )
      
    )
  
}

###############################################################################
# PIPELINE STEP 15: Survival Analysis
###############################################################################

#cox.results

write.csv(
  
  cox.results,
  
  "Cox_results.csv",
  
  row.names=FALSE
  
)
#Kaplan–Meier for significant factors

#Example Factor 1

clinical$Factor1 <- fac[,1]

clinical$Group <-
  
  ifelse(
    
    clinical$Factor1>
      
      median(clinical$Factor1),
    
    "High",
    
    "Low"
    
  )

fit <- survfit(
  
  Surv(time,event)~Group,
  
  data=clinical
  
)

plot(fit,
     col=c("red","blue"),
     lwd=2)

legend(
  
  "bottomleft",
  
  c("High","Low"),
  
  col=c("red","blue"),
  
  lwd=2
  
)

###############################################################################
# PIPELINE STEP 16: Functional Enrichment Analysis
###############################################################################


#Functional enrichment 

#Example for Factor 1:

weights <- get_weights(
  MOFAobject,
  view = "mRNA",
  as.data.frame = TRUE
)

top <- subset(weights, factor == "Factor1")
top <- top[order(abs(top$value), decreasing = TRUE), ]

genes <- sub("\\..*$", "", head(top$feature, 300))

#Convert gene symbols to Entrez IDs:

gene.df <- bitr(
  genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

entrez <- unique(gene.df$ENTREZID)

#GO enrichment:

ego <- enrichGO(
  gene = entrez,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

dotplot(ego)

barplot(ego)

write.csv(
  as.data.frame(ego),
  "Factor1_GO.csv",
  row.names = FALSE
)

#KEGG enrichment:

ekegg <- enrichKEGG(
  gene = entrez,
  organism = "hsa"
)

dotplot(ekegg)
# for entire data
weights <- get_weights(
  MOFAobject,
  view = "mRNA",
  as.data.frame = TRUE
)

GO.results <- list()
KEGG.results <- list()

for(i in 1:15){
  
  top <- subset(weights, factor == paste0("Factor", i))
  top <- top[order(abs(top$value), decreasing = TRUE), ]
  
  genes <- sub("\\..*$", "", head(top$feature, 300))
  
  gene.df <- bitr(
    genes,
    fromType="SYMBOL",
    toType="ENTREZID",
    OrgDb=org.Hs.eg.db
  )
  
  entrez <- unique(gene.df$ENTREZID)
  
  GO.results[[i]] <- enrichGO(
    gene=entrez,
    OrgDb=org.Hs.eg.db,
    ont="BP",
    pAdjustMethod="BH",
    readable=TRUE
  )
  
  KEGG.results[[i]] <- enrichKEGG(
    gene=entrez,
    organism="hsa"
  )
}

dotplot(GO.results[[1]])
dotplot(GO.results[[2]])
dotplot(GO.results[[3]])
dotplot(GO.results[[4]])
dotplot(GO.results[[5]])
dotplot(GO.results[[6]])
dotplot(GO.results[[7]])
dotplot(GO.results[[8]])
dotplot(GO.results[[9]])
dotplot(GO.results[[10]])
dotplot(GO.results[[11]])
dotplot(GO.results[[12]])
dotplot(GO.results[[13]])
dotplot(GO.results[[14]])
dotplot(GO.results[[15]])


#library(enrichplot)
# Compute a similarity matrix first
ego_sim <- pairwise_termsim(GO.results[[1]])
emapplot(ego_sim, showCategory = 20)

library(DOSE)
edo <- enrichDO(
  gene = entrez, # Your converted Entrez IDs
  pvalueCutoff = 0.05,
  readable = TRUE
)
dotplot(edo)

library(MOFA2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggplot2)

dir.create("Enrichment", showWarnings = FALSE)
dir.create("Enrichment/GO", showWarnings = FALSE)
dir.create("Enrichment/KEGG", showWarnings = FALSE)
dir.create("Enrichment/Reactome", showWarnings = FALSE)

dir.create("Enrichment/GO/Dotplots", showWarnings = FALSE)
dir.create("Enrichment/GO/Barplots", showWarnings = FALSE)

dir.create("Enrichment/KEGG/Dotplots", showWarnings = FALSE)
dir.create("Enrichment/KEGG/Barplots", showWarnings = FALSE)

dir.create("Enrichment/Reactome/Dotplots", showWarnings = FALSE)
dir.create("Enrichment/Reactome/Barplots", showWarnings = FALSE)

weights <- get_weights(
  MOFAobject,
  view = "mRNA",
  as.data.frame = TRUE
)

#run enrichement for every factor
weights <- get_weights(
  MOFAobject,
  view = "mRNA",
  as.data.frame = TRUE
)

GO.results <- list()
KEGG.results <- list()
Reactome.results <- list()

for(i in 1:15){
  
  cat("Running Factor",i,"\n")
  
  tmp <- subset(
    weights,
    factor == paste0("Factor",i)
  )
  
  tmp <- tmp[
    order(abs(tmp$value), decreasing=TRUE),
  ]
  
  genes <- sub(
    "\\..*$",
    "",
    head(tmp$feature,300)
  )
  
  gene.df <- bitr(
    genes,
    fromType="SYMBOL",
    toType="ENTREZID",
    OrgDb=org.Hs.eg.db
  )
  
  gene.df <- gene.df[!duplicated(gene.df$ENTREZID),]
  
  entrez <- gene.df$ENTREZID
  
  