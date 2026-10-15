############################################################
## 1. Load packages
############################################################

library(MOFA2)
library(readxl)
library(survival)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(dplyr)

############################################################
## 2. Load trained MOFA model
############################################################

MOFAobject <- load_model("TCGA_MOFA_model.hdf5")
#Load clinical data
clinical <- read_excel("Clinical_data.xlsx")
clinical <- as.data.frame(clinical)

clinical$sample <- gsub("-", ".", clinical$sampleID)

rownames(clinical) <- clinical$sample
clinical$sampleID <- NULL
#Match samples
samples <- samples_names(MOFAobject)[[1]]
samples <- gsub("-", ".", samples)

common <- intersect(samples, rownames(clinical))

length(common)

MOFAobject <- subset_samples(MOFAobject,
                             samples=common)

clinical <- clinical[common, ]
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
#Create survival variables
clinical$event <-
  ifelse(clinical$vital_status=="Dead",1,0)

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
#Visualize latent factors
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
#Association with clinical variables
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
#UMAP
set.seed(123)

MOFAobject <- run_umap(MOFAobject)

plot_dimred(
  
  MOFAobject,
  
  method="UMAP",
  
  color_by="pathologic_stage"
  
)
#Top features
#mRNA
plot_top_weights(
  
  MOFAobject,
  
  view="mRNA",
  
  factor=1,
  
  nfeatures=30
  
)
#methylation
plot_top_weights(
  
  MOFAobject,
  
  view="Methylation",
  
  factor=1,
  
  nfeatures=30
  
)
#miRNA
plot_top_weights(
  
  MOFAobject,
  
  view="miRNA",
  
  factor=1,
  
  nfeatures=30
  
)
#Heatmap of feature weights
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
15. Heatmap of original data
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

#Survival analysis

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