# pathway/regulon analysis

library(Seurat)
library(dplyr)
library(tibble)
library(tidyverse)
library(SCENIC)
library(doSNOW)
library(doParallel)
library(doMPI)
library(multcomp)

# load data

load("../../objs/mda_10x_merged.Rdata")

# extract cell info from seurat objects

MA_info <- data.frame(CellType = mda_10x_merged@active.ident %>% as.character(),
                      DMelNV_percentage = mda_10x_merged$percent.nora,
                      InfectionStatus = "uninfected",
                      nGene = mda_10x_merged$nFeature_NOVIR)

row.names(MA_info) <- colnames(mda_10x_merged)

MA_info[mda_10x_merged$is.dmelnv.infected, "InfectionStatus"] <- "DMelNV-infected"

# assign colors

cell.states <- unique(MA_info$CellType)

colVars <- list(CellType = c("pEC" = "green",
                             "mEC" = "darkorange",
                             "LFC" = "magenta",
                             "aEC" = "hotpink",
                             "copper/iron" = "red",
                             "unk" = "skyblue",
                             "cardia" = "blue",
                             "ISC/EB" = "cyan",
                             "EE" = "tomato"),
                
                InfectionStatus = c("DMelNV-infected" = "red"))

dir.create("int")
saveRDS(MA_info, file = "int/cellInfo.Rds")
saveRDS(colVars, file = "int/colVars.Rds")

# count matrix

MA_counts <- mda_10x_merged@assays$CELLBENDER@counts %>% as.matrix()

# convert gene entrezID to gene symbol

id_symbol <- read.table("../../name_id_symbol_list",
                        stringsAsFactors = FALSE)

row.names(MA_counts) <- id_symbol[match(gsub("-", "_", row.names(MA_counts)), id_symbol$V1), "V3"]

# remove rows with NA

MA_counts <- MA_counts[!is.na(row.names(MA_counts)),]

# remove seurat objects to clear space

rm(mda_10x_merged)

# initialize SCENIC on MA dataset

init.scenic <- FALSE

if(init.scenic) {
  
  scenicOptions <- initializeScenic(org = "dmel",
                                    dbDir = "../EE/dbFiles",
                                    dbs = defaultDbNames[["dmel"]],
                                    datasetTitle = "MA")
  
  scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
  scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
  
} else load("int/scenicOptions_genie3.Rdata")

MA_genesKept <- geneFiltering(MA_counts,
                              scenicOptions = scenicOptions,
                              minCountsPerGene = 3*.01*ncol(MA_counts),
                              minSamples = ncol(MA_counts)*.01)

MA_counts <- MA_counts[MA_genesKept, ]

# only needs to be run once

run.genie3 <- FALSE

if(run.genie3) {
  
  # correlation
  
  runCorrelation(MA_counts, scenicOptions)
  
  # Genie3
  
  runGenie3(MA_counts, scenicOptions)
  
  # save files
  
  save(scenicOptions, file = "int/scenicOptions_genie3.Rdata")
  
  # build and score GRN
  
  runSCENIC_1_coexNetwork2modules(scenicOptions)
  
  scenicOptions@settings$nCores <- 1
  runSCENIC_2_createRegulons(scenicOptions)
  
  runSCENIC_3_scoreCells(scenicOptions, log2(MA_counts + 1))
  
}

# get regulon activity for each cell

MA_regulonAUC.df <- getAUC(readRDS("int/3.4_regulonAUC.Rds")) %>% as.data.frame()

# create a Seurat object and find marker regulons

# too many few cells with 0 DMelNV percentage
#MA_regulonAUC <- CreateSeuratObject(EE_regulonAUC.df[colnames(MA_regulonAUC.df) %in% rownames(MA_info[MA_info$DMelNV_percentage == 0,])])
#Idents(MA_regulonAUC) <- MA_info[MA_info$DMelNV_percentage == 0, "CellType"]


MA_regulonAUC <- CreateSeuratObject(MA_regulonAUC.df)
Idents(MA_regulonAUC) <- MA_info$CellType
MA_regulonAUC.markers <- FindAllMarkers(MA_regulonAUC, slot = "counts", logfc.threshold = 0, only.pos = TRUE)

MA_regulonAUC.markers %>% group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

# glm to test for infection-responding regulons
# only use uninfected cells if they have 0 viral reads

rownames(MA_regulonAUC.df) <- gsub(" |-", "_", rownames(MA_regulonAUC.df)) %>% gsub("\\(|\\)", "", .)

MA.infection.regulon <- data.frame(regulon = character(),
                                   `DMelNV-infected` = numeric(),
                                   `DMelNV-infected_p` = numeric())

Hsf.AUC <- list()

for(regulon in rownames(MA_regulonAUC.df)) {
  
  regulon.df <- data.frame(AUC = MA_regulonAUC.df %>% t() %>% as.data.frame() %>% dplyr::select(regulon),
                           cell = MA_info$CellType,
                           infection = MA_info$InfectionStatus)
  
  colnames(regulon.df)[1] <- "AUC"
  
  #regulon.df$infection <- factor(regulon.df$infection, levels = c("uninfected", "DMelNV-infected"))
  
  regulon.df$group <-paste0(regulon.df$infection, regulon.df$cell)
  
  regulon.glm <- glm(AUC ~ 0 + group, data = regulon.df)
  
  infection.main <- glht(regulon.glm, linfct = matrix(c(1/9,1/9,1/9,1/9,1/9,1/9,1/9,1/9,1/9,-1/9,-1/9,-1/9,-1/9,-1/9,-1/9,-1/9,-1/9,-1/9), 1))
  
  MA.infection.regulon[nrow(MA.infection.regulon) + 1,] <- c(regulon,
                                                             summary(infection.main)$test$coefficients,
                                                             summary(infection.main)$test$pvalues) %>% as.list()
  
  if (length(grep("Hsf", regulon)) > 0) Hsf.AUC[[regulon]] <- regulon.glm
}

MA.infection.regulon$DMelNV.infected_p.adj <- p.adjust(MA.infection.regulon$DMelNV.infected_p, method = "BH")

MA.DMelNV.regulons <- MA.infection.regulon %>% filter(DMelNV.infected_p.adj < 0.05) %>% dplyr::select(c("regulon", "DMelNV.infected", "DMelNV.infected_p.adj"))

# Hsf may have higher activity in EEs

MA_regulonAUC.markers %>% filter(gene == "Hsf (53g)")
Hsf.AUC$Hsf_53g %>% summary()

MA_regulonAUC.markers %>% filter(gene == "Hsf-extended (60g)")
Hsf.AUC$Hsf_extended_60g %>% summary()

# Hsf response to DMelNV in EEs

# add hsf AUC to seurat obj

mda_10x_merged <- AddMetaData(mda_10x_merged, metadata = MA_regulonAUC.df["Hsf (53g)",] %>% t() %>% as.data.frame.matrix(), col.name = "Hsf_53g")

gg.hsf.MA <- FeaturePlot(mda_10x_merged, reduction = "tsne", label = TRUE, pt.size = 0.5, repel = TRUE, features = "Hsf_53g") +
  ylab(label = element_blank()) + xlab(label = element_blank()) + ggtitle("A") + theme(plot.title = element_text(hjust = 0))


