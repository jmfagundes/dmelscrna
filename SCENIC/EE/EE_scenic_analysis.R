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

load("../../objs/dmel_filtered.Rdata")

# download fly data

get_fly <- FALSE

if(get_fly) {
  
  dir.create("dbFiles")
  setwd("dbFiles")
  
  dbFiles <- c("https://resources.aertslab.org/cistarget/databases/drosophila_melanogaster/dm6/flybase_r6.02/mc8nr/gene_based/dm6-5kb-upstream-full-tx-11species.mc8nr.feather")
  for(featherURL in dbFiles) {
    download.file(featherURL, destfile = basename(featherURL))
  }
  
  setwd("..")
  
}

# extract cell info from seurat objects

EE_info <- data.frame(CellType = dmel_filtered@active.ident %>% as.character(),
                      DMelNV_percentage = dmel_filtered$percent.nora,
                      TV_percentage = dmel_filtered$percent.thika,
                      InfectionStatus = "uninfected",
                      nGene = dmel_filtered$nFeature_NOVIR)

row.names(EE_info) <- colnames(dmel_filtered)

EE_info[dmel_filtered$is.dmelnv.infected, "InfectionStatus"] <- "DMelNV-infected"
EE_info[dmel_filtered$is.tv.infected, "InfectionStatus"] <- "TV-infected"
EE_info[dmel_filtered$is.dmelnv.infected & dmel_filtered$is.tv.infected, "InfectionStatus"] <- "coinfected"

# assign colors

cell.states <- unique(EE_info$CellType)

colVars <- list(CellType = c("I-a" = "green",
                             "II-a" = "darkorange",
                             "I-ap-a" = "magenta",
                             "I-m" = "hotpink",
                             "II-m1" = "red",
                             "II-m2" = "skyblue",
                             "I-pCCHa1" = "blue",
                             "I-pAstA" = "cyan",
                             "I-ap-p" = "tomato",
                             "II-p" = "yellow",
                             "III" = "khaki",
                             "EEP" = "lightslateblue"),
                
                InfectionStatus = c("DMelNV-infected" = "red",
                                    "TV-infected" = "blue",
                                    "coinfected" = "green"))

dir.create("int")
saveRDS(EE_info, file = "int/cellInfo.Rds")
saveRDS(colVars, file = "int/colVars.Rds")

# count matrix

EE_counts <- dmel_filtered@assays$CELLBENDER@counts %>% as.matrix()

# convert gene entrezID to gene symbol

id_symbol <- read.table("../../name_id_symbol_list",
                        stringsAsFactors = FALSE)

row.names(EE_counts) <- id_symbol[match(gsub("-", "_", row.names(EE_counts)), id_symbol$V1), "V3"]

# remove rows with NA

EE_counts <- EE_counts[!is.na(row.names(EE_counts)),]

# remove seurat objects to clear space

rm(dmel_filtered)

# initialize SCENIC on EE dataset
# must initialize to import getAUC function

init.scenic <- FALSE

if(init.scenic) {
  
  scenicOptions <- initializeScenic(org = "dmel",
                                    dbDir = "dbFiles",
                                    dbs = defaultDbNames[["dmel"]],
                                    datasetTitle = "EE")
  
  scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
  scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
  
} else load("int/scenicOptions_genie3.Rdata")

EE_genesKept <- geneFiltering(EE_counts,
                              scenicOptions = scenicOptions,
                              minCountsPerGene = 3*.01*ncol(EE_counts),
                              minSamples = ncol(EE_counts)*.01)

EE_counts <- EE_counts[EE_genesKept, ]

# only needs to be run once

run.genie3 <- FALSE

if(run.genie3) {
  
  # correlation
  
  runCorrelation(EE_counts, scenicOptions)
  
  # Genie3
  
  runGenie3(EE_counts, scenicOptions)
  
  # save files
  
  save(scenicOptions, file = "int/scenicOptions_genie3.Rdata")
  
  # build and score GRN
  
  runSCENIC_1_coexNetwork2modules(scenicOptions)
  
  scenicOptions@settings$nCores <- 1
  runSCENIC_2_createRegulons(scenicOptions)
  
  runSCENIC_3_scoreCells(scenicOptions, log2(EE_counts + 1))
  
}

# get regulon activity for each cell

EE_regulonAUC.df <- getAUC(readRDS("int/3.4_regulonAUC.Rds")) %>% as.data.frame()

# create a Seurat object and find marker regulons

EE_regulonAUC <- CreateSeuratObject(EE_regulonAUC.df)
Idents(EE_regulonAUC) <- EE_info$CellType
EE_regulonAUC.markers <- FindAllMarkers(EE_regulonAUC, slot = "counts", logfc.threshold = 0, only.pos = TRUE)

EE_regulonAUC.markers %>% group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

# glm to test for infection-responding regulons
# only use uninfected cells if they have 0 viral reads

rownames(EE_regulonAUC.df) <- gsub(" |-", "_", rownames(EE_regulonAUC.df)) %>% gsub("\\(|\\)", "", .)

EE.infection.regulon <- data.frame(regulon = character(),
                                   `DMelNV-infected` = numeric(),
                                   `TV-infected` = numeric(),
                                   coinfected = numeric(),
                                   `DMelNV-infected_p` = numeric(),
                                   `TV-infected_p` = numeric(),
                                   coinfected_p = numeric())

for(regulon in rownames(EE_regulonAUC.df)) {
  
  regulon.df <- data.frame(AUC = EE_regulonAUC.df %>% t() %>% as.data.frame() %>% dplyr::select(regulon),
                           cell = EE_info$CellType,
                           infection = EE_info$InfectionStatus)
  
  colnames(regulon.df)[1] <- "AUC"
  
  #regulon.df$infection <- factor(regulon.df$infection, levels = c("uninfected", "DMelNV-infected", "TV-infected", "coinfected"))
  
  regulon.df$group <-paste0(regulon.df$infection, regulon.df$cell)
  
  #regulon.glm <- glm(AUC ~ infection + cell, data = regulon.df)
  regulon.glm <- glm(AUC ~ 0 + group, data = regulon.df)
  
  tv.infection.main <- glht(regulon.glm, linfct = matrix(c(0,0,0,0,0,
                                                           0,0,0,0,0,0,0,0,0,
                                                           1/11,1/11,1/11,1/11,1/11,1/11,1/11,1/11,1/11,1/11,1/11,
                                                           0,-1/11,-1/11,-1/11,-1/11,-1/11,-1/11,-1/11,-1/11,-1/11,-1/11,-1/11), 1))
  
  dmelnv.infection.main <- glht(regulon.glm, linfct = matrix(c(0,0,0,0,0,
                                                               1/9,1/9,1/9,1/9,1/9,1/9,1/9,1/9,1/9,
                                                               0,0,0,0,0,0,0,0,0,0,0,
                                                               0,0,-1/9,-1/9,-1/9,-1/9,-1/9,-1/9,-1/9,-1/9,-1/9,0), 1))
  
  coinf.infection.main <- glht(regulon.glm, linfct = matrix(c(1/5,1/5,1/5,1/5,1/5,
                                                              0,0,0,0,0,0,0,0,0,
                                                              0,0,0,0,0,0,0,0,0,0,0,
                                                              0,0,0,-1/5,-1/5,-1/5,0,0,-1/5,0,-1/5,0), 1))
  
  EE.infection.regulon[nrow(EE.infection.regulon) + 1,] <- c(regulon,
                                                             summary(tv.infection.main)$test$coefficients,
                                                             summary(dmelnv.infection.main)$test$coefficients,
                                                             summary(coinf.infection.main)$test$coefficients,
                                                             summary(tv.infection.main)$test$pvalues,
                                                             summary(dmelnv.infection.main)$test$pvalues,
                                                             summary(coinf.infection.main)$test$pvalues) %>% as.list()
  
  #EE.infection.regulon[nrow(EE.infection.regulon) + 1,] <- c(regulon,
  #summary(regulon.glm)$coefficients[2:4, 1],
  #summary(regulon.glm)$coefficients[2:4, 4]) %>% as.list()
  
}

EE.infection.regulon$DMelNV.infected_p.adj <- p.adjust(EE.infection.regulon$DMelNV.infected_p, method = "BH")
EE.infection.regulon$TV.infected_p.adj <- p.adjust(EE.infection.regulon$TV.infected_p, method = "BH")
EE.infection.regulon$coinfected_p.adj <- p.adjust(EE.infection.regulon$coinfected_p, method = "BH")

EE.DMelNV.regulons <- EE.infection.regulon %>% filter(DMelNV.infected_p.adj < 0.05) %>% dplyr::select(c("regulon", "DMelNV.infected", "DMelNV.infected_p.adj"))
EE.TV.regulons <- EE.infection.regulon %>% filter(TV.infected_p.adj < 0.05) %>% dplyr::select(c("regulon", "TV.infected", "TV.infected_p.adj"))
EE.coinf.regulons <- EE.infection.regulon %>% filter(coinfected_p.adj < 0.05) %>% dplyr::select(c("regulon", "coinfected", "coinfected_p.adj"))

# EE and MA DMelNV common altered regulons

EE.MA.DMelNV.regulons <- gsub("_[0-9]*g", "", EE.DMelNV.regulons$regulon)[gsub("_[0-9]*g", "", EE.DMelNV.regulons$regulon) %in%
                                                                            gsub("_[0-9]*g", "", MA.DMelNV.regulons$regulon)]

# add hsf AUC to seurat obj

dmel_filtered <- AddMetaData(dmel_filtered, metadata = EE_regulonAUC.df["Hsf (53g)",] %>% t() %>% as.data.frame.matrix(), col.name = "Hsf_53g")

gg.hsf.EE <- FeaturePlot(dmel_filtered, reduction = "tsne", label = TRUE, pt.size = 0.5, repel = TRUE, features = "Hsf_53g") +
  ylab(label = element_blank()) + xlab(label = element_blank()) + ggtitle("B") + theme(plot.title = element_text(hjust = 0))

# save figures

g <- arrangeGrob(gg.hsf.MA, gg.hsf.EE,
                 ncol = 2, left = "tSNE_2", bottom = "tSNE_1")

ggsave("../../figures_out/fig5.pdf", g, width = 174, height = 85, units = "mm")

# save marker regulons

write_xlsx(list(EE_marker_regulons = EE_regulonAUC.markers,
                MA_marker_regulons = MA_regulonAUC.markers), "../../Supplementary 3.xlsx")