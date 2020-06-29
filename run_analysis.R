#!/usr/bin/Rscript

# all scripts in these folder were meant to be run in an environment (e.g. RStudio) so that the data is analyzed interactivelly
# this script will perform:

# cluster generation
# identification of cell subtypes
# virus replication analyses
# gene expression analysis (partial correlation and anova)
# part of the sensitivity curve analysis

# the data generated in this script must be loaded before execution of make_figs.R, GO.R and graph.R

library(Seurat)
library(dplyr)
library(tibble)
library(viridis)
library(ggplot2)
library(sctransform)
library(reshape2)
library(ppcor)
library(ggpubr)
library(gridExtra)
library(agricolae)
library(GO.db)
library(limma)
library(ggforce)
library(igraph)
library(tidyverse)
library(data.table)
library(car)
library(broom)

###############
## functions ##
###############

# percentage of raw reads derived from TV and DMelNV and their log expression
# a seurat obj must be provided

virus_pct_exp_fun <- function(obj) {
  
  df <- obj@assays$RNA@counts
  dexp <- obj@assays$SCT@data
  df_out <- data.frame(cell = character(0),
                       cluster = character(0),
                       total_counts = numeric(0),
                       thika_counts = numeric(0),
                       nora_counts = numeric(0),
                       thika_pct = numeric(0),
                       nora_pct = numeric(0),
                       thika_exp = numeric(0),
                       nora_exp = numeric(0),
                       stringsAsFactors = FALSE)
  
  for (col in 1:ncol(df)) {
    cell <- colnames(df)[col]
    cluster <- as.character(obj@active.ident[cell])
    counts <- df[,col]
    total_counts <- sum(counts)
    thika_counts <- counts["Thika-virus"]
    nora_counts <- counts["Dmel-nora-virus"]
    thika_pct <- thika_counts/total_counts
    nora_pct <- nora_counts/total_counts
    thika_exp <- dexp[,col]["Thika-virus"]
    nora_exp <- dexp[,col]["Dmel-nora-virus"]
    df_out[col,] <- list(cell, cluster, total_counts, thika_counts, nora_counts, thika_pct, nora_pct,
                         thika_exp, nora_exp)
  }
  return(df_out)
}

# table to convert ids <-> symbol <-> name 

gene_symbol <- read.table("name_id_symbol_list",
                          stringsAsFactors = FALSE)

# changing lrRNA symbol to mt:lrRNA

gene_symbol[gene_symbol$V3 == "lrRNA", 3] <- "mt:lrRNA"

# adding some missing genes

gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel_CR45859", 26067187, "28SrRNA-Psi:CR45859")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel_CR40596", 5740740, "28SrRNA-Psi:CR40596")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel_CR34096", 19893562, "srRNA")

# function to make genes list
# returns genes that are expressed (> min.exp) in at least min.pct cells of each cluster
# a filtered seurat obj must be provided

do_gene_list <- function(min.pct = 0.25,
                         min.exp = 0,
                         obj) {
  
  genes <- rownames(obj@assays$SCT@data)
  genes <- genes[genes != "Thika-virus" & genes != "Dmel-nora-virus" & genes != "Dmel-C-virus"]
  genes_out <- vector()
  total_cells <- table(obj@active.ident)
  for (gene in genes) {
    exp_tb <- data.frame(exp = obj@assays$SCT@data[gene,],
                         cell = obj@active.ident)
    exp_tb <- as.data.frame(table(exp_tb))
    exp_tb$exp <- as.numeric(as.character(exp_tb$exp))
    exp_tb[exp_tb$exp <= min.exp, "Freq"] <- 0
    pcts <- aggregate(exp_tb$Freq, exp_tb["cell"], sum)$`x`/total_cells

    if (length(pcts[pcts >= min.pct]) == length(pcts)) {
      genes_out[length(genes_out) + 1] <- gene
    }
  }
  return(genes_out)
}

# partial correlation functions
# also applying GLM with covariate

exp_cor <- function(obj, x, y) {
  
  x_y_exp <- data.frame(x = obj@assays$SCT@data[x, ],
                        y = obj@assays$SCT@data[y, ],
                        cluster = factor(obj@active.ident))
  x_y_exp$x <- as.numeric(x_y_exp$x)
  x_y_exp$y <- as.numeric(x_y_exp$y)
  
  # make dummy table
  # skip I-ap-p or unk
  # I-ap-p or unk cells will have 0 0 0 ... as dummy holder
  
  subtypes <- unique(x_y_exp$cluster)[!unique(x_y_exp$cluster) %in% list("I-ap-p", "unk")]
  dummy_tb <- data.frame(matrix(ncol = length(subtypes),
                                nrow = nrow(x_y_exp)),
                         stringsAsFactors = FALSE)
  colnames(dummy_tb) <- subtypes
  
  for (subtype in subtypes) {
    
    ones <- as.character(x_y_exp$cluster)
    ones[ones != subtype] <- 0
    ones[ones == subtype] <- 1
    dummy_tb[subtype] <- as.numeric(ones)
  }
  
  x_y_exp <- cbind(x_y_exp, dummy_tb)
  
  # partial correlation
  
  pc <- pcor.test(x_y_exp$x, 
                  x_y_exp$y,
                  dummy_tb)
  
  # GLM with covariate
  
  glm_covariate <- summary(lm(data = x_y_exp, y ~ x + cluster))
  p.value_glm <- glm_covariate$coefficients[2, 4]
  
  return(list(x_y_exp, p.value_glm, pc))
}

# TV and DMelNV partial correlations

do_correlations <- function(dataset = "EE") {
  
  if (dataset == "EE") {
    
    t <- data.frame(estimate = numeric(0),
                    p.value = numeric(0),
                    statistic = numeric(0),
                    n = numeric(0),
                    gp = numeric(0),
                    Method = character(0),
                    GLM_p.value = numeric(0),
                    stringsAsFactors = FALSE)
    
    for (gene in genes) {
      p.cor <- exp_cor(only_thika, "Thika-virus", gene)
      tc <- p.cor[[3]]
      tc[length(tc) + 1] <- p.cor[[2]]
      rownames(tc) <- gene
      t[nrow(t) + 1,] <- tc
    }
    
    n <- data.frame(estimate = numeric(0),
                    p.value = numeric(0),
                    statistic = numeric(0),
                    n = numeric(0),
                    gp = numeric(0),
                    Method = character(0),
                    GLM_p.value = numeric(0),
                    stringsAsFactors = FALSE)
    
    for (gene in genes) {
      p.cor <- exp_cor(only_nora, "Dmel-nora-virus", gene)
      nc <- p.cor[[3]]
      nc[length(nc) + 1] <- p.cor[[2]]
      rownames(nc) <- gene
      n[nrow(n) + 1,] <- nc
    }
    return(list(t, n))
    
    # only analyze DMelNV if analyzing the midgut dataset
    
  } else if (dataset == "midgut") {
    
    n <- data.frame(estimate = numeric(0),
                    p.value = numeric(0),
                    statistic = numeric(0),
                    n = numeric(0),
                    gp = numeric(0),
                    Method = character(0),
                    GLM_p.value = numeric(0),
                    stringsAsFactors = FALSE)
    
    for (gene in mda_genes) {
      
      p.cor <- exp_cor(mda_10x_nora, "Dmel-nora-virus", gene)
      nc <- p.cor[[3]]
      nc[length(nc) + 1] <- p.cor[[2]]
      rownames(nc) <- gene
      n[nrow(n) + 1,] <- nc
    }
    return(n)
  }
}

# two-way anova of every gene
# aov(exp ~ cell + infection + cell:infection)

do_aov <- function(dataset = "EE",
                   obj = NULL) {
  
  if (dataset == "EE") {
    
    thika_aov <- data.frame(gene = character(0),
                            group = character(0),
                            Df = numeric(0),
                            `Sum Sq` = numeric(0),
                            `Mean Sq` = numeric(0),
                            `F value` = numeric(0),
                            `Pr(>F)` = numeric(0),
                            stringsAsFactors = FALSE)
    
    nora_aov <- data.frame(gene = character(0),
                           group = character(0),
                           Df = numeric(0),
                           `Sum Sq` = numeric(0),
                           `Mean Sq` = numeric(0),
                           `F value` = numeric(0),
                           `Pr(>F)` = numeric(0),
                           stringsAsFactors = FALSE)
    
    for (gene in genes) {
      
      # log expression in uninfected cells
      
      gene_exp_uninfected <- data.frame(exp = dmel_no_vir@assays$SCT@data[gene,],
                                        cell = dmel_no_vir@active.ident)
      gene_exp_uninfected$infection <- replicate(nrow(gene_exp_uninfected), "uninfected")
      
      # log expression in TV-infected cells
      
      gene_exp_infected <- data.frame(exp = only_thika@assays$SCT@data[gene,],
                                      cell = only_thika@active.ident)
      gene_exp_infected$infection <- replicate(nrow(gene_exp_infected), "infected")
      
      exp_tb <- rbind(gene_exp_infected, gene_exp_uninfected)
      aov_res <- anova(lm(
        exp ~ cell + infection + cell:infection, exp_tb
      ))
      tmp_tb <- data.frame(gene = replicate(nrow(aov_res), gene),
                           group = rownames(aov_res),
                           Df = aov_res$`Df`,
                           `Sum Sq` = aov_res$`Df`,
                           `Mean Sq` = aov_res$`Mean Sq`,
                           `F value` = aov_res$`F value`,
                           `Pr(>F)` = aov_res$`Pr(>F)`,
                           stringsAsFactors = FALSE)
      thika_aov <- rbind(thika_aov, tmp_tb)
      
      # log expression in DMelNV-infected cells
      
      gene_exp_infected <- data.frame(exp = only_nora@assays$SCT@data[gene,],
                                      cell = only_nora@active.ident)
      gene_exp_infected$infection <- replicate(nrow(gene_exp_infected), "infected")
      
      exp_tb <- rbind(gene_exp_infected, gene_exp_uninfected)
      aov_res <- anova(lm(
        exp ~ cell + infection + cell:infection, exp_tb
      ))
      tmp_tb <- data.frame(gene = replicate(nrow(aov_res), gene),
                           group = rownames(aov_res),
                           Df = aov_res$`Df`,
                           `Sum Sq` = aov_res$`Df`,
                           `Mean Sq` = aov_res$`Mean Sq`,
                           `F value` = aov_res$`F value`,
                           `Pr(>F)` = aov_res$`Pr(>F)`,
                           stringsAsFactors = FALSE)
      nora_aov <- rbind(nora_aov, tmp_tb)
    }
    return(list(thika_aov, nora_aov))
    
    # only analyze DMelNV if analyzing the midgut dataset
    
  } else if (dataset == "midgut") {
    
    if (is.null(obj)) {
      obj <- mda_10x_nora
    }
    
    nora_aov <- data.frame(gene = character(0),
                           group = character(0),
                           Df = numeric(0),
                           `Sum Sq` = numeric(0),
                           `Mean Sq` = numeric(0),
                           `F value` = numeric(0),
                           `Pr(>F)` = numeric(0),
                           stringsAsFactors = FALSE)
    
    for (gene in mda_genes) {
      
      # log expression in uninfected cells
      
      gene_exp_uninfected <- data.frame(exp = mda_10x_uninfected@assays$SCT@data[gene,],
                                        cell = mda_10x_uninfected@active.ident)
      gene_exp_uninfected$infection <- replicate(nrow(gene_exp_uninfected), "uninfected")
      
      # log expression in DMelNV-infected cells
      
      gene_exp_infected <- data.frame(exp = obj@assays$SCT@data[gene,],
                                      cell = obj@active.ident)
      gene_exp_infected$infection <- replicate(nrow(gene_exp_infected), "infected")
      
      exp_tb <- rbind(gene_exp_infected, gene_exp_uninfected)
      aov_res <- anova(lm(
        exp ~ cell + infection + cell:infection, exp_tb
      ))
      tmp_tb <- data.frame(gene = replicate(nrow(aov_res), gene),
                           group = rownames(aov_res),
                           Df = aov_res$`Df`,
                           `Sum Sq` = aov_res$`Df`,
                           `Mean Sq` = aov_res$`Mean Sq`,
                           `F value` = aov_res$`F value`,
                           `Pr(>F)` = aov_res$`Pr(>F)`,
                           stringsAsFactors = FALSE)
      nora_aov <- rbind(nora_aov, tmp_tb)
    }
    return(nora_aov)
  }
}

# function to calculate whether DEGs are up or downregulated
# for each cluster, the mean log fold change is first calculated
# then, they are used to calculate the global mean log fold change
# calculating mean expression out of log space

up_down_fun <- function(obj_vir,
                        gene_lst,
                        obj_no_vir,
                        cluster_names.list) {
  
  out_tb <- data.frame(gene = character(0),
                       symbol = character(0),
                       mean_logFC = character(0),
                       stringsAsFactors = FALSE)
  
  for (i in 1:length(gene_lst)) {
    
    gene <- gsub("_", "-", gene_lst[i])
    mean_expressions <- data.frame(infected = numeric(0),
                                   uninfected = numeric(0),
                                   logFC = numeric(0),
                                   cell = character(0),
                                   stringsAsFactors = FALSE)
    
    for (i in 1:length(cluster_names.list)) {
      
      infected_exp <- obj_vir@assays$SCT@data[gene, obj_vir@active.ident == cluster_names.list[i]]
      uninfected_exp <- obj_no_vir@assays$SCT@data[gene, obj_no_vir@active.ident == cluster_names.list[i]]
      
      if(length(infected_exp) == 0) {
        infected_exp <- NA
      }
      infected_exp <- log1p(mean(expm1(infected_exp)))
      uninfected_exp <- log1p(mean(expm1(uninfected_exp)))
      logFC <- infected_exp - uninfected_exp
      
      mean_expressions[nrow(mean_expressions) + 1,] <- list(infected_exp,
                                                            uninfected_exp,
                                                            logFC,
                                                            cluster_names.list[i])
    }
    
    mean_expressions <- mean_expressions[complete.cases(mean_expressions[1]),]
    mean_logFC <- mean(mean_expressions$logFC)
    symbol <- gene_symbol[gene_symbol$V1 == gsub("-", "_", gene), 3]
    symbol <- as.character(symbol)
    
    out_tb[nrow(out_tb) + 1,] <- list(gsub("-", "_", gene),
                                      symbol,
                                      mean_logFC)
  }
  return(out_tb)
}

##################
## EEs analysis ##
##################

# load raw matrix from the original study for comparison

GSE_raw.data <- Read10X(data.dir = "gene_matrices/GSE132274_RAW/")
GSE_raw <- CreateSeuratObject(counts = GSE_raw.data, project = "GSE_raw", min.cells = 3, min.features = 200)

# load EE data

dmel.data <- Read10X(data.dir = "gene_matrices/EE_cellranger_out/")
dmel <- CreateSeuratObject(counts = dmel.data, project = "EE", min.cells = 3, min.features = 200)

# count mitochondrial and virus genes

dmel[["percent.mito"]] <- PercentageFeatureSet(dmel, pattern = "Dmel-CG34063|Dmel-CG34067|Dmel-CG34069|Dmel-CG34072|Dmel-CG34073|Dmel-CG34074|Dmel-CG34076|Dmel-CG34083|Dmel-CG34085|Dmel-CG34086|Dmel-CG34089|Dmel-CG34090|Dmel-CG34092|Dmel-CR34060|Dmel-CR34061|Dmel-CR34062|Dmel-CR34064|Dmel-CR34065|Dmel-CR34066|Dmel-CR34068|Dmel-CR34070|Dmel-CR34071|Dmel-CR34075|Dmel-CR34077|Dmel-CR34078|Dmel-CR34079|Dmel-CR34080|Dmel-CR34081|Dmel-CR34082|Dmel-CR34084|Dmel-CR34087|Dmel-CR34088|Dmel-CR34091|Dmel-CR34093|Dmel-CR34094|Dmel-CR34095|Dmel-CR34096")
dmel[["percent.vir"]] <- PercentageFeatureSet(dmel, pattern = "Dmel-nora-virus|Dmel-C-virus|Thika-virus")
dmel[["percent.nora"]] <- PercentageFeatureSet(dmel, pattern = "Dmel-nora-virus")
dmel[["percent.dmelc"]] <- PercentageFeatureSet(dmel, pattern = "Dmel-C-virus")
dmel[["percent.thika"]] <- PercentageFeatureSet(dmel, pattern = "Thika-virus")

# filter, scale and normalize
# mitochondrial and viruses genes are regressed out

dmel_filtered <- subset(dmel, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mito < 5)
dmel_filtered <- SCTransform(dmel_filtered, vars.to.regress = c("percent.mito", "percent.vir"), return.only.var.genes = FALSE)

# principal component analysis

dmel_filtered <- RunPCA(dmel_filtered, features = VariableFeatures(object = dmel_filtered))

# cluster generation

dmel_filtered <- FindNeighbors(dmel_filtered, dims = 1:20)
dmel_filtered <- FindClusters(dmel_filtered, resolution = 0.4)

# non-linear dimensional reduction
# subset cells based on https://doi.org/10.1016/j.celrep.2019.11.048
# cluster 9 was originaly annotated as part of I-ap-a + I-ap-p, but since CCHa2 is not highly expressed in these cells,
# here it was renamed as I-ap-a (CCHa2-), and I-ap-a + I-ap-p was reannotated as I-ap-p (CCHa2+)
# merge clusters 1 + 8

dmel_filtered <- RunTSNE(dmel_filtered, dims = 1:20)

# cluster               0       1       2         3       4         5       6           7         8         9      10     11     12
new.cluster.ids <- c("I-m", "II-m1", "I-ap-p", "II-a", "II-m2", "I-pAstA", "II-p", "I-pCCHa1", "II-m1", "I-ap-a", "I-a", "III", "EEP")
names(new.cluster.ids) <- levels(dmel_filtered)
dmel_filtered <- RenameIdents(dmel_filtered, new.cluster.ids)

cluster_names = list("I-m", "II-m1", "I-ap-p", "II-a", "II-m2", "I-pCCHa1", "I-pAstA", "II-p", "I-a", "III", "EEP", "I-ap-a")

cell_order <- c("EEP", # EE progenitor cell
                "I-a", "III", "I-ap-a", "II-a", # anterior
                "I-m", "II-m1",  "II-m2", # medium
                "II-p", "I-ap-p", "I-pCCHa1", "I-pAstA") # posterior

# subset cells infected with viruses

dmel_vir <- subset(dmel_filtered, subset = `Thika-virus` > 1 | `Dmel-nora-virus` > 1 | `Dmel-C-virus` > 1, slot = "counts")

# subset cells uninfected with viruses

dmel_no_vir <- subset(dmel_filtered, subset = `Thika-virus` <= 1 & `Dmel-nora-virus` <= 1 & `Dmel-C-virus` <= 1, slot = "counts")

# subset cells infected with nora virus

nora <- subset(dmel_filtered, subset = `Dmel-nora-virus` > 1, slot = "counts")

# subset cells infected ONLY with nora

only_nora <- subset(dmel_filtered, subset = `Dmel-nora-virus` > 1 & `Thika-virus` <= 1 & `Dmel-C-virus` <= 1, slot = "counts")

# subset cells infected with dmel c virus

dmelc <- subset(dmel_filtered, subset = `Dmel-C-virus` > 1, slot = "counts")

## subset cells infected with thika virus

thika <- subset(dmel_filtered, subset = `Thika-virus` > 1, slot = "counts")

# subset cells infected ONLY with thika

only_thika <- subset(dmel_filtered, subset = `Thika-virus` > 1 & `Dmel-nora-virus` <= 1 & `Dmel-C-virus` <= 1, slot = "counts")

# subset cells infected with nora AND thika virus

nora_thika <- subset(dmel_filtered, subset = `Dmel-nora-virus` > 1 & `Thika-virus` > 1 & `Dmel-C-virus` <= 1, slot = "counts")

# calculate percentage and log expression of TV and DMelNV for each cell

raw_reads_virus_pct <- virus_pct_exp_fun(dmel_filtered)
raw_thika_pct <- sum(raw_reads_virus_pct$thika_counts)/sum(raw_reads_virus_pct$total_counts)
raw_nora_pct <- sum(raw_reads_virus_pct$nora_counts)/sum(raw_reads_virus_pct$total_counts)

thika_pct <- virus_pct_exp_fun(only_thika)
thika_pct$cluster <- factor(thika_pct$cluster, levels = cell_order)

nora_pct <- virus_pct_exp_fun(only_nora)
nora_pct$cluster <- factor(nora_pct$cluster, levels = cell_order)

# one-way anova virus percentage and log expression between cell types

thika_cluster_anova <- list(anova(lm(thika_exp ~ cluster, thika_pct)),
                            anova(lm(thika_pct ~ cluster, thika_pct)))

nora_cluster_anova <- list(anova(lm(nora_exp ~ cluster, nora_pct)),
                           anova(lm(nora_pct ~ cluster, nora_pct)))

# post hoc tukey-HSD test

thika_tukey <- list()
thika_tukey[1] <- TukeyHSD(aov(lm(thika_pct$thika_exp ~ thika_pct$cluster)), 'thika_pct$cluster', conf.level = 0.95)
thika_tukey[2:6] <- HSD.test(lm(thika_pct$thika_exp ~ thika_pct$cluster), 'thika_pct$cluster')

nora_tukey <- list()
nora_tukey[1] <- TukeyHSD(aov(lm(nora_pct$nora_exp ~ nora_pct$cluster)), 'nora_pct$cluster', conf.level = 0.95)
nora_tukey[2:6] <- HSD.test(lm(nora_pct$nora_exp ~ nora_pct$cluster), 'nora_pct$cluster')

# proportion of TV- and DMelNV-infected cells

prc <- data.frame(cell = c("I-m", "II-m1", "I-ap-p", "II-a", "II-m2", "I-pCCHa1", "I-pAstA", "II-p", "I-a", "III", "EEP", "I-ap-a"),
                  thika = numeric(12), nora = numeric(12), coinfected = numeric(12), uninfected = numeric(12),
                  stringsAsFactors = FALSE)

for (i in 1:12) {
  t <- try(length(WhichCells(only_thika, idents = cluster_names[[i]])), silent = TRUE)
  n <- try(length(WhichCells(only_nora, idents = cluster_names[[i]])), silent = TRUE)
  co <- try(length(WhichCells(nora_thika, idents = cluster_names[[i]])), silent = TRUE)
  if("try-error" %in% class(t)) {
    t <- 0
  }
  if("try-error" %in% class(n)) {
    n <- 0
  }
  if("try-error" %in% class(co)) {
    co <- 0
  }
  prc[i, 2] <- t
  prc[i, 3] <- n
  prc[i, 4] <- co
  prc[i, 5] <- length(WhichCells(dmel_no_vir, idents = cluster_names[[i]]))
}

# calculate probability of coinfection

prc_co <- prc
prc_co[nrow(prc_co) + 1,] <- list("total", as.numeric(sum(prc_co$thika)), as.numeric(sum(prc_co$nora)), 
                                  as.numeric(sum(prc_co$coinfected)), as.numeric(sum(prc_co$uninfected)))
p_co <- prc_co
p_co$thika <- (prc_co$thika + prc_co$coinfected)/(prc_co$thika + prc_co$nora + prc_co$coinfected + prc_co$uninfected) 
p_co$nora <- (prc_co$nora + prc_co$coinfected)/(prc_co$thika + prc_co$nora + prc_co$coinfected + prc_co$uninfected) 
p_co$coinfected <- (prc_co$coinfected)/(prc_co$thika + prc_co$nora + prc_co$coinfected + prc_co$uninfected) 
p_co$expected_coinfected <- p_co$thika * p_co$nora
p_co <- p_co[c(1:4,6)]
p_co$"obs/exp" <- p_co$coinfected/p_co$expected_coinfected

# make gene list
# this genes will be used for partial correlation and two-way anova tests

genes <- do_gene_list(min.pct = 0.1,
                      obj = dmel_filtered)

# TV and DmelNV partial correlations

virus_pcorrelations <- do_correlations()
thika_pcor <- virus_pcorrelations[[1]][complete.cases(virus_pcorrelations[[1]]),]
nora_pcor <- virus_pcorrelations[[2]][complete.cases(virus_pcorrelations[[2]]),]
rm(virus_pcorrelations)

# bonferroni correction

thika_pcor$p.adj <- p.adjust(thika_pcor$p.value, "bonferroni")
nora_pcor$p.adj <- p.adjust(nora_pcor$p.value, "bonferroni")

# write files for bingo analysis

thika_pcor_names <- gsub("-", "_", rownames(thika_pcor[thika_pcor$p.adj < 0.05,]))
thika_pcor_names_pos <- gsub("-", "_", rownames(thika_pcor[thika_pcor$p.adj < 0.05 & thika_pcor$estimate > 0,]))
thika_pcor_names_neg <- gsub("-", "_", rownames(thika_pcor[thika_pcor$p.adj < 0.05 & thika_pcor$estimate < 0,]))

nora_pcor_names <- gsub("-", "_", rownames(nora_pcor[nora_pcor$p.adj < 0.05,]))
nora_pcor_names_pos <- gsub("-", "_", rownames(nora_pcor[nora_pcor$p.adj < 0.05 & nora_pcor$estimate > 0,]))
nora_pcor_names_neg <- gsub("-", "_", rownames(nora_pcor[nora_pcor$p.adj < 0.05 & nora_pcor$estimate < 0,]))

thika_pcor_symbol <- gene_symbol[gene_symbol$V1 %in% thika_pcor_names, 3]
thika_pcor_symbol_pos <- gene_symbol[gene_symbol$V1 %in% thika_pcor_names_pos, 3]
thika_pcor_symbol_neg <- gene_symbol[gene_symbol$V1 %in% thika_pcor_names_neg, 3]

nora_pcor_symbol <- gene_symbol[gene_symbol$V1 %in% nora_pcor_names, 3]
nora_pcor_symbol_pos <- gene_symbol[gene_symbol$V1 %in% nora_pcor_names_pos, 3]
nora_pcor_symbol_neg <- gene_symbol[gene_symbol$V1 %in% nora_pcor_names_neg, 3]

write.table(thika_pcor_symbol, file = "GO/partial_correlation/thika_pcor_symbol.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(thika_pcor_symbol_pos, file = "GO/partial_correlation/thika_pcor_pos_symbol.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(thika_pcor_symbol_neg, file = "GO/partial_correlation/thika_pcor_neg_symbol.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(nora_pcor_symbol, file = "GO/partial_correlation/nora_pcor_symbol.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(nora_pcor_symbol_pos, file = "GO/partial_correlation/nora_pcor_pos_symbol.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(nora_pcor_symbol_neg, file = "GO/partial_correlation/nora_pcor_neg_symbol.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# DEGs analysis via two-way anova

infection_aov <- do_aov()
thika_aov <- as.data.frame(infection_aov[1])[complete.cases(as.data.frame(infection_aov[1])),]
nora_aov <- as.data.frame(infection_aov[2])[complete.cases(as.data.frame(infection_aov[2])),]

# bonferroni correction

thika_aov$p.adj <- numeric(nrow(thika_aov))
thika_aov[thika_aov$group == "cell",]$p.adj <- p.adjust(thika_aov[thika_aov$group == "cell",]$Pr..F., "bonferroni" )
thika_aov[thika_aov$group == "infection",]$p.adj <- p.adjust(thika_aov[thika_aov$group == "infection",]$Pr..F., "bonferroni" )
thika_aov[thika_aov$group == "cell:infection",]$p.adj <- p.adjust(thika_aov[thika_aov$group == "cell:infection",]$Pr..F., "bonferroni" )

nora_aov$p.adj <- numeric(nrow(nora_aov))
nora_aov[nora_aov$group == "cell",]$p.adj <- p.adjust(nora_aov[nora_aov$group == "cell",]$Pr..F., "bonferroni" )
nora_aov[nora_aov$group == "infection",]$p.adj <- p.adjust(nora_aov[nora_aov$group == "infection",]$Pr..F., "bonferroni" )
nora_aov[nora_aov$group == "cell:infection",]$p.adj <- p.adjust(nora_aov[nora_aov$group == "cell:infection",]$Pr..F., "bonferroni" )

# make lists of DEGs

thika_infection_genes <- gsub("-", "_", 
                              thika_aov[thika_aov$group == "infection" & thika_aov$p.adj < 0.05,]$gene)
nora_infection_genes <- gsub("-", "_",
                             nora_aov[nora_aov$group == "infection" & nora_aov$p.adj < 0.05,]$gene)

thika_cell.infection_genes <- gsub("-", "_",
                                   thika_aov[thika_aov$group == "cell:infection" & thika_aov$p.adj < 0.05,]$gene)
nora_cell.infection_genes <- gsub("-", "_",
                                  nora_aov[nora_aov$group == "cell:infection" & nora_aov$p.adj < 0.05,]$gene)

thika_exclusively_infection_genes <- thika_infection_genes[! thika_infection_genes %in% thika_cell.infection_genes]
nora_exclusively_infection_genes <- nora_infection_genes[! nora_infection_genes %in% nora_cell.infection_genes]

thika_exclusively_cell.infection_genes <- thika_cell.infection_genes[! thika_cell.infection_genes %in% thika_infection_genes]
nora_exclusively_cell.infection_genes <- nora_cell.infection_genes[! nora_cell.infection_genes %in% nora_infection_genes]

# make list of top 12 DEGs for supplementary figures

top12_thika_infection_genes <- thika_aov[thika_aov$group == "infection" & thika_aov$p.adj < 0.05,] %>% top_n(n = -12, wt = p.adj)
top12_thika_infection_genes <- gsub("-", "_", top12_thika_infection_genes$gene)

top12_nora_infection_genes <- nora_aov[nora_aov$group == "infection" & nora_aov$p.adj < 0.05,] %>% top_n(n = -12, wt = p.adj)
top12_nora_infection_genes <- gsub("-", "_", top12_nora_infection_genes$gene)

top12_thika_cell.infection_genes <- thika_aov[thika_aov$group == "cell:infection" & thika_aov$p.adj < 0.05,] %>% top_n(n = -12, wt = p.adj)
top12_thika_cell.infection_genes <- gsub("-", "_", top12_thika_cell.infection_genes$gene)

top12_nora_cell.infection_genes <- nora_aov[nora_aov$group == "cell:infection" & nora_aov$p.adj < 0.05,] %>% top_n(n = -12, wt = p.adj)
top12_nora_cell.infection_genes <- gsub("-", "_", top12_nora_cell.infection_genes$gene)

# calculate global mean log fold change

# TV infection

mean_logFC_thika_infection <- up_down_fun(obj_vir = only_thika,
                                          gene_lst = thika_infection_genes,
                                          obj_no_vir = dmel_no_vir,
                                          cluster_names.list = cluster_names)

# DMelNV infection

mean_logFC_nora_infection <- up_down_fun(obj_vir = only_nora,
                                         gene_lst = nora_infection_genes,
                                         obj_no_vir = dmel_no_vir,
                                         cluster_names.list = cluster_names)

# write files for bingo analysis

# TV

thika_exclusively_cell.infection_symbol <- gene_symbol[gene_symbol$V1 %in% thika_exclusively_cell.infection_genes, 3]
thika_infection_symbol_pos <- mean_logFC_thika_infection[mean_logFC_thika_infection$mean_logFC > 0, ]$symbol
thika_infection_symbol_neg <- mean_logFC_thika_infection[mean_logFC_thika_infection$mean_logFC < 0, ]$symbol

write.table(thika_exclusively_cell.infection_symbol,
            file = "GO/anova/thika_exclusively_cell.infection_genes_symbol", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(thika_infection_symbol_pos,
            file = "GO/anova/thika_infection_genes_symbol_pos", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(thika_infection_symbol_neg,
            file = "GO/anova/thika_infection_genes_symbol_neg", row.names = FALSE, col.names = FALSE, quote = FALSE)

# DMelNV

nora_exclusively_cell.infection_symbol <- gene_symbol[gene_symbol$V1 %in% nora_exclusively_cell.infection_genes, 3]
nora_infection_symbol_pos <- mean_logFC_nora_infection[mean_logFC_nora_infection$mean_logFC > 0, ]$symbol
nora_infection_symbol_neg <- mean_logFC_nora_infection[mean_logFC_nora_infection$mean_logFC < 0, ]$symbol

write.table(nora_exclusively_cell.infection_symbol,
            file = "GO/anova/nora_exclusively_cell.infection_genes_symbol", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(nora_infection_symbol_pos,
            file = "GO/anova/nora_infection_genes_symbol_pos", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(nora_infection_symbol_neg,
            file = "GO/anova/nora_infection_genes_symbol_neg", row.names = FALSE, col.names = FALSE, quote = FALSE)

#######################
## bulk RNA-seq data ##
#######################

# common cellular and systemic (bulk RNA-seq) response genes to DMelNV

day2 <- fread("gene_matrices/dmelnv_bulk_RNA-seq_data/NV_infected_day2_vs_uninfected.tsv",
              sep = "\t", dec = ",", stringsAsFactors = FALSE)
day2$day <- "2"
day10 <- fread("gene_matrices/dmelnv_bulk_RNA-seq_data/NV_infected_day10_vs_uninfected.tsv",
               sep = "\t", dec = ",", stringsAsFactors = FALSE)
day10$day <- "10"
day20 <- fread("gene_matrices/dmelnv_bulk_RNA-seq_data/NV_infected_day20_vs_uninfected.tsv",
               sep = "\t", dec = ",", stringsAsFactors = FALSE)
day20$day <- "20"
day30 <- fread("gene_matrices/dmelnv_bulk_RNA-seq_data/NV_infected_day30_vs_uninfected.tsv",
               sep = "\t", dec = ",", stringsAsFactors = FALSE)
day30$day <- "30"

bulk_nora <- rbind(day2,
                   day10,
                   day20,
                   day30)

bulk_nora_pos <- rbind(day2[day2$`log2(fold_change)` > 0, ],
                       day10[day10$`log2(fold_change)` > 0, ],
                       day20[day20$`log2(fold_change)` > 0, ],
                       day30[day30$`log2(fold_change)` > 0, ])

bulk_nora_neg <- rbind(day2[day2$`log2(fold_change)` < 0, ],
                       day10[day10$`log2(fold_change)` < 0, ],
                       day20[day20$`log2(fold_change)` < 0, ],
                       day30[day30$`log2(fold_change)` < 0, ])

bulk_nora_neg_common <- bulk_nora_neg[bulk_nora_neg$gene_name %in% mean_logFC_nora_infection[mean_logFC_nora_infection$mean_logFC < 0, ]$symbol, ]
common_tmp <- mean_logFC_nora_infection[mean_logFC_nora_infection$symbol %in% bulk_nora_neg_common$gene_name, ]
common_tmp <- common_tmp[match(bulk_nora_neg_common$gene_name, common_tmp$symbol), ]
bulk_nora_neg_common$`cellular_mean_log_fold_change` <- common_tmp$mean_logFC

bulk_nora_pos_common <- bulk_nora_pos[bulk_nora_pos$gene_name %in% mean_logFC_nora_infection[mean_logFC_nora_infection$mean_logFC > 0, ]$symbol, ]
common_tmp <- mean_logFC_nora_infection[mean_logFC_nora_infection$symbol %in% bulk_nora_pos_common$gene_name, ]
common_tmp <- common_tmp[match(bulk_nora_pos_common$gene_name, common_tmp$symbol), ]
bulk_nora_pos_common$`cellular_mean_log_fold_change` <- common_tmp$mean_logFC

bulk_nora_neg_cell_pos <- bulk_nora_neg[bulk_nora_neg$gene_name %in% mean_logFC_nora_infection[mean_logFC_nora_infection$mean_logFC > 0, ]$symbol, ]
common_tmp <- mean_logFC_nora_infection[mean_logFC_nora_infection$symbol %in% bulk_nora_neg_cell_pos$gene_name, ]
common_tmp <- common_tmp[match(bulk_nora_neg_cell_pos$gene_name, common_tmp$symbol), ]
bulk_nora_neg_cell_pos$`cellular_mean_log_fold_change` <- common_tmp$mean_logFC

bulk_nora_pos_cell_neg <- bulk_nora_pos[bulk_nora_pos$gene_name %in% mean_logFC_nora_infection[mean_logFC_nora_infection$mean_logFC < 0, ]$symbol, ]
common_tmp <- mean_logFC_nora_infection[mean_logFC_nora_infection$symbol %in% bulk_nora_pos_cell_neg$gene_name, ]
common_tmp <- common_tmp[match(bulk_nora_pos_cell_neg$gene_name, common_tmp$symbol), ]
bulk_nora_pos_cell_neg$`cellular_mean_log_fold_change` <- common_tmp$mean_logFC

nora_common <- rbind(bulk_nora_neg_common,
                     bulk_nora_pos_common,
                     bulk_nora_neg_cell_pos,
                     bulk_nora_pos_cell_neg)

###########################
## midgut atlas analysis ##
###########################

# InDrop data has only R1 available

# load expression matrices

# 10x

mda_s1_r1_raw.data <- Read10X_h5("gene_matrices/midgut_atlas/10x_S1_rep1_filtered_feature_bc_matrix.h5")
mda_s1_r1_raw <- CreateSeuratObject(counts = mda_s1_r1_raw.data, project = "mda_s1_r1", min.cells = 3, min.features = 200)

mda_s1_r2_raw.data <- Read10X_h5("gene_matrices/midgut_atlas/10x_S1_rep2_filtered_feature_bc_matrix.h5")
mda_s1_r2_raw <- CreateSeuratObject(counts = mda_s1_r2_raw.data, project = "mda_s1_r2", min.cells = 3, min.features = 200)

mda_s2_r1_raw.data <- Read10X_h5("gene_matrices/midgut_atlas/10x_S2_rep1_filtered_feature_bc_matrix.h5")
mda_s2_r1_raw <- CreateSeuratObject(counts = mda_s2_r1_raw.data, project = "mda_s2_r1", min.cells = 3, min.features = 200)

mda_s2_r2_raw.data <- Read10X_h5("gene_matrices/midgut_atlas/10x_S2_rep2_filtered_feature_bc_matrix.h5")
mda_s2_r2_raw <- CreateSeuratObject(counts = mda_s2_r2_raw.data, project = "mda_s2_r2", min.cells = 3, min.features = 200)

# count mitochondrial and virus genes

mda_s1_r1_raw[["percent.mito"]] <- PercentageFeatureSet(mda_s1_r1_raw, pattern = "Dmel-CG34063|Dmel-CG34067|Dmel-CG34069|Dmel-CG34072|Dmel-CG34073|Dmel-CG34074|Dmel-CG34076|Dmel-CG34083|Dmel-CG34085|Dmel-CG34086|Dmel-CG34089|Dmel-CG34090|Dmel-CG34092|Dmel-CR34060|Dmel-CR34061|Dmel-CR34062|Dmel-CR34064|Dmel-CR34065|Dmel-CR34066|Dmel-CR34068|Dmel-CR34070|Dmel-CR34071|Dmel-CR34075|Dmel-CR34077|Dmel-CR34078|Dmel-CR34079|Dmel-CR34080|Dmel-CR34081|Dmel-CR34082|Dmel-CR34084|Dmel-CR34087|Dmel-CR34088|Dmel-CR34091|Dmel-CR34093|Dmel-CR34094|Dmel-CR34095|Dmel-CR34096")
mda_s1_r1_raw[["percent.nora"]] <- PercentageFeatureSet(mda_s1_r1_raw, pattern = "Dmel-nora-virus")

mda_s1_r2_raw[["percent.mito"]] <- PercentageFeatureSet(mda_s1_r2_raw, pattern = "Dmel-CG34063|Dmel-CG34067|Dmel-CG34069|Dmel-CG34072|Dmel-CG34073|Dmel-CG34074|Dmel-CG34076|Dmel-CG34083|Dmel-CG34085|Dmel-CG34086|Dmel-CG34089|Dmel-CG34090|Dmel-CG34092|Dmel-CR34060|Dmel-CR34061|Dmel-CR34062|Dmel-CR34064|Dmel-CR34065|Dmel-CR34066|Dmel-CR34068|Dmel-CR34070|Dmel-CR34071|Dmel-CR34075|Dmel-CR34077|Dmel-CR34078|Dmel-CR34079|Dmel-CR34080|Dmel-CR34081|Dmel-CR34082|Dmel-CR34084|Dmel-CR34087|Dmel-CR34088|Dmel-CR34091|Dmel-CR34093|Dmel-CR34094|Dmel-CR34095|Dmel-CR34096")
mda_s1_r2_raw[["percent.nora"]] <- PercentageFeatureSet(mda_s1_r2_raw, pattern = "Dmel-nora-virus")

mda_s2_r1_raw[["percent.mito"]] <- PercentageFeatureSet(mda_s2_r1_raw, pattern = "Dmel-CG34063|Dmel-CG34067|Dmel-CG34069|Dmel-CG34072|Dmel-CG34073|Dmel-CG34074|Dmel-CG34076|Dmel-CG34083|Dmel-CG34085|Dmel-CG34086|Dmel-CG34089|Dmel-CG34090|Dmel-CG34092|Dmel-CR34060|Dmel-CR34061|Dmel-CR34062|Dmel-CR34064|Dmel-CR34065|Dmel-CR34066|Dmel-CR34068|Dmel-CR34070|Dmel-CR34071|Dmel-CR34075|Dmel-CR34077|Dmel-CR34078|Dmel-CR34079|Dmel-CR34080|Dmel-CR34081|Dmel-CR34082|Dmel-CR34084|Dmel-CR34087|Dmel-CR34088|Dmel-CR34091|Dmel-CR34093|Dmel-CR34094|Dmel-CR34095|Dmel-CR34096")
mda_s2_r1_raw[["percent.nora"]] <- PercentageFeatureSet(mda_s2_r1_raw, pattern = "Dmel-nora-virus")

mda_s2_r2_raw[["percent.mito"]] <- PercentageFeatureSet(mda_s2_r2_raw, pattern = "Dmel-CG34063|Dmel-CG34067|Dmel-CG34069|Dmel-CG34072|Dmel-CG34073|Dmel-CG34074|Dmel-CG34076|Dmel-CG34083|Dmel-CG34085|Dmel-CG34086|Dmel-CG34089|Dmel-CG34090|Dmel-CG34092|Dmel-CR34060|Dmel-CR34061|Dmel-CR34062|Dmel-CR34064|Dmel-CR34065|Dmel-CR34066|Dmel-CR34068|Dmel-CR34070|Dmel-CR34071|Dmel-CR34075|Dmel-CR34077|Dmel-CR34078|Dmel-CR34079|Dmel-CR34080|Dmel-CR34081|Dmel-CR34082|Dmel-CR34084|Dmel-CR34087|Dmel-CR34088|Dmel-CR34091|Dmel-CR34093|Dmel-CR34094|Dmel-CR34095|Dmel-CR34096")
mda_s2_r2_raw[["percent.nora"]] <- PercentageFeatureSet(mda_s2_r2_raw, pattern = "Dmel-nora-virus")

# combine samples
# raw counts from technical replicates will be merged before integration

mda_s1_r1_filtered <- subset(mda_s1_r1_raw, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mito < 25)
mda_s1_r2_filtered <- subset(mda_s1_r2_raw, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mito < 25)

mda_s1_filtered <- merge(mda_s1_r1_filtered,
                         y = mda_s1_r2_filtered, add.cell.ids = c("S1R1", "S1R2"), project = "10x_S1")

mda_s1_filtered <- SCTransform(object = mda_s1_filtered,
                               vars.to.regress = c("percent.mito", "percent.nora"),
                               return.only.var.genes = FALSE)

mda_s2_r1_filtered <- subset(mda_s2_r1_raw, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mito < 25)
mda_s2_r2_filtered <- subset(mda_s2_r2_raw, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mito < 25)

mda_s2_filtered <- merge(mda_s2_r1_filtered,
                         y = mda_s2_r2_filtered, add.cell.ids = c("S2R1", "S2R2"), project = "10x_S2")

mda_s2_filtered <- SCTransform(object = mda_s2_filtered,
                               vars.to.regress = c("percent.mito", "percent.nora"),
                               return.only.var.genes = FALSE)

mda.lst <- list(mda_s1_filtered, mda_s2_filtered)

# integrate datasets

mda.features <- SelectIntegrationFeatures(object.list = mda.lst, nfeatures = 3000)

mda.lst <- PrepSCTIntegration(object.list = mda.lst,
                              anchor.features = mda.features)

mda_10x_anchors <- FindIntegrationAnchors(object.list = mda.lst,
                                          dims = 1:50,
                                          normalization.method = "SCT",
                                          anchor.features = mda.features)

mda_10x_merged <- IntegrateData(mda_10x_anchors, dims = 1:50,
                                normalization.method = "SCT")

# cluster generation

mda_10x_merged <- RunPCA(mda_10x_merged, features = VariableFeatures(object = mda_10x_merged))
mda_10x_merged <- FindNeighbors(mda_10x_merged, dims = 1:50)
mda_10x_merged <- FindClusters(mda_10x_merged, resolution = 1.0)

mda_10x_merged <- RunTSNE(mda_10x_merged, dims = 1:50)
mda_10x_merged <- RunUMAP(mda_10x_merged, dims = 1:50)

# set SCT normalized values obtained before integration as default assay

DefaultAssay(mda_10x_merged) <- "SCT"

# rename major subtypes

mda_new.cluster.ids <- c("pEC", # 0
                         "aEC", # 1
                         "LFC", # 2
                         "ISC/EB", # 3
                         "aEC", # 4
                         "aEC", # 5
                         "mEC", # 6
                         "pEC", # 7
                         "LFC", # 8
                         "copper/iron", # 9
                         "pEC", # 10
                         "EE", # 11
                         "ISC/EB", # 12
                         "aEC", # 13
                         "cardia", # 14
                         "aEC", # 15
                         "pEC", # 16
                         "pEC", # 17
                         "unk") # 18

mda_cluster_names <- unique(mda_new.cluster.ids)
mda_cluster_order <- c("ISC/EB", # progenitor
                       "EE", # scattered
                       "unk", # unknown location
                       "cardia", # proventriculus
                       "aEC", # anterior
                       "mEC", "copper/iron", # medium
                       "LFC", # posterior medium region
                       "pEC") # posterior 

names(mda_new.cluster.ids) <- levels(mda_10x_merged)
mda_10x_merged <- RenameIdents(mda_10x_merged, mda_new.cluster.ids)

# subset infected cells

mda_10x_uninfected <- subset(mda_10x_merged, subset = `Dmel-nora-virus` <= 1)
mda_10x_nora <- subset(mda_10x_merged, subset = `Dmel-nora-virus` > 1)

# calculate percentage and log expression of DMelNV for each cell

mda_raw_reads_virus_pct <- virus_pct_exp_fun(mda_10x_merged)
mda_raw_nora_pct <- sum(mda_raw_reads_virus_pct$nora_counts)/sum(mda_raw_reads_virus_pct$total_counts)

mda_nora_pct <- virus_pct_exp_fun(mda_10x_nora)
mda_nora_pct$cluster <- factor(mda_nora_pct$cluster, levels = mda_cluster_order)

# one-way anova virus percentage and log expression between cell types

mda_nora_cluster_anova <- list(anova(lm(nora_exp ~ cluster, mda_nora_pct)),
                               anova(lm(nora_pct ~ cluster, mda_nora_pct)))

# post hoc tukey-HSD test

mda_nora_tukey <- list()
mda_nora_tukey[1] <- TukeyHSD(aov(lm(mda_nora_pct$nora_exp ~ mda_nora_pct$cluster)), 'mda_nora_pct$cluster', conf.level = 0.95)
mda_nora_tukey[2:6] <- HSD.test(lm(mda_nora_pct$nora_exp ~ mda_nora_pct$cluster), 'mda_nora_pct$cluster')

# proportion of DMelNV-infected midgut cells

prn <- data.frame(cell = mda_cluster_names,
                  nora = numeric(9),
                  uninfected = numeric(9),
                  stringsAsFactors = FALSE)

for (i in 1:9) {
  n <- try(length(WhichCells(mda_10x_nora, idents = mda_cluster_names[[i]])), silent = TRUE)
  
  if("try-error" %in% class(n)) {
    n <- 0
  }
  
  prn[i, 2] <- n
  prn[i, 3] <- length(WhichCells(mda_10x_uninfected, idents = mda_cluster_names[[i]]))
}

prn_inf <- prn
prn_inf[nrow(prn_inf) + 1,] <- list("total", as.numeric(sum(prn_inf$nora)), 
                                  as.numeric(sum(prn_inf$uninfected)))

# do a different gene list for midgut data

mda_genes <- do_gene_list(min.pct = 0.1,
                          obj = mda_10x_merged)

# partial correlation analysis

mda_pcor <- do_correlations(dataset = "midgut")

# bonferroni correction

mda_pcor <- mda_pcor[complete.cases(mda_pcor),]
mda_pcor$p.adj <- p.adjust(mda_pcor$p.value, "bonferroni")

mda_pcor_names <- gsub("-", "_", rownames(mda_pcor[mda_pcor$p.adj < 0.05,]))
mda_pcor_names_pos <- gsub("-", "_", rownames(mda_pcor[mda_pcor$p.adj < 0.05 & mda_pcor$estimate > 0,]))
mda_pcor_names_neg <- gsub("-", "_", rownames(mda_pcor[mda_pcor$p.adj < 0.05 & mda_pcor$estimate < 0,]))

# aov

mda_aov <- do_aov("midgut")

mda_aov <- mda_aov[complete.cases(mda_aov), ]

mda_aov$p.adj <- numeric(nrow(mda_aov))
mda_aov[mda_aov$group == "cell",]$p.adj <- p.adjust(mda_aov[mda_aov$group == "cell",]$Pr..F., "bonferroni" )
mda_aov[mda_aov$group == "infection",]$p.adj <- p.adjust(mda_aov[mda_aov$group == "infection",]$Pr..F., "bonferroni" )
mda_aov[mda_aov$group == "cell:infection",]$p.adj <- p.adjust(mda_aov[mda_aov$group == "cell:infection",]$Pr..F., "bonferroni" )

# make lists of DEGs

mda_infection_genes <- gsub("-", "_",
                             mda_aov[mda_aov$group == "infection" & mda_aov$p.adj < 0.05,]$gene)

mda_cell.infection_genes <- gsub("-", "_",
                                  mda_aov[mda_aov$group == "cell:infection" & mda_aov$p.adj < 0.05,]$gene)

mda_exclusively_cell.infection_genes <- mda_cell.infection_genes[!mda_cell.infection_genes %in% mda_infection_genes]

# mean logFC of midgut DEGs

mean_logFC_mda_infection <- up_down_fun(obj_vir = mda_10x_nora,
                                        gene_lst = mda_infection_genes,
                                        obj_no_vir = mda_10x_uninfected,
                                        cluster_names.list = mda_cluster_names)

# create symbol lists

mda_pcor_names_pos_symbol <- gene_symbol[gene_symbol$V1 %in% mda_pcor_names_pos, 3]

mda_pcor_names_neg_symbol <- gene_symbol[gene_symbol$V1 %in% mda_pcor_names_neg, 3]

mda_infection_genes_symbol <- gene_symbol[gene_symbol$V1 %in% mda_infection_genes, 3]

mda_infection_genes_symbol_pos <- mean_logFC_mda_infection[mean_logFC_mda_infection$mean_logFC > 0, "symbol"]

mda_infection_genes_symbol_neg <- mean_logFC_mda_infection[mean_logFC_mda_infection$mean_logFC < 0, "symbol"]

mda_cell.infection_genes_symbol <- gene_symbol[gene_symbol$V1 %in% mda_cell.infection_genes, 3]

mda_exclusively_cell.infection_genes_symbol <- gene_symbol[gene_symbol$V1 %in% mda_exclusively_cell.infection_genes, 3]

# save symbol lists for GO enrichment analysis

write.table(mda_pcor_names_pos_symbol, file = "GO/midgut_atlas/pcor_pos_symbol.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(mda_pcor_names_neg_symbol, file = "GO/midgut_atlas/pcor_neg_symbol.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(mda_exclusively_cell.infection_genes_symbol,
            file = "GO/midgut_atlas/exclusively_cell.infection_genes_symbol", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(mda_infection_genes_symbol_pos,
            file = "GO/midgut_atlas/infection_genes_symbol_pos", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(mda_infection_genes_symbol_neg,
            file = "GO/midgut_atlas/infection_genes_symbol_neg", row.names = FALSE, col.names = FALSE, quote = FALSE)

# DEGs that are DE only in a subtype specific manner in the midgut atlas dataset AND DE only in a generic manner in the EE dataset
# those should be EE specific transcriptional DEGs

EE_specific_DEGs <- mda_exclusively_cell.infection_genes[mda_exclusively_cell.infection_genes %in% nora_exclusively_infection_genes]
EE_specific_DEGs_symbol <- gene_symbol[gene_symbol$V1 %in% EE_specific_DEGs, 3]

# common DMelNV DEGs between datsets

nora_mda_common_pcor_names <- mda_pcor_names[mda_pcor_names %in% nora_pcor_names]

nora_mda_common_cell.infection <- mda_cell.infection_genes[mda_cell.infection_genes %in% nora_cell.infection_genes]

nora_mda_common_infection <- mda_infection_genes[mda_infection_genes %in% nora_infection_genes]

#####################
## quality control ##
#####################

# sensitivity curve
# this analysis was performed on a cluster since it takes quite some time (sensitivity_curve.R and mda_sensitivity_curve.R)
# save files to run the script

save(dmel_no_vir, file = "sensitivity_curve/dmel_no_vir.Rdata")
save(genes, file = "sensitivity_curve/genes.Rdata")

save(mda_10x_uninfected, file = "sensitivity_curve/mda_10x_uninfected.Rdata")
save(mda_genes, file = "sensitivity_curve/mda_genes.Rdata")
