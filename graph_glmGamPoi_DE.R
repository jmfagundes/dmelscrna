source("source_me.R")

load("objs/glm.vir_glmGamPoi.Rdata")

# this script will perform all gene network analysis

# load predicted dmel interactome
# this file uses flybase ID, which are converted to symbol with fbgn_annotation_ID_fb_2020_02.tsv

dro_gene_interaction_raw <- read.table("network/networks1.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
fly_symbol <- read.table("network/fbgn_annotation_ID_fb_2020_02.tsv", sep = '\t', quote = "", comment.char = "#")
dro_edges <- tibble(from = dro_gene_interaction_raw$gene1, to = dro_gene_interaction_raw$gene2)

sec_ids <- split_sec_ids()

dro_edges[is.na(fly_symbol[match(dro_edges$from, fly_symbol$V3), 1]), ]$from <- sec_ids[match(
  dro_edges[is.na(fly_symbol[match(dro_edges$from, fly_symbol$V3), 1]), ]$from, sec_ids$sec), 1]

dro_edges[is.na(fly_symbol[match(dro_edges$to, fly_symbol$V3), 1]), ]$to <- sec_ids[match(
  dro_edges[is.na(fly_symbol[match(dro_edges$to, fly_symbol$V3), 1]), ]$to, sec_ids$sec), 1]

# convert flybase ID to symbol for primary flybase IDs

dro_edges[!is.na(fly_symbol[match(dro_edges$from, fly_symbol$V3), 1]), ]$from <- fly_symbol[match(
  dro_edges[!is.na(fly_symbol[match(dro_edges$from, fly_symbol$V3), 1]), ]$from, fly_symbol$V3), 1]

dro_edges[!is.na(fly_symbol[match(dro_edges$to, fly_symbol$V3), 1]), ]$to <- fly_symbol[match(
  dro_edges[!is.na(fly_symbol[match(dro_edges$to, fly_symbol$V3), 1]), ]$to, fly_symbol$V3), 1]

# remove duplicated rows

dro_edges <- unique(dro_edges)

# create dmel interactome network

dro_net <- graph_from_edgelist(as.matrix(dro_edges), directed = FALSE)
dro_net_betweenness <- betweenness(dro_net)

# create subnetworks
# merge generic response genes to correlated genes

thika_glmGamPoi_infection_genes <- thika_gGP_infection[thika_gGP_infection$adj_pval < 0.05,]
thika_glmGamPoi_infection_genes$symbol <- gene_symbol[match(gsub("-", "_", thika_glmGamPoi_infection_genes$name), gene_symbol$V1), "V3"]
nora_EE_glmGamPoi_infection_genes <- nora_EE_gGP_infection[nora_EE_gGP_infection$adj_pval < 0.05,]
nora_EE_glmGamPoi_infection_genes$symbol <- gene_symbol[match(gsub("-", "_", nora_EE_glmGamPoi_infection_genes$name), gene_symbol$V1), "V3"]
nora_MA_glmGamPoi_infection_genes <- nora_MA_gGP_infection[nora_MA_gGP_infection$adj_pval < 0.05,]
nora_MA_glmGamPoi_infection_genes$symbol <- gene_symbol[match(gsub("-", "_", nora_MA_glmGamPoi_infection_genes$name), gene_symbol$V1), "V3"]

thika.cor.gampoi_significant <- thika.cor.gampoi[thika.cor.gampoi$adj_pval < 0.05,]
thika.cor.gampoi_significant$symbol <- gene_symbol[match(gsub("-", "_", thika.cor.gampoi_significant$name), gene_symbol$V1), "V3"]
nora_EE.cor.gampoi_significant <- nora_EE.cor.gampoi[nora_EE.cor.gampoi$adj_pval < 0.05,]
nora_EE.cor.gampoi_significant$symbol <- gene_symbol[match(gsub("-", "_", nora_EE.cor.gampoi_significant$name), gene_symbol$V1), "V3"]
nora_MA.cor.gampoi_significant <- nora_MA.cor.gampoi[nora_MA.cor.gampoi $adj_pval < 0.05,]
nora_MA.cor.gampoi_significant$symbol <- gene_symbol[match(gsub("-", "_", nora_MA.cor.gampoi_significant$name), gene_symbol$V1), "V3"]

thika_posterior <- thika_region_DE$thika_gGP_p[thika_region_DE$thika_gGP_p$adj_pval < 0.05,]

# create lists

thika.generic.response <- list(`TV generic up` = thika_glmGamPoi_infection_genes[thika_glmGamPoi_infection_genes$lfc > 0, "symbol"],
                               `TV generic down` = thika_glmGamPoi_infection_genes[thika_glmGamPoi_infection_genes$lfc < 0, "symbol"])

nora.generic.response <- list(`DMelNV (EE) up` = nora_EE_glmGamPoi_infection_genes[nora_EE_glmGamPoi_infection_genes$lfc > 0, "symbol"],
                              `DMelNV (EE) down` = nora_EE_glmGamPoi_infection_genes[nora_EE_glmGamPoi_infection_genes$lfc < 0, "symbol"],
                              `DMelNV (EE) positive correlation` = nora_EE.cor.gampoi_significant[nora_EE.cor.gampoi_significant$lfc > 0, "symbol"],
                              `DMelNV (EE) negative correlation` = nora_EE.cor.gampoi_significant[nora_EE.cor.gampoi_significant$lfc < 0, "symbol"],
                              `DMelNV (MA) up` = nora_MA_glmGamPoi_infection_genes[nora_MA_glmGamPoi_infection_genes$lfc > 0, "symbol"],
                              `DMelNV (MA) down` = nora_MA_glmGamPoi_infection_genes[nora_MA_glmGamPoi_infection_genes$lfc < 0, "symbol"],
                              `DMelNV (MA) positive correlation` = nora_MA.cor.gampoi_significant[nora_MA.cor.gampoi_significant$lfc > 0, "symbol"],
                              `DMelNV (MA) negative correlation` = nora_MA.cor.gampoi_significant[nora_MA.cor.gampoi_significant$lfc < 0, "symbol"])

# cell-type-specific DEGs

thika.celltype.response <- list(`TV I-m up` = thika_gGP_oneway_DE$I_m[thika_gGP_oneway_DE$I_m$lfc > 0 & thika_gGP_oneway_DE$I_m$adj_pval < 0.05,],
                                `TV II-m1 up` = thika_gGP_oneway_DE$II_m1[thika_gGP_oneway_DE$II_m1$lfc > 0 & thika_gGP_oneway_DE$II_m1$adj_pval < 0.05,],
                                `TV I-ap-p up` = thika_gGP_oneway_DE$I_ap_p[thika_gGP_oneway_DE$I_ap_p$lfc > 0 & thika_gGP_oneway_DE$I_ap_p$adj_pval < 0.05,],
                                `TV II-a up` = thika_gGP_oneway_DE$II_a[thika_gGP_oneway_DE$II_a$lfc > 0 & thika_gGP_oneway_DE$II_a$adj_pval < 0.05,],
                                `TV II-m2 up` = thika_gGP_oneway_DE$II_m2[thika_gGP_oneway_DE$II_m2$lfc > 0 & thika_gGP_oneway_DE$II_m2$adj_pval < 0.05,],
                                `TV I-pAstA up` = thika_gGP_oneway_DE$I_pAstA[thika_gGP_oneway_DE$I_pAstA$lfc > 0 & thika_gGP_oneway_DE$I_pAstA$adj_pval < 0.05,],
                                `TV I-pCCHa1 up` = thika_gGP_oneway_DE$I_pCCHa1[thika_gGP_oneway_DE$I_pCCHa1$lfc > 0 & thika_gGP_oneway_DE$I_pCCHa1$adj_pval < 0.05,],
                                `TV II-p up` = thika_gGP_oneway_DE$II_p[thika_gGP_oneway_DE$II_p$lfc > 0 & thika_gGP_oneway_DE$II_p$adj_pval < 0.05,],
                                `TV I-ap-a up` = thika_gGP_oneway_DE$I_ap_a[thika_gGP_oneway_DE$I_ap_a$lfc > 0 & thika_gGP_oneway_DE$I_ap_a$adj_pval < 0.05,],
                                `TV I-a up` = thika_gGP_oneway_DE$I_a[thika_gGP_oneway_DE$I_a$lfc > 0 & thika_gGP_oneway_DE$I_a$adj_pval < 0.05,],
                                `TV III up` = thika_gGP_oneway_DE$III[thika_gGP_oneway_DE$III$lfc > 0 & thika_gGP_oneway_DE$III$adj_pval < 0.05,],
                                `TV posterior up` = thika_posterior[thika_posterior$lfc > 0 & thika_posterior$adj_pval < 0.05,],
                                 
                                `TV I-m down` = thika_gGP_oneway_DE$I_m[thika_gGP_oneway_DE$I_m$lfc < 0 & thika_gGP_oneway_DE$I_m$adj_pval < 0.05,],
                                `TV II-m1 down` = thika_gGP_oneway_DE$II_m1[thika_gGP_oneway_DE$II_m1$lfc < 0 & thika_gGP_oneway_DE$II_m1$adj_pval < 0.05,],
                                `TV I-ap-p down` = thika_gGP_oneway_DE$I_ap_p[thika_gGP_oneway_DE$I_ap_p$lfc < 0 & thika_gGP_oneway_DE$I_ap_p$adj_pval < 0.05,],
                                `TV II-a down` = thika_gGP_oneway_DE$II_a[thika_gGP_oneway_DE$II_a$lfc < 0 & thika_gGP_oneway_DE$II_a$adj_pval < 0.05,],
                                `TV II-m2 down` = thika_gGP_oneway_DE$II_m2[thika_gGP_oneway_DE$II_m2$lfc < 0 & thika_gGP_oneway_DE$II_m2$adj_pval < 0.05,],
                                `TV I-pAstA down` = thika_gGP_oneway_DE$I_pAstA[thika_gGP_oneway_DE$I_pAstA$lfc < 0 & thika_gGP_oneway_DE$I_pAstA$adj_pval < 0.05,],
                                `TV I-pCCHa1 down` = thika_gGP_oneway_DE$I_pCCHa1[thika_gGP_oneway_DE$I_pCCHa1$lfc < 0 & thika_gGP_oneway_DE$I_pCCHa1$adj_pval < 0.05,],
                                `TV II-p down` = thika_gGP_oneway_DE$II_p[thika_gGP_oneway_DE$II_p$lfc < 0 & thika_gGP_oneway_DE$II_p$adj_pval < 0.05,],
                                `TV I-ap-a down` = thika_gGP_oneway_DE$I_ap_a[thika_gGP_oneway_DE$I_ap_a$lfc < 0 & thika_gGP_oneway_DE$I_ap_a$adj_pval < 0.05,],
                                `TV I-a down` = thika_gGP_oneway_DE$I_a[thika_gGP_oneway_DE$I_a$lfc < 0 & thika_gGP_oneway_DE$I_a$adj_pval < 0.05,],
                                `TV III down` = thika_gGP_oneway_DE$III[thika_gGP_oneway_DE$III$lfc < 0 & thika_gGP_oneway_DE$III$adj_pval < 0.05,],
                                `TV posterior down` = thika_posterior[thika_posterior$lfc < 0 & thika_posterior$adj_pval < 0.05,])

nora.celltype.response <- list(`DMelNV I-m up` = nora_gGP_oneway_DE$I_m[nora_gGP_oneway_DE$I_m$lfc > 0 & nora_gGP_oneway_DE$I_m$adj_pval < 0.05,],
                               `DMelNV II-m1 up` = nora_gGP_oneway_DE$II_m1[nora_gGP_oneway_DE$II_m1$lfc > 0 & nora_gGP_oneway_DE$II_m1$adj_pval < 0.05,],
                               `DMelNV I-ap-p up` = nora_gGP_oneway_DE$I_ap_p[nora_gGP_oneway_DE$I_ap_p$lfc > 0 & nora_gGP_oneway_DE$I_ap_p$adj_pval < 0.05,],
                               `DMelNV II-a up` = nora_gGP_oneway_DE$II_a[nora_gGP_oneway_DE$II_a$lfc > 0 & nora_gGP_oneway_DE$II_a$adj_pval < 0.05,],
                               `DMelNV II-m2 up` = nora_gGP_oneway_DE$II_m2[nora_gGP_oneway_DE$II_m2$lfc > 0 & nora_gGP_oneway_DE$II_m2$adj_pval < 0.05,],
                               `DMelNV I-pAstA up` = nora_gGP_oneway_DE$I_pAstA[nora_gGP_oneway_DE$I_pAstA$lfc > 0 & nora_gGP_oneway_DE$I_pAstA$adj_pval < 0.05,],
                               `DMelNV I-pCCHa1 up` = nora_gGP_oneway_DE$I_pCCHa1[nora_gGP_oneway_DE$I_pCCHa1$lfc > 0 & nora_gGP_oneway_DE$I_pCCHa1$adj_pval < 0.05,],
                               `DMelNV II-p up` = nora_gGP_oneway_DE$II_p[nora_gGP_oneway_DE$II_p$lfc > 0 &  nora_gGP_oneway_DE$II_p$adj_pval < 0.05,],
                               `DMelNV I-ap-a up` = nora_gGP_oneway_DE$I_ap_a[nora_gGP_oneway_DE$I_ap_a$lfc > 0 & nora_gGP_oneway_DE$I_ap_a$adj_pval < 0.05,],
                               `DMelNV EE (MA) up` = nora_gGP_oneway_DE$EE[nora_gGP_oneway_DE$EE$lfc > 0 & nora_gGP_oneway_DE$EE$adj_pval < 0.05,],
                               `DMelNV pEC up` = nora_gGP_oneway_DE$pEC[nora_gGP_oneway_DE$pEC$lfc > 0 & nora_gGP_oneway_DE$pEC$adj_pval < 0.05,],
                               `DMelNV aEC up` = nora_gGP_oneway_DE$aEC[nora_gGP_oneway_DE$aEC$lfc > 0 & nora_gGP_oneway_DE$aEC$adj_pval < 0.05,],
                               `DMelNV cardia up` = nora_gGP_oneway_DE$cardia[nora_gGP_oneway_DE$cardia$lfc > 0 & nora_gGP_oneway_DE$cardia$adj_pval < 0.05,],
                               `DMelNV copper/iron up` = nora_gGP_oneway_DE$`copper/iron`[nora_gGP_oneway_DE$`copper/iron`$lfc > 0 & nora_gGP_oneway_DE$`copper/iron`$adj_pval < 0.05,],
                               `DMelNV ISC/EB up` = nora_gGP_oneway_DE$`ISC/EB`[nora_gGP_oneway_DE$`ISC/EB`$lfc > 0 & nora_gGP_oneway_DE$`ISC/EB`$adj_pval < 0.05,],
                               `DMelNV LFC up` = nora_gGP_oneway_DE$LFC[nora_gGP_oneway_DE$LFC$lfc > 0 & nora_gGP_oneway_DE$LFC$adj_pval < 0.05,],
                               `DMelNV mEC up` = nora_gGP_oneway_DE$mEC[nora_gGP_oneway_DE$mEC$lfc > 0 & nora_gGP_oneway_DE$mEC$adj_pval < 0.05,],
                               
                               `DMelNV I-m down` = nora_gGP_oneway_DE$I_m[nora_gGP_oneway_DE$I_m$lfc < 0 & nora_gGP_oneway_DE$I_m$adj_pval < 0.05,],
                               `DMelNV II-m1 down` = nora_gGP_oneway_DE$II_m1[nora_gGP_oneway_DE$II_m1$lfc < 0 & nora_gGP_oneway_DE$II_m1$adj_pval < 0.05,],
                               `DMelNV I-ap-p down` = nora_gGP_oneway_DE$I_ap_p[nora_gGP_oneway_DE$I_ap_p$lfc < 0 & nora_gGP_oneway_DE$I_ap_p$adj_pval < 0.05,],
                               `DMelNV II-a down` = nora_gGP_oneway_DE$II_a[nora_gGP_oneway_DE$II_a$lfc < 0 & nora_gGP_oneway_DE$II_a$adj_pval < 0.05,],
                               `DMelNV II-m2 down` = nora_gGP_oneway_DE$II_m2[nora_gGP_oneway_DE$II_m2$lfc < 0 & nora_gGP_oneway_DE$II_m2$adj_pval < 0.05,],
                               `DMelNV I-pAstA down` = nora_gGP_oneway_DE$I_pAstA[nora_gGP_oneway_DE$I_pAstA$lfc < 0 & nora_gGP_oneway_DE$I_pAstA$adj_pval < 0.05,],
                               `DMelNV I-pCCHa1 down` = nora_gGP_oneway_DE$I_pCCHa1[nora_gGP_oneway_DE$I_pCCHa1$lfc < 0 & nora_gGP_oneway_DE$I_pCCHa1$adj_pval < 0.05,],
                               `DMelNV II-p down` = nora_gGP_oneway_DE$II_p[nora_gGP_oneway_DE$II_p$lfc < 0 & nora_gGP_oneway_DE$II_p$adj_pval < 0.05,],
                               `DMelNV I-ap-a down` = nora_gGP_oneway_DE$I_ap_a[nora_gGP_oneway_DE$I_ap_a$lfc < 0 & nora_gGP_oneway_DE$I_ap_a$adj_pval < 0.05,],
                               `DMelNV EE (MA) down` = nora_gGP_oneway_DE$EE[nora_gGP_oneway_DE$EE$lfc < 0 & nora_gGP_oneway_DE$EE$adj_pval < 0.05,],
                               `DMelNV pEC down` = nora_gGP_oneway_DE$pEC[nora_gGP_oneway_DE$pEC$lfc < 0 & nora_gGP_oneway_DE$pEC$adj_pval < 0.05,],
                               `DMelNV aEC down` = nora_gGP_oneway_DE$aEC[nora_gGP_oneway_DE$aEC$lfc < 0 & nora_gGP_oneway_DE$aEC$adj_pval < 0.05,],
                               `DMelNV cardia down` = nora_gGP_oneway_DE$cardia[nora_gGP_oneway_DE$cardia$lfc < 0 & nora_gGP_oneway_DE$cardia$adj_pval < 0.05,],
                               `DMelNV copper/iron down` = nora_gGP_oneway_DE$`copper/iron`[nora_gGP_oneway_DE$`copper/iron`$lfc < 0 & nora_gGP_oneway_DE$`copper/iron`$adj_pval < 0.05,],
                               `DMelNV ISC/EB down` = nora_gGP_oneway_DE$`ISC/EB`[nora_gGP_oneway_DE$`ISC/EB`$lfc < 0 & nora_gGP_oneway_DE$`ISC/EB`$adj_pval < 0.05,],
                               `DMelNV LFC down` = nora_gGP_oneway_DE$LFC[nora_gGP_oneway_DE$LFC$lfc < 0 & nora_gGP_oneway_DE$LFC$adj_pval < 0.05,],
                               `DMelNV mEC down` = nora_gGP_oneway_DE$mEC[nora_gGP_oneway_DE$mEC$lfc < 0 & nora_gGP_oneway_DE$mEC$adj_pval < 0.05,])

# add symbols

for (i in 1:length(thika.celltype.response)) {
  if (is.null(thika.celltype.response[[i]])) next
  thika.celltype.response[[i]]$symbol <- gene_symbol[match(gsub("-", "_", thika.celltype.response[[i]]$name), gene_symbol$V1), "V3"]
}

for (i in 1:length(nora.celltype.response)) {
  if (is.null(nora.celltype.response[[i]])) next
  nora.celltype.response[[i]]$symbol <- gene_symbol[match(gsub("-", "_", nora.celltype.response[[i]]$name), gene_symbol$V1), "V3"]
}

# summarize all cell-type-specific response genes and add to lists

thika.celltype.response$`TV all cell-type-specific DEGs up` <- bind_rows(thika.celltype.response[1:11])
thika.celltype.response$`TV all cell-type-specific DEGs down` <- bind_rows(thika.celltype.response[13:23])

nora.celltype.response$`DMelNV all cell-type-specific DEGs up` <- bind_rows(nora.celltype.response[1:17])
nora.celltype.response$`DMelNV all cell-type-specific DEGs down` <- bind_rows(nora.celltype.response[18:34])

# iterate through these lists and for every list, do reactome analysis, calculate betweenness and compare hubness to complete interactome

thika.graph.results <- network.pipe(thika.generic.response)
nora.graph.results <- network.pipe(nora.generic.response)

thika.celltype.graph.results <- network.pipe(thika.celltype.response, celltype = TRUE)
nora.celltype.graph.results <- network.pipe(nora.celltype.response, celltype = TRUE)

# compare cellular response to systemic response to DMelNV

# create subgraphs for bulk RNA-seq DEGs

nora.bulk.day2 <- fread("gene_matrices/dmelnv_bulk_RNA-seq_data/NV_infected_day2_vs_uninfected.tsv",
              sep = "\t", dec = ",", stringsAsFactors = FALSE)
nora.bulk.day2$day <- "2"
nora.bulk.day10 <- fread("gene_matrices/dmelnv_bulk_RNA-seq_data/NV_infected_day10_vs_uninfected.tsv",
               sep = "\t", dec = ",", stringsAsFactors = FALSE)
nora.bulk.day10$day <- "10"
nora.bulk.day20 <- fread("gene_matrices/dmelnv_bulk_RNA-seq_data/NV_infected_day20_vs_uninfected.tsv",
               sep = "\t", dec = ",", stringsAsFactors = FALSE)
nora.bulk.day20$day <- "20"
nora.bulk.day30 <- fread("gene_matrices/dmelnv_bulk_RNA-seq_data/NV_infected_day30_vs_uninfected.tsv",
               sep = "\t", dec = ",", stringsAsFactors = FALSE)
nora.bulk.day30$day <- "30"

nora.bulk.response <- list(`DMelNV systemic day 2 up` = (nora.bulk.day2[nora.bulk.day2$`log2(fold_change)` > 0 & nora.bulk.day2$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 2 down` = (nora.bulk.day2[nora.bulk.day2$`log2(fold_change)` < 0 & nora.bulk.day2$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 10 up` = (nora.bulk.day10[nora.bulk.day10$`log2(fold_change)` > 0 & nora.bulk.day10$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 10 down` = (nora.bulk.day10[nora.bulk.day10$`log2(fold_change)` < 0 & nora.bulk.day10$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 20 up` = (nora.bulk.day20[nora.bulk.day20$`log2(fold_change)` > 0 & nora.bulk.day20$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 20 down` = (nora.bulk.day20[nora.bulk.day20$`log2(fold_change)` < 0 & nora.bulk.day20$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 30 up` = (nora.bulk.day30[nora.bulk.day30$`log2(fold_change)` > 0 & nora.bulk.day30$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 30 down` = (nora.bulk.day30[nora.bulk.day30$`log2(fold_change)` < 0 & nora.bulk.day30$p_value < 0.05, "gene_name"] %>% as.list())$gene_name)

nora.bulk.graph.results <- network.pipe(nora.bulk.response)

# investigate translation related genes in "DMelNV (MA) negative correlation", "DMelNV I-pCCHa1 down" and "DMelNV pEC down"

nora.translation_genes <- list(`DMelNV I-pCCHa1 down-regulated\ntranslation-related genes` = strsplit(nora.celltype.graph.results$`DMelNV I-pCCHa1 down`$reactome[[2]]["R-DME-72766", "geneID"], "/")[[1]],
                               `DMelNV pEC down-regulated\ntranslation-related genes` = strsplit(nora.celltype.graph.results$`DMelNV pEC down`$reactome[[2]]["R-DME-72766", "geneID"], "/")[[1]])

# save figs then rename list to save tables

pdf("figures_out/fig4.pdf", width = 6.85, height = 3.5)
layout(rbind(c(1,2)))
network.pipe(nora.translation_genes)
dev.off()

nora.translation_genes <- list(`DMelNV I-pCCHa1 down-regulated translation-related genes` = strsplit(nora.celltype.graph.results$`DMelNV I-pCCHa1 down`$reactome[[2]]["R-DME-72766", "geneID"], "/")[[1]],
                               `DMelNV pEC down-regulated translation-related genes` = strsplit(nora.celltype.graph.results$`DMelNV pEC down`$reactome[[2]]["R-DME-72766", "geneID"], "/")[[1]])

nora.translation_genes.results <- network.pipe(nora.translation_genes)

# map DMelNV cellular response DEGs to systemic response

nora_dmel_sc2bulk <- list()

for (lst.index in 1:length(nora.celltype.response)) {
  
  lst <- nora.celltype.response[[lst.index]]$symbol
  
  if (length(lst) == 0) next
  
  sc2bulk <- data.frame(systemic.lst = character(0),
                        symbol = character(0))
  
  for (sys.degs in 1:length(nora.bulk.response)) {
    
    common_symbols <- nora.bulk.response[[sys.degs]][nora.bulk.response[[sys.degs]] %in% lst]
    
    if (length(common_symbols) == 0) next
    
    sc2bulk <- rbind(sc2bulk,
                     data.frame(systemic.lst = names(nora.bulk.response[sys.degs]),
                                symbol = common_symbols))
    
  }
  nora_dmel_sc2bulk[[names(nora.celltype.response[lst.index])]] <- sc2bulk
}

for (lst.index in 1:length(nora.generic.response)) {
  
  lst <- nora.generic.response[[lst.index]]
  
  if (length(lst) == 0) next
  
  sc2bulk <- data.frame(symbol = character(0),
                        systemic.lst = character(0))
  
  for (sys.degs in 1:length(nora.bulk.response)) {
    
    common_symbols <- nora.bulk.response[[sys.degs]][nora.bulk.response[[sys.degs]] %in% lst]
    
    if (length(common_symbols) == 0) next
    
    sc2bulk <- rbind(sc2bulk,
                     data.frame(symbol = common_symbols,
                                systemic.lst = names(nora.bulk.response[sys.degs])))
    
  }
  nora_dmel_sc2bulk[[names(nora.generic.response[lst.index])]] <- sc2bulk
}

nora_dmel_sc2bulk.tb <- bind_rows(nora_dmel_sc2bulk, .id = "cellular.lst")

# DEGs in both systemic and cellular responses

nora_common.degs <- c(nora_dmel_sc2bulk$`DMelNV all cell-type-specific DEGs up`$symbol,
                      nora_dmel_sc2bulk$`DMelNV all cell-type-specific DEGs down`$symbol,
                      nora_dmel_sc2bulk$`DMelNV (EE) up`$symbol,
                      nora_dmel_sc2bulk$`DMelNV (EE) down`$symbol,
                      nora_dmel_sc2bulk$`DMelNV (MA) up`$symbol,
                      nora_dmel_sc2bulk$`DMelNV (MA) down`$symbol) %>% unique()

nora_common.cor.degs <- c(nora_dmel_sc2bulk$`DMelNV (EE) positive correlation`$symbol,
                          nora_dmel_sc2bulk$`DMelNV (EE) negative correlation`$symbol,
                          nora_dmel_sc2bulk$`DMelNV (MA) positive correlation`$symbol,
                          nora_dmel_sc2bulk$`DMelNV (MA) negative correlation`$symbol) %>% unique()

cellular.up_systemic.down <- nora_dmel_sc2bulk.tb[grep("up", nora_dmel_sc2bulk.tb$cellular.lst),]
cellular.up_systemic.down <- cellular.up_systemic.down[grep("down", cellular.up_systemic.down$systemic.lst),]

cellular.pos.cor_systemic.down <- nora_dmel_sc2bulk.tb[grep("positive", nora_dmel_sc2bulk.tb$cellular.lst),]
cellular.pos.cor_systemic.down <- cellular.pos.cor_systemic.down[grep("down", cellular.pos.cor_systemic.down$systemic.lst),]

cellular.down_systemic.up <- nora_dmel_sc2bulk.tb[grep("down", nora_dmel_sc2bulk.tb$cellular.lst),]
cellular.down_systemic.up <- cellular.down_systemic.up[grep("up", cellular.down_systemic.up$systemic.lst),]

cellular.neg.cor_systemic.up <- nora_dmel_sc2bulk.tb[grep("negative", nora_dmel_sc2bulk.tb$cellular.lst),]
cellular.neg.cor_systemic.up <- cellular.neg.cor_systemic.up[grep("up", cellular.neg.cor_systemic.up$systemic.lst),]

nora.cellular_bulk_diff <- list(`DMelNV cellular up- systemic down-regulated (20-30 days)` = cellular.up_systemic.down[grep("(20|30)", cellular.up_systemic.down$systemic.lst), "symbol"] %>% unique,
                                `DMelNV cellular positive correlation systemic down-regulated (20-30 days)` = cellular.pos.cor_systemic.down[grep("(20|30)", cellular.pos.cor_systemic.down$systemic.lst), "symbol"] %>% unique)

nora.cellular_bulk_diff.results <- network.pipe(nora.cellular_bulk_diff)

# summarize all graph and reactome results in a table

all.graph.results <- data.frame(`gene_list` = character(0),
                                `subgraph_slope` = numeric(0),
                                `slope_diff` = numeric(0),
                                `slope_diff_p.value` = numeric(0),
                                `mean_betweenness` = numeric(0),
                                `mean_betweenness_diff` = numeric(0),
                                `betweenness_diff_p.value` = numeric(0))

all.react.results <- data.frame(`gene list` = character(0),
                                description = character(0),
                                p.adjust = numeric(0),
                                q.value = numeric(0),
                                count = numeric(0),
                                reactome = character(0))

for (lst in c(thika.graph.results,
              nora.graph.results,
              thika.celltype.graph.results,
              nora.celltype.graph.results,
              nora.cellular_bulk_diff.results,
              nora.translation_genes.results,
              nora.bulk.graph.results)) {
  
  all.graph.results[nrow(all.graph.results) + 1,] <- lst$summary.row
  
  if (nrow(lst$reactome[[2]]) == 0) next # skip if no enrichment found
  
  selected.reactome.cols <- lst$reactome[[2]][lst$reactome[[2]]$p.adjust < 0.05, c("Description", "p.adjust", "qvalue", "Count")]
  
  if (nrow(selected.reactome.cols) == 0) next # skip if no enriched pathway with p.adj < 0.05
  
  selected.reactome.cols$reactome <- rownames(selected.reactome.cols)
  row.names(selected.reactome.cols) <- NULL
  selected.reactome.cols$`gene list` <- lst$summary.row$title
  
  all.react.results <- rbind(all.react.results,
                             selected.reactome.cols)
}

all.graph.results$slope_diff_adj_p.value <- p.adjust(all.graph.results$slope_diff_p.value, method = "BH")
all.graph.results$betweenness_diff_adj_p.value <- p.adjust(all.graph.results$betweenness_diff_p.value, method = "BH")

write.csv(all.graph.results, file = "objs/all.graph.results.csv", row.names = FALSE, quote = FALSE)
write.csv(all.react.results, file = "objs/all.react.results.csv", row.names = FALSE, quote = FALSE)

# also export as xlsx
# DEGs -> supplementary 1
# graph results will be presented as a table and reactome results as supplementary 2
# common cellular and systemic response to DMelNV genes will appended to supplementary 1

# change names and add symbols to keep xlsx clean and tidy

thika_cell_tmp <- thika_gGP_oneway_DE
names(thika_cell_tmp) <- gsub("_", "-", paste0("TV ", names(thika_cell_tmp)))

nora_cell_tmp <- nora_gGP_oneway_DE
names(nora_cell_tmp) <- gsub("/", "_", gsub("_", "-", names(nora_cell_tmp)))
names(nora_cell_tmp) <- paste0("DMelNV ", names(nora_cell_tmp))

for (lst in 1:length(thika_cell_tmp)) {
  thika_cell_tmp[[lst]]$symbol <- gene_symbol[match(gsub("-", "_", thika_cell_tmp[[lst]]$name), gene_symbol$V1), "V3"]
}

for (lst in 1:length(nora_cell_tmp)) {
  nora_cell_tmp[[lst]]$symbol <- gene_symbol[match(gsub("-", "_", nora_cell_tmp[[lst]]$name), gene_symbol$V1), "V3"]
}

thika_posterior_tmp <- thika_region_DE$thika_gGP_p
thika_infection_tmp <- thika_gGP_infection
nora_EE_infection_tmp <- nora_EE_gGP_infection
nora_MA_infection_tmp <- nora_MA_gGP_infection

thika_posterior_tmp$symbol <- gene_symbol[match(gsub("-", "_", thika_posterior_tmp$name), gene_symbol$V1), "V3"]
thika_infection_tmp$symbol <- gene_symbol[match(gsub("-", "_", thika_infection_tmp$name), gene_symbol$V1), "V3"]
nora_EE_infection_tmp$symbol <- gene_symbol[match(gsub("-", "_", nora_EE_infection_tmp$name), gene_symbol$V1), "V3"]
nora_MA_infection_tmp$symbol <- gene_symbol[match(gsub("-", "_", nora_MA_infection_tmp$name), gene_symbol$V1), "V3"]

write_xlsx(c(list(`TV generic response` = thika_infection_tmp,
                  `TV posterior EEs response` = thika_posterior_tmp,
                  `DMelNV (EE) generic response` = nora_EE_infection_tmp,
                  `DMelNV (MA) generic response` = nora_MA_infection_tmp),
             thika_cell_tmp, nora_cell_tmp,
             list(`common cellular systemic DMelNV` = nora_dmel_sc2bulk.tb)), "Supplementary 1.xlsx")

# summarize DEGs into one table

all_DE <- bind_rows(c(list(`TV generic response` = thika_infection_tmp,
                           `TV posterior EEs response` = thika_posterior_tmp,
                           `DMelNV (EE) generic response` = nora_EE_infection_tmp,
                           `DMelNV (MA) generic response` = nora_MA_infection_tmp),
                      thika_cell_tmp, nora_cell_tmp), .id = "response")

all_DE$symbol <- gene_symbol[match(gsub("-", "_", all_DE$name), gene_symbol$V1), "V3"]

all_DE[grep("_", all_DE$response), "response"] <- gsub("_", "/", all_DE[grep("_", all_DE$response), "response"])
all_DE$response <- factor(all_DE$response, levels = all_DE$response %>% unique())

all_DEGs <- all_DE[all_DE$adj_pval < 0.05, "symbol"] %>% unique()

# to better visualize, lfc > 5 will be capped to 5 and lfc < -5 will be capped to -5

all_DE_capped <- all_DE
all_DE_capped[all_DE$lfc < -5, "lfc"] <- -5
all_DE_capped[all_DE$lfc > 5, "lfc"] <- 5

ggvolcano <- ggplot(all_DE_capped[all_DE_capped$symbol %in% all_DEGs,], aes(x = lfc, y = -log10(adj_pval))) +
  geom_point(size = 0.4, aes(color = adj_pval < 0.05)) +
  xlim(-5, 5) +
  facet_wrap(~response, ncol = 4) +
  ylab(label = "adjusted P-values (-log10)") + xlab(label = "fold change (log2)") +
  theme(legend.position = "none", strip.text.x = element_text(size = 7.25))

ggsave("figures_out/fig3.pdf", ggvolcano, width = 174, height = 200, units = "mm")

# number of DEGs for each test

n.DEGs <- data.frame(response = all_DE$response %>% levels(),
                     up = all_DE[all_DE$adj_pval < 0.05 & all_DE$lfc > 0, "response"] %>% table() %>% as.integer(),
                     down = all_DE[all_DE$adj_pval < 0.05 & all_DE$lfc < 0, "response"] %>% table() %>% as.integer())

# clear space

rm(thika_infection_tmp, nora_EE_infection_tmp, nora_MA_infection_tmp,
   thika_cell_tmp, nora_cell_tmp)

# reactome and network analyses results

write_xlsx(list(`reactome results` = all.react.results,
                `network analysis results` = all.graph.results), "Supplementary 2.xlsx")
