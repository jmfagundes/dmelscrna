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

# cuffdiff cannot differentiate the expression of overlapping genes
# split merged genes

nora.bulk.day2 <- nora.bulk.day2 %>% separate_rows(1, sep = ",")
nora.bulk.day10 <- nora.bulk.day10 %>% separate_rows(1, sep = ",")
nora.bulk.day20 <- nora.bulk.day20 %>% separate_rows(1, sep = ",")
nora.bulk.day30 <- nora.bulk.day30 %>% separate_rows(1, sep = ",")

nora.bulk.response <- list(`DMelNV systemic day 2 up` = (nora.bulk.day2[nora.bulk.day2$`log2(fold_change)` > 0 & nora.bulk.day2$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 2 down` = (nora.bulk.day2[nora.bulk.day2$`log2(fold_change)` < 0 & nora.bulk.day2$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 10 up` = (nora.bulk.day10[nora.bulk.day10$`log2(fold_change)` > 0 & nora.bulk.day10$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 10 down` = (nora.bulk.day10[nora.bulk.day10$`log2(fold_change)` < 0 & nora.bulk.day10$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 20 up` = (nora.bulk.day20[nora.bulk.day20$`log2(fold_change)` > 0 & nora.bulk.day20$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 20 down` = (nora.bulk.day20[nora.bulk.day20$`log2(fold_change)` < 0 & nora.bulk.day20$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 30 up` = (nora.bulk.day30[nora.bulk.day30$`log2(fold_change)` > 0 & nora.bulk.day30$p_value < 0.05, "gene_name"] %>% as.list())$gene_name,
                           `DMelNV systemic day 30 down` = (nora.bulk.day30[nora.bulk.day30$`log2(fold_change)` < 0 & nora.bulk.day30$p_value < 0.05, "gene_name"] %>% as.list())$gene_name)

nora.bulk.graph.results <- network.pipe(nora.bulk.response)

# investigate translation related genes in "DMelNV I-pCCHa1 down" and "DMelNV pEC down"

nora.translation_genes <- list(`DMelNV I-pCCHa1 down-regulated translation-related genes` = strsplit(nora.celltype.graph.results$`DMelNV I-pCCHa1 down`$reactome[[2]]["R-DME-72766", "geneID"], "/")[[1]],
                               `DMelNV pEC down-regulated translation-related genes` = strsplit(nora.celltype.graph.results$`DMelNV pEC down`$reactome[[2]]["R-DME-72766", "geneID"], "/")[[1]])

nora.translation_genes.results <- network.pipe(nora.translation_genes)

# map DMelNV cellular response DEGs to systemic response

nora_dmel_sc2bulk <- list()

for (lst.index in 1:length(nora.celltype.response)) {
  
  lst <- nora.celltype.response[[lst.index]]$symbol
  
  if (length(lst) == 0) next
  
  sc2bulk <- data.frame(systemic.lst = character(),
                        symbol = character())
  
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
  
  sc2bulk <- data.frame(symbol = character(),
                        systemic.lst = character())
  
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

all.graph.results <- data.frame(`gene_list` = character(),
                                `subgraph_slope` = numeric(),
                                `slope_diff` = numeric(),
                                `slope_diff_p.value` = numeric(),
                                `mean_betweenness` = numeric(),
                                `mean_betweenness_diff` = numeric(),
                                `betweenness_diff_p.value` = numeric())

all.react.results <- data.frame(`gene list` = character(),
                                description = character(),
                                p.adjust = numeric(),
                                q.value = numeric(),
                                count = numeric(),
                                geneID = character(),
                                reactome = character())

for (lst in c(thika.graph.results,
              nora.graph.results,
              thika.celltype.graph.results,
              nora.celltype.graph.results,
              nora.cellular_bulk_diff.results,
              nora.translation_genes.results,
              nora.bulk.graph.results)) {
  
  all.graph.results[nrow(all.graph.results) + 1,] <- lst$summary.row
  
  if (nrow(lst$reactome[[2]]) == 0) next # skip if no enrichment found
  
  selected.reactome.cols <- lst$reactome[[2]][lst$reactome[[2]]$p.adjust < 0.05, c("Description", "p.adjust", "qvalue", "Count", "geneID")]
  
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

thika_cor.tmp <- thika.cor.gampoi
nora_EE_cor.tmp <- nora_EE.cor.gampoi
nora_MA_cor.tmp <- nora_MA.cor.gampoi

thika_cor.tmp$symbol <- gene_symbol[match(gsub("-", "_", thika_cor.tmp$name), gene_symbol$V1), "V3"]
nora_EE_cor.tmp$symbol <- gene_symbol[match(gsub("-", "_", nora_EE_cor.tmp$name), gene_symbol$V1), "V3"]
nora_MA_cor.tmp$symbol <- gene_symbol[match(gsub("-", "_", nora_MA_cor.tmp$name), gene_symbol$V1), "V3"]

write_xlsx(c(list(`TV generic response` = thika_infection_tmp,
                  `TV posterior EEs response` = thika_posterior_tmp,
                  `DMelNV (EE) generic response` = nora_EE_infection_tmp,
                  `DMelNV (MA) generic response` = nora_MA_infection_tmp,
                  `TV correlation` = thika_cor.tmp,
                  `DMelNV (EE) correlation` = nora_EE_cor.tmp,
                  `DMelNV (MA) correlation` = nora_MA_cor.tmp,
                  coinfection = coinfection.DE.tmp),
             thika_cell_tmp, nora_cell_tmp,
             list(`common cellular systemic DMelNV` = nora_dmel_sc2bulk.tb)), "Supplementary 1.xlsx")

# summarize DEGs into one table

all_DE <- bind_rows(c(list(`TV generic response` = thika_infection_tmp,
                           `TV posterior EEs response` = thika_posterior_tmp,
                           `DMelNV (EE) generic response` = nora_EE_infection_tmp,
                           `DMelNV (MA) generic response` = nora_MA_infection_tmp,
                           `TV correlation` = thika_cor.tmp,
                           `DMelNV (EE) correlation` = nora_EE_cor.tmp,
                           `DMelNV (MA) correlation` = nora_MA_cor.tmp),
                      thika_cell_tmp, nora_cell_tmp), .id = "response")

all_DE$symbol <- gene_symbol[match(gsub("-", "_", all_DE$name), gene_symbol$V1), "V3"]

all_DE[grep("_", all_DE$response), "response"] <- gsub("_", "/", all_DE[grep("_", all_DE$response), "response"])
all_DE$response <- factor(all_DE$response, levels = all_DE$response %>% unique())

all_DEGs <- all_DE[all_DE$adj_pval < 0.05, "symbol"] %>% unique()

# to better visualize, lfc > 5 will be capped to 5 and lfc < -5 will be capped to -5

all_DE_capped <- all_DE
all_DE_capped[all_DE$lfc < -5, "lfc"] <- -5
all_DE_capped[all_DE$lfc > 5, "lfc"] <- 5

# complete volcano plots without labels

ggvolcano <- ggplot(all_DE_capped, aes(x = lfc, y = -log10(adj_pval))) +
  geom_point(size = 0.4, aes(color = adj_pval < 0.05)) +
  facet_wrap(~response, ncol = 4) +
  ylab(label = "adjusted P-values (-log10)") + xlab(label = "fold change (log2)") +
  theme(legend.position = "none", strip.text.x = element_text(size = 7.25))

ggsave("figures_out/complete_volcano_plots.pdf", ggvolcano, width = 174, height = 200, units = "mm")

# volcano plots with heat shock-related labels
# get list of proteins from Cellular response to heat stress

R_DME_3371556 <- read.table("Participating Molecules [R-DME-3371556].tsv", sep = "\t", header = TRUE)[1:3]

heatshock.proteins <- select(org.Dm.eg.db, R_DME_3371556$Identifier, "SYMBOL", "UNIPROT")

all_DE_capped_heat.shock <- all_DE_capped %>%
  filter(symbol %in% heatshock.proteins$SYMBOL)

all_DE_capped_heat.shock <- all_DE_capped_heat.shock %>%
  filter(response %in% unique(all_DE_capped_heat.shock[all_DE_capped_heat.shock$adj_pval < 0.05, "response"]))

ggvolcano.heat <- ggplot(all_DE_capped_heat.shock,
                         aes(x = lfc, y = -log10(adj_pval))) +
  geom_point(size = 0.4, aes(color = adj_pval < 0.05)) +
  geom_text_repel(data = all_DE_capped_heat.shock %>% filter(adj_pval < 0.05),
                  aes(x = lfc, y = -log10(adj_pval),
                  label = symbol), size = 2, segment.size = .2, seed = 255) +
  facet_wrap(~response, ncol = 4) +
  ylab(label = "adjusted P-values (-log10)") + xlab(label = "fold change (log2)") +
  theme(legend.position = "none", strip.text.x = element_text(size = 7.25))

ggsave("figures_out/fig4.pdf", ggvolcano.heat, width = 174, height = 150, units = "mm")

# number of DEGs for each test

n.DEGs <- data.frame(response = all_DE$response %>% levels(),
                     up = all_DE[all_DE$adj_pval < 0.05 & all_DE$lfc > 0, "response"] %>% table() %>% as.integer(),
                     down = all_DE[all_DE$adj_pval < 0.05 & all_DE$lfc < 0, "response"] %>% table() %>% as.integer())

n.DEGs$virus <- n.DEGs$response %>% gsub(" .*", "", .)
n.DEGs$response <- n.DEGs$response %>% gsub("(TV|DMelNV|DMelNV \\(EE\\)) ", "", .)
n.DEGs$response <- factor(n.DEGs$response, levels = c("generic response", "I-a", "III", "I-ap-a", "II-a",
                                                      "I-m", "II-m1", "II-m2",
                                                      "II-p", "I-ap-p", "I-pCCHa1", "I-pAstA", "posterior EEs response",
                                                      "(MA) generic response", "ISC/EB",
                                                      "EE", "unk", "cardia", "aEC",
                                                      "mEC", "copper/iron", "LFC",
                                                      "pEC", "correlation", "(MA) correlation"))

# barplot of DEGs

n.DEGs$down <- -n.DEGs$down
n.DEGs[n.DEGs$up == 0, "up"] <- NA
n.DEGs[n.DEGs$down == 0, "down"] <- NA

# reactome figure

cell.react.results <- list(DMelNV = all.react.results[!grepl("systemic", all.react.results$`gene list`) &
                                                        grepl("DMelNV", all.react.results$`gene list`),],
                           TV = all.react.results[grepl("TV", all.react.results$`gene list`),])

# plot reactome pathways as networks

# TV

TV.reactome.edges <- cell.react.results$TV %>%
  dplyr::select(c("Description", "geneID", "gene list")) %>%
  separate_rows(2, sep = "/")

TV.reactome.net <- graph_from_edgelist(as.matrix(TV.reactome.edges[1:2]))

V(TV.reactome.net)$color <- ifelse(names(V(TV.reactome.net)) %in% TV.reactome.edges$geneID, "#00BFC4", "grey")
V(TV.reactome.net)$pathway <- TV.reactome.edges$Description[match(names(V(TV.reactome.net)), TV.reactome.edges$Description)]
V(TV.reactome.net)$gene <- TV.reactome.edges$geneID[match(names(V(TV.reactome.net)), TV.reactome.edges$geneID)]

# set seed to layout

set.seed(255)
TV.reactome.net.dh.lay <- create_layout(TV.reactome.net, layout = "dh")
set.seed(2)
TV.reactome.net.fr.lay <- create_layout(TV.reactome.net, layout = "fr")

# DMelNV 

DMelNV.reactome.edges.down <- cell.react.results$DMelNV %>%
  filter(grepl("down", cell.react.results$DMelNV$`gene list`) & !grepl("translation", cell.react.results$DMelNV$`gene list`)) %>%
  dplyr::select(c("Description", "geneID", "gene list")) %>%
  separate_rows(2, sep = "/")

DMelNV.reactome.net.down <- graph_from_edgelist(as.matrix(DMelNV.reactome.edges.down[1:2]))

V(DMelNV.reactome.net.down)$color <- ifelse(names(V(DMelNV.reactome.net.down)) %in% DMelNV.reactome.edges.down$geneID, "#F8766D", "grey")
V(DMelNV.reactome.net.down)$pathway <- DMelNV.reactome.edges.down$Description[match(names(V(DMelNV.reactome.net.down)), DMelNV.reactome.edges.down$Description)]
V(DMelNV.reactome.net.down)$gene <- DMelNV.reactome.edges.down$geneID[match(names(V(DMelNV.reactome.net.down)), DMelNV.reactome.edges.down$geneID)]

# set seed to dh layout

set.seed(100)
DMelNV.reactome.net.down.dh.lay <- create_layout(DMelNV.reactome.net.down, layout = "dh")

DMelNV.reactome.edges.up <- cell.react.results$DMelNV %>%
  filter(grepl("up", cell.react.results$DMelNV$`gene list`)) %>%
  dplyr::select(c("Description", "geneID", "gene list")) %>%
  separate_rows(2, sep = "/")

DMelNV.reactome.net.up <- graph_from_edgelist(as.matrix(DMelNV.reactome.edges.up[1:2]))

V(DMelNV.reactome.net.up)$color <- ifelse(names(V(DMelNV.reactome.net.up)) %in% DMelNV.reactome.edges.up$geneID, "#00BFC4", "grey")
V(DMelNV.reactome.net.up)$pathway <- DMelNV.reactome.edges.up$Description[match(names(V(DMelNV.reactome.net.up)), DMelNV.reactome.edges.up$Description)]
V(DMelNV.reactome.net.up)$gene <- DMelNV.reactome.edges.up$geneID[match(names(V(DMelNV.reactome.net.up)), DMelNV.reactome.edges.up$geneID)]

# set seed to dh layout

set.seed(255)
DMelNV.reactome.net.up.dh.lay <- create_layout(DMelNV.reactome.net.up, layout = "dh")

# make figure 3

ggbarplot <- ggplot(n.DEGs %>% 
                      pivot_longer(!response & !virus, names_to = "regulation"),
                    aes(x = response, y = value, fill = regulation)) +
  geom_col() +
  ggtitle("A") +
  facet_wrap(~virus, ncol = 1) +
  theme(strip.background = element_blank(),
        strip.text = element_text(angle = 0, hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "top",
        legend.text = element_text(size = 8),
        plot.title = element_text(margin = margin(t = 0, b = -15), hjust = -.2),
        legend.key.size = unit(.5, "cm")) +
  scale_y_continuous(breaks = c(-100, 0 , 100, 200), labels = c(100, 0 , 100, 200)) +
  ylab(label = "number of DEGs") +
  scale_fill_discrete(name = "", labels = c("down-regulated", "up-regulated")) +
  xlab("")

gg.TV.reactome.net <- ggraph(TV.reactome.net.fr.lay) +
  geom_edge_link(color = "grey85") + 
  geom_node_point(color = V(TV.reactome.net)$color) +
  geom_node_text(aes(label = V(TV.reactome.net)$pathway), size = 1.25, color = "black", repel = TRUE, seed = 10, segment.size = .2) +
  geom_node_text(aes(label = V(TV.reactome.net)$gene), size = .75, color = "grey10") +
  theme_void() +
  ggtitle("C    TV")

gg.DMelNV.reactome.net.down <- ggraph(DMelNV.reactome.net.down.dh.lay) +
  geom_edge_link(color = "grey85") + 
  geom_node_point(color = V(DMelNV.reactome.net.down)$color) +
  geom_node_text(aes(label = V(DMelNV.reactome.net.down)$pathway), size = 1.25, color = "black", repel = TRUE, force = 100, seed = 2, segment.size = .2) +
  geom_node_text(aes(label = V(DMelNV.reactome.net.down)$gene), size = .75, color = "grey10") +
  theme_void() +
  ggtitle("B                                                          DMelNV")

gg.DMelNV.reactome.net.up <- ggraph(DMelNV.reactome.net.up.dh.lay) +
  geom_edge_link(color = "grey85") + 
  geom_node_point(color = V(DMelNV.reactome.net.up)$color) +
  geom_node_text(aes(label = V(DMelNV.reactome.net.up)$pathway), size = 1.25, color = "black", repel = TRUE, force = 100, seed = 16, segment.size = .2, max.time = 5) +
  geom_node_text(aes(label = V(DMelNV.reactome.net.up)$gene), size = .75, color = "grey10") +
  theme_void()

g <- ggarrange(ggbarplot, gg.TV.reactome.net,
               gg.DMelNV.reactome.net.down, gg.DMelNV.reactome.net.up,
               common.legend = TRUE, legend = "top",
               heights = c(1, 1.5))

ggsave("figures_out/fig3.pdf", g, width = 174, height = 234, units = "mm")

# clear space

rm(thika_infection_tmp, nora_EE_infection_tmp, nora_MA_infection_tmp,
   thika_cell_tmp, nora_cell_tmp)

# test if any individual gene in the cellular response is an articulation point

dro_net.articulation.points <- articulation_points(dro_net) %>% names()

DE_articulation.point <- all_DE %>% filter(symbol %in% dro_net.articulation.points & adj_pval < 0.05)

# articulation points in late systemic response

bulk_articulation.point <- list(down20 = nora.bulk.response$`DMelNV systemic day 20 down`[nora.bulk.response$`DMelNV systemic day 20 down` %in% dro_net.articulation.points],
                                down30 = nora.bulk.response$`DMelNV systemic day 30 down`[nora.bulk.response$`DMelNV systemic day 30 down` %in% dro_net.articulation.points])

# test if correlated genes are articulation points

cor_articulation.point <- rbind(thika.cor.gampoi_significant %>% filter(symbol %in% dro_net.articulation.points & adj_pval < 0.05),
                                nora_EE.cor.gampoi_significant %>% filter(symbol %in% dro_net.articulation.points & adj_pval < 0.05),
                                nora_MA.cor.gampoi_significant %>% filter(symbol %in% dro_net.articulation.points & adj_pval < 0.05))

cor_articulation.point$response <- c(rep("TV", nrow(thika.cor.gampoi_significant %>% filter(symbol %in% dro_net.articulation.points & adj_pval < 0.05))),
                                     rep("DMelNV EE", nrow(nora_EE.cor.gampoi_significant %>% filter(symbol %in% dro_net.articulation.points & adj_pval < 0.05))),
                                     rep("DMelNV MA", nrow(nora_MA.cor.gampoi_significant %>% filter(symbol %in% dro_net.articulation.points & adj_pval < 0.05))))

# plot number of DEGs in bulk RNA-seq and pathway enrichment

n.DEGs_nora.bulk.response <- data.frame(`up-regulated` = c(length(nora.bulk.response$`DMelNV systemic day 2 up`),
                                                           length(nora.bulk.response$`DMelNV systemic day 10 up`),
                                                           length(nora.bulk.response$`DMelNV systemic day 20 up`),
                                                           length(nora.bulk.response$`DMelNV systemic day 30 up`)),
                                        `down-regulated` = -c(length(nora.bulk.response$`DMelNV systemic day 2 down`),
                                                              length(nora.bulk.response$`DMelNV systemic day 10 down`),
                                                              length(nora.bulk.response$`DMelNV systemic day 20 down`),
                                                              length(nora.bulk.response$`DMelNV systemic day 30 down`)),
                                        `days post-eclosion` = c("day 2", "day 10", "day 20", "day 30"))

n.DEGs_nora.bulk.response$days.post.eclosion <- factor(n.DEGs_nora.bulk.response$days.post.eclosion, levels = n.DEGs_nora.bulk.response$days.post.eclosion)

ggbarplot.bulk <- ggplot(n.DEGs_nora.bulk.response %>% 
                           pivot_longer(!`days.post.eclosion`, names_to = "regulation"),
                    aes(x = `days.post.eclosion`, y = value, fill = regulation)) +
  geom_col() +
  ggtitle("A") +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "top",
        legend.text = element_text(size = 6),
        plot.title = element_text(margin = margin(t = 0, b = -15), hjust = -.4),
        legend.key.size = unit(.5, "cm")) +
  scale_y_continuous(breaks = c(-2300, 0 , 1350, 2700), labels = c(2300, 0 , 1350, 2700)) +
  ylab(label = "number of DEGs") +
  scale_fill_discrete(name = "", labels = c("down-regulated", "up-regulated")) +
  xlab("")

# plot complement cascade network from bulk RNA-seq and cellular positive correlation down-regulated at 20-30 days post-eclosion in systemic response

systemic.complement <- all.react.results %>%
  filter(grepl("Complement cascade", Description) & grepl("DMelNV systemic", `gene list`)) %>%
  dplyr::select(c("Description", "geneID", "gene list")) %>%
  separate_rows(2, sep = "/")

systemic.complement.net <- graph_from_edgelist(as.matrix(systemic.complement[1:2]))

V(systemic.complement.net)$color <- ifelse(names(V(systemic.complement.net)) %in% systemic.complement$geneID, "#00BFC4", "grey")
V(systemic.complement.net)$pathway <- systemic.complement$Description[match(names(V(systemic.complement.net)), systemic.complement$Description)]
V(systemic.complement.net)$gene <- systemic.complement$geneID[match(names(V(systemic.complement.net)), systemic.complement$geneID)]

set.seed(100)
systemic.complement.net.dh.lay <- create_layout(systemic.complement.net, layout = "dh")

gg.systemic.complement.net <- ggraph(systemic.complement.net.dh.lay) +
  geom_edge_link(color = "grey85") + 
  geom_node_point(color = V(systemic.complement.net)$color) +
  geom_node_text(aes(label = V(systemic.complement.net)$pathway), size = 1.25, color = "black", repel = TRUE, force = 100, seed = 2, segment.size = .2) +
  geom_node_text(aes(label = V(systemic.complement.net)$gene), size = .75, color = "grey10") +
  theme_void() +
  ggtitle(" B")

positive.cor_systemic.down.edges <- all.react.results %>%
  filter(grepl("positive.*systemic down", `gene list`)) %>%
  dplyr::select(c("Description", "geneID", "gene list")) %>%
  separate_rows(2, sep = "/")

positive.cor_systemic.down.net <- graph_from_edgelist(as.matrix(positive.cor_systemic.down.edges[1:2]))

V(positive.cor_systemic.down.net)$color <- ifelse(names(V(positive.cor_systemic.down.net)) %in% positive.cor_systemic.down.edges$geneID, "#C77CFF", "grey")
V(positive.cor_systemic.down.net)$pathway <- positive.cor_systemic.down.edges$Description[match(names(V(positive.cor_systemic.down.net)), positive.cor_systemic.down.edges$Description)]
V(positive.cor_systemic.down.net)$gene <- positive.cor_systemic.down.edges$geneID[match(names(V(positive.cor_systemic.down.net)), positive.cor_systemic.down.edges$geneID)]

set.seed(100)
positive.cor_systemic.down.net.dh.lay <- create_layout(positive.cor_systemic.down.net, layout = "dh")

gg.positive.cor_systemic.down.net <- ggraph(positive.cor_systemic.down.net.dh.lay) +
  geom_edge_link(color = "grey85") + 
  geom_node_point(color = V(positive.cor_systemic.down.net)$color) +
  geom_node_text(aes(label = V(positive.cor_systemic.down.net)$pathway), size = 1.25, color = "black", repel = TRUE, force = 100, seed = 2, segment.size = .2) +
  geom_node_text(aes(label = V(positive.cor_systemic.down.net)$gene), size = .75, color = "grey10") +
  theme_void() +
  ggtitle("C")

g <- ggarrange(ggbarplot.bulk, gg.systemic.complement.net,
               gg.positive.cor_systemic.down.net,
               legend = "top", nrow = 1)

ggsave("figures_out/fig6.pdf", g, width = 174, height = 70, units = "mm")

# articulation points and translation-related hubs

# number of articulation points in each DEG list

n.articulation_nora.bulk.response <- data.frame(response = c("day 2", "day 10", "day 20", "day 30"),
                                                up = c(length(nora.bulk.response$`DMelNV systemic day 2 up`[nora.bulk.response$`DMelNV systemic day 2 up` %in% dro_net.articulation.points]),
                                                       length(nora.bulk.response$`DMelNV systemic day 10 up`[nora.bulk.response$`DMelNV systemic day 10 up` %in% dro_net.articulation.points]),
                                                       length(nora.bulk.response$`DMelNV systemic day 20 up`[nora.bulk.response$`DMelNV systemic day 20 up` %in% dro_net.articulation.points]),
                                                       length(nora.bulk.response$`DMelNV systemic day 30 up`[nora.bulk.response$`DMelNV systemic day 30 up` %in% dro_net.articulation.points])),
                                                down = -c(length(nora.bulk.response$`DMelNV systemic day 2 down`[nora.bulk.response$`DMelNV systemic day 2 down` %in% dro_net.articulation.points]),
                                                          length(nora.bulk.response$`DMelNV systemic day 10 down`[nora.bulk.response$`DMelNV systemic day 10 down` %in% dro_net.articulation.points]),
                                                          length(nora.bulk.response$`DMelNV systemic day 20 down`[nora.bulk.response$`DMelNV systemic day 20 down` %in% dro_net.articulation.points]),
                                                          length(nora.bulk.response$`DMelNV systemic day 30 down`[nora.bulk.response$`DMelNV systemic day 30 down` %in% dro_net.articulation.points])),
                                                virus ="DMelNV systemic")

n.articulation_nora.bulk.response$response <- factor(n.articulation_nora.bulk.response$response, levels = c("day 2", "day 10", "day 20", "day 30"))

n.articulation_nora.bulk.response[n.articulation_nora.bulk.response$up == 0, "up"] <- NA
n.articulation_nora.bulk.response[n.articulation_nora.bulk.response$down == 0, "down"] <- NA

gg.n.articulation_nora.bulk.response <- ggplot(n.articulation_nora.bulk.response %>% 
                                                 pivot_longer(!response & !virus, names_to = "regulation"),
                                               aes(x = response, y = value, fill = regulation)) +
  geom_col() +
  facet_wrap(~virus, ncol = 1) +
  theme(strip.background = element_blank(),
        strip.text = element_text(angle = 0, hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "top",
        legend.text = element_text(size = 8),
        plot.title = element_text(margin = margin(t = 0, b = -15), hjust = -.2, vjust = .1),
        legend.key.size = unit(.5, "cm")) +
  scale_y_continuous(breaks = c(-100, -50, 0, 50, 100), labels = c(100, 50, 0, 50, 100)) +
  ylab(label = "") +
  scale_fill_discrete(name = "", labels = c("down-regulated", "up-regulated")) +
  xlab("")

n.articulation <- data.frame(response = all_DE$response %>% levels(),
                             up = all_DE[all_DE$adj_pval < 0.05 & all_DE$lfc > 0 & all_DE$symbol %in% dro_net.articulation.points, "response"] %>% table() %>% as.integer(),
                             down = all_DE[all_DE$adj_pval < 0.05 & all_DE$lfc < 0 & all_DE$symbol %in% dro_net.articulation.points, "response"] %>% table() %>% as.integer())

n.articulation$virus <- n.articulation$response %>% gsub(" .*", "", .)
n.articulation$response <- n.articulation$response %>% gsub("(TV|DMelNV|DMelNV \\(EE\\)) ", "", .)
n.articulation$response <- factor(n.articulation$response, levels = c("generic response", "I-a", "III", "I-ap-a", "II-a",
                                                      "I-m", "II-m1", "II-m2",
                                                      "II-p", "I-ap-p", "I-pCCHa1", "I-pAstA", "posterior EEs response",
                                                      "(MA) generic response", "ISC/EB",
                                                      "EE", "unk", "cardia", "aEC",
                                                      "mEC", "copper/iron", "LFC",
                                                      "pEC", "correlation", "(MA) correlation"))

n.articulation$down <- -n.articulation$down
n.articulation[n.articulation$up == 0, "up"] <- NA
n.articulation[n.articulation$down == 0, "down"] <- NA

all_n.articulation <- rbind(n.articulation_nora.bulk.response,
                            n.articulation)

gg.n.articulation <- ggplot(n.articulation %>% 
         pivot_longer(!response & !virus, names_to = "regulation"),
       aes(x = response, y = value, fill = regulation)) +
  geom_col() +
  ggtitle("A") +
  facet_wrap(~virus, ncol = 1) +
  theme(strip.background = element_blank(),
        strip.text = element_text(angle = 0, hjust = 0),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "top",
        legend.text = element_text(size = 8),
        plot.title = element_text(margin = margin(t = 0, b = -15), hjust = -.1),
        legend.key.size = unit(.5, "cm")) +
  scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), labels = c(10, 5, 0, 5, 10)) +
  ylab(label = "number of articulation points") +
  scale_fill_discrete(name = "", labels = c("down-regulated", "up-regulated")) +
  xlab("")

g.articulation <- ggarrange(gg.n.articulation, ggarrange(gg.n.articulation_nora.bulk.response, "", ncol = 1, legend = "none", heights = c(1, .8)),
                            common.legend = TRUE, legend = "top")

# hubs

nora.translation_genes.newline <- list(`DMelNV I-pCCHa1 down-regulated\ntranslation-related genes` = strsplit(nora.celltype.graph.results$`DMelNV I-pCCHa1 down`$reactome[[2]]["R-DME-72766", "geneID"], "/")[[1]],
                                       `DMelNV pEC down-regulated\ntranslation-related genes` = strsplit(nora.celltype.graph.results$`DMelNV pEC down`$reactome[[2]]["R-DME-72766", "geneID"], "/")[[1]])

# reactome and network analyses results

write_xlsx(list(`reactome results` = all.react.results,
                `network analysis results` = all.graph.results,
                `articulation points` = dro_net.articulation.points %>% as.data.frame()), "Supplementary 2.xlsx")

# save figs

ggsave("figures_out/fig7A.pdf", g.articulation, width = 6.85, height = 3.5)

pdf("figures_out/fig7B.pdf", width = 6.85, height = 3.5)
layout(rbind(c(1,2)))
network.pipe(nora.translation_genes.newline)
dev.off()

