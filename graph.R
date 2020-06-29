#!/usr/bin/Rscript

# this script will perform all gene network analysis

# load predicted dmel interactome
# this file uses flybase ID, which are converted to symbol with fbgn_annotation_ID_fb_2020_02.tsv

dro_gene_interaction_raw <- read.table("network/networks1.csv", sep = ',', header = TRUE, stringsAsFactors = FALSE)
fly_symbol <- read.table("network/fbgn_annotation_ID_fb_2020_02.tsv", sep = '\t', quote = "", comment.char = "#")
dro_edges <- tibble(from = dro_gene_interaction_raw$gene1, to = dro_gene_interaction_raw$gene2)

# function to create a table to covert secondary flybase IDs to symbol first since they return NA when trying to fetch primary IDs

split_sec_ids <- function() {
  sec_ids_tmp <- data.frame(symbol = fly_symbol[fly_symbol$V4 != "", 1],
                            sec = fly_symbol[fly_symbol$V4 != "", 4])
  sec_ids <- data.frame(symbol = character(0),
                        sec = character(0),
                        stringsAsFactors = FALSE)
  for (i in sec_ids_tmp$symbol) {
    sec_split <- as.vector(str_split(sec_ids_tmp[sec_ids_tmp$symbol == i, 2], ",", simplify = TRUE))
    df <- data.frame(symbol = replicate(length(sec_split), i),
                     sec = sec_split,
                     stringsAsFactors = FALSE)
    sec_ids <- rbind(sec_ids, df)
  }
  return(sec_ids)
}
sec_ids <- split_sec_ids()

dro_edges[is.na(fly_symbol[match(dro_edges$from, fly_symbol$V3), 1]), ]$from <- sec_ids[match(
  dro_edges[is.na(fly_symbol[match(dro_edges$from, fly_symbol$V3), 1]), ]$from, sec_ids$sec), 1]

dro_edges[is.na(fly_symbol[match(dro_edges$to, fly_symbol$V3), 1]), ]$to <- sec_ids[match(
  dro_edges[is.na(fly_symbol[match(dro_edges$to, fly_symbol$V3), 1]), ]$to, sec_ids$sec), 1]

# function to compare subgraphs to interactome

hub_fun <- function(subgraph, title) {
  
  # calculate degree probability distribution and fit power law to subgraphs
  
  subgraph_d <- degree(subgraph, mode = "total", loops = FALSE)
  subgraph_dd <- degree.distribution(subgraph, cumulative = FALSE, mode = "total", loops = FALSE)
  subgraph_degree <- 1:max(subgraph_d)
  subgraph_prob_d <- subgraph_dd[-1]
  subgraph_prob_d <- subgraph_prob_d[which(subgraph_prob_d != 0)]
  subgraph_degree <- subgraph_degree[which(subgraph_prob_d != 0)]
  
  # linear regression in the log-log space
  
  subgraph_reg <- lm(log(subgraph_prob_d) ~ log(subgraph_degree))
  subgraph_cozf <- coef(subgraph_reg)
  subgraph_power_law <- function(x) exp(subgraph_cozf[[1]] + subgraph_cozf[[2]] * log(x))
  subgraph_slope <- -subgraph_cozf[[2]]
  
  # make plot
  # will produce some figures but they will not be saved
  # xlim and ylim were precalculated to avoid padding axes
  
  xrange <- c(0.7142857 /2, 142.8000000 *2)
  yrange <- c(0.0002763522 /7, 0.6494275563 *2)
  par(xaxs = "i", yaxs = "i")
  plot(dro_net_prob_d ~ dro_net_degree, log = "xy", xlab = "degree (log)", ylab = "probability (log)",
       ylim = yrange,
       xlim = xrange,
       col = 1, main = title, pch = 15)
  curve(dro_net_power_law, col = "black", add = TRUE, n = length(dro_net_d))
  
  # compute ranges for xlim and ylim
  # xrange <- c(10^par("usr")[1], 10^par("usr")[2])
  # print(xrange)
  # yrange <- c(10^par("usr")[3], 10^par("usr")[4])
  # print(yrange)
  
  par(new = TRUE)
  plot(subgraph_prob_d ~ subgraph_degree, log = "xy", xlab = "", ylab = "",
       ylim = yrange,
       xlim = xrange,
       col = "red", main = "", pch = 17)
  curve(subgraph_power_law, col = "red", add = TRUE, n = length(subgraph_d))
  
  # compare regression coefficients
  
  n_dro_net <- length(dro_net_prob_d)
  n_subgraph <- length(subgraph_prob_d)
  
  dro_net_slope_se <- summary(dro_net_reg)$coefficients[2, 2]
  subgraph_slope_se <- summary(subgraph_reg)$coefficients[2, 2]
  
  # test if the slope difference is positive or negative to set upper or lower tail
  
  if (-dro_net_cozf[[2]] + subgraph_cozf[[2]] < 0) lower.tail <- FALSE
  else lower.tail <- TRUE
  
  # slope difference
  
  slope_diff <- dro_net_cozf[[2]] - subgraph_cozf[[2]]
  
  # standard error difference
  
  se_diff <- sqrt(dro_net_slope_se^2 + subgraph_slope_se^2)
  
  # t statistic
  
  t.stat <- slope_diff/se_diff
  p <- 2*pt(t.stat, df = n_dro_net + n_subgraph - 4, lower.tail = lower.tail)
  
  return(list(title, dro_net_slope, subgraph_slope, p))
}

# convert flybase ID to symbol for primary flybase IDs

dro_edges[!is.na(fly_symbol[match(dro_edges$from, fly_symbol$V3), 1]), ]$from <- fly_symbol[match(
  dro_edges[!is.na(fly_symbol[match(dro_edges$from, fly_symbol$V3), 1]), ]$from, fly_symbol$V3), 1]

dro_edges[!is.na(fly_symbol[match(dro_edges$to, fly_symbol$V3), 1]), ]$to <- fly_symbol[match(
  dro_edges[!is.na(fly_symbol[match(dro_edges$to, fly_symbol$V3), 1]), ]$to, fly_symbol$V3), 1]

# remove duplicated rows

dro_edges <- unique(dro_edges)

# create dmel interactome network

dro_net <- graph_from_edgelist(as.matrix(dro_edges), directed = FALSE)

# create subnetworks

# DMelNV partial correlations

subgraph_nora_pcor_pos <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% nora_pcor_symbol_pos |
                                                                    dro_edges$to %in% nora_pcor_symbol_pos, ]),
                                              directed = FALSE)

subgraph_nora_pcor_neg <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% nora_pcor_symbol_neg |
                                                                    dro_edges$to %in% nora_pcor_symbol_neg, ]),
                                              directed = FALSE)

# TV DEGs

subgraph_thika_aov_pos <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% thika_infection_symbol_pos |
                                                                    dro_edges$to %in% thika_infection_symbol_pos, ]),
                                              directed = FALSE)

subgraph_thika_aov_neg <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% thika_infection_symbol_neg |
                                                                    dro_edges$to %in% thika_infection_symbol_neg, ]),
                                              directed = FALSE)

subgraph_thika_exclusively_cell.infection <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% thika_exclusively_cell.infection_symbol |
                                                                                       dro_edges$to %in% thika_exclusively_cell.infection_symbol, ]),
                                                                 directed = FALSE)

# DMelNV DEGs

subgraph_nora_aov_pos <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% nora_infection_symbol_pos |
                                                                    dro_edges$to %in% nora_infection_symbol_pos, ]),
                                              directed = FALSE)

subgraph_nora_aov_neg <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% nora_infection_symbol_neg |
                                                                    dro_edges$to %in% nora_infection_symbol_neg, ]),
                                              directed = FALSE)

subgraph_nora_exclusively_cell.infection <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% nora_exclusively_cell.infection_symbol |
                                                                                       dro_edges$to %in% nora_exclusively_cell.infection_symbol, ]),
                                                                 directed = FALSE)

# calculate degree probability distribution and fit power law to the interactome

dro_net_d <- degree(dro_net, mode = "total", loops = FALSE)
dro_net_dd <- degree.distribution(dro_net, cumulative = FALSE, mode = "total", loops = FALSE)

dro_net_degree <- 1:max(dro_net_d)
dro_net_prob_d <- dro_net_dd[-1]
dro_net_prob_d <- dro_net_prob_d[which(dro_net_prob_d != 0)]
dro_net_degree <- dro_net_degree[which(dro_net_prob_d != 0)]

# linear regression in the log-log space

dro_net_reg <- lm(log(dro_net_prob_d) ~ log(dro_net_degree))
dro_net_cozf <- coef(dro_net_reg)
dro_net_power_law <- function(x) exp(dro_net_cozf[[1]] + dro_net_cozf[[2]] * log(x))
dro_net_slope <- -dro_net_cozf[[2]]

# compare cellular response to systemic response to DMelNV

# create subgraphs for bulk RNA-seq DEGs

subgraph_nora_bulk_neg <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% bulk_nora_neg$gene_name |
                                                                    dro_edges$to %in% bulk_nora_neg$gene_name, ]),
                                              directed = FALSE)
subgraph_nora_bulk_pos <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% bulk_nora_pos$gene_name |
                                                                    dro_edges$to %in% bulk_nora_pos$gene_name, ]),
                                              directed = FALSE)

# create subgraphs for midgut atlas dataset


subgraph_mda_pcor_pos <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% mda_pcor_names_pos_symbol |
                                                                    dro_edges$to %in% mda_pcor_names_pos_symbol, ]),
                                              directed = FALSE)

subgraph_mda_pcor_neg <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% mda_pcor_names_neg_symbol |
                                                                    dro_edges$to %in% mda_pcor_names_neg_symbol, ]),
                                              directed = FALSE)


subgraph_mda_aov_pos <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% mda_infection_genes_symbol_pos |
                                                                   dro_edges$to %in% mda_infection_genes_symbol_pos, ]),
                                             directed = FALSE)

subgraph_mda_aov_neg <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% mda_infection_genes_symbol_neg |
                                                                   dro_edges$to %in% mda_infection_genes_symbol_neg, ]),
                                             directed = FALSE)

subgraph_mda_exclusively_cell.infection <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% mda_exclusively_cell.infection_genes_symbol |
                                                                                      dro_edges$to %in% mda_exclusively_cell.infection_genes_symbol, ]),
                                                                directed = FALSE)

# calculate betweenness of all graphs

dro_net_betweenness <- betweenness(dro_net)

# correlated genes

thika_pcor_pos_betweenness <- betweenness(dro_net, v = thika_pcor_symbol_pos[thika_pcor_symbol_pos %in% V(dro_net)$name])
thika_pcor_neg_betweenness <- betweenness(dro_net, v = thika_pcor_symbol_neg[thika_pcor_symbol_neg %in% V(dro_net)$name])

nora_pcor_pos_betweenness <- betweenness(dro_net, v = nora_pcor_symbol_pos[nora_pcor_symbol_pos %in% V(dro_net)$name])
nora_pcor_neg_betweenness <- betweenness(dro_net, v = nora_pcor_symbol_neg[nora_pcor_symbol_neg %in% V(dro_net)$name])

# DEGs genes

thika_infecion_pos_betweenness <- betweenness(dro_net, 
                                              v = thika_infection_symbol_pos[thika_infection_symbol_pos %in% V(dro_net)$name])
thika_infecion_neg_betweenness <- betweenness(dro_net, 
                                              v = thika_infection_symbol_neg[thika_infection_symbol_neg %in% V(dro_net)$name])
thika_exclusively_cell.infecion_betweenness <- betweenness(dro_net, 
                                                           v = thika_exclusively_cell.infection_symbol[thika_exclusively_cell.infection_symbol %in% V(dro_net)$name])

nora_infecion_pos_betweenness <- betweenness(dro_net, 
                                             v = nora_infection_symbol_pos[nora_infection_symbol_pos %in% V(dro_net)$name])
nora_infecion_neg_betweenness <- betweenness(dro_net, 
                                             v = nora_infection_symbol_neg[nora_infection_symbol_neg %in% V(dro_net)$name])
nora_exclusively_cell.infecion_betweenness <- betweenness(dro_net, 
                                                          v = nora_exclusively_cell.infection_symbol[nora_exclusively_cell.infection_symbol %in% V(dro_net)$name])

# DMelNV bulk RNA-seq

bulk_nora_neg_betweenness <- betweenness(dro_net, 
                                         v = bulk_nora_neg[bulk_nora_neg$gene_name %in% V(dro_net)$name, ]$gene_name)
bulk_nora_pos_betweenness <- betweenness(dro_net, 
                                         v = bulk_nora_pos[bulk_nora_pos$gene_name %in% V(dro_net)$name, ]$gene_name)

# midgut atlas dataset

mda_pcor_pos_betweenness <- betweenness(dro_net, v = mda_pcor_names_pos_symbol[mda_pcor_names_pos_symbol %in% V(dro_net)$name])
mda_pcor_neg_betweenness <- betweenness(dro_net, v = mda_pcor_names_neg_symbol[mda_pcor_names_neg_symbol %in% V(dro_net)$name])

mda_infecion_pos_betweenness <- betweenness(dro_net, 
                                             v = mda_infection_genes_symbol_pos[mda_infection_genes_symbol_pos %in% V(dro_net)$name])
mda_infecion_neg_betweenness <- betweenness(dro_net, 
                                             v = mda_infection_genes_symbol_neg[mda_infection_genes_symbol_neg %in% V(dro_net)$name])
mda_exclusively_cell.infecion_betweenness <- betweenness(dro_net, 
                                                          v = mda_cell.infection_genes_symbol[mda_cell.infection_genes_symbol %in% V(dro_net)$name])


# one-tailed (upper tail) mann-whitney u-test between betweenness
# skipping thika_pcor

betweenness_u.test <- data.frame(DE = character(0),
                                 p.value = numeric(0),
                                 mean.diff = numeric(0),
                                 mean.betweenness = numeric(0),
                                 stringsAsFactors = FALSE)

betweenness_u.test[1,] <- list("nora_pcor_pos",
                               wilcox.test(nora_pcor_pos_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(nora_pcor_pos_betweenness) - mean(dro_net_betweenness),
                               mean(nora_pcor_pos_betweenness))
betweenness_u.test[2,] <- list("nora_pcor_neg",
                               wilcox.test(nora_pcor_neg_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(nora_pcor_neg_betweenness) - mean(dro_net_betweenness),
                               mean(nora_pcor_neg_betweenness))

betweenness_u.test[3,] <- list("nora_aov_pos",
                               wilcox.test(nora_infecion_pos_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(nora_infecion_pos_betweenness) - mean(dro_net_betweenness),
                               mean(nora_infecion_pos_betweenness))
betweenness_u.test[4,] <- list("nora_aov_neg",
                               wilcox.test(nora_infecion_neg_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(nora_infecion_neg_betweenness) - mean(dro_net_betweenness),
                               mean(nora_infecion_neg_betweenness))
betweenness_u.test[5,] <- list("nora_only_cell.infection",
                               wilcox.test(nora_exclusively_cell.infecion_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(nora_exclusively_cell.infecion_betweenness) - mean(dro_net_betweenness),
                               mean(nora_exclusively_cell.infecion_betweenness))

betweenness_u.test[6,] <- list("thika_aov_pos",
                               wilcox.test(thika_infecion_pos_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(thika_infecion_pos_betweenness) - mean(dro_net_betweenness),
                               mean(thika_infecion_pos_betweenness))
betweenness_u.test[7,] <- list("thika_aov_neg",
                               wilcox.test(thika_infecion_neg_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(thika_infecion_neg_betweenness) - mean(dro_net_betweenness),
                               mean(thika_infecion_neg_betweenness))
betweenness_u.test[8,] <- list("thika_only_cell.infection",
                               wilcox.test(thika_exclusively_cell.infecion_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(thika_exclusively_cell.infecion_betweenness) - mean(dro_net_betweenness),
                               mean(thika_exclusively_cell.infecion_betweenness))

betweenness_u.test[9,] <- list("bulk_nora_neg",
                               wilcox.test(bulk_nora_neg_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(bulk_nora_neg_betweenness) - mean(dro_net_betweenness),
                               mean(bulk_nora_neg_betweenness))
betweenness_u.test[10,] <- list("bulk_nora_pos",
                                wilcox.test(bulk_nora_pos_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                                mean(bulk_nora_pos_betweenness) - mean(dro_net_betweenness),
                                mean(bulk_nora_pos_betweenness))

betweenness_u.test[11,] <- list("mda_pcor_pos",
                               wilcox.test(mda_pcor_pos_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(mda_pcor_pos_betweenness) - mean(dro_net_betweenness),
                               mean(mda_pcor_pos_betweenness))
betweenness_u.test[12,] <- list("mda_pcor_neg",
                               wilcox.test(mda_pcor_neg_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(mda_pcor_neg_betweenness) - mean(dro_net_betweenness),
                               mean(mda_pcor_neg_betweenness))

betweenness_u.test[13,] <- list("mda_aov_pos",
                               wilcox.test(mda_infecion_pos_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(mda_infecion_pos_betweenness) - mean(dro_net_betweenness),
                               mean(mda_infecion_pos_betweenness))
betweenness_u.test[14,] <- list("mda_aov_neg",
                               wilcox.test(mda_infecion_neg_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(mda_infecion_neg_betweenness) - mean(dro_net_betweenness),
                               mean(mda_infecion_neg_betweenness))
betweenness_u.test[15,] <- list("mda_only_cell.infection",
                               wilcox.test(mda_exclusively_cell.infecion_betweenness, dro_net_betweenness, alternative = "greater")$p.value,
                               mean(mda_exclusively_cell.infecion_betweenness) - mean(dro_net_betweenness),
                               mean(mda_exclusively_cell.infecion_betweenness))
