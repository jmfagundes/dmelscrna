#!/usr/bin/Rscript

# script to generate supplementary files

# make supplementary file lists

sup_tb <- data.frame(gene = character(0),
                     symbol = character(0),
                     
                     TV_pcor_estimate = numeric(0),
                     TV_pcor_p.adj = numeric(0),
                     TV_aov_cell_p.adj = numeric(0),
                     TV_aov_infection_p.adj = numeric(0),
                     TV_aov_infection_mean_logFC = numeric(0),
                     TV_aov_cell.infection_p.adj = numeric(0),
                     
                     DMelNV_pcor_estimate = numeric(0),
                     DMelNV_pcor_p.adj = numeric(0),
                     DMelNV_aov_cell_p.adj = numeric(0),
                     DMelNV_aov_infection_p.adj = numeric(0),
                     DMelNV_aov_infection_mean_logFC = numeric(0),
                     DMelNV_aov_cell.infection_p.adj = numeric(0),
                     stringsAsFactors = FALSE)

for (i in 1:1460) {
  gene <- row.names(thika_pcor[i, ])
  gene_ <- gsub("-", "_", row.names(thika_pcor[i, ]))
  symbol <- gene_symbol[gene_symbol$V1 == gene_, 3]
  
  TV_pcor_estimate <- thika_pcor[gene, ]$estimate
  TV_pcor_p.adj <- thika_pcor[gene, ]$p.adj
  TV_aov_cell_p.adj <- thika_aov[thika_aov$group == "cell" & thika_aov$gene == gene, ]$p.adj
  TV_aov_infection_p.adj <- thika_aov[thika_aov$group == "infection" & thika_aov$gene == gene, ]$p.adj
  TV_aov_cell.infection_p.adj <- thika_aov[thika_aov$group == "cell:infection" & thika_aov$gene == gene, ]$p.adj
  
  TV_aov_infection_mean_logFC <- mean_logFC_thika_infection[mean_logFC_thika_infection$gene == gene_, ]$mean_logFC
  if (length(TV_aov_infection_mean_logFC) == 0) {
    TV_aov_infection_mean_logFC <- NA
  }
  
  DMelNV_pcor_estimate <- nora_pcor[gene, ]$estimate
  DMelNV_pcor_p.adj <- nora_pcor[gene, ]$p.adj
  DMelNV_aov_cell_p.adj <- nora_aov[nora_aov$group == "cell" & nora_aov$gene == gene, ]$p.adj
  DMelNV_aov_infection_p.adj <- nora_aov[nora_aov$group == "infection" & nora_aov$gene == gene, ]$p.adj
  DMelNV_aov_cell.infection_p.adj <- nora_aov[nora_aov$group == "cell:infection" & nora_aov$gene == gene, ]$p.adj
  
  DMelNV_aov_infection_mean_logFC <- mean_logFC_nora_infection[mean_logFC_nora_infection$gene == gene_, ]$mean_logFC
  if (length(DMelNV_aov_infection_mean_logFC) == 0) {
    DMelNV_aov_infection_mean_logFC <- NA
  }
  sup_tb[i, ] <- list(gene_,
                      symbol,
                      
                      TV_pcor_estimate,
                      TV_pcor_p.adj,
                      TV_aov_cell_p.adj,
                      TV_aov_infection_p.adj,
                      TV_aov_infection_mean_logFC,
                      TV_aov_cell.infection_p.adj,
                      
                      DMelNV_pcor_estimate,
                      DMelNV_pcor_p.adj,
                      DMelNV_aov_cell_p.adj,
                      DMelNV_aov_infection_p.adj,
                      DMelNV_aov_infection_mean_logFC,
                      DMelNV_aov_cell.infection_p.adj)
}

mda_tb <- data.frame(gene = character(0),
                     symbol = character(0),
                     
                     DMelNV_pcor_estimate = numeric(0),
                     DMelNV_pcor_p.adj = numeric(0),
                     DMelNV_aov_cell_p.adj = numeric(0),
                     DMelNV_aov_infection_p.adj = numeric(0),
                     DMelNV_aov_infection_mean_logFC = numeric(0),
                     DMelNV_aov_cell.infection_p.adj = numeric(0),
                     stringsAsFactors = FALSE)

for (i in 1:271) {
  gene <- row.names(mda_pcor[i, ])
  gene_ <- gsub("-", "_", row.names(mda_pcor[i, ]))
  symbol <- gene_symbol[gene_symbol$V1 == gene_, 3]
  
  DMelNV_pcor_estimate <- mda_pcor[gene, ]$estimate
  DMelNV_pcor_p.adj <- mda_pcor[gene, ]$p.adj
  DMelNV_aov_cell_p.adj <- mda_aov[mda_aov$group == "cell" & mda_aov$gene == gene, ]$p.adj
  DMelNV_aov_infection_p.adj <- mda_aov[mda_aov$group == "infection" & mda_aov$gene == gene, ]$p.adj
  DMelNV_aov_cell.infection_p.adj <- mda_aov[mda_aov$group == "cell:infection" & mda_aov$gene == gene, ]$p.adj
  
  DMelNV_aov_infection_mean_logFC <- mean_logFC_mda_infection[mean_logFC_mda_infection$gene == gene_, ]$mean_logFC
  if (length(DMelNV_aov_infection_mean_logFC) == 0) {
    DMelNV_aov_infection_mean_logFC <- NA
  }
  
  mda_tb[i, ] <- list(gene_,
                      symbol,
                      
                      DMelNV_pcor_estimate,
                      DMelNV_pcor_p.adj,
                      DMelNV_aov_cell_p.adj,
                      DMelNV_aov_infection_p.adj,
                      DMelNV_aov_infection_mean_logFC,
                      DMelNV_aov_cell.infection_p.adj)
}

write.table(file = "Supplementary_file_1.tsv", nora_common, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(file = "Supplementary_file_2.tsv", sup_tb, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(file = "Supplementary_file_3.tsv", mda_tb, row.names = FALSE, quote = FALSE, sep = "\t")

# Supplementary figs

# marker genes

# EE

EE_markers <- c("Dmel-CG1147", "Dmel-CG4587", "Dmel-CG43207",
                "Dmel-CG10342",
                "Dmel-CG14375", "Dmel-CG5399", "Dmel-CG10621",
                "Dmel-CG15169", "Dmel-CG7191",
                "Dmel-CG4190", "Dmel-CG6489",
                "Dmel-CG42826", "Dmel-CG13229",
                "Dmel-CG13633",
                "Dmel-CG13094", "Dmel-CG8774",
                "Dmel-CG30340", "Dmel-CG33639",
                "Dmel-CG7266",
                "Dmel-CG12249")

EE_markers_symbol <- c()

for (i in EE_markers) {
  EE_markers_symbol[length(EE_markers_symbol) + 1] <- gene_symbol[gene_symbol$V1 == gsub("-", "_", i), 3]
}

EE_markers_symbol <- rev(EE_markers_symbol)

EE_markers_heatmap <- DoHeatmap(dmel_filtered, features = EE_markers, assay = "SCT", slot = "data", size = 4, angle = 90) +
  scale_y_discrete(labels = EE_markers_symbol) + 
  theme(plot.margin = unit(c(1.2,0.5,0.5,0.5), units = "cm"), legend.position = "bottom", text = element_text(size = 9.2)) + 
  guides(color = FALSE) + scale_fill_viridis()

# midgut atlas

mda_markers <- c("Dmel-CG31956", # cardia
                 "Dmel-CG3758", # ISC\EB
                 
                 "Dmel-CG7939", # some ribosomal proteins for ISC/EB and dEC
                 "Dmel-CG6779",
                 "Dmel-CG10944",
                 "Dmel-CG17596", 
                 
                 "Dmel-CG17228", # EEs
                 "Dmel-CG13633", # AstA-EE    
                 "Dmel-CG10342", # NPF-EE     
                 "Dmel-CG14919", # AstC-EE
                 
                 "Dmel-CG13315", # dEC        
                 
                 "Dmel-CG6295", # aEC 1
                 "Dmel-CG18211", # aEC 2 (1 and 3)
                 "Dmel-CG6164", # aEC 3       
                 "Dmel-CG18730", # aEC 4 
                 
                 "Dmel-CG3161", # mEC
                 
                 "Dmel-CG5097", # copper/iron
                 
                 "Dmel-CG5932", # LFC
                 
                 "Dmel-CG12350", # pEC 1, 2 and 3
                 "Dmel-CG9468", # pEC 1
                 "Dmel-CG1743", # pEC 2
                 "Dmel-CG31901", # pEC 3
                 
                 "Dmel-CG1462", # EC-like 1
                 "Dmel-CG12374", # EC-like 2 - low UMI
                 "Dmel-CR45846", # EC-like 3 - low UMI
                 "Dmel-CG11390") # others

mda_markers_symbol <- c()

for (i in mda_markers) {
  mda_markers_symbol[length(mda_markers_symbol) + 1] <- gene_symbol[gene_symbol$V1 == gsub("-", "_", i), 3]
}

mda_markers_symbol <- rev(mda_markers_symbol)

mda_markers_heatmap <- DoHeatmap(mda_10x_merged, features = mda_markers, assay = "SCT", slot = "data", size = 4, angle = 90) +
  scale_y_discrete(labels = mda_markers_symbol) + 
  theme(plot.margin = unit(c(0.5,0.5,0,0.5), units = "cm"), legend.position = "bottom", text = element_text(size = 9.2)) + 
  guides(color = FALSE) + scale_fill_viridis()

g <- ggarrange(EE_markers_heatmap, mda_markers_heatmap, ncol = 1, common.legend = TRUE, legend = "bottom")

ggsave("figures_out/figS1_marker_genes.eps", g, width = 6.85, height = 8.21)

# replication and infection dynamics for the midgut atlas dataset

mda_clusters_plot <- DimPlot(mda_10x_merged, reduction = "tsne", label = TRUE, pt.size = 0.5, repel = TRUE) + NoLegend() +
  ylab(label = element_blank()) + xlab(label = element_blank())

mda_nora_featplot <- FeaturePlot(mda_10x_nora, features = c("Dmel-nora-virus"), label = TRUE, reduction = "tsne", repel = TRUE) + 
  labs(title = "DMelNV") +
  ylab(label = element_blank()) + xlab(label = element_blank())

mda_pct_plot <- ggplot(mda_nora_pct, aes(x=cluster, y=nora_pct)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16) + 
  xlab(label = "") + ylab(label = "percentage of viral RNA") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

mda_exp_plot <- ggplot(mda_nora_pct, aes(x=cluster, y=nora_exp)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(shape=16) + 
  xlab(label = "") + ylab(label = "log-expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")


prn_plot <- prn %>% gather(variable, value, -cell)
prn_plot$cell <- factor(prn_plot$cell, levels = mda_cluster_order)
prn_plot$variable <- factor(prn_plot$variable, levels = c("nora", "uninfected"))
mda_cells_barplot <- ggplot(prn_plot, aes(fill = variable, y = value, x = cell)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis(discrete = TRUE, name = "", 
                     labels = c("DMelNV infected", "uninfected")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "top") + 
  ylab(label = "number of cells") + xlab(label = element_blank())


prn_pct_plot <- data.frame(cell = prn_inf$cell, value = prn_inf$nora/(prn_inf$nora + prn_inf$uninfected))
prn_pct_plot$cell <- factor(prn_pct_plot$cell, levels = c(mda_cluster_order, "total"))
mda_prob_inf <- ggplot(prn_pct_plot, aes(y = value, x = cell, group = 1)) +
  geom_line() +
  geom_point() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_blank()) +
  ylab(label = "probability of infection") + xlab(label = element_blank())


lay <- rbind(c(1, 1),
             c(1, 1),
             c(1, 1),
             c(1, 1),
             c(2, 3),
             c(2, 3),
             c(2, 3),
             c(4, 5),
             c(4, 5),
             c(4, 5))

g <- arrangeGrob(mda_clusters_plot, mda_nora_featplot,
                 ncol = 2, left = "tSNE_2", bottom = "tSNE_1")

g <- arrangeGrob(g, mda_pct_plot, mda_exp_plot,
                 mda_cells_barplot, mda_prob_inf,
                 layout_matrix = lay)

ggsave("figures_out/figS2_midgut_atlas_feat_plots.eps", g, width = 174, height = 234, units = "mm")

# midgut atlas GO

g <- arrangeGrob(bingo_bubble_pcor_mda, bingo_bubble_aov_mda, nrow = 2)

ggsave("figures_out/figS3_midgut_atlas_GO.eps", g, width = 6.85, height = 9.21)

# DEGs

do_aov_graph_loop <- function(obj_vir,
                              gene_lst,
                              vir_acro,
                              nrow = 6) {
  
  plist <- list()
  for (i in 1:length(gene_lst)) {
    gene <- gene_lst[i]
    gene_exp_infected <- data.frame(exp = obj_vir@assays$SCT@data[gsub("_", "-", gene),],
                                    cell = obj_vir@active.ident)
    gene_exp_infected$infection <- replicate(nrow(gene_exp_infected), paste(vir_acro, "infected"))
    
    gene_exp_uninfected <- data.frame(exp = dmel_no_vir@assays$SCT@data[gsub("_", "-", gene),],
                                      cell = dmel_no_vir@active.ident)
    gene_exp_uninfected$infection <- replicate(nrow(gene_exp_uninfected), "uninfected")
    
    exp_tb <- rbind(gene_exp_infected, gene_exp_uninfected)
    exp_tb$cell <- factor(exp_tb$cell, levels = cell_order)
    
    symbol <- gene_symbol[gene_symbol$V1 == gene, 3]
    p <- ggline(exp_tb, x = "cell", y = "exp", color = "infection", add = c("mean_se"))
    p <- ggpar(p, title = symbol, x.text.angle = 45, 
               ylab = FALSE, xlab = FALSE, legend.title = "", 
               font.main = c(8, "plain", "black"), font.xtickslab = c(8, "plain", "black"), font.ytickslab = c(8, "plain", "black"))
    plist[[length(plist) + 1]] <- p
  }
  g <- do.call("ggarrange", c(plist, ncol = 4, nrow = nrow, common.legend = TRUE, legend = "top"))
  return(g)
}

g <- do_aov_graph_loop(only_thika, thika_cell.infection_genes, "TV")
ggexport(filename = "figures_out/thika_cell:infection.pdf", g, width = 6.85, height = 9.21)

g <- do_aov_graph_loop(only_nora, nora_cell.infection_genes, "DMelNV")
ggexport(filename = "figures_out/nora_cell:infection.pdf", g, width = 6.85, height = 9.21)

g <- do_aov_graph_loop(only_thika, thika_infection_genes, "TV")
ggexport(filename = "figures_out/thika_infection.pdf", g, width = 6.85, height = 9.21)

g <- do_aov_graph_loop(only_nora, nora_infection_genes, "DMelNV")
ggexport(filename = "figures_out/nora_infection.pdf", g, width = 6.85, height = 9.21)

# top 12

tg <- do_aov_graph_loop(only_thika, top12_thika_cell.infection_genes, "TV", 3)
ng <- do_aov_graph_loop(only_nora, top12_nora_cell.infection_genes, "DMelNV", 3)

g <- arrangeGrob(tg, ng, nrow = 2)
ggsave(filename = "figures_out/figS4_top12_cell:infection.eps", g, width = 6.85, height = 9.21)

tg <- do_aov_graph_loop(only_thika, top12_thika_infection_genes, "TV", 3)
ng <- do_aov_graph_loop(only_nora, top12_nora_infection_genes, "DMelNV", 3)

g <- arrangeGrob(ng, tg, nrow = 2)
ggsave(filename = "figures_out/figS5_top12_infection.eps", g, width = 6.85, height = 9.21)

# midgut atlas + bulk RNA-seq network figure

lay <- rbind(c(1,1,1,2,2,2),
             c(3,3,4,4,5,5),
             c(6,6,6,7,7,7))

pdf("figures_out/figS6_network_analysis.pdf", width = 6.85, height = 9.21)
layout(lay)

# partial correlations

hub_fun(subgraph_mda_pcor_pos, "DMelNV midgut atlas\npositively correlated")
hub_fun(subgraph_mda_pcor_neg, "DMelNV midgut atlas\nnegatively correlated")

# DEGs

hub_fun(subgraph_mda_aov_pos, "DMelNV midgut atlas\nup-regulated")
hub_fun(subgraph_mda_aov_neg, "DMelNV midgut atlas\ndown-regulated")
hub_fun(subgraph_mda_exclusively_cell.infection, "DMelNV midgut atlas\ncell-subtype-specific")

# bulk RNA-seq

hub_fun(subgraph_nora_bulk_neg, "DMelNV bulk RNA-seq\ndown-regulated")
hub_fun(subgraph_nora_bulk_pos, "DMelNV bulk RNA-seq\nup-regulated")

dev.off()
