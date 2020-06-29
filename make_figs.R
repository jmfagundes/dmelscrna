#!/usr/bin/Rscript

# this script will generate figures
# it must only be run after run_analysis.R, GO.R and graph.R

# 0.3.2
# merging some figures

# Fig 1
# this figure was latter edited in order to include a scheme of the fruit fly midgut

clusters_plot <- DimPlot(dmel_filtered, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend() +
  ylab(label = element_blank()) + xlab(label = element_blank())

nora_featplot <- FeaturePlot(nora, features = c("Dmel-nora-virus"), label = TRUE) + labs(title = "DMelNV") +
  ylab(label = element_blank()) + xlab(label = element_blank()) + xlim(-45, 50) + ylim(-45, 50)

dmelc_featplot <- FeaturePlot(dmelc, features = c("Dmel-C-virus"), label = TRUE) + labs(title = "DCV") +
  ylab(label = element_blank()) + xlab(label = element_blank()) + xlim(-45, 55) + ylim(-45, 50)

thika_featplot <- FeaturePlot(thika, features = c("Thika-virus"), label = TRUE) + labs(title = "TV") +
  ylab(label = element_blank()) + xlab(label = element_blank()) + xlim(-45, 50) + ylim(-45, 50)

g <-arrangeGrob(clusters_plot, dmelc_featplot, thika_featplot, nora_featplot,
                ncol = 2, left = "tSNE_2", bottom = "tSNE_1")
ggsave("figures_out/fig1_feat_plots.eps", g, width = 174, height = 174, units = "mm")

# Fig 2

t_pct_plot <- ggplot(thika_pct, aes(x = cluster, y = thika_pct)) +
                 geom_boxplot(outlier.shape = NA) + geom_jitter(shape = 16) + 
                 xlab(label = "") + ylab(label = "percentage of viral RNA") + ggtitle(label = "TV") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

t_exp_plot <- ggplot(thika_pct, aes(x = cluster, y = thika_exp)) +
                 geom_boxplot(outlier.shape = NA) + geom_jitter(shape = 16) + 
                 xlab(label = "") + ylab(label = "log-expression") + ggtitle(label = "TV") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

n_pct_plot <- ggplot(nora_pct, aes(x = cluster, y = nora_pct)) +
                 geom_boxplot(outlier.shape = NA) + geom_jitter(shape = 16) + 
                 xlab(label = "") + ylab(label = "") + ggtitle(label = "DMelNV") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

n_exp_plot <- ggplot(nora_pct, aes(x = cluster, y = nora_exp)) +
                 geom_boxplot(outlier.shape = NA) + geom_jitter(shape = 16) + 
                 xlab(label = "") + ylab(label = "") + ggtitle(label = "DMelNV") +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

prc_plot <- prc %>% gather(variable, value, -cell)
prc_plot$cell <- factor(prc_plot$cell, levels = cell_order)
prc_plot$variable <- factor(prc_plot$variable, levels = c("thika", "nora", "coinfected", "uninfected"))
cells_barplot <- ggplot(prc_plot, aes(fill = variable, y = value, x = cell)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis(discrete = TRUE, name = "", 
                     labels = c("TV infected", "DMelNV infected", "coinfected", "uninfected")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "top") + guides(fill = guide_legend(nrow = 2,byrow = TRUE)) +
  ylab(label = "number of cells") + xlab(label = element_blank())

pco_plot <- p_co[1:3] %>% gather(variable, value, -cell)
pco_plot$cell <- factor(pco_plot$cell, levels = c(cell_order, "total"))
pco_plot$variable <- factor(pco_plot$variable, levels = c("thika", "nora"))
prob_inf <- ggplot(pco_plot, aes(group = variable, y = value, x = cell)) +
  geom_line(aes(color = variable)) +
  geom_point(aes(color = variable)) +
  scale_color_discrete(name = "", labels = c("TV", "DMelNV")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_blank()) +
  ylab(label = "probability of infection") + xlab(label = element_blank())

g <- ggarrange(t_pct_plot, n_pct_plot, 
               t_exp_plot, n_exp_plot,
               cells_barplot, prob_inf,
               ncol = 2, nrow = 3)
ggsave("figures_out/fig2.eps", g, width = 174, height = 210, units = "mm")


# Fig 3

ggsave("figures_out/fig3_GO_bubble_nora_pcor.pdf", bingo_bubble_pcor_nora, width = 6.85, height = 6)


# Fig 4

lay <- rbind(c(1),
             c(1),
             c(1),
             c(1),
             c(1),
             c(2),
             c(2),
             c(2))

g <- arrangeGrob(bingo_bubble_aov_thika + theme(text = element_text(size = 9.2)) + ggtitle("TV"),
                 bingo_bubble_aov_nora + theme(text = element_text(size = 9.2)) + ggtitle("DMelNV"),
                 layout_matrix = lay)

ggsave("figures_out/fig4_GO_bubble_aov.eps", g, width = 6.85, height = 9.21)

# Fig 5

# load results

load("sensitivity_curve/mda_sens_tb.Rdata")
mda_sens.tb <- sens.tb
load("sensitivity_curve/sens_tb.Rdata")

sens.tb_plot <- sens.tb[1:3] %>% gather(variable, value, -step)
mda_sens.tb_plot <- mda_sens.tb[1:3] %>% gather(variable, value, -step)

# dividing by the number of genes on the genes list 

sens.tb_plot$value <- sens.tb_plot$value/length(genes)
mda_sens.tb_plot$value <- mda_sens.tb_plot$value/length(mda_genes)

options(scipen = 100)

s <- ggplot(sens.tb_plot, aes(group = variable, y = value, x = step)) +
  geom_line(aes(color = variable)) +
  geom_point(aes(color = variable)) +
  scale_color_discrete(name = "", labels = c("cell:infection", "infection")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_blank()) +
  ylab(label = "number of false positives (%)") + xlab(label = "number of uninfected cells assigned as infected")

m <- ggplot(mda_sens.tb_plot, aes(group = variable, y = value, x = step)) +
  geom_line(aes(color = variable)) +
  geom_point(aes(color = variable)) +
  scale_color_discrete(name = "", labels = c("cell:infection", "infection")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_blank()) +
  ylab(label = "number of false positives (%)") + xlab(label = "number of uninfected cells assigned as infected")

g <- arrangeGrob(s, m, nrow = 2)

ggsave("figures_out/fig5_sensitivity_curve.eps", g, width = 174, height = 160, units = "mm")

# Fig 6

lay <- rbind(c(1,1,1,2,2,2),
             c(3,3,4,4,5,5),
             c(6,6,7,7,8,8))

pdf("figures_out/fig6_hub_analysis.pdf", width = 6.85, height = 9.21)
layout(lay)

# DMelNV partial correlations

hub_fun(subgraph_nora_pcor_pos, "DMelNV positively correlated genes")
hub_fun(subgraph_nora_pcor_neg, "DMelNV negatively correlated genes")

## TV DEGs

hub_fun(subgraph_thika_aov_pos, "TV generic response\nup-regulated genes")
hub_fun(subgraph_thika_aov_neg, "TV generic response\ndown-regulated genes")
hub_fun(subgraph_thika_exclusively_cell.infection, "TV cell-subtype-specific\nresponse genes")

## DMelNV DEGs

hub_fun(subgraph_nora_aov_pos, "DMelNV generic response\nup-regulated genes")
hub_fun(subgraph_nora_aov_neg, "DMelNV generic response\ndown-regulated genes")
hub_fun(subgraph_nora_exclusively_cell.infection,  "DMelNV cell-subtype-specific\nresponse genes")

dev.off()
