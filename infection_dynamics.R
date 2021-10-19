source("source_me.R")

load("objs/dmel_filtered.Rdata")
load("objs/mda_10x_merged.Rdata")

# percentage and log-percentage per cell

virus_pct_per_cell <- data.frame(dataset = c(rep("EE", ncol(dmel_filtered)), rep("MA", ncol(mda_10x_merged))),
                                 cell = c(dmel_filtered@active.ident %>% as.character(), mda_10x_merged@active.ident %>% as.character()),
                                 tv.infection = c(dmel_filtered$is.tv.infected, rep(FALSE, ncol(mda_10x_merged))),
                                 dmelnv.infection = c(dmel_filtered$is.dmelnv.infected, mda_10x_merged$is.dmelnv.infected),
                                 tv.pct = c(dmel_filtered$percent.thika, rep(0, ncol(mda_10x_merged))),
                                 dmelnv.pct = c(dmel_filtered$percent.nora, mda_10x_merged$percent.nora))

# log1p expression, scale = 10,000

virus_pct_per_cell$tv.exp <- log1p(virus_pct_per_cell$tv.pct*100)
virus_pct_per_cell$dmelnv.exp <- log1p(virus_pct_per_cell$dmelnv.pct*100)

EE_cell_order <- c("EEP", # EE progenitor cell
                   "I-a", "III", "I-ap-a", "II-a", # anterior
                   "I-m", "II-m1",  "II-m2", # medium
                   "II-p", "I-ap-p", "I-pCCHa1", "I-pAstA") # posterior

MA_cell_order <- c("ISC/EB", # progenitor
                   "EE", # scattered
                   "unk", # unknown location
                   "cardia", # proventriculus
                   "aEC", # anterior
                   "mEC", "copper/iron", # medium
                   "LFC", # posterior medium region
                   "pEC") # posterior 

virus_pct_per_cell$cell <- factor(virus_pct_per_cell$cell, levels = c(EE_cell_order, MA_cell_order))

gg_tv.pct <- ggplot(virus_pct_per_cell[virus_pct_per_cell$tv.infection,], aes(x = cell, y = tv.pct)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.3, shape = 16) + 
  xlab(label = "") + ylab(label = "percentage of viral counts") + ggtitle(label = "TV") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

gg_tv.exp <- ggplot(virus_pct_per_cell[virus_pct_per_cell$tv.infection,], aes(x = cell, y = tv.exp)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.3, shape = 16) + 
  xlab(label = "") + ylab(label = "log-normalized viral counts") + ggtitle(label = "TV") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

gg_dv.pct <- ggplot(virus_pct_per_cell[virus_pct_per_cell$dmelnv.infection,], aes(x = cell, y = dmelnv.pct)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.3, shape = 16) + 
  xlab(label = "") + ylab(label = "percentage of viral counts") + ggtitle(label = "DMelNV") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

gg_dv.exp <- ggplot(virus_pct_per_cell[virus_pct_per_cell$dmelnv.infection,], aes(x = cell, y = dmelnv.exp)) +
  geom_boxplot(outlier.shape = NA) + geom_jitter(alpha = 0.3, shape = 16) + 
  xlab(label = "") + ylab(label = "log-normalized viral counts") + ggtitle(label = "DMelNV") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")

# one-way anova virus percentage and log expression between cell types
# do not include unk cells for DMelNV

thika_cluster_anova <- list(log = Anova(lm(tv.exp ~ cell, virus_pct_per_cell[virus_pct_per_cell$tv.infection,]), type = 2),
                            pct = Anova(lm(tv.pct ~ cell, virus_pct_per_cell[virus_pct_per_cell$tv.infection,]), type = 2))

nora_EE_cluster_anova <- list(log = Anova(lm(dmelnv.exp ~ cell, virus_pct_per_cell[virus_pct_per_cell$dmelnv.infection & virus_pct_per_cell$dataset == "EE",]), type = 2),
                              pct = Anova(lm(dmelnv.pct ~ cell, virus_pct_per_cell[virus_pct_per_cell$dmelnv.infection & virus_pct_per_cell$dataset == "EE",]), type = 2))

nora_MA_cluster_anova <- list(log = Anova(lm(dmelnv.exp ~ cell, virus_pct_per_cell[virus_pct_per_cell$dmelnv.infection & virus_pct_per_cell$dataset == "MA" & virus_pct_per_cell$cell != "unk",]), type = 2),
                              pct = Anova(lm(dmelnv.pct ~ cell, virus_pct_per_cell[virus_pct_per_cell$dmelnv.infection & virus_pct_per_cell$dataset == "MA" & virus_pct_per_cell$cell != "unk",]), type = 2))

# Tukey's HSD tests

tukey.tests <- list(tv.exp = HSD.test(lm(tv.exp ~ cell, virus_pct_per_cell[virus_pct_per_cell$tv.infection,]), 'cell', unbalanced = TRUE, alpha = 0.05),
                    dmelnv.EE.exp = HSD.test(lm(dmelnv.exp ~ cell, virus_pct_per_cell[virus_pct_per_cell$dmelnv.infection & virus_pct_per_cell$dataset == "EE",]), 'cell', unbalanced = TRUE, alpha = 0.05),
                    dmelnv.MA.exp = HSD.test(lm(dmelnv.exp ~ cell, virus_pct_per_cell[virus_pct_per_cell$dmelnv.infection & virus_pct_per_cell$dataset == "MA" & virus_pct_per_cell$cell != "unk",]), 'cell', unbalanced = TRUE, alpha = 0.05))

# infected cells

infect.stats <- data.frame(cell = levels(virus_pct_per_cell$cell),
                           tv.pct.infected_cells = ((table(virus_pct_per_cell[virus_pct_per_cell$tv.infection, "cell"]) / table(virus_pct_per_cell$cell))*100) %>% as.numeric(),
                           dmelnv.pct.infected_cells = ((table(virus_pct_per_cell[virus_pct_per_cell$dmelnv.infection, "cell"]) / table(virus_pct_per_cell$cell))*100) %>% as.numeric(),
                           n.uninfected.cells = (table(virus_pct_per_cell[!virus_pct_per_cell$tv.infection & !virus_pct_per_cell$dmelnv.infection, "cell"])) %>% as.numeric(),
                           n.coinfected.cells = (table(virus_pct_per_cell[virus_pct_per_cell$tv.infection & virus_pct_per_cell$dmelnv.infection, "cell"])) %>% as.numeric(),
                           tv.n.infected_cells = table(virus_pct_per_cell[virus_pct_per_cell$tv.infection & !virus_pct_per_cell$dmelnv.infection, "cell"]) %>% as.numeric(),
                           dmelnv.n.infected_cells = table(virus_pct_per_cell[virus_pct_per_cell$dmelnv.infection & !virus_pct_per_cell$tv.infection, "cell"]) %>% as.numeric())

# integrate EEs from both datasets

EE_mda <- subset(mda_10x_merged, idents = "EE")
DefaultAssay(dmel_filtered) <- "NOVIR"
DefaultAssay(EE_mda) <- "NOVIR"

EE.lst <- list(EE = FindVariableFeatures(dmel_filtered, nfeatures = 5000),
               mda = FindVariableFeatures(EE_mda, nfeatures = 5000))

EE_anchors <- FindIntegrationAnchors(object.list = EE.lst,
                                     dims = 1:25, k.filter = 50)

EE_mda_merged <- IntegrateData(EE_anchors, dims = 1:25)
EE_mda_merged$orig.cell.type <- EE_mda_merged@active.ident
EE_mda_merged <- ScaleData(EE_mda_merged)
EE_mda_merged <- RunPCA(EE_mda_merged, features = VariableFeatures(object = EE_mda_merged))
EE_mda_merged <- FindNeighbors(EE_mda_merged, dims = 1:25)
EE_mda_merged <- FindClusters(EE_mda_merged, resolution = 0.4)
EE_mda_merged <- RunTSNE(EE_mda_merged, dims = 1:25)

# reannotate

# cluster                 0       1       2         3        4         5         6         7          8       9      10     11     12
EE_new.cluster.ids <- c("I-m", "II-m1", "II-a", "I-ap-p", "II-m2", "I-pAstA", "II-p", "I-pCCHa1", "I-ap-a", "I-a", "I-m", "III", "EEP")
names(EE_new.cluster.ids) <- levels(EE_mda_merged)
EE_mda_merged <- RenameIdents(EE_mda_merged, EE_new.cluster.ids)
DimPlot(EE_mda_merged, reduction = "tsne", label = TRUE, pt.size = 0.5)

# check EE subtypes infected with DMelNV in MA dataset

EE_mda_merged@active.ident[colnames(EE_mda_merged) %in% colnames(mda_10x_merged)[mda_10x_merged$is.dmelnv.infected]]

## Figs
# number of infected cells barplot

gather_n.vir <- infect.stats[c("cell", "tv.n.infected_cells", "dmelnv.n.infected_cells", "n.uninfected.cells", "n.coinfected.cells")] %>% gather(variable, value, -cell)
gather_n.vir$cell <- factor(gather_n.vir$cell, levels = levels(virus_pct_per_cell$cell))
gather_n.vir$variable <- factor(gather_n.vir$variable, levels = c("tv.n.infected_cells", "dmelnv.n.infected_cells", "n.coinfected.cells", "n.uninfected.cells"))

cells_barplot <- ggplot(gather_n.vir, aes(fill = variable, y = value, x = cell)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis(discrete = TRUE, name = "", 
                     labels = c("TV-infected", "DMelNV-infected", "coinfected", "uninfected")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "top") + guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  ylab(label = "number of cells") + xlab(label = element_blank())

# percentage of infected cells

gather_pct.vir <- infect.stats[c("cell", "tv.pct.infected_cells", "dmelnv.pct.infected_cells")] %>% gather(variable, value, -cell)
gather_pct.vir[gather_pct.vir$value == 0, "value"] <- NA
gather_pct.vir$cell <- factor(gather_pct.vir$cell, levels = levels(virus_pct_per_cell$cell))
gather_pct.vir$variable <- factor(gather_pct.vir$variable, levels = c("tv.pct.infected_cells", "dmelnv.pct.infected_cells"))

prob_inf <- ggplot(gather_pct.vir, aes(group = variable, y = value, x = cell)) +
  geom_line(aes(color = variable)) +
  geom_point(aes(color = variable)) +
  scale_color_discrete(name = "", labels = c("TV", "DMelNV")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_blank()) +
  ylab(label = "percentage of infected cells") + xlab(label = element_blank())

g <- arrangeGrob(gg_tv.pct, gg_dv.pct,
                 gg_tv.exp, gg_dv.exp,
                 cells_barplot, prob_inf, nrow = 3)

ggsave("figures_out/fig2.pdf", g, width = 174, height = 234, units = "mm")

# correlation between mean expression and virus susceptibility
# see if correlated genes are DE between cell types

load("objs/glm.vir_glmGamPoi.Rdata")

EE.uninfected.seu <- subset(dmel_filtered, subset = percent.thika == 0 & percent.nora == 0 & percent.dmelc == 0,
                            idents = c("I-m", "II-m1", "I-ap-p", "II-a", "II-m2", "I-pAstA", "I-pCCHa1", "II-p", "I-ap-a", "I-a", "III"))

# log-normalize based on size factors obtained from glmGamPoi

EE.uninfected.seu[["CELLBENDER"]]@data <- t(glm.vir$EE$fit$size_factors[match(colnames(EE.uninfected.seu), names(glm.vir$EE$fit$size_factors))]*
                                              t(EE.uninfected.seu@assays$CELLBENDER@counts)) %>% log1p()

EE.mean.exp <- AverageExpression(EE.uninfected.seu, assays = "CELLBENDER", slot = "data")

# investigate whether infection can change cell clusters

tv.EE.exp <- virus_pct_per_cell[virus_pct_per_cell$dataset == "EE", "tv.exp"]
tv.EE.exp[!virus_pct_per_cell[virus_pct_per_cell$dataset == "EE", "tv.infection"]] <- 0

dmelnv.EE.exp <- virus_pct_per_cell[virus_pct_per_cell$dataset == "EE", "dmelnv.exp"]
dmelnv.EE.exp[!virus_pct_per_cell[virus_pct_per_cell$dataset == "EE", "dmelnv.infection"]] <- 0

dmelnv.MA.exp <- virus_pct_per_cell[virus_pct_per_cell$dataset == "MA", "dmelnv.exp"]
dmelnv.MA.exp[!virus_pct_per_cell[virus_pct_per_cell$dataset == "MA", "dmelnv.infection"]] <- 0

dmel_filtered <- AddMetaData(dmel_filtered, metadata = tv.EE.exp, col.name = "TV.log.normalized.counts")
dmel_filtered <- AddMetaData(dmel_filtered, metadata = dmelnv.EE.exp, col.name = "DMelNV.log.normalized.counts")
mda_10x_merged <- AddMetaData(mda_10x_merged, metadata = dmelnv.MA.exp, col.name = "DMelNV.log.normalized.counts")

gg.tv.exp <- FeaturePlot(dmel_filtered, features = "TV.log.normalized.counts", label = TRUE, label.size = 2) + ggtitle("A                       TV") +
  theme(plot.title = element_text(hjust = -3)) +
  ylab(label = element_blank()) + xlab(label = element_blank())

gg.dmelnv.EE.exp <- FeaturePlot(dmel_filtered, features = "DMelNV.log.normalized.counts", label = TRUE, label.size = 2) + ggtitle("DMelNV (EE)") +
  ylab(label = element_blank()) + xlab(label = element_blank())

gg.dmelnv.MA.exp <- FeaturePlot(mda_10x_merged, features = "DMelNV.log.normalized.counts", reduction = "tsne", label = TRUE, label.size = 2) + ggtitle("DMelNV (MA)") +
  ylab(label = element_blank()) + xlab(label = element_blank())

# test I-pAstA cluster

ipasta <- subset(dmel_filtered, ident = "I-pAstA")
ipasta@active.assay <- "CELLBENDER"

ipasta <- NormalizeData(ipasta)
ipasta <- ScaleData(ipasta, vars.to.regress = c("percent.mito", "percent.vir"))
ipasta <- FindVariableFeatures(ipasta)

ipasta <- RunPCA(ipasta, features = VariableFeatures(object = ipasta))

# cluster generation

ipasta <- FindNeighbors(ipasta, dims = 1:20)
ipasta <- FindClusters(ipasta, resolution = .5)

ipasta <- RunTSNE(ipasta, dims = 1:20)

gg.tv.exp.ipasta <- FeaturePlot(ipasta, features = "TV.log.normalized.counts", label = TRUE, pt.size = 0.5) + ggtitle("B                       TV") +
  theme(plot.title = element_text(hjust = -3)) +
  ylab(label = element_blank()) + xlab(label = element_blank())

gg.dmelnv.EE.exp.ipasta <- FeaturePlot(ipasta, features = "DMelNV.log.normalized.counts", label = TRUE, pt.size = 0.5) + ggtitle("DMelNV") +
  ylab(label = element_blank()) + xlab(label = element_blank())

ipasta.markers <- FindAllMarkers(ipasta, only.pos = TRUE)
ipasta.markers$symbol <- gene_symbol[match(gsub("-", "_", ipasta.markers$gene), gene_symbol$V1), "V3"]

ipasta.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC) %>% as.data.frame()

# save fig

gg.fig8a <- arrangeGrob(gg.tv.exp, gg.dmelnv.EE.exp, gg.dmelnv.MA.exp,
                        ncol = 2, left = "tSNE_2", bottom = "tSNE_1")

gg.fig8b <- arrangeGrob(gg.tv.exp.ipasta, gg.dmelnv.EE.exp.ipasta,
                        ncol = 2, left = "tSNE_2", bottom = "tSNE_1")

ggsave("figures_out/fig8.pdf", arrangeGrob(gg.fig8a, gg.fig8b,
                                           layout_matrix = rbind(1,
                                                                 1,
                                                                 2)),
                                           width = 174, height = 234, units = "mm")
