source("source_me.R")

########
## EE ##
########

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

# cluster               0       1       2         3       4       5        6           7         8       9        10      11     12      13
new.cluster.ids <- c("I-m", "II-m1", "I-ap-p", "II-a", "I-m", "II-m2", "I-pAstA", "I-pCCHa1", "II-p", "II-m1", "I-ap-a", "I-a", "III", "EEP")
names(new.cluster.ids) <- levels(dmel_filtered)
dmel_filtered <- RenameIdents(dmel_filtered, new.cluster.ids)

cluster_names = list("I-m", "II-m1", "I-ap-p", "II-a", "II-m2", "I-pCCHa1", "I-pAstA", "II-p", "I-a", "III", "EEP", "I-ap-a")

cell_order <- c("EEP", # EE progenitor cell
                "I-a", "III", "I-ap-a", "II-a", # anterior
                "I-m", "II-m1",  "II-m2", # medium
                "II-p", "I-ap-p", "I-pCCHa1", "I-pAstA") # posterior

# add assay without virus counts and log-normalize it by cell size

dmel_filtered[["NOVIR"]] <- CreateAssayObject(counts = dmel_filtered@assays$RNA@counts[!rownames(dmel_filtered@assays$RNA@counts) %in% c("Dmel-nora-virus", "Thika-virus", "Dmel-C-virus"),])
dmel_filtered <- NormalizeData(dmel_filtered, assay = "NOVIR")

##################
## midgut atlas ##
##################

# InDrop data has only R1 available

# load expression matrices

# sctransforms fails with seurat 3.2.2
# load mda_10x_merged data

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
# without integration cells did not cluster based on original sample, but integrate datasets nonetheless

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
                         "ISC/EB", # 11
                         "EE", # 12
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

# add assay without virus counts and log-normalize it by cell size

mda_10x_merged[["NOVIR"]] <- CreateAssayObject(counts = mda_10x_merged@assays$RNA@counts[!rownames(mda_10x_merged@assays$RNA@counts) %in% c("Dmel-nora-virus", "Thika-virus", "Dmel-C-virus"),])
mda_10x_merged <- NormalizeData(mda_10x_merged, assay = "NOVIR")

# estimate probability of cell being infected

# get raw cellranger output and subset cells that were filtered out

dmel_raw <- CreateSeuratObject(counts = Read10X("gene_matrices/raw_matrices/EE_raw_feature_bc_matrix/"), project = "EE_raw")

MA_raw_s1r1 <- CreateSeuratObject(counts = Read10X_h5("gene_matrices/raw_matrices/S1_rep1_10x_raw_feature_bc_matrix.h5"), project = "S1R1")
MA_raw_s1r2 <- CreateSeuratObject(counts = Read10X_h5("gene_matrices/raw_matrices/S1_rep2_10x_raw_feature_bc_matrix.h5"), project = "S1R2")
MA_raw_s2r1 <- CreateSeuratObject(counts = Read10X_h5("gene_matrices/raw_matrices/S2_rep1_10x_raw_feature_bc_matrix.h5"), project = "S2R1")
MA_raw_s2r2 <- CreateSeuratObject(counts = Read10X_h5("gene_matrices/raw_matrices/S2_rep2_10x_raw_feature_bc_matrix.h5"), project = "S2R2")

# add assay without virus counts

dmel_raw[["NOVIR"]] <- CreateAssayObject(counts = dmel_raw@assays$RNA@counts[!rownames(dmel_raw@assays$RNA@counts) %in% c("Dmel-nora-virus", "Thika-virus", "Dmel-C-virus"),])
MA_raw_s1r1[["NOVIR"]] <- CreateAssayObject(counts = MA_raw_s1r1@assays$RNA@counts[rownames(MA_raw_s1r1@assays$RNA@counts) != "Dmel-nora-virus",])
MA_raw_s1r2[["NOVIR"]] <- CreateAssayObject(counts = MA_raw_s1r2@assays$RNA@counts[rownames(MA_raw_s1r2@assays$RNA@counts) != "Dmel-nora-virus",])
MA_raw_s2r1[["NOVIR"]] <- CreateAssayObject(counts = MA_raw_s2r1@assays$RNA@counts[rownames(MA_raw_s2r1@assays$RNA@counts) != "Dmel-nora-virus",])
MA_raw_s2r2[["NOVIR"]] <- CreateAssayObject(counts = MA_raw_s2r2@assays$RNA@counts[rownames(MA_raw_s2r2@assays$RNA@counts) != "Dmel-nora-virus",])

dmel_EE_empty_cells <- subset(dmel_raw, cells = colnames(dmel_raw)[!colnames(dmel_raw) %in% colnames(dmel_filtered)], subset = nCount_NOVIR < 100)
dmel_MA_s1r1_empty_cells <- subset(MA_raw_s1r1, cells = colnames(MA_raw_s1r1)[!paste0("S1R1_", colnames(MA_raw_s1r1)) %in% colnames(dmel_filtered)], subset = nCount_NOVIR < 100)
dmel_MA_s1r2_empty_cells <- subset(MA_raw_s1r2, cells = colnames(MA_raw_s1r2)[!paste0("S1R2_", colnames(MA_raw_s1r2)) %in% colnames(dmel_filtered)], subset = nCount_NOVIR < 100)
dmel_MA_s2r1_empty_cells <- subset(MA_raw_s2r1, cells = colnames(MA_raw_s2r1)[!paste0("S2R1_", colnames(MA_raw_s2r1)) %in% colnames(dmel_filtered)], subset = nCount_NOVIR < 100)
dmel_MA_s2r2_empty_cells <- subset(MA_raw_s2r2, cells = colnames(MA_raw_s2r2)[!paste0("S2R2_", colnames(MA_raw_s2r2)) %in% colnames(dmel_filtered)], subset = nCount_NOVIR < 100)

# write h5 files without virus reads to estimate the proportion of genes in each cell that are due to contamination

write10xCounts("cellbender/dmel_raw_novir.h5", dmel_raw@assays$NOVIR@counts, barcodes = colnames(dmel_raw),
               gene.id = rownames(dmel_raw@assays$NOVIR@counts), gene.symbol = rownames(dmel_raw@assays$NOVIR@counts),
               gene.type = "Gene Expression", genome = "dmel_nv_dav")
write10xCounts("cellbender/MA_raw_s1r1_novir.h5", MA_raw_s1r1@assays$NOVIR@counts, barcodes = colnames(MA_raw_s1r1),
               gene.id = rownames(MA_raw_s1r1@assays$NOVIR@counts), gene.symbol = rownames(MA_raw_s1r1@assays$NOVIR@counts),
               gene.type = "Gene Expression", genome = "dmel_nv_dav")
write10xCounts("cellbender/MA_raw_s1r2_novir.h5", MA_raw_s1r2@assays$NOVIR@counts, barcodes = colnames(MA_raw_s1r2),
               gene.id = rownames(MA_raw_s1r2@assays$NOVIR@counts), gene.symbol = rownames(MA_raw_s1r2@assays$NOVIR@counts),
               gene.type = "Gene Expression", genome = "dmel_nv_dav")
write10xCounts("cellbender/MA_raw_s2r1_novir.h5", MA_raw_s2r1@assays$NOVIR@counts, barcodes = colnames(MA_raw_s2r1),
               gene.id = rownames(MA_raw_s2r1@assays$NOVIR@counts), gene.symbol = rownames(MA_raw_s2r1@assays$NOVIR@counts),
               gene.type = "Gene Expression", genome = "dmel_nv_dav")
write10xCounts("cellbender/MA_raw_s2r2_novir.h5", MA_raw_s2r2@assays$NOVIR@counts, barcodes = colnames(MA_raw_s2r2),
               gene.id = rownames(MA_raw_s2r2@assays$NOVIR@counts), gene.symbol = rownames(MA_raw_s2r2@assays$NOVIR@counts),
               gene.type = "Gene Expression", genome = "dmel_nv_dav")

# read cellbender .h5 output and add matrix to the corresponding object

# arguments used for cellbender

# cellbender remove-background --input dmel_raw_novir.h5 --output dmel_novir_cellbender.h5 --expected-cells 5000 --total-droplets-included 20000 --fpr 0.01 --epochs 150
# cellbender remove-background --input MA_raw_s1r1_novir.h5 --output MA_raw_s1r1_novir_cellbender.h5 --expected-cells 200 --total-droplets-included 5000 --fpr 0.01 --epochs 150 --low-count-threshold 5
# cellbender remove-background --input MA_raw_s1r2_novir.h5 --output MA_raw_s1r2_novir_cellbender.h5 --expected-cells 200 --total-droplets-included 5000 --fpr 0.01 --epochs 150 --low-count-threshold 5
# cellbender remove-background --input MA_raw_s2r1_novir.h5 --output MA_raw_s2r1_novir_cellbender.h5 --expected-cells 1500 --total-droplets-included 15000 --fpr 0.01 --epochs 160 --low-count-threshold 15 --learning-rate 0.5e-4
# cellbender remove-background --input MA_raw_s2r2_novir.h5 --output MA_raw_s2r2_novir_cellbender.h5 --expected-cells 1500 --total-droplets-included 15000 --fpr 0.01 --epochs 150 --low-count-threshold 15 --learning-rate 1e-5

# subset same cells called by cellranger and rearrange cell order

EE_cellbender <- Read10X_h5("cellbender/dmel_novir_cellbender_filtered.h5")
EE_cellbender <- EE_cellbender[,match(colnames(dmel_filtered), colnames(EE_cellbender))]
dmel_filtered[["CELLBENDER"]] <- CreateAssayObject(counts = EE_cellbender)

MA_s1r1_cellbender <- Read10X_h5("cellbender/MA_raw_s1r1_novir_cellbender_filtered.h5")
MA_s1r2_cellbender <- Read10X_h5("cellbender/MA_raw_s1r2_novir_cellbender_filtered.h5")
MA_s2r1_cellbender <- Read10X_h5("cellbender/MA_raw_s2r1_novir_cellbender_filtered.h5")
MA_s2r2_cellbender <- Read10X_h5("cellbender/MA_raw_s2r2_novir_cellbender_filtered.h5")

colnames(MA_s1r1_cellbender) <- paste0("S1R1_", colnames(MA_s1r1_cellbender))
colnames(MA_s1r2_cellbender) <- paste0("S1R2_", colnames(MA_s1r2_cellbender))
colnames(MA_s2r1_cellbender) <- paste0("S2R1_", colnames(MA_s2r1_cellbender))
colnames(MA_s2r2_cellbender) <- paste0("S2R2_", colnames(MA_s2r2_cellbender))

MA_cellbender <- cbind(MA_s1r1_cellbender,
                       MA_s1r2_cellbender,
                       MA_s2r1_cellbender,
                       MA_s2r2_cellbender)

MA_cellbender <- MA_cellbender[,match(colnames(mda_10x_merged), colnames(MA_cellbender))]
mda_10x_merged[["CELLBENDER"]] <- CreateAssayObject(counts = MA_cellbender)

# percentage of viral reads in ambient RNA

pct.nora_EE_ambient <- sum(dmel_EE_empty_cells@assays$RNA@counts["Dmel-nora-virus",])/sum(dmel_EE_empty_cells@assays$RNA@counts)

pct.nora_MA_ambient_s1r1 <- sum(dmel_MA_s1r1_empty_cells@assays$RNA@counts["Dmel-nora-virus",])/sum(dmel_MA_s1r1_empty_cells@assays$RNA@counts)
pct.nora_MA_ambient_s1r2 <- sum(dmel_MA_s1r2_empty_cells@assays$RNA@counts["Dmel-nora-virus",])/sum(dmel_MA_s1r2_empty_cells@assays$RNA@counts)
pct.nora_MA_ambient_s2r1 <- sum(dmel_MA_s2r1_empty_cells@assays$RNA@counts["Dmel-nora-virus",])/sum(dmel_MA_s2r1_empty_cells@assays$RNA@counts)
pct.nora_MA_ambient_s2r2 <- sum(dmel_MA_s2r2_empty_cells@assays$RNA@counts["Dmel-nora-virus",])/sum(dmel_MA_s2r2_empty_cells@assays$RNA@counts)

# no TV reads found in ambient RNA, infected cells will be called if TV reads > 1

dmel_filtered$is.tv.infected <- dmel_filtered@assays$RNA@counts["Thika-virus",] > 1
dmel_filtered$is.dmelnv.infected <- is.infected(dmel_filtered, pct.nora_EE_ambient)

MA_is.dmelnv.infected_s1r1 <- is.infected(subset(mda_10x_merged, cells = colnames(mda_10x_merged)[mda_10x_merged$orig.ident == "mda_s1_r1"]),
                                          pct.nora_MA_ambient_s1r1)
MA_is.dmelnv.infected_s1r2 <- is.infected(subset(mda_10x_merged, cells = colnames(mda_10x_merged)[mda_10x_merged$orig.ident == "mda_s1_r2"]),
                                          pct.nora_MA_ambient_s1r2)
MA_is.dmelnv.infected_s2r1 <- is.infected(subset(mda_10x_merged, cells = colnames(mda_10x_merged)[mda_10x_merged$orig.ident == "mda_s2_r1"]),
                                          pct.nora_MA_ambient_s2r1)
MA_is.dmelnv.infected_s2r2 <- is.infected(subset(mda_10x_merged, cells = colnames(mda_10x_merged)[mda_10x_merged$orig.ident == "mda_s2_r2"]),
                                          pct.nora_MA_ambient_s2r2)

MA_is.dmelnv.infected <- c(MA_is.dmelnv.infected_s1r1,
                           MA_is.dmelnv.infected_s1r2,
                           MA_is.dmelnv.infected_s2r1,
                           MA_is.dmelnv.infected_s2r2)

MA_is.dmelnv.infected <- MA_is.dmelnv.infected[match(colnames(mda_10x_merged), names(MA_is.dmelnv.infected))]

mda_10x_merged$is.dmelnv.infected <- MA_is.dmelnv.infected

# save objs

save(dmel_filtered, file = "objs/dmel_filtered.Rdata")
save(mda_10x_merged, file = "objs/mda_10x_merged.Rdata")

# plot clusters

clusters_plot <- DimPlot(dmel_filtered, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend() +
  ylab(label = element_blank()) + xlab(label = element_blank())

MA_clusters_plot <- DimPlot(mda_10x_merged, reduction = "tsne", label = TRUE, pt.size = 0.5, repel = TRUE) + NoLegend() +
  ylab(label = element_blank()) + xlab(label = element_blank())

cell_clusters <- arrangeGrob(clusters_plot, MA_clusters_plot,
                             ncol = 2, left = "tSNE_2", bottom = "tSNE_1")

ggsave("figures_out/fig1_unedited.pdf", cell_clusters, width = 174, height = 85, units = "mm")

# plot markers

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

ggsave("figures_out/figS1.pdf", g, width = 6.85, height = 8.21)
