source("source_me.R")

load("objs/dmel_filtered.Rdata")
load("objs/mda_10x_merged.Rdata")

save_steps <- FALSE

# GLM fit

glm.vir <- DE_gGP(dmel_filtered, mda_10x_merged,
                  batch = FALSE)

# set contrasts

thika_contrasts <- makeContrasts(infection = (groupingEETV_infected.I_m +
                                                groupingEETV_infected.II_m1 +
                                                groupingEETV_infected.I_ap_p +
                                                groupingEETV_infected.II_a +
                                                groupingEETV_infected.II_m2 +
                                                groupingEETV_infected.I_pAstA +
                                                groupingEETV_infected.I_pCCHa1 +
                                                groupingEETV_infected.II_p +
                                                groupingEETV_infected.I_ap_a +
                                                groupingEETV_infected.I_a +
                                                groupingEETV_infected.III)/11 -
                                   (groupingEEuninfected.I_m +
                                      groupingEEuninfected.II_m1 +
                                      groupingEEuninfected.I_ap_p +
                                      groupingEEuninfected.II_a +
                                      groupingEEuninfected.II_m2 +
                                      groupingEEuninfected.I_pAstA +
                                      groupingEEuninfected.I_pCCHa1 +
                                      groupingEEuninfected.II_p +
                                      groupingEEuninfected.I_ap_a +
                                      groupingEEuninfected.I_a +
                                      groupingEEuninfected.III)/11,
                                 a_p = (groupingEETV_infected.I_ap_p +
                                          groupingEETV_infected.II_a +
                                          groupingEETV_infected.I_pAstA +
                                          groupingEETV_infected.I_pCCHa1 +
                                          groupingEETV_infected.II_p +
                                          groupingEETV_infected.I_ap_a +
                                          groupingEETV_infected.I_a +
                                          groupingEETV_infected.III)/8 -
                                   (groupingEEuninfected.I_ap_p +
                                      groupingEEuninfected.II_a +
                                      groupingEEuninfected.I_pAstA +
                                      groupingEEuninfected.I_pCCHa1 +
                                      groupingEEuninfected.II_p +
                                      groupingEEuninfected.I_ap_a +
                                      groupingEEuninfected.I_a +
                                      groupingEEuninfected.III)/8,
                                 a = (groupingEETV_infected.II_a +
                                        groupingEETV_infected.I_ap_a +
                                        groupingEETV_infected.I_a)/3 -
                                   (groupingEEuninfected.II_a +
                                      groupingEEuninfected.I_ap_a +
                                      groupingEEuninfected.I_a)/3,
                                 p = (groupingEETV_infected.I_ap_p +
                                        groupingEETV_infected.I_pAstA +
                                        groupingEETV_infected.I_pCCHa1 +
                                        groupingEETV_infected.II_p)/4 -
                                   (groupingEEuninfected.I_ap_p +
                                      groupingEEuninfected.I_pAstA +
                                      groupingEEuninfected.I_pCCHa1 +
                                      groupingEEuninfected.II_p)/4,
                                 I_m = groupingEETV_infected.I_m - groupingEEuninfected.I_m,
                                 II_m1 = groupingEETV_infected.II_m1 - groupingEEuninfected.II_m1,
                                 I_ap_p = groupingEETV_infected.I_ap_p - groupingEEuninfected.I_ap_p,
                                 II_a = groupingEETV_infected.II_a - groupingEEuninfected.II_a,
                                 II_m2 = groupingEETV_infected.II_m2 - groupingEEuninfected.II_m2,
                                 I_pAstA = groupingEETV_infected.I_pAstA - groupingEEuninfected.I_pAstA,
                                 I_pCCHa1 = groupingEETV_infected.I_pCCHa1 - groupingEEuninfected.I_pCCHa1,
                                 II_p = groupingEETV_infected.II_p - groupingEEuninfected.II_p,
                                 I_ap_a = groupingEETV_infected.I_ap_a - groupingEEuninfected.I_ap_a,
                                 I_a = groupingEETV_infected.I_a - groupingEEuninfected.I_a,
                                 III = groupingEETV_infected.III - groupingEEuninfected.III,
                                 levels = glm.vir$EE$design)

nora_EE_contrasts <- makeContrasts(I_m = groupingEEDMelNV_infected.I_m - groupingEEuninfected.I_m,
                                   II_m1 = groupingEEDMelNV_infected.II_m1 - groupingEEuninfected.II_m1,
                                   I_ap_p = groupingEEDMelNV_infected.I_ap_p - groupingEEuninfected.I_ap_p,
                                   II_a = groupingEEDMelNV_infected.II_a - groupingEEuninfected.II_a,
                                   II_m2 = groupingEEDMelNV_infected.II_m2 - groupingEEuninfected.II_m2,
                                   I_pAstA = groupingEEDMelNV_infected.I_pAstA - groupingEEuninfected.I_pAstA,
                                   I_pCCHa1 = groupingEEDMelNV_infected.I_pCCHa1 - groupingEEuninfected.I_pCCHa1,
                                   II_p = groupingEEDMelNV_infected.II_p - groupingEEuninfected.II_p,
                                   I_ap_a = groupingEEDMelNV_infected.I_ap_a - groupingEEuninfected.I_ap_a,
                                   EE = (groupingEEDMelNV_infected.I_m +
                                           groupingEEDMelNV_infected.II_m1 +
                                           groupingEEDMelNV_infected.I_ap_p +
                                           groupingEEDMelNV_infected.II_a +
                                           groupingEEDMelNV_infected.II_m2 +
                                           groupingEEDMelNV_infected.I_pAstA +
                                           groupingEEDMelNV_infected.I_pCCHa1 +
                                           groupingEEDMelNV_infected.II_p +
                                           groupingEEDMelNV_infected.I_ap_a)/9 -
                                     (groupingEEuninfected.I_m +
                                        groupingEEuninfected.II_m1 +
                                        groupingEEuninfected.I_ap_p +
                                        groupingEEuninfected.II_a +
                                        groupingEEuninfected.II_m2 +
                                        groupingEEuninfected.I_pAstA +
                                        groupingEEuninfected.I_pCCHa1 +
                                        groupingEEuninfected.II_p +
                                        groupingEEuninfected.I_ap_a)/9,
                                   levels = glm.vir$EE$design)

# do not include unk cells

nora_MA_contrasts <- makeContrasts(MA = (groupingMADMelNV_infected.EE +
                                           groupingMADMelNV_infected.pEC +
                                           groupingMADMelNV_infected.aEC +
                                           groupingMADMelNV_infected.cardia +
                                           groupingMADMelNV_infected.copperiron +
                                           groupingMADMelNV_infected.ISCEB +
                                           groupingMADMelNV_infected.LFC +
                                           groupingMADMelNV_infected.mEC)/8 -
                                     (groupingMAuninfected.EE +
                                        groupingMAuninfected.pEC +
                                        groupingMAuninfected.aEC +
                                        groupingMAuninfected.cardia +
                                        groupingMAuninfected.copperiron +
                                        groupingMAuninfected.ISCEB +
                                        groupingMAuninfected.LFC +
                                        groupingMAuninfected.mEC)/8,
                                   EE = groupingMADMelNV_infected.EE - groupingMAuninfected.EE,
                                   pEC = groupingMADMelNV_infected.pEC - groupingMAuninfected.pEC,
                                   aEC = groupingMADMelNV_infected.aEC - groupingMAuninfected.aEC,
                                   cardia = groupingMADMelNV_infected.cardia - groupingMAuninfected.cardia,
                                   `copper/iron` = groupingMADMelNV_infected.copperiron - groupingMAuninfected.copperiron,
                                   `ISC/EB` = groupingMADMelNV_infected.ISCEB - groupingMAuninfected.ISCEB,
                                   LFC = groupingMADMelNV_infected.LFC - groupingMAuninfected.LFC,
                                   mEC = groupingMADMelNV_infected.mEC - groupingMAuninfected.mEC,
                                   ECs = (groupingMADMelNV_infected.pEC +
                                            groupingMADMelNV_infected.aEC +
                                            groupingMADMelNV_infected.mEC)/3 -
                                     (groupingMAuninfected.pEC +
                                        groupingMAuninfected.aEC +
                                        groupingMAuninfected.mEC)/3,
                                   levels = glm.vir$MA$design)

# make tests

thika_gGP_infection <- test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"infection"],
                               sort_by = pval)

nora_EE_gGP_infection <- test_de(glm.vir$EE$fit, contrast = nora_EE_contrasts[,"EE"],
                                 sort_by = pval)

nora_MA_gGP_infection <- test_de(glm.vir$MA$fit, contrast = nora_MA_contrasts[,"MA"],
                                 sort_by = pval)

# generic response genes

thika_glmGamPoi_infection_genes <- thika_gGP_infection[thika_gGP_infection$adj_pval < 0.05,]
thika_glmGamPoi_infection_genes$symbol <- gene_symbol[match(gsub("-", "_", thika_glmGamPoi_infection_genes$name), gene_symbol$V1), "V3"]
thika_glmGamPoi_infection_genes_pos <- thika_glmGamPoi_infection_genes[thika_glmGamPoi_infection_genes$adj_pval < 0.05 & thika_glmGamPoi_infection_genes$lfc > 0,]
thika_glmGamPoi_infection_genes_neg <- thika_glmGamPoi_infection_genes[thika_glmGamPoi_infection_genes$adj_pval < 0.05 & thika_glmGamPoi_infection_genes$lfc < 0,]

nora_EE_glmGamPoi_infection_genes <- nora_EE_gGP_infection[nora_EE_gGP_infection$adj_pval < 0.05,]
nora_EE_glmGamPoi_infection_genes$symbol <- gene_symbol[match(gsub("-", "_", nora_EE_glmGamPoi_infection_genes$name), gene_symbol$V1), "V3"]
nora_EE_glmGamPoi_infection_genes_pos <- nora_EE_glmGamPoi_infection_genes[nora_EE_glmGamPoi_infection_genes$adj_pval < 0.05 & nora_EE_glmGamPoi_infection_genes$lfc > 0,]
nora_EE_glmGamPoi_infection_genes_neg <- nora_EE_glmGamPoi_infection_genes[nora_EE_glmGamPoi_infection_genes$adj_pval < 0.05 & nora_EE_glmGamPoi_infection_genes$lfc < 0,]

nora_MA_glmGamPoi_infection_genes <- nora_MA_gGP_infection[nora_MA_gGP_infection$adj_pval < 0.05,]
nora_MA_glmGamPoi_infection_genes$symbol <- gene_symbol[match(gsub("-", "_", nora_MA_glmGamPoi_infection_genes$name), gene_symbol$V1), "V3"]
nora_MA_glmGamPoi_infection_genes_pos <- nora_MA_glmGamPoi_infection_genes[nora_MA_glmGamPoi_infection_genes$adj_pval < 0.05 & nora_MA_glmGamPoi_infection_genes$lfc > 0,]
nora_MA_glmGamPoi_infection_genes_neg <- nora_MA_glmGamPoi_infection_genes[nora_MA_glmGamPoi_infection_genes$adj_pval < 0.05 & nora_MA_glmGamPoi_infection_genes$lfc < 0,]

if(save_steps) {
  save(glm.vir, thika_gGP_infection, nora_EE_gGP_infection, nora_MA_gGP_infection, file = "objs/glm.vir_glmGamPoi.Rdata")
}

# cell-subtype-specific response genes

thika_gGP_oneway_DE <- list(I_m = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"I_m"]),
                            II_m1 = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"II_m1"]),
                            I_ap_p = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"I_ap_p"]),
                            II_a = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"II_a"]),
                            II_m2 = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"II_m2"]),
                            I_pAstA = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"I_pAstA"]),
                            I_pCCHa1 = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"I_pCCHa1"]),
                            II_p = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"II_p"]),
                            I_ap_a = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"I_ap_a"]),
                            I_a = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"I_a"]),
                            III = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"III"]))

nora_gGP_oneway_DE <- list(I_m = test_de(glm.vir$EE$fit, contrast = nora_EE_contrasts[,"I_m"]),
                           II_m1 = test_de(glm.vir$EE$fit, contrast = nora_EE_contrasts[,"II_m1"]),
                           I_ap_p = test_de(glm.vir$EE$fit, contrast = nora_EE_contrasts[,"I_ap_p"]),
                           II_a = test_de(glm.vir$EE$fit, contrast = nora_EE_contrasts[,"II_a"]),
                           II_m2 = test_de(glm.vir$EE$fit, contrast = nora_EE_contrasts[,"II_m2"]),
                           I_pAstA = test_de(glm.vir$EE$fit, contrast = nora_EE_contrasts[,"I_pAstA"]),
                           I_pCCHa1 = test_de(glm.vir$EE$fit, contrast = nora_EE_contrasts[,"I_pCCHa1"]),
                           II_p = test_de(glm.vir$EE$fit, contrast = nora_EE_contrasts[,"II_p"]),
                           I_ap_a = test_de(glm.vir$EE$fit, contrast = nora_EE_contrasts[,"I_ap_a"]),
                           EE = test_de(glm.vir$MA$fit, contrast = nora_MA_contrasts[,"EE"]),
                           pEC = test_de(glm.vir$MA$fit, contrast = nora_MA_contrasts[,"pEC"]),
                           aEC = test_de(glm.vir$MA$fit, contrast = nora_MA_contrasts[,"aEC"]),
                           cardia = test_de(glm.vir$MA$fit, contrast = nora_MA_contrasts[,"cardia"]),
                           `copper/iron` = test_de(glm.vir$MA$fit, contrast = nora_MA_contrasts[,"copper/iron"]),
                           `ISC/EB` = test_de(glm.vir$MA$fit, contrast = nora_MA_contrasts[,"ISC/EB"]),
                           LFC = test_de(glm.vir$MA$fit, contrast = nora_MA_contrasts[,"LFC"]),
                           mEC = test_de(glm.vir$MA$fit, contrast = nora_MA_contrasts[,"mEC"]))

if(save_steps) {
  save(glm.vir, thika_gGP_infection, nora_EE_gGP_infection, nora_MA_gGP_infection, 
     thika_gGP_oneway_DE, nora_gGP_oneway_DE,
     file = "objs/glm.vir_glmGamPoi.Rdata")
}


# make lists

thika_gGP_celltype_genes <- list(I_m_pos = thika_gGP_oneway_DE$I_m[thika_gGP_oneway_DE$I_m$lfc > 0 & thika_gGP_oneway_DE$I_m$adj_pval < 0.05,],
                                 II_m1_pos = thika_gGP_oneway_DE$II_m1[thika_gGP_oneway_DE$II_m1$lfc > 0 & thika_gGP_oneway_DE$II_m1$adj_pval < 0.05,],
                                 I_ap_p_pos = thika_gGP_oneway_DE$I_ap_p[thika_gGP_oneway_DE$I_ap_p$lfc > 0 & thika_gGP_oneway_DE$I_ap_p$adj_pval < 0.05,],
                                 II_a_pos = thika_gGP_oneway_DE$II_a[thika_gGP_oneway_DE$II_a$lfc > 0 & thika_gGP_oneway_DE$II_a$adj_pval < 0.05,],
                                 II_m2_pos = thika_gGP_oneway_DE$II_m2[thika_gGP_oneway_DE$II_m2$lfc > 0 & thika_gGP_oneway_DE$II_m2$adj_pval < 0.05,],
                                 I_pAstA_pos = thika_gGP_oneway_DE$I_pAstA[thika_gGP_oneway_DE$I_pAstA$lfc > 0 & thika_gGP_oneway_DE$I_pAstA$adj_pval < 0.05,],
                                 I_pCCHa1_pos = thika_gGP_oneway_DE$I_pCCHa1[thika_gGP_oneway_DE$I_pCCHa1$lfc > 0 & thika_gGP_oneway_DE$I_pCCHa1$adj_pval < 0.05,],
                                 II_p_pos = thika_gGP_oneway_DE$II_p[thika_gGP_oneway_DE$II_p$lfc > 0 & thika_gGP_oneway_DE$II_p$adj_pval < 0.05,],
                                 I_ap_a_pos = thika_gGP_oneway_DE$I_ap_a[thika_gGP_oneway_DE$I_ap_a$lfc > 0 & thika_gGP_oneway_DE$I_ap_a$adj_pval < 0.05,],
                                 I_a_pos = thika_gGP_oneway_DE$I_a[thika_gGP_oneway_DE$I_a$lfc > 0 & thika_gGP_oneway_DE$I_a$adj_pval < 0.05,],
                                 III_pos = thika_gGP_oneway_DE$III[thika_gGP_oneway_DE$III$lfc > 0 & thika_gGP_oneway_DE$III$adj_pval < 0.05,],
                                   
                                 I_m_neg = thika_gGP_oneway_DE$I_m[thika_gGP_oneway_DE$I_m$lfc < 0 & thika_gGP_oneway_DE$I_m$adj_pval < 0.05,],
                                 II_m1_neg = thika_gGP_oneway_DE$II_m1[thika_gGP_oneway_DE$II_m1$lfc < 0 & thika_gGP_oneway_DE$II_m1$adj_pval < 0.05,],
                                 I_ap_p_neg = thika_gGP_oneway_DE$I_ap_p[thika_gGP_oneway_DE$I_ap_p$lfc < 0 & thika_gGP_oneway_DE$I_ap_p$adj_pval < 0.05,],
                                 II_a_neg = thika_gGP_oneway_DE$II_a[thika_gGP_oneway_DE$II_a$lfc < 0 & thika_gGP_oneway_DE$II_a$adj_pval < 0.05,],
                                 II_m2_neg = thika_gGP_oneway_DE$II_m2[thika_gGP_oneway_DE$II_m2$lfc < 0 & thika_gGP_oneway_DE$II_m2$adj_pval < 0.05,],
                                 I_pAstA_neg = thika_gGP_oneway_DE$I_pAstA[thika_gGP_oneway_DE$I_pAstA$lfc < 0 & thika_gGP_oneway_DE$I_pAstA$adj_pval < 0.05,],
                                 I_pCCHa1_neg = thika_gGP_oneway_DE$I_pCCHa1[thika_gGP_oneway_DE$I_pCCHa1$lfc < 0 & thika_gGP_oneway_DE$I_pCCHa1$adj_pval < 0.05,],
                                 II_p_neg = thika_gGP_oneway_DE$II_p[thika_gGP_oneway_DE$II_p$lfc < 0 & thika_gGP_oneway_DE$II_p$adj_pval < 0.05,],
                                 I_ap_a_neg = thika_gGP_oneway_DE$I_ap_a[thika_gGP_oneway_DE$I_ap_a$lfc < 0 & thika_gGP_oneway_DE$I_ap_a$adj_pval < 0.05,],
                                 I_a_neg = thika_gGP_oneway_DE$I_a[thika_gGP_oneway_DE$I_a$lfc < 0 & thika_gGP_oneway_DE$I_a$adj_pval < 0.05,],
                                 III_neg = thika_gGP_oneway_DE$III[thika_gGP_oneway_DE$III$lfc < 0 & thika_gGP_oneway_DE$III$adj_pval < 0.05,])

thika_gGP_celltype_exclusive_genes <- data.frame(name = character(), symbol = character())

for (i in 1:length(thika_gGP_celltype_genes)) {
  if (is.null(thika_gGP_celltype_genes[[i]])) next
  thika_gGP_celltype_genes[[i]]$symbol <- gene_symbol[match(gsub("-", "_", thika_gGP_celltype_genes[[i]]$name), gene_symbol$V1), "V3"]
  tmp_tb <- thika_gGP_celltype_genes[[i]][!thika_gGP_celltype_genes[[i]]$symbol %in% thika_glmGamPoi_infection_genes$symbol,]
  thika_gGP_celltype_exclusive_genes <- rbind(thika_gGP_celltype_exclusive_genes, tmp_tb[c("name", "symbol", "lfc")])
}

thika_unique_celltype_genes_pos <- thika_gGP_celltype_exclusive_genes[thika_gGP_celltype_exclusive_genes$lfc > 0, "symbol"] %>% unique()
thika_unique_celltype_genes_neg <- thika_gGP_celltype_exclusive_genes[thika_gGP_celltype_exclusive_genes$lfc < 0, "symbol"] %>% unique()

nora_gGP_celltype_genes <- list(I_m_pos = nora_gGP_oneway_DE$I_m[nora_gGP_oneway_DE$I_m$lfc > 0 & nora_gGP_oneway_DE$I_m$adj_pval < 0.05,],
                                II_m1_pos = nora_gGP_oneway_DE$II_m1[nora_gGP_oneway_DE$II_m1$lfc > 0 & nora_gGP_oneway_DE$II_m1$adj_pval < 0.05,],
                                I_ap_p_pos = nora_gGP_oneway_DE$I_ap_p[nora_gGP_oneway_DE$I_ap_p$lfc > 0 & nora_gGP_oneway_DE$I_ap_p$adj_pval < 0.05,],
                                II_a_pos = nora_gGP_oneway_DE$II_a[nora_gGP_oneway_DE$II_a$lfc > 0 & nora_gGP_oneway_DE$II_a$adj_pval < 0.05,],
                                II_m2_pos = nora_gGP_oneway_DE$II_m2[nora_gGP_oneway_DE$II_m2$lfc > 0 & nora_gGP_oneway_DE$II_m2$adj_pval < 0.05,],
                                I_pAstA_pos = nora_gGP_oneway_DE$I_pAstA[nora_gGP_oneway_DE$I_pAstA$lfc > 0 & nora_gGP_oneway_DE$I_pAstA$adj_pval < 0.05,],
                                I_pCCHa1_pos = nora_gGP_oneway_DE$I_pCCHa1[nora_gGP_oneway_DE$I_pCCHa1$lfc > 0 & nora_gGP_oneway_DE$I_pCCHa1$adj_pval < 0.05,],
                                II_p_pos = nora_gGP_oneway_DE$II_p[nora_gGP_oneway_DE$II_p$lfc > 0 &  nora_gGP_oneway_DE$II_p$adj_pval < 0.05,],
                                I_ap_a_pos = nora_gGP_oneway_DE$I_ap_a[nora_gGP_oneway_DE$I_ap_a$lfc > 0 & nora_gGP_oneway_DE$I_ap_a$adj_pval < 0.05,],
                                EE_pos = nora_gGP_oneway_DE$EE[nora_gGP_oneway_DE$EE$lfc > 0 & nora_gGP_oneway_DE$EE$adj_pval < 0.05,],
                                pEC_pos = nora_gGP_oneway_DE$pEC[nora_gGP_oneway_DE$pEC$lfc > 0 & nora_gGP_oneway_DE$pEC$adj_pval < 0.05,],
                                aEC_pos = nora_gGP_oneway_DE$aEC[nora_gGP_oneway_DE$aEC$lfc > 0 & nora_gGP_oneway_DE$aEC$adj_pval < 0.05,],
                                cardia_pos = nora_gGP_oneway_DE$cardia[nora_gGP_oneway_DE$cardia$lfc > 0 & nora_gGP_oneway_DE$cardia$adj_pval < 0.05,],
                                `copper/iron_pos` = nora_gGP_oneway_DE$`copper/iron`[nora_gGP_oneway_DE$`copper/iron`$lfc > 0 & nora_gGP_oneway_DE$`copper/iron`$adj_pval < 0.05,],
                                `ISC/EB_pos` = nora_gGP_oneway_DE$`ISC/EB`[nora_gGP_oneway_DE$`ISC/EB`$lfc > 0 & nora_gGP_oneway_DE$`ISC/EB`$adj_pval < 0.05,],
                                LFC_pos = nora_gGP_oneway_DE$LFC[nora_gGP_oneway_DE$LFC$lfc > 0 & nora_gGP_oneway_DE$LFC$adj_pval < 0.05,],
                                mEC_pos = nora_gGP_oneway_DE$mEC[nora_gGP_oneway_DE$mEC$lfc > 0 & nora_gGP_oneway_DE$mEC$adj_pval < 0.05,],

                                I_m_neg = nora_gGP_oneway_DE$I_m[nora_gGP_oneway_DE$I_m$lfc < 0 & nora_gGP_oneway_DE$I_m$adj_pval < 0.05,],
                                II_m1_neg = nora_gGP_oneway_DE$II_m1[nora_gGP_oneway_DE$II_m1$lfc < 0 & nora_gGP_oneway_DE$II_m1$adj_pval < 0.05,],
                                I_ap_p_neg = nora_gGP_oneway_DE$I_ap_p[nora_gGP_oneway_DE$I_ap_p$lfc < 0 & nora_gGP_oneway_DE$I_ap_p$adj_pval < 0.05,],
                                II_a_neg = nora_gGP_oneway_DE$II_a[nora_gGP_oneway_DE$II_a$lfc < 0 & nora_gGP_oneway_DE$II_a$adj_pval < 0.05,],
                                II_m2_neg = nora_gGP_oneway_DE$II_m2[nora_gGP_oneway_DE$II_m2$lfc < 0 & nora_gGP_oneway_DE$II_m2$adj_pval < 0.05,],
                                I_pAstA_neg = nora_gGP_oneway_DE$I_pAstA[nora_gGP_oneway_DE$I_pAstA$lfc < 0 & nora_gGP_oneway_DE$I_pAstA$adj_pval < 0.05,],
                                I_pCCHa1_neg = nora_gGP_oneway_DE$I_pCCHa1[nora_gGP_oneway_DE$I_pCCHa1$lfc < 0 & nora_gGP_oneway_DE$I_pCCHa1$adj_pval < 0.05,],
                                II_p_neg = nora_gGP_oneway_DE$II_p[nora_gGP_oneway_DE$II_p$lfc < 0 & nora_gGP_oneway_DE$II_p$adj_pval < 0.05,],
                                I_ap_a_neg = nora_gGP_oneway_DE$I_ap_a[nora_gGP_oneway_DE$I_ap_a$lfc < 0 & nora_gGP_oneway_DE$I_ap_a$adj_pval < 0.05,],
                                EE_neg = nora_gGP_oneway_DE$EE[nora_gGP_oneway_DE$EE$lfc < 0 & nora_gGP_oneway_DE$EE$adj_pval < 0.05,],
                                pEC_neg = nora_gGP_oneway_DE$pEC[nora_gGP_oneway_DE$pEC$lfc < 0 & nora_gGP_oneway_DE$pEC$adj_pval < 0.05,],
                                aEC_neg = nora_gGP_oneway_DE$aEC[nora_gGP_oneway_DE$aEC$lfc < 0 & nora_gGP_oneway_DE$aEC$adj_pval < 0.05,],
                                cardia_neg = nora_gGP_oneway_DE$cardia[nora_gGP_oneway_DE$cardia$lfc < 0 & nora_gGP_oneway_DE$cardia$adj_pval < 0.05,],
                                `copper/iron_neg` = nora_gGP_oneway_DE$`copper/iron`[nora_gGP_oneway_DE$`copper/iron`$lfc < 0 & nora_gGP_oneway_DE$`copper/iron`$adj_pval < 0.05,],
                                `ISC/EB_neg` = nora_gGP_oneway_DE$`ISC/EB`[nora_gGP_oneway_DE$`ISC/EB`$lfc < 0 & nora_gGP_oneway_DE$`ISC/EB`$adj_pval < 0.05,],
                                LFC_neg = nora_gGP_oneway_DE$LFC[nora_gGP_oneway_DE$LFC$lfc < 0 & nora_gGP_oneway_DE$LFC$adj_pval < 0.05,],
                                mEC_neg = nora_gGP_oneway_DE$mEC[nora_gGP_oneway_DE$mEC$lfc < 0 & nora_gGP_oneway_DE$mEC$adj_pval < 0.05,])

nora_gGP_celltype_exclusive_genes <- data.frame(name = character(), symbol = character())

for (i in 1:length(nora_gGP_celltype_genes)) {
  if (is.null(nora_gGP_celltype_genes[[i]])) next
  nora_gGP_celltype_genes[[i]]$symbol <- gene_symbol[match(gsub("-", "_", nora_gGP_celltype_genes[[i]]$name), gene_symbol$V1), "V3"]
  tmp_tb <- nora_gGP_celltype_genes[[i]][!nora_gGP_celltype_genes[[i]]$symbol %in% c(nora_EE_glmGamPoi_infection_genes$symbol, nora_MA_glmGamPoi_infection_genes$symbol),]
  nora_gGP_celltype_exclusive_genes <- rbind(nora_gGP_celltype_exclusive_genes, tmp_tb[c("name", "symbol", "lfc")])
}

rm(tmp_tb)

nora_unique_celltype_genes_pos <- nora_gGP_celltype_exclusive_genes[nora_gGP_celltype_exclusive_genes$lfc > 0, "symbol"] %>% unique()
nora_unique_celltype_genes_neg <- nora_gGP_celltype_exclusive_genes[nora_gGP_celltype_exclusive_genes$lfc < 0, "symbol"] %>% unique()

common_unique_celltype_genes_pos <- nora_unique_celltype_genes_pos[nora_unique_celltype_genes_pos %in% thika_unique_celltype_genes_pos]
common_unique_celltype_genes_neg <- nora_unique_celltype_genes_neg[nora_unique_celltype_genes_neg %in% thika_unique_celltype_genes_neg]

# additional tests

thika_region_DE <- list(thika_gGP_a_p = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"a_p"]),
                        thika_gGP_a = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"a"]),
                        thika_gGP_p = test_de(glm.vir$EE$fit, contrast = thika_contrasts[,"p"]))

nora_major_DE <- list(EC = test_de(glm.vir$MA$fit, contrast = nora_MA_contrasts[,"ECs"]))

for (i in 1:length(thika_region_DE)) {
  if (is.null(thika_region_DE[[i]])) next
  thika_region_DE[[i]]$symbol <- gene_symbol[match(gsub("-", "_", thika_region_DE[[i]]$name), gene_symbol$V1), "V3"]
}

for (i in 1:length(nora_major_DE)) {
  if (is.null(nora_major_DE[[i]])) next
  nora_major_DE[[i]]$symbol <- gene_symbol[match(gsub("-", "_", nora_major_DE[[i]]$name), gene_symbol$V1), "V3"]
}

if(save_steps) {
  save(glm.vir, thika_gGP_infection, nora_EE_gGP_infection, nora_MA_gGP_infection, 
       thika_gGP_oneway_DE, nora_gGP_oneway_DE,
       thika_region_DE, nora_major_DE,
       file = "objs/glm.vir_glmGamPoi.Rdata")
}

thika_p_pos <- thika_region_DE$thika_gGP_p[thika_region_DE$thika_gGP_p$adj_pval < 0.05 & thika_region_DE$thika_gGP_p$lfc > 0,]
thika_p_neg <- thika_region_DE$thika_gGP_p[thika_region_DE$thika_gGP_p$adj_pval < 0.05 & thika_region_DE$thika_gGP_p$lfc < 0,]

nora_EC_pos <- nora_major_DE$EC[nora_major_DE$EC$adj_pval < 0.05 & nora_major_DE$EC$lfc > 0,]
nora_EC_neg <- nora_major_DE$EC[nora_major_DE$EC$adj_pval < 0.05 & nora_major_DE$EC$lfc < 0,]

# correlation glm

cor.glm.GamPoi <- cor_gGP(dmel_filtered, mda_10x_merged)

thika.cor.gampoi <- test_de(cor.glm.GamPoi$EE$TV$fit, contrast = "tv.pct")
nora_EE.cor.gampoi <- test_de(cor.glm.GamPoi$EE$DMelNV$fit, contrast = "dmelnvEE.pct")
nora_MA.cor.gampoi <- test_de(cor.glm.GamPoi$MA$fit, contrast = "dmelnvMA.pct")

thika.cor.gampoi_significant <- thika.cor.gampoi[thika.cor.gampoi$adj_pval < 0.05,]
thika.cor.gampoi_significant$symbol <- gene_symbol[match(gsub("-", "_", thika.cor.gampoi_significant$name), gene_symbol$V1), "V3"]

nora_EE.cor.gampoi_significant <- nora_EE.cor.gampoi[nora_EE.cor.gampoi$adj_pval < 0.05,]
nora_EE.cor.gampoi_significant$symbol <- gene_symbol[match(gsub("-", "_", nora_EE.cor.gampoi_significant$name), gene_symbol$V1), "V3"]

nora_MA.cor.gampoi_significant <- nora_MA.cor.gampoi[nora_MA.cor.gampoi $adj_pval < 0.05,]
nora_MA.cor.gampoi_significant$symbol <- gene_symbol[match(gsub("-", "_", nora_MA.cor.gampoi_significant$name), gene_symbol$V1), "V3"]

if(save_steps) {
  save(glm.vir, thika_gGP_infection, nora_EE_gGP_infection, nora_MA_gGP_infection, 
       thika_gGP_oneway_DE, nora_gGP_oneway_DE,
       thika_region_DE, nora_major_DE,
       cor.glm.GamPoi,
       thika.cor.gampoi, nora_EE.cor.gampoi, nora_MA.cor.gampoi,
       file = "objs/glm.vir_glmGamPoi.Rdata")
}
