# use cellbender counts for DE

library(Seurat)
library(dplyr)
library(tibble)
library(viridis)
library(ggplot2)
library(sctransform)
library(limma)
library(reshape2)
library(ggpubr)
library(gridExtra)
library(agricolae)
library(ggforce)
library(igraph)
library(tidyverse)
library(data.table)
library(car)
library(broom)
library(glmGamPoi)
library(ReactomePA)
library(DropletUtils)
library(ggrepel)
library(writexl)

# create folders that will contain output files

dir.create("objs") # precessed Seurat objects and DE analysis results
dir.create("cellbender") # CellBender analysis should be performed on this directory
dir.create("figure_out") # output figures

# table to convert ids <-> symbol <-> name 

gene_symbol <- read.table("name_id_symbol_list",
                          stringsAsFactors = FALSE)

# changing lrRNA symbol to mt:lrRNA

gene_symbol[gene_symbol$V3 == "lrRNA", 3] <- "mt:lrRNA"

# adding some missing genes

gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel_CR45859", 26067187, "28SrRNA-Psi:CR45859")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel_CR40596", 5740740, "28SrRNA-Psi:CR40596")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel_CR34096", 19893562, "srRNA")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel_CR40741", 5740241, "28SrRNA-Psi:CR40741")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel_CR41609", 5740375, "28SrRNA-Psi:CR41609")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel_CR45853", 26067181, "28SrRNA-Psi:CR45853")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel_CR41602", 5740420, "28SrRNA-Psi:CR41602")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel_CG7298", 40210, "CG7298")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR45855", 26067183, "28SrRNA-Psi:CR45855")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR45857", 26067185, "5.8SrRNA-Psi:CR45857")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR45860", 26067188, "28SrRNA-Psi:CR45860")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG4184", 33223, "MED15")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG14031", 33756, "Cyp4ac3")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG31713", 318908, "Apf")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG5603", 34380, "CYLD")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG10846", 34873, "DCTN5-p25")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG10195", 35228, "CG10195")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG34392", 35588, "Epac")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG5498", 40250, "CG5498")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG33217", 2768951, "CG33217")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG11984", 41082, "Kcmf1")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG43143", 41256, "Nuak1")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG6195", 42343, "CG6195")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG10618", 43067, "CHKov1")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR45848", 26067176, "28SrRNA-Psi:CR4584")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR45849", 26067177, "5.8SrRNA-Psi:CR45849")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34077", 19893543, "trnA")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34078", 19893544, "trnR")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34095", 19893561, "trnV")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR45863", 26067191, "5.8SrRNA-Psi:CR45863")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG6711", 39164, "Taf2")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG18596", 42622, "CG18596")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG11839", 43006, "CG11839")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34075", 19893541, "trnG")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34079", 19893545, "trnN")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34080", 19893546, "trnS2")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34093", 19893559, "trnL2")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG3206", 31207, "Or2a")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG15899", 31550, "Ca-alpha1T")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG42395", 7354436, "CG42395")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG7068", 34045, "Tep3")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG33267", 2768969, "CG33267")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG6284", 41254, "Sirt6")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34064", 19893530, "trnW")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34066", 19893532, "trnY")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34068", 19893534, "trnL1")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34070", 19893536, "trnK")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34071", 19893537, "trnD")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR45851", 26067179, "28SrRNA-Psi:CR45851")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34065", 19893531, "trnC")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CR34082", 19893548, "trnF")
gene_symbol[nrow(gene_symbol) + 1, ] <- list("Dmel-CG32261", 117481, "Gr64a")

# calculate a threshold for each cell to call them infected with DMelNV based on the percentage of ambient RNA and percentage of viral RNA in the ambient RNA
# no TV reads were found in empty cells

is.infected <- function(obj,
                        pct.nora.ambientRNA) {
  
  is.inf <- rep(NA, ncol(obj))
  names(is.inf) <- colnames(obj)
  
  for (cell in names(is.inf)) {
    
    counts.vir <- obj@assays$RNA@counts["Dmel-nora-virus", cell]
    
    if (counts.vir == 0) {
      
      is.inf[cell] <- FALSE 
      next
    }
    
    nUMI_novir <- obj$nCount_NOVIR[cell]
    nUMI_cellbender <- obj$nCount_CELLBENDER[cell]
    nUMI_RNA <- obj$nCount_RNA[cell]
    
    # some cells have more UMIs after cellbender correction, these will be considered to contain no ambient RNA
    
    if (nUMI_cellbender >= nUMI_novir) {
      
      is.inf[cell] <- FALSE 
      next
    }
    
    pct.ambientRNA_in_cell <- 1 - (nUMI_cellbender/nUMI_novir)
    prob_ambient.vir <- pct.ambientRNA_in_cell*pct.nora.ambientRNA
    
    # calculate probability of having at least counts.vir reads due to ambient viral RNA
    # if counts.vir or more (survival function > 0.01) is due to ambient viral RNA, cell is not infected
    # infected cells need to have > 1 viral reads
    
    prob_infected <- 1 - pbinom(q = counts.vir, size = nUMI_RNA, prob = prob_ambient.vir)
    
    is.inf[[cell]] <- prob_infected < 0.01 & counts.vir > 1
  }
  return(is.inf)
}

DE_gGP <- function(EE.seu, MA.seu,
                   batch = FALSE, intercept = FALSE,
                   oneway = TRUE,
                   verbose = TRUE) {
  
  # uninfected cells only will be selected if they have 0% viral reads
  
  obj_tv <- subset(EE.seu, cells = colnames(EE.seu)[EE.seu$is.tv.infected])
  obj_dmelnv_EE <- subset(EE.seu, cells = colnames(EE.seu)[EE.seu$is.dmelnv.infected])
  obj_EE_no_vir <- subset(EE.seu, subset = percent.thika == 0 & percent.nora == 0 & percent.dmelc == 0,
                          idents = c(obj_tv@active.ident %>% as.character(),
                                     obj_dmelnv_EE@active.ident %>% as.character()) %>% unique())
  
  obj_dmelnv_MA <- subset(MA.seu, cells = colnames(MA.seu)[MA.seu$is.dmelnv.infected])
  obj_MA_no_vir <- subset(MA.seu, subset = percent.nora == 0,
                          idents = obj_dmelnv_MA@active.ident %>% unique())

  count_tv <- obj_tv@assays$CELLBENDER@counts
  status_tv_EE <- rep("TV_infected", ncol(count_tv))
  cells_tv_EE <- obj_tv@active.ident %>% as.character()

  count_dmelnv_EE <- obj_dmelnv_EE@assays$CELLBENDER@counts
  status_dmelnv_EE <- rep("DMelNV_infected", ncol(count_dmelnv_EE))
  cells_dmelnv_EE <- obj_dmelnv_EE@active.ident %>% as.character()
  
  uninfected_count_EE <- obj_EE_no_vir@assays$CELLBENDER@counts
  uninfected_status_EE <- rep("uninfected", ncol(uninfected_count_EE))
  uninfected_cells_EE <- obj_EE_no_vir@active.ident %>% as.character()
  
  infected_count_dmelnv_MA <- obj_dmelnv_MA@assays$CELLBENDER@counts
  infected_status_dmelnv_MA <- rep("DMelNV_infected", ncol(infected_count_dmelnv_MA))
  infected_cells_dmelnv_MA <- obj_dmelnv_MA@active.ident %>% as.character()
  infected_batch_dmelnv_MA <- obj_dmelnv_MA$orig.ident
  
  uninfected_count_MA <- obj_MA_no_vir@assays$CELLBENDER@counts
  uninfected_status_MA <- rep("uninfected", ncol(uninfected_count_MA))
  uninfected_cells_MA <- obj_MA_no_vir@active.ident %>% as.character()
  uninfected_batch_MA <- obj_MA_no_vir$orig.ident
  
  # one count matrix for each virus dataset
  
  count.mtx_EE <- cbind(count_tv, count_dmelnv_EE, uninfected_count_EE)
  infection_EE <- c(status_tv_EE, status_dmelnv_EE, uninfected_status_EE) %>% factor()
  cell_EE <- c(cells_tv_EE, cells_dmelnv_EE, uninfected_cells_EE) %>% factor()

  count.mtx_dmelnv_MA <- cbind(infected_count_dmelnv_MA, uninfected_count_MA)
  infection_dmelnv_MA <- c(infected_status_dmelnv_MA, uninfected_status_MA) %>% factor()
  cell_dmelnv_MA <- c(infected_cells_dmelnv_MA, uninfected_cells_MA) %>% factor()
  sample.names_dmelnv_MA <- c(infected_batch_dmelnv_MA, uninfected_batch_MA) %>% factor()
  
  if (oneway) {
    
    groupingEE <- paste0(infection_EE, ".", cell_EE)
    groupingMA <- paste0(infection_dmelnv_MA, ".", cell_dmelnv_MA)
    
    if (intercept) {
      
      if (batch) {
        
        design_EE <- model.matrix(~groupingEE)
        design_MA <- model.matrix(~groupingMA + sample.names_dmelnv_MA)
        
      } else {
        
        design_EE <- model.matrix(~groupingEE)
        design_MA <- model.matrix(~groupingMA)
        
      }
      
    } else {
      
      if (batch) {
        
        design_EE <- model.matrix(~0 + groupingEE)
        design_MA <- model.matrix(~0 + groupingMA + sample.names_dmelnv_MA)
        
      } else {
        
        design_EE <- model.matrix(~0 + groupingEE)
        design_MA <- model.matrix(~0 + groupingMA)
        
      }
      
    }
    
  } else {
    
    if (intercept) {
      
      if (batch) {
        
        design_EE <- model.matrix(~infection_EE * cell_EE)
        design_MA <- model.matrix(~infection_dmelnv_MA * cell_dmelnv_MA + sample.names_dmelnv_MA)

      } else {
        
        design_EE <- model.matrix(~infection_EE * cell_EE)
        design_MA <- model.matrix(~infection_dmelnv_MA * cell_dmelnv_MA)
        
      }

    } else {
      
      if (batch) {
        
        design_EE <- model.matrix(~0 + infection_EE * cell_EE + sample.names_thika)
        design_MA <- model.matrix(~0 + infection_dmelnv_MA * cell_dmelnv_MA + sample.names_dmelnv_MA)
        
      } else {
        
        design_EE <- model.matrix(~0 + infection_EE * cell_EE)
        design_MA <- model.matrix(~0 + infection_dmelnv_MA * cell_dmelnv_MA)
        
      }
      
    } 
  }
  
  colnames(design_EE) <- gsub("/", "", gsub("-", "_", gsub(":", ".", colnames(design_EE))))
  colnames(design_MA) <- gsub("/", "", gsub("-", "_", gsub(":", ".", colnames(design_MA))))
  
  if (verbose) print("fitting to EE dataset")
  fit_EE <- glm_gp(count.mtx_EE %>% as.matrix(), design_EE, on_disk = FALSE)
  
  if (verbose) print("fitting to MA dataset")
  fit_MA <- glm_gp(count.mtx_dmelnv_MA %>% as.matrix(), design_MA, on_disk = FALSE)
  
  return(list(EE = list(fit = fit_EE,
                           design = design_EE),
              MA = list(fit = fit_MA,
                          design = design_MA,
                          samples = sample.names_dmelnv_MA)))
}

# function to do DE tests on uninfected cells
# a fraction of the uninfected cells will be marked as infected

robustness_DE_gGP <- function(EE.seu, MA.seu,
                              batch = FALSE, intercept = FALSE,
                              oneway = TRUE,
                              seed,
                              verbose = FALSE) {
  
  # first, elimanate EEP cells from EE dataset
  
  EE.seu <- subset(EE.seu, idents = unique(EE.seu@active.ident)[unique(EE.seu@active.ident) != "EEP"])
  
  obj_tv <- subset(EE.seu, cells = colnames(EE.seu)[EE.seu$is.tv.infected])
  obj_dmelnv_EE <- subset(EE.seu, cells = colnames(EE.seu)[EE.seu$is.dmelnv.infected])
  obj_EE_no_vir <- subset(EE.seu, subset = percent.thika == 0 & percent.nora == 0 & percent.dmelc == 0,
                          idents = c(obj_tv@active.ident %>% as.character(),
                                     obj_dmelnv_EE@active.ident %>% as.character()) %>% unique())
  
  # uninfected cells only will be selected if they have 0% viral reads
  
  obj_dmelnv_MA <- subset(MA.seu, cells = colnames(MA.seu)[MA.seu$is.dmelnv.infected])
  obj_MA_no_vir <- subset(MA.seu, subset = percent.nora == 0,
                          idents = obj_dmelnv_MA@active.ident %>% unique())
  
  # calculate percentage of infected cells for each subtype
  
  tv.pct <- (obj_tv@active.ident %>% table())/
    (subset(EE.seu, idents = obj_tv@active.ident %>% unique())@active.ident %>% table())
  
  dmelnv_EE.pct <- (obj_dmelnv_EE@active.ident %>% table())/
    (subset(EE.seu, idents = obj_dmelnv_EE@active.ident %>% unique())@active.ident %>% table())
  
  dmelnv_MA.pct <- (obj_dmelnv_MA@active.ident %>% table())/
    (subset(MA.seu, idents = obj_dmelnv_MA@active.ident %>% unique())@active.ident %>% table())
  
  # total uninfected cells
  
  EE.nCell <- obj_EE_no_vir@active.ident %>% table()
  MA.nCell <- obj_MA_no_vir@active.ident %>% table()
  
  if (verbose) {
    
    print(tv.pct)
    print(dmelnv_EE.pct)
    print(dmelnv_MA.pct)
    
  }
  
  # get cells to resubset objs of infected cells using non-infected cells based on the percentage calculated above
  
  tv_cells_to_subset <- c()
  dmelnv_EE_cells_to_subset <- c()
  
  for (subtype in obj_EE_no_vir@active.ident %>% unique()) {
    
    set.seed(seed)
    
    tv_cells_to_subset <- c(tv_cells_to_subset,
                            sample(obj_EE_no_vir@active.ident[obj_EE_no_vir@active.ident == subtype] %>% names(),
                                   ceiling(tv.pct[names(tv.pct) == subtype]*EE.nCell[names(EE.nCell) == subtype])))
    
    set.seed(seed)
    
    dmelnv_EE_cells_to_subset <- c(dmelnv_EE_cells_to_subset,
                                   sample(obj_EE_no_vir@active.ident[obj_EE_no_vir@active.ident == subtype] %>% names(),
                                          ceiling(dmelnv_EE.pct[names(dmelnv_EE.pct) == subtype]*EE.nCell[names(EE.nCell) == subtype])))
    
  }
  
  dmelnv_MA_cells_to_subset <- c()
  
  for (subtype in obj_MA_no_vir@active.ident %>% unique()) {
    
    set.seed(seed)
    
    dmelnv_MA_cells_to_subset <- c(dmelnv_MA_cells_to_subset,
                                   sample(obj_MA_no_vir@active.ident[obj_MA_no_vir@active.ident == subtype] %>% names(),
                                          ceiling(dmelnv_MA.pct[names(dmelnv_MA.pct) == subtype]*MA.nCell[names(MA.nCell) == subtype])))
    
  }
  
  # resubset
  
  obj_tv <- subset(EE.seu, cells = tv_cells_to_subset)
  obj_dmelnv_EE <- subset(EE.seu, cells = dmelnv_EE_cells_to_subset)
  obj_EE_no_vir <- subset(EE.seu, cells = colnames(obj_EE_no_vir)[!colnames(obj_EE_no_vir) %in% c(tv_cells_to_subset, dmelnv_EE_cells_to_subset)])
  
  obj_dmelnv_MA <- subset(MA.seu, cells = dmelnv_MA_cells_to_subset)
  obj_MA_no_vir <- subset(MA.seu, cells = colnames(obj_MA_no_vir)[!colnames(obj_MA_no_vir) %in% dmelnv_MA_cells_to_subset])
  
  
  # DE just like the DE_gGP function
  
  count_tv <- obj_tv@assays$CELLBENDER@counts
  status_tv_EE <- rep("TV_infected", ncol(count_tv))
  cells_tv_EE <- obj_tv@active.ident %>% as.character()
  
  count_dmelnv_EE <- obj_dmelnv_EE@assays$CELLBENDER@counts
  status_dmelnv_EE <- rep("DMelNV_infected", ncol(count_dmelnv_EE))
  cells_dmelnv_EE <- obj_dmelnv_EE@active.ident %>% as.character()
  
  uninfected_count_EE <- obj_EE_no_vir@assays$CELLBENDER@counts
  uninfected_status_EE <- rep("uninfected", ncol(uninfected_count_EE))
  uninfected_cells_EE <- obj_EE_no_vir@active.ident %>% as.character()
  
  infected_count_dmelnv_MA <- obj_dmelnv_MA@assays$CELLBENDER@counts
  infected_status_dmelnv_MA <- rep("DMelNV_infected", ncol(infected_count_dmelnv_MA))
  infected_cells_dmelnv_MA <- obj_dmelnv_MA@active.ident %>% as.character()
  infected_batch_dmelnv_MA <- obj_dmelnv_MA$orig.ident
  
  uninfected_count_MA <- obj_MA_no_vir@assays$CELLBENDER@counts
  uninfected_status_MA <- rep("uninfected", ncol(uninfected_count_MA))
  uninfected_cells_MA <- obj_MA_no_vir@active.ident %>% as.character()
  uninfected_batch_MA <- obj_MA_no_vir$orig.ident
  
  # one count matrix for each virus dataset
  
  count.mtx_EE <- cbind(count_tv, count_dmelnv_EE, uninfected_count_EE)
  infection_EE <- c(status_tv_EE, status_dmelnv_EE, uninfected_status_EE) %>% factor()
  cell_EE <- c(cells_tv_EE, cells_dmelnv_EE, uninfected_cells_EE) %>% factor()
  
  count.mtx_dmelnv_MA <- cbind(infected_count_dmelnv_MA, uninfected_count_MA)
  infection_dmelnv_MA <- c(infected_status_dmelnv_MA, uninfected_status_MA) %>% factor()
  cell_dmelnv_MA <- c(infected_cells_dmelnv_MA, uninfected_cells_MA) %>% factor()
  sample.names_dmelnv_MA <- c(infected_batch_dmelnv_MA, uninfected_batch_MA) %>% factor()
  
  if (oneway) {
    
    groupingEE <- paste0(infection_EE, ".", cell_EE)
    groupingMA <- paste0(infection_dmelnv_MA, ".", cell_dmelnv_MA)
    
    if (intercept) {
      
      if (batch) {
        
        design_EE <- model.matrix(~groupingEE)
        design_MA <- model.matrix(~groupingMA + sample.names_dmelnv_MA)
        
      } else {
        
        design_EE <- model.matrix(~groupingEE)
        design_MA <- model.matrix(~groupingMA)
        
      }
      
    } else {
      
      if (batch) {
        
        design_EE <- model.matrix(~0 + groupingEE)
        design_MA <- model.matrix(~0 + groupingMA + sample.names_dmelnv_MA)
        
      } else {
        
        design_EE <- model.matrix(~0 + groupingEE)
        design_MA <- model.matrix(~0 + groupingMA)
        
      }
      
    }
    
  } else {
    
    if (intercept) {
      
      if (batch) {
        
        design_EE <- model.matrix(~infection_EE * cell_EE)
        design_MA <- model.matrix(~infection_dmelnv_MA * cell_dmelnv_MA + sample.names_dmelnv_MA)
        
      } else {
        
        design_EE <- model.matrix(~infection_EE * cell_EE)
        design_MA <- model.matrix(~infection_dmelnv_MA * cell_dmelnv_MA)
        
      }
      
    } else {
      
      if (batch) {
        
        design_EE <- model.matrix(~0 + infection_EE * cell_EE + sample.names_thika)
        design_MA <- model.matrix(~0 + infection_dmelnv_MA * cell_dmelnv_MA + sample.names_dmelnv_MA)
        
      } else {
        
        design_EE <- model.matrix(~0 + infection_EE * cell_EE)
        design_MA <- model.matrix(~0 + infection_dmelnv_MA * cell_dmelnv_MA)
        
      }
      
    } 
  }
  
  colnames(design_EE) <- gsub("/", "", gsub("-", "_", gsub(":", ".", colnames(design_EE))))
  colnames(design_MA) <- gsub("/", "", gsub("-", "_", gsub(":", ".", colnames(design_MA))))
  
  fit_EE <- glm_gp(count.mtx_EE %>% as.matrix(), design_EE, on_disk = FALSE)
  fit_MA <- glm_gp(count.mtx_dmelnv_MA %>% as.matrix(), design_MA, on_disk = FALSE)
  
  return(list(EE = list(fit = fit_EE,
                        design = design_EE),
              MA = list(fit = fit_MA,
                        design = design_MA,
                        samples = sample.names_dmelnv_MA)))
}

# viral load correlation with glmGamPoi

cor_gGP <- function(EE.seu, MA.seu,
                    verbose = TRUE) {
  
  obj_tv <- subset(EE.seu, cells = colnames(EE.seu)[EE.seu$is.tv.infected])
  obj_dmelnv_EE <- subset(EE.seu, cells = colnames(EE.seu)[EE.seu$is.dmelnv.infected])
  obj_dmelnv_MA <- subset(MA.seu, cells = colnames(MA.seu)[MA.seu$is.dmelnv.infected])
  
  count_tv <- obj_tv@assays$CELLBENDER@counts
  cells_tv_EE <- obj_tv@active.ident %>% as.character()
  tv.pct <- obj_tv$percent.thika
  
  count_dmelnv_EE <- obj_dmelnv_EE@assays$CELLBENDER@counts
  cells_dmelnv_EE <- obj_dmelnv_EE@active.ident %>% as.character()
  dmelnvEE.pct <- obj_dmelnv_EE$percent.nora
  
  count_dmelnv_MA <- obj_dmelnv_MA@assays$CELLBENDER@counts
  cells_dmelnv_MA <- obj_dmelnv_MA@active.ident %>% as.character()
  batch_dmelnv_MA <- obj_dmelnv_MA$orig.ident
  dmelnvMA.pct <- obj_dmelnv_MA$percent.nora
  
  design_tv <- model.matrix(~tv.pct + cells_tv_EE)
  design_EE_dmelnv <- model.matrix(~dmelnvEE.pct + cells_dmelnv_EE)
  design_MA_dmelnv <- model.matrix(~dmelnvMA.pct + cells_dmelnv_MA)
  
  colnames(design_tv) <- gsub("/", "", gsub("-", "_", gsub(":", ".", colnames(design_tv))))
  colnames(design_EE_dmelnv) <- gsub("/", "", gsub("-", "_", gsub(":", ".", colnames(design_EE_dmelnv))))
  colnames(design_MA_dmelnv) <- gsub("/", "", gsub("-", "_", gsub(":", ".", colnames(design_MA_dmelnv))))
  
  if (verbose) print("fitting TV EE dataset")
  fit_tv <- glm_gp(count_tv %>% as.matrix(), design_tv, on_disk = FALSE)
  
  if (verbose) print("fitting DMelNV EE dataset")
  fit_dmelnv_EE <- glm_gp(count_dmelnv_EE %>% as.matrix(), design_EE_dmelnv, on_disk = FALSE)
  
  if (verbose) print("fitting DMelNV MA dataset")
  fit_dmelnv_MA <- glm_gp(count_dmelnv_MA %>% as.matrix(), design_MA_dmelnv, on_disk = FALSE)
  
  return(list(EE = list(TV = list(fit = fit_tv,
                                  design = design_tv),
                        DMelNV = list(fit = fit_dmelnv_EE,
                                      design = design_EE_dmelnv)),
              MA = list(fit = fit_dmelnv_MA,
                        design = design_MA_dmelnv,
                        samples = batch_dmelnv_MA)))
}

# function to do correlation tests but randomly sorting the viral load covariate

cor_scramble_gGP <- function(EE.seu, MA.seu,
                             seed,
                             verbose = TRUE) {
  
  obj_tv <- subset(EE.seu, cells = colnames(EE.seu)[EE.seu$is.tv.infected])
  obj_dmelnv_EE <- subset(EE.seu, cells = colnames(EE.seu)[EE.seu$is.dmelnv.infected])
  obj_dmelnv_MA <- subset(MA.seu, cells = colnames(MA.seu)[MA.seu$is.dmelnv.infected])
  
  count_tv <- obj_tv@assays$CELLBENDER@counts
  cells_tv_EE <- obj_tv@active.ident %>% as.character()
  tv.pct <- obj_tv$percent.thika
  
  count_dmelnv_EE <- obj_dmelnv_EE@assays$CELLBENDER@counts
  cells_dmelnv_EE <- obj_dmelnv_EE@active.ident %>% as.character()
  dmelnvEE.pct <- obj_dmelnv_EE$percent.nora
  
  count_dmelnv_MA <- obj_dmelnv_MA@assays$CELLBENDER@counts
  cells_dmelnv_MA <- obj_dmelnv_MA@active.ident %>% as.character()
  batch_dmelnv_MA <- obj_dmelnv_MA$orig.ident
  dmelnvMA.pct <- obj_dmelnv_MA$percent.nora
  
  # scramble pct vectors
  
  set.seed(seed)
  tv.pct <- sample(tv.pct)
  set.seed(seed)
  dmelnvEE.pct <- sample(dmelnvEE.pct)
  set.seed(seed)
  dmelnvMA.pct <- sample(dmelnvMA.pct)
  
  design_tv <- model.matrix(~tv.pct + cells_tv_EE)
  design_EE_dmelnv <- model.matrix(~dmelnvEE.pct + cells_dmelnv_EE)
  design_MA_dmelnv <- model.matrix(~dmelnvMA.pct + cells_dmelnv_MA)
  
  colnames(design_tv) <- gsub("/", "", gsub("-", "_", gsub(":", ".", colnames(design_tv))))
  colnames(design_EE_dmelnv) <- gsub("/", "", gsub("-", "_", gsub(":", ".", colnames(design_EE_dmelnv))))
  colnames(design_MA_dmelnv) <- gsub("/", "", gsub("-", "_", gsub(":", ".", colnames(design_MA_dmelnv))))
  
  if (verbose) print("fitting TV EE dataset")
  fit_tv <- glm_gp(count_tv %>% as.matrix(), design_tv, on_disk = FALSE)
  
  if (verbose) print("fitting DMelNV EE dataset")
  fit_dmelnv_EE <- glm_gp(count_dmelnv_EE %>% as.matrix(), design_EE_dmelnv, on_disk = FALSE)
  
  if (verbose) print("fitting DMelNV MA dataset")
  fit_dmelnv_MA <- glm_gp(count_dmelnv_MA %>% as.matrix(), design_MA_dmelnv, on_disk = FALSE)
  
  return(list(EE = list(TV = list(fit = fit_tv,
                                  design = design_tv),
                        DMelNV = list(fit = fit_dmelnv_EE,
                                      design = design_EE_dmelnv)),
              MA = list(fit = fit_dmelnv_MA,
                        design = design_MA_dmelnv,
                        samples = batch_dmelnv_MA)))
}

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

# Reactome function

reactome_pipe <- function(symbols) {
  
  entrezIDs <- gene_symbol[gene_symbol$V3 %in% symbols, "V2"]
  reactome_out <- enrichPathway(entrezIDs, pvalueCutoff = 1, readable = TRUE, organism = "fly")
  reactome_df <- as.data.frame(reactome_out)
  
  return(list(reactome_out, reactome_df))
}

# function to compare subgraphs to interactome

hub_fun <- function(graph,
                    subgraph,
                    title,
                    print.P = FALSE) {
  
  # calculate degree probability distribution and fit power law to graphs
  
  graph_d <- degree(graph, mode = "total", loops = FALSE)
  graph_dd <- degree.distribution(graph, cumulative = FALSE, mode = "total", loops = FALSE)
  graph_degree <- 1:max(graph_d)
  graph_prob_d <- graph_dd[-1]
  graph_prob_d <- graph_prob_d[which(graph_prob_d != 0)]
  graph_degree <- graph_degree[which(graph_prob_d != 0)]
  
  subgraph_d <- degree(subgraph, mode = "total", loops = FALSE)
  subgraph_dd <- degree.distribution(subgraph, cumulative = FALSE, mode = "total", loops = FALSE)
  subgraph_degree <- 1:max(subgraph_d)
  subgraph_prob_d <- subgraph_dd[-1]
  subgraph_prob_d <- subgraph_prob_d[which(subgraph_prob_d != 0)]
  subgraph_degree <- subgraph_degree[which(subgraph_prob_d != 0)]
  
  # linear regression in the log-log space
  
  graph_reg <- lm(log(graph_prob_d) ~ log(graph_degree))
  graph_cozf <- coef(graph_reg)
  graph_power_law <- function(x) exp(graph_cozf[[1]] + graph_cozf[[2]] * log(x))
  graph_slope <- -graph_cozf[[2]]
  
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
  plot(graph_prob_d ~ graph_degree, log = "xy", xlab = "degree (log)", ylab = "probability (log)",
       ylim = yrange,
       xlim = xrange,
       cex.main = 0.8,
       col = 1, main = title, pch = 15)
  curve(graph_power_law, col = "black", add = TRUE, n = length(graph_d))
  
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
  
  n_graph <- length(graph_prob_d)
  n_subgraph <- length(subgraph_prob_d)
  
  graph_slope_se <- summary(graph_reg)$coefficients[2, 2]
  subgraph_slope_se <- summary(subgraph_reg)$coefficients[2, 2]
  
  # test if the slope difference is positive or negative to set upper or lower tail
  
  if (-graph_cozf[[2]] + subgraph_cozf[[2]] < 0) lower.tail <- FALSE
  else lower.tail <- TRUE
  
  # slope difference
  
  slope_diff <- graph_cozf[[2]] - subgraph_cozf[[2]]
  
  # standard error difference
  
  se_diff <- sqrt(graph_slope_se^2 + subgraph_slope_se^2)
  
  # t statistic
  
  t.stat <- slope_diff/se_diff
  p <- 2*pt(t.stat, df = n_graph + n_subgraph - 4, lower.tail = lower.tail)
  
  # plot text
  
  if (print.P) {
    
    text(2, 2e-4, paste("diff = ", formatC(slope_diff, format = "e", digits = 2),
                        "\nP = ", formatC(p, format = "e", digits = 2), sep = ""))
  } else {
    
    text(2, 2e-4, paste("diff = ", formatC(slope_diff, digits = 2), sep = ""))
  }
  
  return(list(title = title, 
              subgraph_beta = subgraph_slope,
              slope_diff = slope_diff,
              slope_diff_p.value = p))
}

# do reactome analysis, compute betweenness and compare hubness for a gene symbols list list

network.pipe <- function(symbols.lst,
                         celltype = FALSE) {
  
  out.lst <- list()
  
  for (lst.index in 1:length(symbols.lst)) {
    
    if (celltype) lst <- symbols.lst[[lst.index]]$symbol
    else lst <- symbols.lst[[lst.index]]
    
    lst <- unique(lst)
    
    # check if symbols can be mapped to the interactome
    # skip if mapped genes < 3
    
    mapped.genes <- V(dro_net)$name[V(dro_net)$name %in% lst]
    
    if (length(mapped.genes) < 3) next
    
    reactome_results <- reactome_pipe(lst)
    subgraph.betweenness <- betweenness(dro_net, v = lst[lst %in% V(dro_net)$name])
    subgraph_mean_betweenness <- mean(subgraph.betweenness)
    betweenness_diff <- mean(subgraph.betweenness) - mean(dro_net_betweenness)
    betweenness_diff_p.value <- wilcox.test(subgraph.betweenness, dro_net_betweenness, alternative = "greater")$p.value
    
    subgraph <- graph_from_edgelist(as.matrix(dro_edges[dro_edges$from %in% lst |
                                                          dro_edges$to %in% lst, ]),
                                    directed = FALSE)
    
    hub <- hub_fun(dro_net, subgraph, names(symbols.lst[lst.index]))
    summary.row <- c(hub, list(subgraph_mean_betweenness = subgraph_mean_betweenness,
                                   mean_betweenness_diff = betweenness_diff,
                                   betweenness_diff_p.value = betweenness_diff_p.value))
    
    
    out.res <- list(reactome = reactome_results,
                    summary.row = summary.row)

    out.lst[[names(symbols.lst[lst.index])]] <- out.res
  }
  
  return(out.lst)
}

