#!/usr/bin/Rscript

library(Seurat)

#######################
## sensitivity curve ##
#######################

# ANOVA function to perform the analysis based on a random subset of cells from uninfected cells

random_aov_fun <- function(n.replicates = 10,
                           n.infected = 200,
                           return.table = FALSE) {
  
  # output will be a dataframe with the ANOVA results if return.table = true
  
  aov.list <- data.frame(gene = character(0),
                         group = character(0),
                         Df = numeric(0),
                         `Sum Sq` = numeric(0),
                         `Mean Sq` = numeric(0),
                         `F value` = numeric(0),
                         `Pr(>F)` = numeric(0),
                         p.adj = numeric(0),
                         replicate = numeric(0),
                         stringsAsFactors = FALSE)
  cell_names <- colnames(dmel_no_vir@assays$SCT)
  
  for (replicata in 1:n.replicates) {
    
    set.seed(replicata)
    
    # create "infected" object
    
    infected_cells <- sample(x = cell_names,
                             size = n.infected,
                             replace = FALSE)
    
    random_vir <- subset(dmel_no_vir, cells = infected_cells)
    
    # create "uninfected" object
    
    uninfected_cells <- cell_names[!cell_names %in% infected_cells]
    random_no_vir <- subset(dmel_no_vir, cells = uninfected_cells)
    
    # aov
    
    random_aov <- data.frame(gene = character(0),
                             group = character(0),
                             Df = numeric(0),
                             `Sum Sq` = numeric(0),
                             `Mean Sq` = numeric(0),
                             `F value` = numeric(0),
                             `Pr(>F)` = numeric(0),
                             stringsAsFactors = FALSE)
    
    for (gene in genes) {
      
      # log expression in uninfected cells
      
      gene_exp_uninfected <- data.frame(exp = random_no_vir@assays$SCT@data[gene,],
                                        cell = random_no_vir@active.ident)
      gene_exp_uninfected$infection <- replicate(nrow(gene_exp_uninfected), "uninfected")
      
      # log expression in infected cells
      
      gene_exp_infected <- data.frame(exp = random_vir@assays$SCT@data[gene,],
                                      cell = random_vir@active.ident)
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
      
      random_aov <- rbind(random_aov, tmp_tb)
    }
    # bonferroni correction
    
    random_aov$p.adj <- numeric(nrow(random_aov))
    random_aov[random_aov$group == "cell",]$p.adj <- p.adjust(random_aov[random_aov$group == "cell",]$Pr..F., "bonferroni" )
    random_aov[random_aov$group == "infection",]$p.adj <- p.adjust(random_aov[random_aov$group == "infection",]$Pr..F., "bonferroni" )
    random_aov[random_aov$group == "cell:infection",]$p.adj <- p.adjust(random_aov[random_aov$group == "cell:infection",]$Pr..F., "bonferroni" )
    
    random_aov$replicate <- replicate(nrow(random_aov), replicata)
    
    aov.list <- rbind(aov.list, random_aov)
  }
  if (return.table) return(aov.list)
  if (!return.table) {
    
    n.false_pos_infection <- nrow(aov.list[aov.list$p.adj < 0.05 & aov.list$group == "infection",])
    n.false_pos_cell.infection <- nrow(aov.list[aov.list$p.adj < 0.05 & aov.list$group == "cell:infection",])
    n.cell <- nrow(aov.list[aov.list$p.adj < 0.05 & aov.list$group == "cell",])
    
    return(list(n.false_pos_infection,
                n.false_pos_cell.infection,
                n.cell))
  }
}

load("mda_10x_uninfected.Rdata")
load("mda_genes.Rdata")

dmel_no_vir <- mda_10x_uninfected
genes <- mda_genes

sens.tb <- data.frame(step = numeric(0),
                      infection = numeric(0),
                      cell.infection = numeric(0),
                      cell = numeric(0))

for (i in seq(10, 500, by = 10)) {
  
  print(paste("step", i))
  
  rep.res <- random_aov_fun(n.replicates = 100,
                            n.infected = i)
  sens.tb[nrow(sens.tb) + 1,] <- c(i, rep.res)
  
  print(rep.res)
}

save(sens.tb, file = "mda_sens_tb.Rdata")


