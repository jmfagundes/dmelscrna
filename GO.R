#!/usr/bin/Rscript

# script to perform GO enrichment analysis

###############
## functions ##
###############

# simplify GO enrichment output
# exclude child terms if parent is also present
# exception: translation (GO:0006412) and neuropeptide signaling (GO:0007218)

simplify_go_terms <- function(tb) {
  goids <- tb$GO.ID
  out_tb <- tb
  child <- list()
  for (i in 1:length(goids)) {
    nzeros <- 7 - nchar(goids[[i]])
    goids[[i]] <- paste("GO:", paste(replicate(nzeros, "0"), collapse = ""), goids[[i]], sep = "")
    child <- c(child, GOBPCHILDREN[[goids[[i]]]])
  }
  out_tb$GO.ID <- goids
  child <- unique(child)
  exc <- out_tb[out_tb$GO.ID %in% list("GO:0006412", "GO:0007218"), ]
  out_tb <- out_tb[! out_tb$GO.ID %in% child, ]
  out_tb <- rbind(out_tb, exc)
  return(out_tb)
}

# function to draw GO enrichment results
# ggplot objects are created here but saved to a figure with make_figs.R

bingo_bubble <-function(pos = NULL,
                        neg,
                        cell.inf = NULL,
                        simplify = FALSE,
                        label = "pcor") {
  
  if (simplify) {
    
    if (! is.null(pos)) {
      pos <- simplify_go_terms(pos)
    }
    
    neg <- simplify_go_terms(neg)
    
    if (! is.null(cell.inf)) {
      cell.inf <- simplify_go_terms(cell.inf)
    }
    
  }
  
  if (is.null(cell.inf) & !is.null(pos)) {
    
    if (label == "pcor") {
      
      poslabel <- "positive correlation"
      neglabel <- "negative correlation"
      
    } else if (label == "aov") {
      
      poslabel <- "up-regulated"
      neglabel <- "down-regulated"
    }
    
    pos$sample <- replicate(nrow(pos), poslabel)
    neg$sample <- replicate(nrow(neg), neglabel)
    tb <- rbind(pos, neg)
    
  } else if (is.null(cell.inf) & is.null(pos)) {
    
    neg$sample <- replicate(nrow(neg), "negative correlation")
    tb <- neg
    
  } else if (is.null(pos) & !is.null(cell.inf)) {
    
    neg$sample <- replicate(nrow(neg), "down-regulated")
    cell.inf$sample <- replicate(nrow(cell.inf), "cell-subtype-dependent")
    tb <- rbind(neg, cell.inf)
    
  } else {
    
    pos$sample <- replicate(nrow(pos), "up-regulated")
    neg$sample <- replicate(nrow(neg), "down-regulated")
    cell.inf$sample <- replicate(nrow(cell.inf), "cell-subtype-dependent")
    tb <- rbind(pos, neg, cell.inf)
    
  }
  
  gg <- ggplot(tb, aes(x = sample, y = Description)) +
    geom_point(aes(color = corr.p.value, size = x)) +
    scale_color_viridis(direction = -1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab(label = "GO terms") + xlab(label = element_blank()) +
    guides(size = "none") + guides(colour = guide_colorbar("adjusted\nP-value"))
  return(gg)
}

##################
## EEs analysis ##
##################

# load enrichment results

# partial correlations

bingo_nora_pcor_neg <- read.table("GO/partial_correlation/GO_terms/nora_pcor_neg_go-basic.bgo",
                                  skip = 19, sep = "\t", header = TRUE)

# anova

bingo_nora_infection_pos <- read.table("GO/anova/GO_terms/nora_infection_pos_go-basic.bgo",
                                       skip = 19, sep = "\t", header = TRUE)
bingo_nora_infection_neg <- read.table("GO/anova/GO_terms/nora_infection_neg_go-basic.bgo",
                                       skip = 19, sep = "\t", header = TRUE)
bingo_nora_only_cell.infection <- read.table("GO/anova/GO_terms/nora_excusively_cell.infection_go-basic.bgo",
                                             skip = 19, sep = "\t", header = TRUE)

bingo_thika_infection_pos <- read.table("GO/anova/GO_terms/thika_infection_pos_go-basic.bgo",
                                        skip = 19, sep = "\t", header = TRUE)
bingo_thika_infection_neg <- read.table("GO/anova/GO_terms/thika_infection_neg_go-basic.bgo",
                                        skip = 19, sep = "\t", header = TRUE)
bingo_thika_only_cell.infection <- read.table("GO/anova/GO_terms/thika_excusively_cell.infection_go-basic.bgo",
                                              skip = 19, sep = "\t", header = TRUE)

# partial correlation GO enrichment

bingo_bubble_pcor_nora <- bingo_bubble(neg = bingo_nora_pcor_neg,
                                       simplify = TRUE)

# DEGs GO enrichment

bingo_bubble_aov_thika <- bingo_bubble(pos = bingo_thika_infection_pos, 
                                       neg = bingo_thika_infection_neg,
                                       cell.inf = bingo_thika_only_cell.infection,
                                       simplify = TRUE)

bingo_bubble_aov_nora <- bingo_bubble(pos = bingo_nora_infection_pos, 
                                      neg = bingo_nora_infection_neg,
                                      cell.inf = bingo_nora_only_cell.infection,
                                      simplify = TRUE)

# midgut atlas dataset

bingo_mda_pcor_neg <- read.table("GO/midgut_atlas/pcor_neg.bgo",
                                  skip = 19, sep = "\t", header = TRUE)
bingo_mda_pcor_pos <- read.table("GO/midgut_atlas/pcor_pos.bgo",
                                 skip = 19, sep = "\t", header = TRUE)

bingo_mda_infection_neg <- read.table("GO/midgut_atlas/aov_neg.bgo",
                                 skip = 19, sep = "\t", header = TRUE)
bingo_mda_cell.infection <- read.table("GO/midgut_atlas/aov_cell.infection.bgo",
                                      skip = 19, sep = "\t", header = TRUE)

bingo_bubble_pcor_mda <- bingo_bubble(neg = bingo_mda_pcor_neg,
                                      pos = bingo_mda_pcor_pos,
                                      simplify = TRUE)

bingo_bubble_aov_mda <- bingo_bubble(neg = bingo_mda_pcor_neg,
                                     cell.inf = bingo_mda_cell.infection,
                                     simplify = TRUE)


