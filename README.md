# Heterogeneity in the response of different subtypes of Drosophila melanogaster enteroendocrine cells to viral infections

### Description
This page contains all R scripts used for the analysis of viral infections on public single-cell RNA sequencing (scRNA-seq) data from Drosophila melanogaster enteroendocrine cells (EEs) and midgut cells.

### Overview
Raw data was downloaded from the sequence read archives (SRA) portal from NCBI (https://www.ncbi.nlm.nih.gov/sra; BioProject accession numbers: PRJNA547484 and PRJNA493298). Gene count matrices were generated with cellranger v3.1.0 using a concatenated reference genome containing the D. melanogaster, Thika virus (TV) and Drosophila melanogaster Nora virus (DMelNV) genomic sequences. Downstream analyses were performed with Seurat v3.1.4 with the R scripts present in the current project, with the exception of GO terms enrichment analysis which was conducted using BiNGO v3.0.4 plugin for Cytoscape v3.7.2.

#### Contents
##### run_analysis.R
Main R script in the project. This script will load all necessary libraries, do the preprocessing of the data and produce the lists of differentially expressed genes (DEGs) for the GO enrichment and network analyses.

##### GO.R
Load and processes the results of the GO enrichment analysis.

##### graph.R
Performs all network analyses.

##### make_figs.R and make_supplementary.R
Scripts to produce figures and supplementary material. For the manuscript, some figures were latter edited with Adobe Illustrator.

##### sensitivity_curve/
Contains all necessary data and scripts to perform the robustness test of DEGs analyses.

##### gene_matrices/
Contains count matrices, as well as a list of DEGs from bulk RNA-seq data from DMelNV-infected flies.

##### figures_out/
Output folder for all figures produced by make_figs.R and make_supplementary.R files.

##### GO/
Folder that contains all files necessary for the GO terms enrichment analysis. DEGs lists were produced by the run_analysis.R script.

##### network/
Contains the D. melanogaster interactome files.

##### name_id_symbol_list
File to convert gene IDs to symbol.

#### References
##### cellranger
10X Genomics, CA, USA
https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome

##### Seurat
Stuart, T., Butler, A., Hoffman, P., Hafemeister, C., Papalexi, E., Mauck III, W.M., Hao, Y., Stoeckius, M., Smibert, P. and Satija, R. (2019). Comprehensive integration of single-cell data. Cell 177, 1888-1902.

##### BiNGO
Maere, S., Heymans, K. and Kuiper, M. (2005). BiNGO: a Cytoscape plugin to assess overrepresentation of gene ontology categories in biological networks. Bioinformatics 21, 3448-3449.

##### Cytoscape
Shannon, P., Markiel, A., Ozier, O., Baliga, N.S., Wang, J.T., Ramage, D., Amin, N., Schwikowski, B. and Ideker, T. (2003). Cytoscape: a software environment for integrated models of biomolecular interaction networks. Genome Res. 13, 2498-2504.

##### scRNA-seq data
Guo, X., Yin, C., Yang, F., Zhang, Y., Huang, H., Wang, J., Deng, B., Cai T., Rao, Y., and Xi, R. (2019). The cellular diversity and transcription factor code of Drosophila enteroendocrine cells. Cell Rep. 29, 4172-4185.
Hung, R.J., Hu, Y., Kirchner, R., Liu, Y., Xu, C., Comjean, A., Tattikota, S.G., Li, F., Song, W., Sui, S.H. and Perrimon, N. (2020). A cell atlas of the adult Drosophila midgut. Proc. Natl. Acad. Sci. USA. 21, 1514-1523.

##### bulk RNA-seq data
Lopez, W., Page, A.M., Carlson, D.J., Ericson, B.L., Cserhati, M.F., Guda, C. and Carlson, K.A. (2018). Analysis of immune-related genes during Nora virus infection of Drosophila melanogaster using next generation sequencing. AIMS Microbiol. 4, 123.

##### D. melanogaster interactome
Ding, X.B., Jin, J., Tao, Y.T., Guo, W.P., Ruan, L., Yang, Q.L., Chen, P.C., Yao, H., Zhang, H.B. and Chen, X. (2020). Predicted Drosophila Interactome Resource and web tool for functional interpretation of differentially expressed genes. Database 2020, baaa005.
