# INIT --------------------------------------------------------------------

## Installation of Bioconductor package management system
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# # install packages
# BiocManager::install(c("diffuStats","STRINGdb","supraHex","Rgraphviz","graph"))
# BiocManager::install("Category")
# 
# # 'igraph' is a CRAN Paket and can be installed via the integrated package management system
# install.packages("igraph")
# # 'dnet'
# install.packages("dnet", repos="http://R-Forge.R-project.org")
# install.packages("dnet")
# # 'openxlsx'
# install.packages("openxlsx", dependencies=TRUE)

source("scripts/_network_diffusion.R")
source("scripts/_annovar_filtering.R")
set.seed(16341)
data.dir <- "data"

# PREPARE GENE-SETS -------------------------------------------------------

# Set-up local storage for gene-sets
eh = ExperimentHub()
query(eh , 'msigdb')

# Download molecular signatures database (MSigDB) hosted on the ExperimentHub 
msigdb.hs = getMsigdb(org = 'hs', id = 'SYM', version = '7.5.1')

# Include KEGG-DB (required workaround due to licensing issues)
msigdb.hs = appendKEGG(msigdb.hs)

# The output of this command makes more sense when you take a look at https://www.gsea-msigdb.org/gsea/msigdb
listCollections(msigdb.hs)

# List the names of the sub-collections, i.e., show all databases available
listSubCollections(msigdb.hs)

# select C2 (among others REACTOME) and C5 (GO terms)
Hs.c2 <- subsetCollection(msigdb.hs, 'c2')

# create subsets for REACTOME and KEGG gene-sets
Hs.c2.REACTOME <- geneIds(Hs.c2[grep("REACTOME_", names(Hs.c2))])
Hs.c2.KEGG <- geneIds(Hs.c2[grep("KEGG_", names(Hs.c2))])
Hs.c2.CP <- geneIds(Hs.c2[grep("CP_", names(Hs.c2))])

## aaand oncosets 
Hs.c5 <- subsetCollection(msigdb.hs, 'c5')
Hs.c5.GO <- geneIds(Hs.c5[grep("GOBP_", names(Hs.c5))])
Hs.c5.HP <- geneIds(Hs.c5[grep("HP_", names(Hs.c5))])

names(Hs.c5)



# PREPARE NETWORK ---------------------------------------------------------

# This will create a folder in the work-directory, download the stringDB and
# prepare the graph object. Takes a while when executed for the first time.
prepareSTRINGdb()

# The previous function created the stringDB object. It provies access to stringDB
# functionality with 'stringDB$string.DB' and contains a graph of the PPI
#names(stringDB)

# The previous call the 'stringDB' object including the PPI database and the 
# corresponding graph object. We will limit the vertex-degree to 900 in order
# to facilitate faster computation and reduce noise.
#gi <- getSubDegree(stringDB$string.g, 900)

# And save the subnetwork for easy access. 
#saveRDS(gi, file = file.path(FILESOURCE, "RDS", "network_d900.rds"), compress = TRUE)


# PREPARE KERNEL ----------------------------------------------------------

# load subnetwork
#  gi <- readRDS(file=file.path(FILESOURCE, "RDS", "network_d900.rds"), refhook = NULL)
# #
# # Take all the nodes in the network, or take a subset, if desired
# nodes_mapped <- V(gi)$name
# # Take the largest connected component of the network and remove self loops
# network <- dnet::dNetInduce(g=gi, nodes_query=nodes_mapped, knn=0, remove.loops=T,
#                             largest.comp=T)
# # save-it
# saveRDS(network, file = file.path(FILESOURCE, "RDS", "network.rds"), compress = TRUE)
# 
# #network <- readRDS(file = file.path(FILESOURCE, "RDS", "network.rds"), refhook = NULL)
# 
# ### BUILD AND STORE KERNEL
# # Compute the Laplacian Kernel, takes long, but do this only once
# K_rl <- diffuStats::regularisedLaplacianKernel(network)
# # Save the Kernel
# saveRDS(K_rl, file = file.path(FILESOURCE, "RDS", "network_kernel.rds"),
#         compress = TRUE)
### BUILD AND STORE KERNEL

# SINGLE SAMPLE + FILTERING FOR CONTROL VARIANTS --------------------------

## For every single sample filter against the variants present matched controls.
## Then repeat diffusion for single samples and for cohort gene overlap and compare
## to the original output.

meta <- readRDS(file.path("data/RDS", "VRI_Pilot_meta_v3.RDS"))

## 0. load data
## 1. get a sample
## 2. find all controls for that sample
## 3. retrieve controls from control cohort and make sure that only shared variants are retained

# full ctrl and case cohort
df.cohort <- readRDS(file.path("data/RDS", "VRI_Pilot_Cohort_annotated.filtered.26122024.RDS"))
# single samples, frequency filtered, but not background removed
sample.list <- readRDS(file = file.path("data", "RDS", "VRI_Pilot_SingleSample_annotated.filtered.26122024.RDS"))
# meta-data
meta <- readRDS(file.path("data/RDS", "VRI_Pilot_meta_v3.RDS"))

df.noBG <- list()
df.noBG[["breast"]] <- removeCtrlBackgroundMatched(sample.list$breast, df.cohort$breast_ctrl, df.meta=meta, match.col="m_breast_4")
df.noBG[["colon"]] <- removeCtrlBackgroundMatched(sample.list$colon, df.cohort$colon_ctrl, df.meta=meta, match.col="m_colon_4")
df.noBG[["leukemia"]] <- removeCtrlBackgroundMatched(sample.list$leukemia, df.cohort$leukemia_ctrl, df.meta=meta, match.col="m_leukemia_4")

## sanity
dim(df.noBG[["breast"]]$DXB_HLA_20_102468281_00426$potential)
dim(sample.list$breast$DXB_HLA_20_102468281_00426$potential)

saveRDS(df.noBG, file = file.path("data", "RDS", "VRI_Pilot_SingleSample_annotated.filtered.noBG.26122024.RDS"), compress = TRUE)

# SINGLE SAMPLE DIFFUSION -------------------------------------------------


## diffusion network and kernel
network <- readRDS(file = file.path(FILESOURCE, "RDS","network.rds"), refhook = NULL)
network_kernel <- readRDS(file = file.path(FILESOURCE,"RDS", "network_kernel.rds"), 
                          refhook = NULL)

## load the filtered sample list that only contains breast, colon, and leukemia single sample dfs
sample.list <- readRDS(file = file.path("data", "RDS", "VRI_Pilot_SingleSample_annotated.filtered.26122024.RDS"))

### RUN DIffusion
reactome.df <- createMergeDF()
kegg.df <- createMergeDF()
go.df <- createMergeDF()
hp.df <- createMergeDF()

for( leCohort in names(sample.list) ) {
  
  cohort.list <- sample.list[[eval(leCohort)]]
  
  for( leSample in names(cohort.list) ) {
    
    m = unlist(strsplit(leSample, '_'))
    df.list <- cohort.list[[eval(leSample)]]
    
    # score variants and prepare diffusion input
    input_oncogenes <- prepDiffSet(df.list)
    
    # take remaining Genes and create a binary input_vector for the diffusion
    tmp_input_vec <- rep(0,length(igraph::V(network)$name))
    names(tmp_input_vec) <- igraph::V(network)$name
    idx <- intersect( names(tmp_input_vec), as.character(input_oncogenes$node_id))
    
    # in case the remaining genes in input_oncogenes cannot be mapped to the network 
    # remove unmapped genes from input_oncogenes
    input_oncogenes <- input_oncogenes[which(input_oncogenes$node_id %in% idx),]
    
    ## draw network of input genes for diffusion
    input_oncogenes.mapping <- stringDB$string.DB$map(input_oncogenes, "node_id")
    
    # display a STRING network png with the "halo"
    pdf(file = file.path(FILESOURCE, "output", paste(leSample,"_input_diffnet.pdf", sep="")), 
        height = 16, width = 16, useDingbats = F, onefile=TRUE)
    
    # finally, use stringDB functionality to produce the plot
    stringDB$string.DB$plot_network( input_oncogenes.mapping$STRING_id )
    
    dev.off()
    
    # and set the input value for the oncogenes
    tmp_input_vec[match(input_oncogenes$node_id, names(tmp_input_vec))] <- 
      input_oncogenes$node_score
    
    # Kernel based diffusion with z-method
    df_z <- diffuStats::diffuse_grid(K = network_kernel, scores=tmp_input_vec,
                                     grid_param = expand.grid(method = "z"), 
                                     n.perm = 1)
    
    s.genes.z <- df_z[abs(df_z$node_score) > 0.5,]
    
    ## stringDB provides 'map()' to match gene symbols to Ensembl-Protein-IDs
    s.genes.z.mapping <- stringDB$string.DB$map(s.genes.z, "node_id")
    
    ## Export results for later use.
    openxlsx::write.xlsx(s.genes.z.mapping,
                         file = file.path(FILESOURCE,"output",paste(leSample,"_diffusion_zscore.xlsx", sep="")),
                         rowNames = TRUE,colNames=TRUE)
    
    gage.input <- s.genes.z$node_score
    names(gage.input) <- s.genes.z$node_id
    
    ## REACTOME
    gage.REACTOME <- gage(gage.input,gsets = Hs.c2.REACTOME,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
    gage.REACTOME <- sigGeneSet(gage.REACTOME, qpval=c("q.val"), cutoff=0.05)
    ## KEGG
    gage.KEGG <- gage(gage.input,gsets = Hs.c2.KEGG,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
    gage.KEGG <- sigGeneSet(gage.KEGG, qpval=c("q.val"), cutoff=0.05)
    ## GO
    gage.GO <- gage(gage.input,gsets = Hs.c5.GO,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
    gage.GO <- sigGeneSet(gage.GO, qpval=c("q.val"), cutoff=0.05)
    ## HPO Human Phenotype Ontology
    gage.HP <- gage(gage.input,gsets = Hs.c5.HP,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
    gage.HP <- sigGeneSet(gage.HP, qpval=c("q.val"), cutoff=0.05)
    
    reactome.df <- mergeGage( reactome.df, gage.REACTOME, s.idx=leSample, c.idx=leCohort )
    kegg.df <- mergeGage( kegg.df, gage.KEGG, s.idx=leSample, c.idx=leCohort )
    go.df <- mergeGage( go.df, gage.GO, s.idx=leSample, c.idx=leCohort )
    hp.df <- mergeGage( hp.df, gage.HP, s.idx=leSample, c.idx=leCohort )
    
    
  }
  
}


reactome.df$greater$pathway <- gsub("REACTOME_", "",reactome.df$greater$pathway)
reactome.df$less$pathway <- gsub("REACTOME_", "",reactome.df$less$pathway)
kegg.df$greater$pathway <- gsub("KEGG_", "",kegg.df$greater$pathway)
kegg.df$less$pathway <- gsub("KEGG_", "",kegg.df$less$pathway)
go.df$greater$pathway <- gsub("GOBP_", "",go.df$greater$pathway)
go.df$less$pathway <- gsub("GOBP_", "",go.df$less$pathway)
hp.df$greater$pathway <- gsub("HP_", "",hp.df$greater$pathway)
hp.df$less$pathway <- gsub("HP_", "",hp.df$less$pathway)


saveRDS(reactome.df, file = file.path("data", "RDS", "singlesample_withBG_reactome.df.RDS"), compress = TRUE)
saveRDS(kegg.df, file = file.path("data", "RDS", "singlesample_withBG_kegg.df.RDS"), compress = TRUE)
saveRDS(go.df, file = file.path("data", "RDS", "singlesample_withBG_go.df.RDS"), compress = TRUE)
saveRDS(hp.df, file = file.path("data", "RDS", "singlesample_withBG_hp.df.RDS"), compress = TRUE)

## write to file
openxlsx::write.xlsx(reactome.df, file = file.path("data", "output", "SingleSample_cancers_withBG_GAGE_REACTOME.xlsx"))
openxlsx::write.xlsx(kegg.df, file = file.path("data", "output", "SingleSample_cancers_withBG_GAGE_KEGG.xlsx"))
openxlsx::write.xlsx(go.df, file = file.path("data", "output", "SingleSample_cancers_withBG_GAGE_GO.xlsx"))
openxlsx::write.xlsx(hp.df, file = file.path("data", "output", "SingleSample_cancers_withBG_GAGE_HP.xlsx"))


## REACTOME
reactome.df$greater <- reactome.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

reactome.df$less <- reactome.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(reactome.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_withBG_GAGE_REACTOME.xlsx"))

## KEGG
kegg.df$greater <- kegg.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

kegg.df$less <- kegg.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(kegg.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_withBG_GAGE_KEGG.xlsx"))

## GO
go.df$greater <- go.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

go.df$less <- go.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(go.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_withBG_GAGE_GO.xlsx"))

## HP
hp.df$greater <- hp.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

hp.df$less <- hp.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(hp.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_withBG_GAGE_HP.xlsx"))










# BACKGROUND FILTERED SINGLE SAMPLE DIFFUSION -----------------------------

## diffusion network and kernel
network <- readRDS(file = file.path(FILESOURCE, "RDS","network.rds"), refhook = NULL)
network_kernel <- readRDS(file = file.path(FILESOURCE,"RDS", "network_kernel.rds"), 
                          refhook = NULL)

## load the filtered sample list that only contains breast, colon, and leukemia single sample dfs
sample.list <- readRDS(file = file.path("data", "RDS", "VRI_Pilot_SingleSample_annotated.filtered.noBG.26122024.RDS"))

### RUN DIffusion
reactome.df <- createMergeDF()
kegg.df <- createMergeDF()
go.df <- createMergeDF()
hp.df <- createMergeDF()

for( leCohort in names(sample.list) ) {
  
  cohort.list <- sample.list[[eval(leCohort)]]
  
  for( leSample in names(cohort.list) ) {
    
    m = unlist(strsplit(leSample, '_'))
    df.list <- cohort.list[[eval(leSample)]]
    
    # score variants and prepare diffusion input
    input_oncogenes <- prepDiffSet(df.list)
    
    # take remaining Genes and create a binary input_vector for the diffusion
    tmp_input_vec <- rep(0,length(igraph::V(network)$name))
    names(tmp_input_vec) <- igraph::V(network)$name
    idx <- intersect( names(tmp_input_vec), as.character(input_oncogenes$node_id))
    
    # in case the remaining genes in input_oncogenes cannot be mapped to the network 
    # remove unmapped genes from input_oncogenes
    input_oncogenes <- input_oncogenes[which(input_oncogenes$node_id %in% idx),]
    
    ## draw network of input genes for diffusion
    input_oncogenes.mapping <- stringDB$string.DB$map(input_oncogenes, "node_id")
    
    # display a STRING network png with the "halo"
    pdf(file = file.path(FILESOURCE, "output", paste(leSample,"_backgroundRemoved_input_diffnet.pdf", sep="")), 
        height = 16, width = 16, useDingbats = F, onefile=TRUE)
    
    # finally, use stringDB functionality to produce the plot
    stringDB$string.DB$plot_network( input_oncogenes.mapping$STRING_id )
    
    dev.off()
    
    # and set the input value for the oncogenes
    tmp_input_vec[match(input_oncogenes$node_id, names(tmp_input_vec))] <- 
      input_oncogenes$node_score
    
    # Kernel based diffusion with z-method
    df_z <- diffuStats::diffuse_grid(K = network_kernel, scores=tmp_input_vec,
                                     grid_param = expand.grid(method = "z"), 
                                     n.perm = 1)
    
    s.genes.z <- df_z[abs(df_z$node_score) > 0.5,]
    
    ## stringDB provides 'map()' to match gene symbols to Ensembl-Protein-IDs
    s.genes.z.mapping <- stringDB$string.DB$map(s.genes.z, "node_id")
    
    ## Export results for later use.
    openxlsx::write.xlsx(s.genes.z.mapping,
                         file = file.path(FILESOURCE,"output",paste(leSample,"_backgroundRemoved_diffusion_zscore.xlsx", sep="")),
                         rowNames = TRUE,colNames=TRUE)
    
    gage.input <- s.genes.z$node_score
    names(gage.input) <- s.genes.z$node_id
    
    ## REACTOME
    gage.REACTOME <- gage(gage.input,gsets = Hs.c2.REACTOME,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
    gage.REACTOME <- sigGeneSet(gage.REACTOME, qpval=c("q.val"), cutoff=0.05)
    ## KEGG
    gage.KEGG <- gage(gage.input,gsets = Hs.c2.KEGG,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
    gage.KEGG <- sigGeneSet(gage.KEGG, qpval=c("q.val"), cutoff=0.05)
    ## GO
    gage.GO <- gage(gage.input,gsets = Hs.c5.GO,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
    gage.GO <- sigGeneSet(gage.GO, qpval=c("q.val"), cutoff=0.05)
    ## HPO Human Phenotype Ontology
    gage.HP <- gage(gage.input,gsets = Hs.c5.HP,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
    gage.HP <- sigGeneSet(gage.HP, qpval=c("q.val"), cutoff=0.05)
    
    reactome.df <- mergeGage( reactome.df, gage.REACTOME, s.idx=leSample, c.idx=leCohort )
    kegg.df <- mergeGage( kegg.df, gage.KEGG, s.idx=leSample, c.idx=leCohort )
    go.df <- mergeGage( go.df, gage.GO, s.idx=leSample, c.idx=leCohort )
    hp.df <- mergeGage( hp.df, gage.HP, s.idx=leSample, c.idx=leCohort )
    
    
  }
  
}


reactome.df$greater$pathway <- gsub("REACTOME_", "",reactome.df$greater$pathway)
reactome.df$less$pathway <- gsub("REACTOME_", "",reactome.df$less$pathway)
kegg.df$greater$pathway <- gsub("KEGG_", "",kegg.df$greater$pathway)
kegg.df$less$pathway <- gsub("KEGG_", "",kegg.df$less$pathway)
go.df$greater$pathway <- gsub("GOBP_", "",go.df$greater$pathway)
go.df$less$pathway <- gsub("GOBP_", "",go.df$less$pathway)
hp.df$greater$pathway <- gsub("HP_", "",hp.df$greater$pathway)
hp.df$less$pathway <- gsub("HP_", "",hp.df$less$pathway)


saveRDS(reactome.df, file = file.path("data", "RDS", "backgroundRemoved_reactome.df.RDS"), compress = TRUE)
saveRDS(kegg.df, file = file.path("data", "RDS", "backgroundRemoved_kegg.df.RDS"), compress = TRUE)
saveRDS(go.df, file = file.path("data", "RDS", "backgroundRemoved_go.df.RDS"), compress = TRUE)
saveRDS(hp.df, file = file.path("data", "RDS", "backgroundRemoved_hp.df.RDS"), compress = TRUE)

## write to file
openxlsx::write.xlsx(reactome.df, file = file.path("data", "output", "SingleSample_cancers_noBG_GAGE_REACTOME.xlsx"))
openxlsx::write.xlsx(kegg.df, file = file.path("data", "output", "SingleSample_cancers_noBG_GAGE_KEGG.xlsx"))
openxlsx::write.xlsx(go.df, file = file.path("data", "output", "SingleSample_cancers_noBG_GAGE_GO.xlsx"))
openxlsx::write.xlsx(hp.df, file = file.path("data", "output", "SingleSample_cancers_noBG_GAGE_HP.xlsx"))


## REACTOME
reactome.df$greater <- reactome.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

reactome.df$less <- reactome.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(reactome.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_noBG_GAGE_REACTOME.xlsx"))

## KEGG
kegg.df$greater <- kegg.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

kegg.df$less <- kegg.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(kegg.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_noBG_GAGE_KEGG.xlsx"))

## GO
go.df$greater <- go.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

go.df$less <- go.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(go.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_noBG_GAGE_GO.xlsx"))

## HP
hp.df$greater <- hp.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

hp.df$less <- hp.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(hp.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_noBG_GAGE_HP.xlsx"))



# PLOT RESULTS - BACKGROUND REMOVED ---------------------------------------

reactome.df <- readRDS(file = file.path("data", "RDS", "backgroundRemoved_reactome.df.RDS"))
kegg.df <- readRDS(file = file.path("data", "RDS", "backgroundRemoved_kegg.df.RDS"))
go.df <- readRDS(file = file.path("data", "RDS", "backgroundRemoved_go.df.RDS"))
hp.df <- readRDS(file = file.path("data", "RDS", "backgroundRemoved_hp.df.RDS"))

reactome.df$greater <- reactome.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia"))
reactome.df$less <- reactome.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia"))

reactome.up <- reactome.df$greater %>%
  filter(pathway %in% names(table(reactome.df$greater$pathway))[which(table(reactome.df$greater$pathway) > 6)])
reactome.down <- reactome.df$less %>%
  filter(pathway %in% names(table(reactome.df$less$pathway))[which(table(reactome.df$less$pathway) > 5)])

## REACTOME
pdf(file.path("data", "output", "SingleSample_cancers_filtered_GAGE_REACTOME.pdf"), width = 10, height = 16)
plotPathwayBar(reactome.up, xlabel="Number of UP-regulated pathways", ylabel="REACTOME pathways" , w.label = 55 )
plotPathwayBar(reactome.down, xlabel="Number of DOWN-regulated pathways", ylabel="REACTOME pathways", w.label = 55  )
dev.off()


kegg.up <- kegg.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(kegg.df$greater$pathway))[which(table(kegg.df$greater$pathway) > 6)])
kegg.down <- kegg.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(kegg.df$less$pathway))[which(table(kegg.df$less$pathway) > 4)])


## KEGG
pdf(file.path("data", "output", "SingleSample_cancers_filtered_GAGE_KEGG.pdf"), width = 10, height = 10)
plotPathwayBar(kegg.up, xlabel="Number of UP-regulated pathways", ylabel="KEGG pathways" )
plotPathwayBar(kegg.down, xlabel="Number of DOWN-regulated pathways", ylabel="KEGG pathways" )
dev.off()

## GO
## subset to reasonable numbers
go.up <- go.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(go.df$greater$pathway))[which(table(go.df$greater$pathway) > 10)])
go.down <- go.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(go.df$less$pathway))[which(table(go.df$less$pathway) > 10)])

pdf(file.path("data", "output", "SingleSample_cancers_filtered_GAGE_GO.pdf"), width = 10, height = 15)
plotPathwayBar(go.up, xlabel="Number of UP-regulated pathways", ylabel="GO pathways", w.label = 55 )
plotPathwayBar(go.down, xlabel="Number of DOWN-regulated pathways", ylabel="GO pathways", w.label = 55 )
dev.off()

## HP
## subset to reasonable numbers
hp.up <- hp.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(hp.df$greater$pathway))[which(table(hp.df$greater$pathway) > 1)])
hp.down <- hp.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(hp.df$less$pathway))[which(table(hp.df$less$pathway) > 6)])

pdf(file.path("data", "output", "SingleSample_cancers_filtered_GAGE_HP.pdf"), width = 10, height = 14)
plotPathwayBar(hp.up, xlabel="Number of UP-regulated pathways", ylabel="HP pathways", w.label = 55 )
plotPathwayBar(hp.down, xlabel="Number of DOWN-regulated pathways", ylabel="HP pathways", w.label = 55 )
dev.off()







































# SINGLE-SAMPLE DIFFUSION -------------------------------------------------

### LOAD the data
## diffusion network and kernel
network <- readRDS(file = file.path(FILESOURCE, "RDS","network.rds"), refhook = NULL)
network_kernel <- readRDS(file = file.path(FILESOURCE,"RDS", "network_kernel.rds"), 
                          refhook = NULL)

### ToDo: load le samples and do diffusion things
samples <- list.files(path = file.path(FILESOURCE, "RAW","single_sample_annotation"), pattern = "multianno.txt")
tbl.mutations <- readRDS(file = file.path("data", "RDS", "VRI_mutations_curated.RDS"))




sample.list <- list()

for( leSample in 1:length(samples) ) {
  
  m = unlist(strsplit(gsub(".hg38_multianno.txt","",samples[leSample]), '_'))
  df <- data.table::fread(file.path(FILESOURCE, "RAW","single_sample_annotation", samples[leSample]), header = TRUE, 
                          na.strings = ".")
  ## filter for exac frequencies smaller than 0.001
  if( m[1] == "breast" ) {
    ass.table <- tbl.mutations$breast$Gene
  } else if( m[1] == "colon" ) {
    ass.table <- tbl.mutations$colon$Gene
  } else if( m[1] == "leukemia" ) {
    ass.table <- tbl.mutations$leukemia$Gene
  } else {
    ass.table = NULL
  }
  
  df.list <- prepATable(df, 0.001, ass.table)
  
  sample.list[[eval(gsub(".hg38_multianno.txt","",samples[leSample]))]] <- df.list
  
  #names(df.list); dim(df.list$pathogenic); dim(df.list$potential); dim(df.list$associated); dim(df.list$unscored); dim(df.list$benign)
  
  #openxlsx::write.xlsx(df.list, file = file.path(FILESOURCE, "output", paste(m[1],"_",m[2],"_",m[3],"_variants.fitlered.xlsx", sep="")))
  
}

saveRDS(sample.list, file = file.path("data", "RDS", "samples_annotated.RDS"), compress = TRUE)


## !!!!! THIS NEEDS ATTENtiON BEForE RUNNINg BECAUSE OF SWITCH IN NAMING !!!!!
### RUN DIffusion
reactome.df <- createMergeDF()
kegg.df <- createMergeDF()
go.df <- createMergeDF()
hp.df <- createMergeDF()

for( leSample in names(sample.list) ) {
  
  
  df.list <- sample.list[[eval(leSample)]]
  
  # score variants and prepare diffusion input
  input_oncogenes <- prepDiffSet(df.list)
  
  # take remaining Genes and create a binary input_vector for the diffusion
  tmp_input_vec <- rep(0,length(igraph::V(network)$name))
  names(tmp_input_vec) <- igraph::V(network)$name
  idx <- intersect( names(tmp_input_vec), as.character(input_oncogenes$node_id))
  
  # in case the remaining genes in input_oncogenes cannot be mapped to the network 
  # remove unmapped genes from input_oncogenes
  input_oncogenes <- input_oncogenes[which(input_oncogenes$node_id %in% idx),]
  
  ## draw network of input genes for diffusion
  input_oncogenes.mapping <- stringDB$string.DB$map(input_oncogenes, "node_id")
  
  # display a STRING network png with the "halo"
  pdf(file = file.path(FILESOURCE, "output", paste(m[1],"_",m[2],"_",m[3],"_input_diffnet.pdf", sep="")), 
      height = 16, width = 16, useDingbats = F, onefile=TRUE)
  
  # finally, use stringDB functionality to produce the plot
  stringDB$string.DB$plot_network( input_oncogenes.mapping$STRING_id )
  
  dev.off()
  
  # and set the input value for the oncogenes
  tmp_input_vec[match(input_oncogenes$node_id, names(tmp_input_vec))] <- 
    input_oncogenes$node_score
  
  # Kernel based diffusion with z-method
  df_z <- diffuStats::diffuse_grid(K = network_kernel, scores=tmp_input_vec,
                                   grid_param = expand.grid(method = "z"), 
                                   n.perm = 1)
  
  s.genes.z <- df_z[abs(df_z$node_score) > 0.5,]
  
  ## stringDB provides 'map()' to match gene symbols to Ensembl-Protein-IDs
  s.genes.z.mapping <- stringDB$string.DB$map(s.genes.z, "node_id")
  
  ## Export results for later use.
  openxlsx::write.xlsx(s.genes.z.mapping,
                       file = file.path(FILESOURCE,"output",paste(m[1],"_",m[2],"_",m[3],"_diffusion_zscore.xlsx", sep="")),
                       rowNames = TRUE,colNames=TRUE)
  
  gage.input <- s.genes.z$node_score
  names(gage.input) <- s.genes.z$node_id
  
  ## REACTOME
  gage.REACTOME <- gage(gage.input,gsets = Hs.c2.REACTOME,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
  gage.REACTOME <- sigGeneSet(gage.REACTOME, qpval=c("q.val"), cutoff=0.05)
  ## KEGG
  gage.KEGG <- gage(gage.input,gsets = Hs.c2.KEGG,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
  gage.KEGG <- sigGeneSet(gage.KEGG, qpval=c("q.val"), cutoff=0.05)
  ## GO
  gage.GO <- gage(gage.input,gsets = Hs.c5.GO,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
  gage.GO <- sigGeneSet(gage.GO, qpval=c("q.val"), cutoff=0.05)
  ## HPO Human Phenotype Ontology
  gage.HP <- gage(gage.input,gsets = Hs.c5.HP,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
  gage.HP <- sigGeneSet(gage.HP, qpval=c("q.val"), cutoff=0.05)
  
  reactome.df <- mergeGage( reactome.df, gage.REACTOME, s.idx=m[3], c.idx=m[1] )
  kegg.df <- mergeGage( kegg.df, gage.KEGG, s.idx=m[3], c.idx=m[1] )
  go.df <- mergeGage( go.df, gage.GO, s.idx=m[3], c.idx=m[1] )
  hp.df <- mergeGage( hp.df, gage.HP, s.idx=m[3], c.idx=m[1] )
  
  
}

reactome.df$greater$pathway <- gsub("REACTOME_", "",reactome.df$greater$pathway)
reactome.df$less$pathway <- gsub("REACTOME_", "",reactome.df$less$pathway)
kegg.df$greater$pathway <- gsub("KEGG_", "",kegg.df$greater$pathway)
kegg.df$less$pathway <- gsub("KEGG_", "",kegg.df$less$pathway)
go.df$greater$pathway <- gsub("GOBP_", "",go.df$greater$pathway)
go.df$less$pathway <- gsub("GOBP_", "",go.df$less$pathway)
hp.df$greater$pathway <- gsub("HP_", "",hp.df$greater$pathway)
hp.df$less$pathway <- gsub("HP_", "",hp.df$less$pathway)


saveRDS(reactome.df, file = file.path("data", "RDS", "reactome.df.RDS"), compress = TRUE)
saveRDS(kegg.df, file = file.path("data", "RDS", "kegg.df.RDS"), compress = TRUE)
saveRDS(go.df, file = file.path("data", "RDS", "go.df.RDS"), compress = TRUE)
saveRDS(hp.df, file = file.path("data", "RDS", "hp.df.RDS"), compress = TRUE)


# SAVE RESULTS ------------------------------------------------------------

reactome.df <- readRDS(file = file.path("data", "RDS", "reactome.df.RDS"))
kegg.df <- readRDS(file = file.path("data", "RDS", "kegg.df.RDS"))
go.df <- readRDS(file = file.path("data", "RDS", "go.df.RDS"))
hp.df <- readRDS(file = file.path("data", "RDS", "hp.df.RDS"))

## wider .. if this is needed
# test <- reshape(reactome.df.breast, idvar = "pathway", timevar = "sID", direction = "wide")

## write to file
openxlsx::write.xlsx(reactome.df, file = file.path("data", "output", "SingleSample_cancers_GAGE_REACTOME.xlsx"))
openxlsx::write.xlsx(kegg.df, file = file.path("data", "output", "SingleSample_cancers_GAGE_KEGG.xlsx"))
openxlsx::write.xlsx(go.df, file = file.path("data", "output", "SingleSample_cancers_GAGE_GO.xlsx"))
openxlsx::write.xlsx(hp.df, file = file.path("data", "output", "SingleSample_cancers_GAGE_HP.xlsx"))


# TABULARIZE AND SAVE -----------------------------------------------------

reactome.df <- readRDS(file = file.path("data", "RDS", "reactome.df.RDS"))
kegg.df <- readRDS(file = file.path("data", "RDS", "kegg.df.RDS"))
go.df <- readRDS(file = file.path("data", "RDS", "go.df.RDS"))
hp.df <- readRDS(file = file.path("data", "RDS", "hp.df.RDS"))

## REACTOME
reactome.df$greater <- reactome.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  dplyr::count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

reactome.df$less <- reactome.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(reactome.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_GAGE_REACTOME.xlsx"))

## KEGG
kegg.df$greater <- kegg.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

kegg.df$less <- kegg.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(kegg.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_GAGE_KEGG.xlsx"))

## GO
go.df$greater <- go.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

go.df$less <- go.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(go.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_GAGE_GO.xlsx"))

## HP
hp.df$greater <- hp.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

hp.df$less <- hp.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(hp.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_GAGE_HP.xlsx"))


# PLOT RESULTS ------------------------------------------------------------

reactome.df <- readRDS(file = file.path("data", "RDS", "reactome.df.RDS"))
kegg.df <- readRDS(file = file.path("data", "RDS", "kegg.df.RDS"))
go.df <- readRDS(file = file.path("data", "RDS", "go.df.RDS"))
hp.df <- readRDS(file = file.path("data", "RDS", "hp.df.RDS"))

reactome.df$greater <- reactome.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia"))
reactome.df$less <- reactome.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia"))

reactome.up <- reactome.df$greater %>%
  filter(pathway %in% names(table(reactome.df$greater$pathway))[which(table(reactome.df$greater$pathway) > 6)])
reactome.down <- reactome.df$less %>%
  filter(pathway %in% names(table(reactome.df$less$pathway))[which(table(reactome.df$less$pathway) > 5)])

## REACTOME
pdf(file.path("data", "output", "SingleSample_cancers_GAGE_REACTOME.pdf"), width = 10, height = 16)
plotPathwayBar(reactome.up, xlabel="Number of UP-regulated pathways", ylabel="REACTOME pathways" , w.label = 55 )
plotPathwayBar(reactome.down, xlabel="Number of DOWN-regulated pathways", ylabel="REACTOME pathways", w.label = 55  )
dev.off()


kegg.up <- kegg.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(kegg.df$greater$pathway))[which(table(kegg.df$greater$pathway) > 6)])
kegg.down <- kegg.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(kegg.df$less$pathway))[which(table(kegg.df$less$pathway) > 4)])


## KEGG
pdf(file.path("data", "output", "SingleSample_cancers_GAGE_KEGG.pdf"), width = 10, height = 10)
plotPathwayBar(kegg.up, xlabel="Number of UP-regulated pathways", ylabel="KEGG pathways" )
plotPathwayBar(kegg.down, xlabel="Number of DOWN-regulated pathways", ylabel="KEGG pathways" )
dev.off()

## GO
## subset to reasonable numbers
go.up <- go.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(go.df$greater$pathway))[which(table(go.df$greater$pathway) > 10)])
go.down <- go.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(go.df$less$pathway))[which(table(go.df$less$pathway) > 10)])

pdf(file.path("data", "output", "SingleSample_cancers_GAGE_GO.pdf"), width = 10, height = 15)
plotPathwayBar(go.up, xlabel="Number of UP-regulated pathways", ylabel="GO pathways", w.label = 55 )
plotPathwayBar(go.down, xlabel="Number of DOWN-regulated pathways", ylabel="GO pathways", w.label = 55 )
dev.off()

## HP
## subset to reasonable numbers
hp.up <- hp.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(hp.df$greater$pathway))[which(table(hp.df$greater$pathway) > 1)])
hp.down <- hp.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(hp.df$less$pathway))[which(table(hp.df$less$pathway) > 6)])

pdf(file.path("data", "output", "SingleSample_cancers_GAGE_HP.pdf"), width = 10, height = 14)
plotPathwayBar(hp.up, xlabel="Number of UP-regulated pathways", ylabel="HP pathways", w.label = 55 )
plotPathwayBar(hp.down, xlabel="Number of DOWN-regulated pathways", ylabel="HP pathways", w.label = 55 )
dev.off()


# COHORT DIFFUSION --------------------------------------------------------

## We take alle the single sample outputs - subset by cohort - 
## merge pathogenic, potential and associated matrices for every sample -->
## extract all genes that overlap between all samples (if there is such a thing)

## diffusion network and kernel
network <- readRDS(file = file.path(FILESOURCE, "RDS","network.rds"), refhook = NULL)
network_kernel <- readRDS(file = file.path(FILESOURCE,"RDS", "network_kernel.rds"), 
                          refhook = NULL)

sample.list <- readRDS(file = file.path("data", "RDS", "samples_annotated.RDS"))


## 1. subset to cohorts of interest
sample.list <- sample.list[grepl("colon", names(sample.list)) | grepl("breast", names(sample.list)) | grepl("leukemia", names(sample.list))]
s.select <- names(sample.list)

## BREAST
tmp.breast <- sample.list[grepl("breast", s.select)]
breast.list <- list()
for( s.idx in names(tmp.breast) ) {
  breast.list[[s.idx]] <- rbind.data.frame(tmp.breast[[s.idx]]$pathogenic, tmp.breast[[s.idx]]$potential, tmp.breast[[s.idx]]$associated)
}
## get intersect
breast.intersect <- Reduce(intersect, lapply(breast.list, function(df) df$Gene.refGene))

## COLON
tmp.colon <- sample.list[grepl("colon", s.select)]
colont.list <- list()
for( s.idx in names(tmp.colon) ) {
  colont.list[[s.idx]] <- rbind.data.frame(tmp.colon[[s.idx]]$pathogenic, tmp.colon[[s.idx]]$potential, tmp.colon[[s.idx]]$associated)
}
## get intersect
colon.intersect <- Reduce(intersect, lapply(colont.list, function(df) df$Gene.refGene))


## LEUKEMIA
tmp.leukemia <- sample.list[grepl("leukemia", s.select)]
leukemia.list <- list()
for( s.idx in names(tmp.leukemia) ) {
  leukemia.list[[s.idx]] <- rbind.data.frame(tmp.leukemia[[s.idx]]$pathogenic, tmp.leukemia[[s.idx]]$potential, tmp.leukemia[[s.idx]]$associated)
}
## get intersect
leukemia.intersect <- Reduce(intersect, lapply(leukemia.list, function(df) df$Gene.refGene))



# node_id - node_score
cohort_label <- "breast_cohort"
# score variants and prepare diffusion input
input_oncogenes <- data.frame("node_id"=breast.intersect, "node_score"=1)

# take remaining Genes and create a binary input_vector for the diffusion
tmp_input_vec <- rep(0,length(igraph::V(network)$name))
names(tmp_input_vec) <- igraph::V(network)$name
idx <- intersect( names(tmp_input_vec), as.character(input_oncogenes$node_id))

# in case the remaining genes in input_oncogenes cannot be mapped to the network 
# remove unmapped genes from input_oncogenes
input_oncogenes <- input_oncogenes[which(input_oncogenes$node_id %in% idx),]

## draw network of input genes for diffusion
input_oncogenes.mapping <- stringDB$string.DB$map(input_oncogenes, "node_id")

# display a STRING network png with the "halo"
pdf(file = file.path(FILESOURCE, "output", paste(eval(cohort_label),"_input_diffnet.pdf", sep="")), 
    height = 16, width = 16, useDingbats = F, onefile=TRUE)

# finally, use stringDB functionality to produce the plot
stringDB$string.DB$plot_network( input_oncogenes.mapping$STRING_id )

dev.off()

# and set the input value for the oncogenes
tmp_input_vec[match(input_oncogenes$node_id, names(tmp_input_vec))] <- 
  input_oncogenes$node_score

# Kernel based diffusion with z-method
df_z <- diffuStats::diffuse_grid(K = network_kernel, scores=tmp_input_vec,
                                 grid_param = expand.grid(method = "z"), 
                                 n.perm = 1)

s.genes.z <- df_z[abs(df_z$node_score) > 0.5,]

## stringDB provides 'map()' to match gene symbols to Ensembl-Protein-IDs
s.genes.z.mapping <- stringDB$string.DB$map(s.genes.z, "node_id")

## Export results for later use.
openxlsx::write.xlsx(s.genes.z.mapping,
                     file = file.path(FILESOURCE,"output",paste(eval(cohort_label),"_diffusion_zscore.xlsx", sep="")),
                     rowNames = TRUE,colNames=TRUE)

gage.input <- s.genes.z$node_score
names(gage.input) <- s.genes.z$node_id

## REACTOME
gage.REACTOME <- gage(gage.input,gsets = Hs.c2.REACTOME,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.REACTOME <- sigGeneSet(gage.REACTOME, qpval=c("q.val"), cutoff=0.05)
## KEGG
gage.KEGG <- gage(gage.input,gsets = Hs.c2.KEGG,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.KEGG <- sigGeneSet(gage.KEGG, qpval=c("q.val"), cutoff=0.05)
## GO
gage.GO <- gage(gage.input,gsets = Hs.c5.GO,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.GO <- sigGeneSet(gage.GO, qpval=c("q.val"), cutoff=0.05)
## HPO Human Phenotype Ontology
gage.HP <- gage(gage.input,gsets = Hs.c5.HP,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.HP <- sigGeneSet(gage.HP, qpval=c("q.val"), cutoff=0.05)


saveRDS(gage.REACTOME, file = file.path("data", "RDS", paste(eval(cohort_label),"_reactome.df.RDS", sep="")), compress = TRUE)
saveRDS(gage.KEGG, file = file.path("data", "RDS", paste(eval(cohort_label),"_kegg.df.RDS", sep="")), compress = TRUE)
saveRDS(gage.GO, file = file.path("data", "RDS", paste(eval(cohort_label),"_go.df.RDS", sep="")), compress = TRUE)
saveRDS(gage.HP, file = file.path("data", "RDS", paste(eval(cohort_label),"_hp.df.RDS", sep="")), compress = TRUE)

openxlsx::write.xlsx(gage.REACTOME, file = file.path("data", "output", paste(eval(cohort_label),"_GAGE_REACTOME.xlsx", sep="")))
openxlsx::write.xlsx(gage.KEGG, file = file.path("data", "output", paste(eval(cohort_label),"_GAGE_KEGG.xlsx", sep="")))
openxlsx::write.xlsx(gage.GO, file = file.path("data", "output", paste(eval(cohort_label),"_GAGE_GO.xlsx", sep="")))
openxlsx::write.xlsx(gage.HP, file = file.path("data", "output", paste(eval(cohort_label),"_GAGE_HP.xlsx", sep="")))






# SINGLE SAMPLE + FILTERING FOR CONTROL VARIANTS --------------------------

## For every single sample filter against the variants present in the control cohort.
## Then repeat diffusion for single samples and for cohort gene overlap and compare
## to the original output.

## diffusion network and kernel
network <- readRDS(file = file.path(FILESOURCE, "RDS","network.rds"), refhook = NULL)
network_kernel <- readRDS(file = file.path(FILESOURCE,"RDS", "network_kernel.rds"), 
                          refhook = NULL)

sample_filtered.list <- list()
sample.list <- readRDS(file = file.path("data", "RDS", "samples_annotated.RDS"))
## I. get sample names

## II. subset to cohorts of interest
sample.list <- sample.list[grepl("colon", names(sample.list)) | grepl("breast", names(sample.list)) | grepl("leukemia", names(sample.list))]

s.select <- names(sample.list)

## BREAST ##
breast.ctrl <- readRDS(file = file.path("data", "RDS", "breast_controls_variants.RDS"))
breast.ctrl <- rbind.data.frame(breast.ctrl$pathogenic, breast.ctrl$potential, breast.ctrl$associated, breast.ctrl$benign,breast.ctrl$likely_benign, breast.ctrl$unscored)

## 1. filter against the cohort variants
tmp.breast <- sample.list[grepl("breast", s.select)]
for( s.idx in names(tmp.breast) ) {
  tmp.breast[[s.idx]] <- removeCtrlBackground( tmp.breast[[s.idx]], breast.ctrl )
}

sample_filtered.list[["breast"]] <- tmp.breast

## 2. get the cohort intersect
breast.list <- list()
for( s.idx in names(tmp.breast) ) {
  breast.list[[s.idx]] <- rbind.data.frame(tmp.breast[[s.idx]]$pathogenic, tmp.breast[[s.idx]]$potential, tmp.breast[[s.idx]]$associated)
}
## get intersect - OLD: "CTBP2"    "TEKT4"    "ANKRD36C" "ANKRD36B" "PABPC1"   "NF1"      "URI1"     "TTN"      "CPS1"     "ERBB4"    "COL18A1"  "FGFR3"    "MET"      "WRN" 
breast.intersect <- Reduce(intersect, lapply(breast.list, function(df) df$Gene.refGene))


## COLON ##
colon.ctrl <- readRDS(file = file.path("data", "RDS", "colon_controls_variants.RDS"))
colon.ctrl <- rbind.data.frame(colon.ctrl$pathogenic, colon.ctrl$potential, colon.ctrl$associated, colon.ctrl$benign,colon.ctrl$likely_benign, colon.ctrl$unscored)

## 1. filter against the cohort variants
tmp.colon <- sample.list[grepl("colon", s.select)]
for( s.idx in names(tmp.colon) ) {
  tmp.colon[[s.idx]] <- removeCtrlBackground( tmp.colon[[s.idx]], colon.ctrl )
}

sample_filtered.list[["colon"]] <- tmp.colon

## 2. get the cohort intersect
colon.list <- list()
for( s.idx in names(tmp.colon) ) {
  colon.list[[s.idx]] <- rbind.data.frame(tmp.colon[[s.idx]]$pathogenic, tmp.colon[[s.idx]]$potential, tmp.colon[[s.idx]]$associated)
}
## get intersect - OLD: "CTBP2"    "ANKRD36C" "ANKRD36B" "POTED"    "COL11A1"  "ATM"      "BRCA1"    "MSH6"     "COL18A1"  "FGFR3"    "TLR4"
colon.intersect <- Reduce(intersect, lapply(colon.list, function(df) df$Gene.refGene))


## LEUKEMIA ##
leukemia.ctrl <- readRDS(file = file.path("data", "RDS", "leukemia_controls_variants.RDS"))
leukemia.ctrl <- rbind.data.frame(leukemia.ctrl$pathogenic, leukemia.ctrl$potential, leukemia.ctrl$associated, leukemia.ctrl$benign,leukemia.ctrl$likely_benign, leukemia.ctrl$unscored)

## 1. filter against the cohort variants
tmp.leukemia <- sample.list[grepl("leukemia", s.select)]
for( s.idx in names(tmp.leukemia) ) {
  tmp.leukemia[[s.idx]] <- removeCtrlBackground( tmp.leukemia[[s.idx]], leukemia.ctrl )
}

sample_filtered.list[["leukemia"]] <- tmp.leukemia

## 2. get the cohort intersect
leukemia.list <- list()
for( s.idx in names(tmp.leukemia) ) {
  leukemia.list[[s.idx]] <- rbind.data.frame(tmp.leukemia[[s.idx]]$pathogenic, tmp.leukemia[[s.idx]]$potential, tmp.leukemia[[s.idx]]$associated)
}
## get intersect - OLD: "CTBP2"    "ANKRD36C" "TTN"      "NBPF1"    "AKNAD1"   "ASTN1"    "SMC3"     "ATM"      "TCTN2"    "SNX29"    "RBBP6"    "NF1"      "MUC16"    "DNMT1"    "SCN1A"    "ZNFX1"   
## "COL18A1"  "MUC4"     "FGFR3"    "CWH43"    "TPMT"     "EYS"      "WRN"
leukemia.intersect <- Reduce(intersect, lapply(leukemia.list, function(df) df$Gene.refGene))


## Save singe sample list 
saveRDS(sample_filtered.list, file = file.path("data", "RDS", "samples_annotated_filtered.RDS"), compress = TRUE)






#### AAAAND DO DIFFUSION AND STUFF

# node_id - node_score
cohort_label <- "breast_cohort_filtered"
# score variants and prepare diffusion input
input_oncogenes <- data.frame("node_id"=breast.intersect, "node_score"=1)

# take remaining Genes and create a binary input_vector for the diffusion
tmp_input_vec <- rep(0,length(igraph::V(network)$name))
names(tmp_input_vec) <- igraph::V(network)$name
idx <- intersect( names(tmp_input_vec), as.character(input_oncogenes$node_id))

# in case the remaining genes in input_oncogenes cannot be mapped to the network 
# remove unmapped genes from input_oncogenes
input_oncogenes <- input_oncogenes[which(input_oncogenes$node_id %in% idx),]

## draw network of input genes for diffusion
input_oncogenes.mapping <- stringDB$string.DB$map(input_oncogenes, "node_id")

# display a STRING network png with the "halo"
pdf(file = file.path(FILESOURCE, "output", paste(eval(cohort_label),"_input_diffnet.pdf", sep="")), 
    height = 16, width = 16, useDingbats = F, onefile=TRUE)

# finally, use stringDB functionality to produce the plot
stringDB$string.DB$plot_network( input_oncogenes.mapping$STRING_id )

dev.off()

# and set the input value for the oncogenes
tmp_input_vec[match(input_oncogenes$node_id, names(tmp_input_vec))] <- 
  input_oncogenes$node_score

# Kernel based diffusion with z-method
df_z <- diffuStats::diffuse_grid(K = network_kernel, scores=tmp_input_vec,
                                 grid_param = expand.grid(method = "z"), 
                                 n.perm = 1)

s.genes.z <- df_z[abs(df_z$node_score) > 0.5,]

## stringDB provides 'map()' to match gene symbols to Ensembl-Protein-IDs
s.genes.z.mapping <- stringDB$string.DB$map(s.genes.z, "node_id")


## Export results for later use.
openxlsx::write.xlsx(s.genes.z.mapping,
                     file = file.path(FILESOURCE,"output",paste(eval(cohort_label),"_diffusion_zscore.xlsx", sep="")),
                     rowNames = TRUE,colNames=TRUE)

gage.input <- s.genes.z$node_score
names(gage.input) <- s.genes.z$node_id

## REACTOME
gage.REACTOME <- gage(gage.input,gsets = Hs.c2.REACTOME,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.REACTOME <- sigGeneSet(gage.REACTOME, qpval=c("q.val"), cutoff=0.05)
## KEGG
gage.KEGG <- gage(gage.input,gsets = Hs.c2.KEGG,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.KEGG <- sigGeneSet(gage.KEGG, qpval=c("q.val"), cutoff=0.05)
## GO
gage.GO <- gage(gage.input,gsets = Hs.c5.GO,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.GO <- sigGeneSet(gage.GO, qpval=c("q.val"), cutoff=0.05)
## HPO Human Phenotype Ontology
gage.HP <- gage(gage.input,gsets = Hs.c5.HP,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.HP <- sigGeneSet(gage.HP, qpval=c("q.val"), cutoff=0.05)


saveRDS(gage.REACTOME, file = file.path("data", "RDS", paste(eval(cohort_label),"_reactome.df.RDS", sep="")), compress = TRUE)
saveRDS(gage.KEGG, file = file.path("data", "RDS", paste(eval(cohort_label),"_kegg.df.RDS", sep="")), compress = TRUE)
saveRDS(gage.GO, file = file.path("data", "RDS", paste(eval(cohort_label),"_go.df.RDS", sep="")), compress = TRUE)
saveRDS(gage.HP, file = file.path("data", "RDS", paste(eval(cohort_label),"_hp.df.RDS", sep="")), compress = TRUE)


openxlsx::write.xlsx(gage.REACTOME, file = file.path("data", "output", paste(eval(cohort_label),"_GAGE_REACTOME.xlsx", sep="")))
openxlsx::write.xlsx(gage.KEGG, file = file.path("data", "output", paste(eval(cohort_label),"_GAGE_KEGG.xlsx", sep="")))
openxlsx::write.xlsx(gage.GO, file = file.path("data", "output", paste(eval(cohort_label),"_GAGE_GO.xlsx", sep="")))
openxlsx::write.xlsx(gage.HP, file = file.path("data", "output", paste(eval(cohort_label),"_GAGE_HP.xlsx", sep="")))




# BACKGROUND FILTERED SINGLE SAMPLE DIFFUSION -----------------------------

## diffusion network and kernel
network <- readRDS(file = file.path(FILESOURCE, "RDS","network.rds"), refhook = NULL)
network_kernel <- readRDS(file = file.path(FILESOURCE,"RDS", "network_kernel.rds"), 
                          refhook = NULL)

## load the filtered sample list that only contains breast, colon, and leukemia single sample dfs
sample.list <- readRDS(file = file.path("data", "RDS", "samples_annotated_filtered.RDS"))



### RUN DIffusion
reactome.df <- createMergeDF()
kegg.df <- createMergeDF()
go.df <- createMergeDF()
hp.df <- createMergeDF()

for( leCohort in names(sample.list) ) {
  
  cohort.list <- sample.list[[eval(leCohort)]]
  
  for( leSample in names(cohort.list) ) {
    
    m = unlist(strsplit(leSample, '_'))
    df.list <- cohort.list[[eval(leSample)]]
    
    # score variants and prepare diffusion input
    input_oncogenes <- prepDiffSet(df.list)
    
    # take remaining Genes and create a binary input_vector for the diffusion
    tmp_input_vec <- rep(0,length(igraph::V(network)$name))
    names(tmp_input_vec) <- igraph::V(network)$name
    idx <- intersect( names(tmp_input_vec), as.character(input_oncogenes$node_id))
    
    # in case the remaining genes in input_oncogenes cannot be mapped to the network 
    # remove unmapped genes from input_oncogenes
    input_oncogenes <- input_oncogenes[which(input_oncogenes$node_id %in% idx),]
    
    ## draw network of input genes for diffusion
    input_oncogenes.mapping <- stringDB$string.DB$map(input_oncogenes, "node_id")
    
    # display a STRING network png with the "halo"
    pdf(file = file.path(FILESOURCE, "output", paste(m[1],"_",m[2],"_",m[3],"_backgroundRemoved_input_diffnet.pdf", sep="")), 
        height = 16, width = 16, useDingbats = F, onefile=TRUE)
    
    # finally, use stringDB functionality to produce the plot
    stringDB$string.DB$plot_network( input_oncogenes.mapping$STRING_id )
    
    dev.off()
    
    # and set the input value for the oncogenes
    tmp_input_vec[match(input_oncogenes$node_id, names(tmp_input_vec))] <- 
      input_oncogenes$node_score
    
    # Kernel based diffusion with z-method
    df_z <- diffuStats::diffuse_grid(K = network_kernel, scores=tmp_input_vec,
                                     grid_param = expand.grid(method = "z"), 
                                     n.perm = 1)
    
    s.genes.z <- df_z[abs(df_z$node_score) > 0.5,]
    
    ## stringDB provides 'map()' to match gene symbols to Ensembl-Protein-IDs
    s.genes.z.mapping <- stringDB$string.DB$map(s.genes.z, "node_id")
    
    ## Export results for later use.
    openxlsx::write.xlsx(s.genes.z.mapping,
                         file = file.path(FILESOURCE,"output",paste(m[1],"_",m[2],"_",m[3],"_backgroundRemoved_diffusion_zscore.xlsx", sep="")),
                         rowNames = TRUE,colNames=TRUE)
    
    gage.input <- s.genes.z$node_score
    names(gage.input) <- s.genes.z$node_id
    
    ## REACTOME
    gage.REACTOME <- gage(gage.input,gsets = Hs.c2.REACTOME,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
    gage.REACTOME <- sigGeneSet(gage.REACTOME, qpval=c("q.val"), cutoff=0.05)
    ## KEGG
    gage.KEGG <- gage(gage.input,gsets = Hs.c2.KEGG,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
    gage.KEGG <- sigGeneSet(gage.KEGG, qpval=c("q.val"), cutoff=0.05)
    ## GO
    gage.GO <- gage(gage.input,gsets = Hs.c5.GO,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
    gage.GO <- sigGeneSet(gage.GO, qpval=c("q.val"), cutoff=0.05)
    ## HPO Human Phenotype Ontology
    gage.HP <- gage(gage.input,gsets = Hs.c5.HP,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
    gage.HP <- sigGeneSet(gage.HP, qpval=c("q.val"), cutoff=0.05)
    
    reactome.df <- mergeGage( reactome.df, gage.REACTOME, s.idx=m[3], c.idx=m[1] )
    kegg.df <- mergeGage( kegg.df, gage.KEGG, s.idx=m[3], c.idx=m[1] )
    go.df <- mergeGage( go.df, gage.GO, s.idx=m[3], c.idx=m[1] )
    hp.df <- mergeGage( hp.df, gage.HP, s.idx=m[3], c.idx=m[1] )
    
    
  }
  
}



reactome.df$greater$pathway <- gsub("REACTOME_", "",reactome.df$greater$pathway)
reactome.df$less$pathway <- gsub("REACTOME_", "",reactome.df$less$pathway)
kegg.df$greater$pathway <- gsub("KEGG_", "",kegg.df$greater$pathway)
kegg.df$less$pathway <- gsub("KEGG_", "",kegg.df$less$pathway)
go.df$greater$pathway <- gsub("GOBP_", "",go.df$greater$pathway)
go.df$less$pathway <- gsub("GOBP_", "",go.df$less$pathway)
hp.df$greater$pathway <- gsub("HP_", "",hp.df$greater$pathway)
hp.df$less$pathway <- gsub("HP_", "",hp.df$less$pathway)


saveRDS(reactome.df, file = file.path("data", "RDS", "backgroundRemoved_reactome.df.RDS"), compress = TRUE)
saveRDS(kegg.df, file = file.path("data", "RDS", "backgroundRemoved_kegg.df.RDS"), compress = TRUE)
saveRDS(go.df, file = file.path("data", "RDS", "backgroundRemoved_go.df.RDS"), compress = TRUE)
saveRDS(hp.df, file = file.path("data", "RDS", "backgroundRemoved_hp.df.RDS"), compress = TRUE)

## write to file
openxlsx::write.xlsx(reactome.df, file = file.path("data", "output", "SingleSample_cancers_filtered_GAGE_REACTOME.xlsx"))
openxlsx::write.xlsx(kegg.df, file = file.path("data", "output", "SingleSample_cancers_filtered_GAGE_KEGG.xlsx"))
openxlsx::write.xlsx(go.df, file = file.path("data", "output", "SingleSample_cancers_filtered_GAGE_GO.xlsx"))
openxlsx::write.xlsx(hp.df, file = file.path("data", "output", "SingleSample_cancers_filtered_GAGE_HP.xlsx"))


## REACTOME
reactome.df$greater <- reactome.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

reactome.df$less <- reactome.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(reactome.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_filtered_GAGE_REACTOME.xlsx"))

## KEGG
kegg.df$greater <- kegg.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

kegg.df$less <- kegg.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(kegg.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_filtered_GAGE_KEGG.xlsx"))

## GO
go.df$greater <- go.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

go.df$less <- go.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(go.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_filtered_GAGE_GO.xlsx"))

## HP
hp.df$greater <- hp.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

hp.df$less <- hp.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  group_by(pathway, cohort) %>%
  count() %>%
  pivot_wider(names_from = cohort, values_from = n, values_fill = 0) %>%
  arrange(desc(breast), desc(colon), desc(leukemia))

openxlsx::write.xlsx(hp.df, file = file.path("data", "output", "COUNTS_SingleSample_cancers_filtered_GAGE_HP.xlsx"))



# PLOT RESULTS - BACKGROUND REMOVED ---------------------------------------

reactome.df <- readRDS(file = file.path("data", "RDS", "backgroundRemoved_reactome.df.RDS"))
kegg.df <- readRDS(file = file.path("data", "RDS", "backgroundRemoved_kegg.df.RDS"))
go.df <- readRDS(file = file.path("data", "RDS", "backgroundRemoved_go.df.RDS"))
hp.df <- readRDS(file = file.path("data", "RDS", "backgroundRemoved_hp.df.RDS"))

reactome.df$greater <- reactome.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia"))
reactome.df$less <- reactome.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia"))

reactome.up <- reactome.df$greater %>%
  filter(pathway %in% names(table(reactome.df$greater$pathway))[which(table(reactome.df$greater$pathway) > 6)])
reactome.down <- reactome.df$less %>%
  filter(pathway %in% names(table(reactome.df$less$pathway))[which(table(reactome.df$less$pathway) > 5)])

## REACTOME
pdf(file.path("data", "output", "SingleSample_cancers_filtered_GAGE_REACTOME.pdf"), width = 10, height = 16)
plotPathwayBar(reactome.up, xlabel="Number of UP-regulated pathways", ylabel="REACTOME pathways" , w.label = 55 )
plotPathwayBar(reactome.down, xlabel="Number of DOWN-regulated pathways", ylabel="REACTOME pathways", w.label = 55  )
dev.off()


kegg.up <- kegg.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(kegg.df$greater$pathway))[which(table(kegg.df$greater$pathway) > 6)])
kegg.down <- kegg.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(kegg.df$less$pathway))[which(table(kegg.df$less$pathway) > 4)])


## KEGG
pdf(file.path("data", "output", "SingleSample_cancers_filtered_GAGE_KEGG.pdf"), width = 10, height = 10)
plotPathwayBar(kegg.up, xlabel="Number of UP-regulated pathways", ylabel="KEGG pathways" )
plotPathwayBar(kegg.down, xlabel="Number of DOWN-regulated pathways", ylabel="KEGG pathways" )
dev.off()

## GO
## subset to reasonable numbers
go.up <- go.df$greater %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(go.df$greater$pathway))[which(table(go.df$greater$pathway) > 10)])
go.down <- go.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(go.df$less$pathway))[which(table(go.df$less$pathway) > 10)])

pdf(file.path("data", "output", "SingleSample_cancers_filtered_GAGE_GO.pdf"), width = 10, height = 15)
plotPathwayBar(go.up, xlabel="Number of UP-regulated pathways", ylabel="GO pathways", w.label = 55 )
plotPathwayBar(go.down, xlabel="Number of DOWN-regulated pathways", ylabel="GO pathways", w.label = 55 )
dev.off()

## HP
## subset to reasonable numbers
hp.up <- hp.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(hp.df$greater$pathway))[which(table(hp.df$greater$pathway) > 1)])
hp.down <- hp.df$less %>%
  dplyr::filter(cohort %in% c("breast","colon","leukemia")) %>%
  filter(pathway %in% names(table(hp.df$less$pathway))[which(table(hp.df$less$pathway) > 6)])

pdf(file.path("data", "output", "SingleSample_cancers_filtered_GAGE_HP.pdf"), width = 10, height = 14)
plotPathwayBar(hp.up, xlabel="Number of UP-regulated pathways", ylabel="HP pathways", w.label = 55 )
plotPathwayBar(hp.down, xlabel="Number of DOWN-regulated pathways", ylabel="HP pathways", w.label = 55 )
dev.off()

































# OLD SHIT ----------------------------------------------------------------



tmp.breast <- sample.list[grepl("breast", s.select)]
breast.list <- list()
for( s.idx in names(tmp.breast) ) {
  breast.list[[s.idx]] <- rbind.data.frame(tmp.breast[[s.idx]]$pathogenic, tmp.breast[[s.idx]]$potential, tmp.breast[[s.idx]]$associated)
}
## get intersect
breast.intersect <- Reduce(intersect, lapply(breast.list, function(df) df$Gene.refGene))

removeCtrlBackground <- function( input.list, ctrl.df ) {
  
  for( type.idx in names(input.list) ) {
    ## get intersect
    tmp <- intersect(ctrl.df$varID, input.list[[eval(type.idx)]]$varID)
    ## keep everything NOT in the intersect
    input.list[[eval(type.idx)]] <- input.list[[eval(type.idx)]] %>%
      dplyr::filter(!varID %in% tmp)
  }
  
  return( input.list )
  
}


tmpCutoff <- quantile(abs(df_z$node_score), 0.9975)
s.genes.z <- df_z[abs(df_z$node_score) > 0.5,]



dummy <- diffuseWrapper(network, K_rl=network_kernel, input=input_oncogenes.mapping, method="z")

#'@param input, output of 'diffuseWrapper' containing input-GOI (x) and diffusion results (df_z)
#'@param network, background-network the kernel is based on
#'@param outdir, path to plot to
#'@param score, String denoting the initial diffusion value
#'@param cutoff, threshold for the the diffusion scores to include in graph (quantile [0,1])
plotDiffNets(input=dummy, network=network, outdir="data", score="leTest", cutoff=0.9975)

stringDB$string.DB$plot_network( s.genes.z.mapping$STRING_id )



/Volumes/WorkDrive/VDIStuff/MAF


df <- data.table::fread(file = file.path("/Volumes/WorkDrive/VDIStuff/MAF", "maf_Illumina_norm_chr19.vcf"))

head(df)

dummy <- df[df$POS > 15879615,]

head(dummy)






###### 26.12.2024

#### AAAAND DO DIFFUSION AND STUFF

# node_id - node_score
cohort_label <- "breast_cohort_filtered"
# score variants and prepare diffusion input
input_oncogenes <- data.frame("node_id"=breast.intersect, "node_score"=1)

# take remaining Genes and create a binary input_vector for the diffusion
tmp_input_vec <- rep(0,length(igraph::V(network)$name))
names(tmp_input_vec) <- igraph::V(network)$name
idx <- intersect( names(tmp_input_vec), as.character(input_oncogenes$node_id))

# in case the remaining genes in input_oncogenes cannot be mapped to the network 
# remove unmapped genes from input_oncogenes
input_oncogenes <- input_oncogenes[which(input_oncogenes$node_id %in% idx),]

## draw network of input genes for diffusion
input_oncogenes.mapping <- stringDB$string.DB$map(input_oncogenes, "node_id")

# display a STRING network png with the "halo"
pdf(file = file.path(FILESOURCE, "output", paste(eval(cohort_label),"_input_diffnet.pdf", sep="")), 
    height = 16, width = 16, useDingbats = F, onefile=TRUE)

# finally, use stringDB functionality to produce the plot
stringDB$string.DB$plot_network( input_oncogenes.mapping$STRING_id )

dev.off()

# and set the input value for the oncogenes
tmp_input_vec[match(input_oncogenes$node_id, names(tmp_input_vec))] <- 
  input_oncogenes$node_score

# Kernel based diffusion with z-method
df_z <- diffuStats::diffuse_grid(K = network_kernel, scores=tmp_input_vec,
                                 grid_param = expand.grid(method = "z"), 
                                 n.perm = 1)

s.genes.z <- df_z[abs(df_z$node_score) > 0.5,]

## stringDB provides 'map()' to match gene symbols to Ensembl-Protein-IDs
s.genes.z.mapping <- stringDB$string.DB$map(s.genes.z, "node_id")


## Export results for later use.
openxlsx::write.xlsx(s.genes.z.mapping,
                     file = file.path(FILESOURCE,"output",paste(eval(cohort_label),"_diffusion_zscore.xlsx", sep="")),
                     rowNames = TRUE,colNames=TRUE)

gage.input <- s.genes.z$node_score
names(gage.input) <- s.genes.z$node_id

## REACTOME
gage.REACTOME <- gage(gage.input,gsets = Hs.c2.REACTOME,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.REACTOME <- sigGeneSet(gage.REACTOME, qpval=c("q.val"), cutoff=0.05)
## KEGG
gage.KEGG <- gage(gage.input,gsets = Hs.c2.KEGG,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.KEGG <- sigGeneSet(gage.KEGG, qpval=c("q.val"), cutoff=0.05)
## GO
gage.GO <- gage(gage.input,gsets = Hs.c5.GO,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.GO <- sigGeneSet(gage.GO, qpval=c("q.val"), cutoff=0.05)
## HPO Human Phenotype Ontology
gage.HP <- gage(gage.input,gsets = Hs.c5.HP,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.HP <- sigGeneSet(gage.HP, qpval=c("q.val"), cutoff=0.05)


saveRDS(gage.REACTOME, file = file.path("data", "RDS", paste(eval(cohort_label),"_reactome.df.RDS", sep="")), compress = TRUE)
saveRDS(gage.KEGG, file = file.path("data", "RDS", paste(eval(cohort_label),"_kegg.df.RDS", sep="")), compress = TRUE)
saveRDS(gage.GO, file = file.path("data", "RDS", paste(eval(cohort_label),"_go.df.RDS", sep="")), compress = TRUE)
saveRDS(gage.HP, file = file.path("data", "RDS", paste(eval(cohort_label),"_hp.df.RDS", sep="")), compress = TRUE)


openxlsx::write.xlsx(gage.REACTOME, file = file.path("data", "output", paste(eval(cohort_label),"_GAGE_REACTOME.xlsx", sep="")))
openxlsx::write.xlsx(gage.KEGG, file = file.path("data", "output", paste(eval(cohort_label),"_GAGE_KEGG.xlsx", sep="")))
openxlsx::write.xlsx(gage.GO, file = file.path("data", "output", paste(eval(cohort_label),"_GAGE_GO.xlsx", sep="")))
openxlsx::write.xlsx(gage.HP, file = file.path("data", "output", paste(eval(cohort_label),"_GAGE_HP.xlsx", sep="")))



