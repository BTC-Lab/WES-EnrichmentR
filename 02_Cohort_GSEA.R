# INIT --------------------------------------------------------------------

# BiocManager::install(c("DESeq2","airway","gage","gageData","org.Hs.eg.db",
#  "pathview","msigdb","ExperimentHub","AnnotationHub", "htmltools"))
# BiocManager::install("GSEABase")

library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
# BiocManager::install("maftools")
library(readr)
library(tidyverse)
data.dir <- "/Users/rmolbrich/Workdir/Projects/VRI/data"
library(pals)
library("MatchIt")
# Load required library
library(VIM)
library(plyr)

new_theme <- theme_bw() + 
  # adjustments for the legend
  theme(legend.position="bottom",
        legend.text = element_text(color = "black", size=12),
        legend.key = element_rect(size=12))


library(DESeq2)
library(gage)
library(plyr)
library(airway)
library(dplyr)

library(igraph)
library(readr)
library(org.Hs.eg.db)
library(pathview)
library(GSEABase)
library(msigdb)
library(ExperimentHub)

# # les bases de data
# Hs.hm <- MSigDB[["HALLMARK"]]
# Hs.c2 <- MSigDB[["C2_CURATED"]]
# Hs.c3 <- MSigDB[["C3_MOTIF"]]
# Hs.c5 <- MSigDB[["C5_GENE_ONTOLOGY"]]
# Hs.c7 <- MSigDB[["C7_IMMUNOLOGIC_SIGNATURES"]]
# # extract what we need
# Hs.c2.KEGG <- Hs.c2[grep("KEGG_", names(Hs.c2))]
# names(Hs.c2.KEGG) <- gsub("^.*?_", "", names(Hs.c2.KEGG))
# Hs.c2.REACTOME <- Hs.c2[grep("REACTOME_", names(Hs.c2))]
# names(Hs.c2.REACTOME) <- gsub("^.*?_", "", names(Hs.c2.REACTOME))
# Hs.c7.IMMUDB <- Hs.c7[names(Hs.c7)]
# names(Hs.c7.IMMUDB) <- gsub("^.*?_", "", names(Hs.c7.IMMUDB))

getPlinkPCAMetrics <- function(val.obj, df.obj, t.var) {
  
  out.df <- data.frame(
    pc.id = paste("PC", seq_len(length(val.obj$V1)), sep = ""),
    var.explained = round(val.obj$V1 / t.var * 100, 2)) %>%
    dplyr::mutate(vr.label = paste0(pc.id, ": ", var.explained, "%")) %>%
    dplyr::mutate(axis.min = (apply(df.obj[,2:21], 2, function(col)
      min(col)))) %>%
    dplyr::mutate(axis.max = (apply(df.obj[,2:21], 2, function(col)
      max(col))))
  
  return(out.df)
}

predictBMI <- function(df) {
  
  # Filter rows with missing BMI values
  missing_bmi <- df %>%
    filter(is.na(BMI.imp))
  # Filter rows with non-missing BMI values
  non_missing_bmi <- df %>%
    filter(!is.na(BMI.imp))
  # Perform linear regression to predict BMI from Weight
  lm_model <- lm(BMI.imp ~ Weight.imp, data = non_missing_bmi)
  # Impute missing BMI values based on the linear regression model
  missing_bmi$BMI.imp <- predict(lm_model, newdata = missing_bmi)
  # Combine the imputed data with the non-missing data
  imputed_df <- bind_rows(missing_bmi, non_missing_bmi)
  # Assign the imputed values back to the original data frame
  df$BMI.imp <- imputed_df$BMI.imp
  
  return(df)
}

p.axes <- c("PC1", "PC2")
cov.col <- "Cancer"

require(pals)
pal.bands(coolwarm, parula, ocean.haline, brewer.blues, cubicl, kovesi.rainbow, ocean.phase, brewer.paired(12), stepped, brewer.seqseq2,
          main="Colormap suggestions")

cols <- pals::tableau20(20)
scales::show_col(cols)

this_theme <- function() {
  return(theme_bw() + 
           theme(#axis.text.x=element_blank(),
             axis.title.y=element_blank(),
             #axis.ticks=element_blank(),
             axis.line=element_blank(),
             panel.border = element_blank(),
             legend.title=element_text(size = 12, face = "bold"),
             legend.text=element_text(size = 10, face = "bold"),
             axis.title = element_text(face = "plain"),
             plot.margin=unit(c(0.05,0.05,0.05,0.05), "cm"),
             legend.position='bottom'))
}

## For the matching matrix of matchIT generate a new column of pairs in meta.df
matchItCol <- function(meta.df, mmat, t.prefix) {
  ## get number of assigned controls
  nm <- dim(mmat)[2]
  ## build col.name
  tmp.name <- paste('m',t.prefix,nm, sep='_')
  ## set-up new column
  meta.df[[tmp.name]] <- ifelse(meta.df$Cancer == "yes", meta.df$S.merge, NA)
  ## for every case aka row
  for(r.idx in 1:dim(mmat)[1]) {
    ## get case-ID
    caseID <- rownames(mmat)[r.idx]
    # for every control aka column
    for(c.idx in 1:dim(mmat)[2]) {
      ## get index of matching sample
      s.idx <- which(rownames(meta.df) %in% mmat[r.idx,c.idx])
      ## assign caseID
      meta.df[[tmp.name]][s.idx] <- caseID
    }
  }
  return(meta.df)
}


colfunc <- colorRampPalette(c("#2CA02C", "#D62728"))
colfunc(10)

p.cols <- c(cols[1],cols[3],cols[5],cols[9],cols[11],cols[13],cols[17],cols[19],cols[6])


# FUNCTION ----------------------------------------------------------------

prepGeneSet <- function( df.list, merge.scheme = 1 ) {
  
  ## merge pathogenic, potential and associated into a single dataframe
  if( merge.scheme == 1 ) {
    df <- rbind.data.frame(df.list$pathogenic, df.list$potential, df.list$unscored)
  } else {
    df <- rbind.data.frame(df.list$pathogenic, df.list$potential, df.list$associated)
  }
  
  df$Polyphen2_HDIV_pred <- factor(df$Polyphen2_HDIV_pred, levels = c("D","P","B"))
  df$Polyphen2_HVAR_pred <- factor(df$Polyphen2_HVAR_pred, levels = c("D","P","B"))
  df$SIFT4G_pred <- factor(df$SIFT4G_pred, levels = c("D","T"))
  
  df$Func.refGene <- factor(df$Func.refGene,
                            levels = c("exonic","exonic;splicing","UTR5","UTR3",
                                       "ncRNA_exonic","splicing","ncRNA_splicing",
                                       "ncRNA_exonic;splicing","intronic",
                                       "ncRNA_intronic","downstream","upstream",
                                       "upstream;downstream","intergenic"))
  df$ExonicFunc.refGene <- factor(df$ExonicFunc.refGene, 
                                  levels = c("stopgain","startloss",
                                             "nonsynonymous SNV","frameshift insertion",
                                             "frameshift deletion","synonymous SNV",
                                             "nonframeshift deletion","nonframeshift insertion",
                                             "unknown" ))
  ## order by function and by deleteriousness scores
  df <- df[order(-df$CADD_phred, df$SIFT4G_score,
                 df$Func.refGene,df$ExonicFunc.refGene,
                 df$Polyphen2_HDIV_pred, df$Polyphen2_HVAR_pred, df$SIFT4G_pred), ]
  
  summarized_df <- df %>%
    dplyr::mutate(gene.importance = (n() - row_number() + 1) / log1p(row_number())) %>%
    #dplyr::mutate(gene.importance = log1p(n() - row_number() + 1)) %>%
    dplyr::group_by(Gene.refGene) %>%
    dplyr::summarize(gsea.rank.score = sum(gene.importance)) %>%
    dplyr::arrange(desc(gsea.rank.score))
  
  vec.genes <- summarized_df$gsea.rank.score
  names(vec.genes) <- summarized_df$Gene.refGene
  
  return(vec.genes)
}

prepPlot <- function( df.gage, n.path=10 ) {
  tmp.g <- as.data.frame(df.gage$greater)
  tmp.l <- as.data.frame(df.gage$less)
  
  if(dim(tmp.g)[1] > 0) tmp.g$class <- "up"
  if(dim(tmp.l)[1] > 0) tmp.l$class <- "down"
  
  tmp.g <- head(tmp.g, n=n.path)
  tmp.l <- head(tmp.l, n=n.path)
  
  tmp <- rbind.data.frame(tmp.g, tmp.l)
  tmp$names <- rownames(tmp)
  tmp$names <- gsub("REACTOME_", "",tmp$names)
  tmp$names <- gsub("KEGG_", "",tmp$names)
  tmp$names <- gsub("CP_", "",tmp$names)
  tmp$names <- gsub("HP_", "",tmp$names)
  tmp$size <- -log2(tmp$p.val)
  
  return(tmp)
}

plotDot <- function( res.df , title = "", n.path = 10 ) {
  
  ggplot(prepPlot(res.df, n.path), aes(x = p.val, y = names, color = class)) +
    geom_point(aes(size=size)) +  # Add dots
    #geom_text(aes(label = p.val), hjust = -0.2) +  # Add p-values as labels
    labs(x = "P-Value", y = "") +  # Label axes
    scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "_" , " "),
                                                   width = 25)) +
    ggtitle(eval(title)) +
    #facet_grid(cols = vars(class)) +
    theme_minimal()
  
}


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



# GSEA - BREAST -----------------------------------------------------------
data.dir = "data"
## LOAD gene-list
df.list <- readRDS(file = file.path(data.dir, "RDS", "breast_cancer_variants.RDS"))

vec.genes <- prepGeneSet(df.list, merge.scheme = 1)

## REACTOME
gage.REACTOME <- gage(vec.genes,gsets = Hs.c2.REACTOME,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.REACTOME <- sigGeneSet(gage.REACTOME, qpval=c("p.val"), cutoff=0.05)
## KEGG
gage.KEGG <- gage(vec.genes,gsets = Hs.c2.KEGG,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.KEGG <- sigGeneSet(gage.KEGG, qpval=c("p.val"), cutoff=0.1)
## GO
gage.GO <- gage(vec.genes,gsets = Hs.c5.GO,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.GO <- sigGeneSet(gage.GO, qpval=c("p.val"), cutoff=0.025)
## HPO Human Phenotype Ontology
gage.HP <- gage(vec.genes,gsets = Hs.c5.HP,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.HP <- sigGeneSet(gage.HP, qpval=c("p.val"), cutoff=0.025)

openxlsx::write.xlsx(list("REACTOME_up"=gage.REACTOME$greater,
                          "REACTOME_down"=gage.REACTOME$less,
                          "KEGG_up"=gage.KEGG$greater,
                          "KEGG_down"=gage.KEGG$less,
                          "GO_up"=gage.GO$greater,
                          "GO_down"=gage.GO$less,
                          "HP_up"=gage.HP$greater,
                          "HP_down"=gage.HP$less),
                     file = file.path(data.dir, "output", "breast_cancer_pathways.xlsx"),
                     rowNames =TRUE)


pdf(file.path(data.dir, "output", "Breast_cancer_GAGE.pdf"), width = 5, height = 8)
plotDot(gage.REACTOME, "REACTOME")
plotDot(gage.KEGG, "KEGG")
plotDot(gage.GO, "Gene Ontology", 20)
plotDot(gage.HP, "Human Phenotype Ontology", 20)
dev.off()



# GSEA - COLON ------------------------------------------------------------

## LOAD gene-list
df.list <- readRDS(file = file.path(data.dir, "RDS", "colon_cancer_variants.RDS"))

vec.genes <- prepGeneSet(df.list)

## REACTOME
gage.REACTOME <- gage(vec.genes,gsets = Hs.c2.REACTOME,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.REACTOME <- sigGeneSet(gage.REACTOME, qpval=c("p.val"), cutoff=0.05)
## KEGG
gage.KEGG <- gage(vec.genes,gsets = Hs.c2.KEGG,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.KEGG <- sigGeneSet(gage.KEGG, qpval=c("p.val"), cutoff=0.05)
## GO
gage.GO <- gage(vec.genes,gsets = Hs.c5.GO,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.GO <- sigGeneSet(gage.GO, qpval=c("p.val"), cutoff=0.025)
## HPO Human Phenotype Ontology
gage.HP <- gage(vec.genes,gsets = Hs.c5.HP,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.HP <- sigGeneSet(gage.HP, qpval=c("p.val"), cutoff=0.025)

openxlsx::write.xlsx(list("REACTOME_up"=gage.REACTOME$greater,
                          "REACTOME_down"=gage.REACTOME$less,
                          "KEGG_up"=gage.KEGG$greater,
                          "KEGG_down"=gage.KEGG$less,
                          "GO_up"=gage.GO$greater,
                          "GO_down"=gage.GO$less,
                          "HP_up"=gage.HP$greater,
                          "HP_down"=gage.HP$less),
                     file = file.path(data.dir, "output", "colon_cancer_pathways.xlsx"),
                     rowNames =TRUE)

pdf(file.path(data.dir, "output", "Colon_cancer_GAGE.pdf"), width = 5, height = 8)
plotDot(gage.REACTOME, "REACTOME")
#plotDot(gage.KEGG, "KEGG")
plotDot(gage.GO, "Gene Ontology", 20)
plotDot(gage.HP, "Human Phenotype Ontology", 20)
dev.off()


# GSEA - LEUKEMIA ---------------------------------------------------------

## LOAD gene-list
df.list <- readRDS(file = file.path(data.dir, "RDS", "leukemia_cancer_variants.RDS"))

vec.genes <- prepGeneSet(df.list)

## REACTOME
gage.REACTOME <- gage(vec.genes,gsets = Hs.c2.REACTOME,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.REACTOME <- sigGeneSet(gage.REACTOME, qpval=c("p.val"), cutoff=0.05)
## KEGG
gage.KEGG <- gage(vec.genes,gsets = Hs.c2.KEGG,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.KEGG <- sigGeneSet(gage.KEGG, qpval=c("p.val"), cutoff=0.05)
## GO
gage.GO <- gage(vec.genes,gsets = Hs.c5.GO,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.GO <- sigGeneSet(gage.GO, qpval=c("p.val"), cutoff=0.05)
## HPO Human Phenotype Ontology
gage.HP <- gage(vec.genes,gsets = Hs.c5.HP,    ref = NULL, samp = NULL,rank.test = TRUE, saaTest = gs.tTest, compare="unpaired")
gage.HP <- sigGeneSet(gage.HP, qpval=c("p.val"), cutoff=0.05)

openxlsx::write.xlsx(list("REACTOME_up"=gage.REACTOME$greater,
                          "REACTOME_down"=gage.REACTOME$less,
                          "KEGG_up"=gage.KEGG$greater,
                          "KEGG_down"=gage.KEGG$less,
                          "GO_up"=gage.GO$greater,
                          "GO_down"=gage.GO$less,
                          "HP_up"=gage.HP$greater,
                          "HP_down"=gage.HP$less),
                     file = file.path(data.dir, "output", "leukemia_cancer_pathways.xlsx"),
                     rowNames =TRUE)

pdf(file.path(data.dir, "output", "Leukemia_cancer_GAGE.pdf"), width = 5, height = 8)
plotDot(gage.REACTOME, "REACTOME")
plotDot(gage.KEGG, "KEGG")
plotDot(gage.GO, "Gene Ontology", 20)
plotDot(gage.HP, "Human Phenotype Ontology", 20)
dev.off()











