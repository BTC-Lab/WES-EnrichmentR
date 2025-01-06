# INIT --------------------------------------------------------------------

source("scripts/_network_diffusion.R")
source("scripts/_annovar_filtering.R")
set.seed(16341)

# BiocManager::install(c("DESeq2","airway","gage","gageData","org.Hs.eg.db",
#  "pathview","msigdb","ExperimentHub","AnnotationHub", "htmltools"))
# BiocManager::install("GSEABase")
library(ggpubr)
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
library(gridExtra)
library(ggpubr)

new_theme <- theme_bw() + 
  # adjustments for the legend
  theme(legend.position="bottom",
        legend.text = element_text(color = "black", size=12),
        legend.key = element_rect(size=12))

library(gridExtra)
library(ggpubr)
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
p.cols <- c(cols[1],cols[3],cols[5],cols[9],cols[11],cols[13],cols[15],cols[17],cols[19],cols[2],cols[4],cols[6],cols[8],cols[10],cols[7])

plt.cols <- c(cols[1],cols[3],cols[6],cols[15],cols[17],cols[19],cols[2],cols[4],cols[6],cols[8],cols[10],cols[7])
names(plt.cols)[1:3] <- c("Breast", "Leukemia", "Colon")

# FUNCTION ----------------------------------------------------------------

## Function to name facet wraps nicely
variable_labeller <- function(variable,value){
  return(variable_names[value])
}

panelPCA <- function(plot.df, metric.df, p.axes, leg.pos = "bottom") {
  
  df.healthy <- plot.df %>%
    dplyr::filter(Cancer == "no")
  
  df.case <- plot.df %>%
    dplyr::filter(Cancer == "yes")
  
  tmp.return <- ggplot(data = plot.df, aes(x=get(p.axes[1]), y=get(p.axes[2]), fill=get(cov.col))) +
    geom_point(data = df.healthy, aes(x=get(p.axes[1]), y=get(p.axes[2]), size = p.size, colour=get(cov.col)), alpha = 0.75, shape = 19) +
    geom_point(data = df.case, aes(x=get(p.axes[1]), y=get(p.axes[2]), size = p.size), alpha = 0.85, shape = 21) + 
    theme_bw() + 
    scale_size_manual(values=c(3,4)) +
    scale_color_manual(values=plt.cols, na.value = "lightgrey") +
    scale_fill_manual(values=plt.cols, na.value = "lightgrey") +
    xlab(metric.df$var.label[which(metric.df$pc.id == p.axes[1])]) +
    ylab(metric.df$var.label[which(metric.df$pc.id == p.axes[2])]) +
    labs(fill = "Cancer Status") +
    guides(fill = guide_legend(override.aes = list(size = 10)), colour="none", size="none") +
    theme(legend.position=eval(leg.pos), 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x = element_text(size =26),
          axis.title.y = element_text(size =26),
          #axis.line=element_blank(),
          panel.border=element_blank(),
          legend.title =  element_text(size =26),
          legend.text =  element_text(size =26))
  
  return(tmp.return)
  
}

prepGeneSet <- function( df.list ) {
  
  ## merge pathogenic, potential and associated into a single dataframe
  df <- rbind.data.frame(df.list$pathogenic, df.list$potential, df.list$associated)
  df <- df %>%
    separate_rows(Gene.refGene, sep = ";")
  ## order by function and by deleteriousness scores
  df <- df[order(-df$CADD_phred, df$SIFT4G_score), ]
  ## create importance score
  summarized_df <- df %>%
    mutate(gene.importance = n() - row_number() + 1) %>%
    group_by(Gene.refGene) %>%
    summarise(gene.importance = sum(gene.importance)) %>%
    arrange(desc(gene.importance)) %>%
    mutate(gsea.rank.score = row_number())
  
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


getPlinkPCAMetrics <- function(val.obj, df.obj, t.var) {
  
  out.df <- data.frame(
    pc.id = paste("PC", seq_len(length(val.obj$V1)), sep = ""),
    var.explained = round(val.obj$V1 / t.var * 100, 2)) %>%
    dplyr::mutate(var.label = paste0(pc.id, ": ", var.explained, "%")) %>%
    dplyr::mutate(axis.min = (apply(df.obj[,2:41], 2, function(col)
      min(col)))) %>%
    dplyr::mutate(axis.max = (apply(df.obj[,2:41], 2, function(col)
      max(col))))
  
  return(out.df)
}



# LEGACY: SINGLE SAMPLE - GO PATHWAY DIFFUSION - BARPLOT ------------------

meta <- readRDS(file.path("data", "RDS", "VRI_Pilot_meta_v3.RDS")) 
go.df <- readRDS(file = file.path("data", "RDS", "backgroundRemoved_go.df.RDS"))

## legacy pathway selection
l.pw.up <- c("DNA_REPAIR","DNA_METABOLIC_PROCESS", "DOUBLE-STRAND_BREAK_REPAIR",
             "DNA_CONFORMATION_CHANGE","CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS",
             "PROTEIN-DNA_COMPLEX_ASSEMBLY","TELOMERE_ORGANIZATION","TELOMERE_MAINTENANCE",
             "PROTEIN-DNA_COMPLEX_SUBUNIT_ORGANIZATION","RDNA_HETEROCHROMATIN_ASSEMBLY")
l.pw.down <- c("MITOCHONDRIAL_TRANSLATION","MITOCHONDRIAL_GENE_EXPRESSION","TRANSLATIONAL_ELONGATION",
               "MITOCHONDRIAL_TRANSLATIONAL_ELONGATION",
               "TRANSLATIONAL_TERMINATION","TRANSLATION","PEPTIDE_BIOSYNTHETIC_PROCESS",
               "AMIDE_BIOSYNTHETIC_PROCESS","PEPTIDE_METABOLIC_PROCESS",
               "CELLULAR_AMIDE_METABOLIC_PROCESS")

test2 <- subsetPWCohort(go.df$greater, n.top = 20, by.pw = l.pw.up)
plt.go_up <- plotPathwayBarCohort(test2, xlabel="Number of UP-regulated pathways", ylabel="GO pathways", w.label = 35 )

test3 <- subsetPWCohort(go.df$less, n.top = 20, by.pw = l.pw.down)
plt.go_down <- plotPathwayBarCohort(test3, xlabel="Number of DOWN-regulated pathways", ylabel="GO pathways", w.label = 35 ) 


mf1 <- ggarrange(plt.go_up, plt.go_down, nrow = 2, labels = c("A)", "B)"), 
                 font.label = list(size = 24, color = "black", 
                                   face = "bold", family = NULL), common.legend = TRUE, legend = "bottom")

ggsave(file.path("data/output", "Pathway_GO_panel_legacy.png"), plot = mf1, width = 13, height = 12, dpi = "retina")
ggsave(file.path("data/output", "Pathway_GO_panel_legacy.pdf"), plot = mf1, width = 13, height = 12)




# REVISION: SINGLE SAMPLE - GO PATHWAY DIFFUSION - BARPLOT ----------------

## load data
meta <- readRDS(file.path("data", "RDS", "VRI_Pilot_meta_v3.RDS")) 
go.df <- readRDS(file = file.path("data", "RDS", "backgroundRemoved_go.df.RDS"))

## revision pathway selection
l.pw.up <- c("DNA_REPAIR",
             "DOUBLE-STRAND_BREAK_REPAIR",
             "CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS",
             "DNA_METABOLIC_PROCESS",
             "TELOMERE_MAINTENANCE",
             "TELOMERE_ORGANIZATION",
             "CIRCULATORY_SYSTEM_DEVELOPMENT",
             "CELL_MIGRATION",
             "TRANSMEMBRANE_RECEPTOR_PROTEIN_TYROSINE_KINASE_SIGNALING_PATHWAY",
             "CELLULAR_COMPONENT_MORPHOGENESIS")
l.pw.down <- c("TRANSLATION",
               "MITOCHONDRIAL_GENE_EXPRESSION",
               "MITOCHONDRIAL_TRANSLATION",
               "PEPTIDE_BIOSYNTHETIC_PROCESS",
               "PEPTIDE_METABOLIC_PROCESS",
               "RESPIRATORY_ELECTRON_TRANSPORT_CHAIN",
               "CELLULAR_RESPIRATION",
               "PURINE-CONTAINING_COMPOUND_BIOSYNTHETIC_PROCESS",
               "SPHINGOLIPID_BIOSYNTHETIC_PROCESS",
               "MITOCHONDRIAL_RESPIRATORY_CHAIN_COMPLEX_I_ASSEMBLY")


## figure out which top ten pathways to select for each type
test2 <- subsetPWCohort(go.df$greater, n.top = 20, by.pw = l.pw.up)
plt.go_up <- plotPathwayBarCohort(test2, xlabel="Number of UP-regulated pathways", ylabel="GO pathways", w.label = 35 )

test3 <- subsetPWCohort(go.df$less, n.top = 20, by.pw = l.pw.down)
plt.go_down <- plotPathwayBarCohort(test3, xlabel="Number of DOWN-regulated pathways", ylabel="GO pathways", w.label = 35 ) 


mf1 <- ggarrange(plt.go_up, plt.go_down, nrow = 2, labels = c("A)", "B)"), 
                 font.label = list(size = 24, color = "black", 
                                   face = "bold", family = NULL), common.legend = TRUE, legend = "none")

ggsave(file.path("data/output", "Pathway_GO_panel_revision.png"), plot = mf1, width = 13, height = 12, dpi = "retina")
ggsave(file.path("data/output", "Pathway_GO_panel_revision.pdf"), plot = mf1, width = 13, height = 12) 





# REVISION: SINGLE SAMPLE - GO PATHWAY DIFFUSION - COHORT-BARPLOT ---------

## load data
meta <- readRDS(file.path("data", "RDS", "VRI_Pilot_meta_v3.RDS")) 
go.df <- readRDS(file = file.path("data", "RDS", "backgroundRemoved_go.df.RDS"))

## figure out which top ten pathways to select for each type
test2 <- subsetPWCohort(go.df$greater, n.top = 10)
plt.go_up <- plotPathwayBarCohort(test2, xlabel="Number of UP-regulated pathways", ylabel="GO pathways", w.label = 35 )

test3 <- subsetPWCohort(go.df$less, n.top = 10)
plt.go_down <- plotPathwayBarCohort(test3, xlabel="Number of DOWN-regulated pathways", ylabel="GO pathways", w.label = 35 ) 

mf1 <- ggarrange(plt.go_up, plt.go_down, nrow = 2, labels = c("A)", "B)"), 
                 font.label = list(size = 24, color = "black", 
                                   face = "bold", family = NULL), common.legend = TRUE, legend = "bottom")

ggsave(file.path("data/output", "Pathway_GO_panel_byCohort.png"), plot = mf1, width = 13, height = 12, dpi = "retina")
ggsave(file.path("data/output", "Pathway_GO_panel_byCohort.pdf"), plot = mf1, width = 13, height = 12)




#


# REVISION: SINGLE SAMPLE - GENE OVERLAP ----------------------------------

meta <- readRDS(file.path("data", "RDS", "VRI_Pilot_meta_v3.RDS")) 
colnames(meta)[1] <- "sampleID"
## get all associated genes
tbl.mutations <- readRDS(file = file.path("data", "RDS", "VRI_mutations_curated.RDS"))
a.genes <- unique(c(unique(tbl.mutations$colon$Gene), unique(tbl.mutations$breast$Gene),unique(tbl.mutations$leukemia$Gene)))
## get curated 'manuscript-list' for gene overlap
df.curated <- openxlsx::read.xlsx(file.path("manuscript", "Supplementary_Table.revision.xlsx"), sheet = 2, startRow = 2)
## get merged tables for gene overlap in cases
save.it <- readRDS(file.path("data/RDS", "VRI_Pilot_SingleSample_Variant_Gene_Overlap.28122024.RDS"))
## get all the single sample controls
out.list <- readRDS(file.path("data/RDS", "VRI_Pilot_CONTROLS_SingleSample_Variant_Gene_Overlap.31122024.RDS"))

#### We want to identify the variants shared across cohorts
## 0. get control cohort
df.ctrl <- rbind(out.list$breast, out.list$colon, out.list$leukemia)
df.ctrl.full <- cleanHomHet(df.ctrl, by.vaf = TRUE, by.dp = 10)  %>%
  distinct(sampleID, varID, .keep_all = TRUE) ## just make sure that the matched controls are not duplicated for different cohorts (inflating numbers)

## 1. get all. pathogenic and potential variants
df.all <- rbind(save.it$pathogenic$full, save.it$potential$full)#, save.it$associated$full)
df.all <- cleanHomHet(df.all, by.vaf = TRUE, by.dp = 10)

## 1. GET TOP 30 from total cohort but only the ones that are also in the controls
## then check for significant differences and keep only the p < 0.05 ones
df.all.onlyBG <- df.all %>%
  dplyr::filter((varID %in% unique(df.ctrl.full$varID)))

c.top.30 <- getGeneOverlap(df.all.onlyBG, n.top = 30, clean.variants = T)
table(c.top.30$Gene.refGene)[order(table(c.top.30$Gene.refGene), decreasing = TRUE)]

plt.heat.ct30 <- plotGeneOverlapHeatmap(c.top.30, meta)
draw(plt.heat.ct30, heatmap_legend_side = "bottom")

## check for p < 0.05 and keep only significant ones
test <- compareVariantCohortDist(c.top.30, df.ctrl.full)

sig.diff <- test$results %>%
  dplyr::filter(p_value < 0.05)

see <- test$contingency_details[sig.diff$varID]

# Variant: chr10_124994489_124994489_T_G is overrepresented in controls

c.top.30.sig <- c.top.30 %>%
  dplyr::filter(varID %in% sig.diff$varID) %>%
  dplyr::filter(varID != "chr10_124994489_124994489_T_G")
table(c.top.30.sig$Gene.refGene)[order(table(c.top.30.sig$Gene.refGene), decreasing = TRUE)]

plt.heat.ct30 <- plotGeneOverlapHeatmap(c.top.30.sig, meta)
draw(plt.heat.ct30, heatmap_legend_side = "bottom")

## generate outputs
plt.vaf <- plotVAF(c.top.30.sig)
ggsave(file.path("data/output", "Revision_VAF_Study_Gene_Overlap_top30_sigBG.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")
openxlsx::write.xlsx(c.top.30.sig, file.path("data/output", "Revision_VAF_Study_Gene_Overlap_top30_sigBG.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(c.top.30.sig, meta)

pdf(file.path("data/output", "Revision_Study_Gene_Overlap_top30_sigBG.pdf"), width = 13, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()




## 2. GET TOP 30 WITHOUT BACKGROUND
df.all.noBG <- df.all %>%
  dplyr::filter(!(varID %in% unique(df.ctrl.full$varID)))

c.top.30.noBG <- getGeneOverlap(df.all.noBG, n.top = 30, clean.variants = T)
table(c.top.30.noBG$Gene.refGene)[order(table(c.top.30.noBG$Gene.refGene), decreasing = TRUE)]
## make heatmap
plt.heat.ct30.noBG <- plotGeneOverlapHeatmap(c.top.30.noBG, meta)
draw(plt.heat.ct30.noBG, heatmap_legend_side = "bottom")

## generate outputs
plt.vaf <- plotVAF(c.top.30.noBG)
ggsave(file.path("data/output", "Revision_VAF_Study_Gene_Overlap_top30_noBG.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")
openxlsx::write.xlsx(c.top.30.noBG, file.path("data/output", "Revision_VAF_Study_Gene_Overlap_top30_noBG.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(c.top.30.noBG, meta)

pdf(file.path("data/output", "Revision_Study_Gene_Overlap_top30_noBG.pdf"), width = 13, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()



## 3. GET TOP 30 from total cohort + associated WITHOUT BACKGROUND
df.all <- rbind(save.it$pathogenic$full, save.it$potential$full, save.it$associated$full)
df.all <- cleanHomHet(df.all, by.vaf = TRUE, by.dp = 10)

df.all.noBG <- df.all %>%
  dplyr::filter(!(varID %in% unique(df.ctrl.full$varID)))

c.top.30 <- getGeneOverlap(df.all.noBG, n.top = 30, clean.variants = T)
table(c.top.30$Gene.refGene)[order(table(c.top.30$Gene.refGene), decreasing = TRUE)]

plt.heat.ct30 <- plotGeneOverlapHeatmap(c.top.30, meta)
draw(plt.heat.ct30, heatmap_legend_side = "bottom")

## generate outputs
plt.vaf <- plotVAF(c.top.30)
ggsave(file.path("data/output", "Revision_VAF_Study_Gene_Overlap_top30_associated_noBG.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")
openxlsx::write.xlsx(c.top.30, file.path("data/output", "Revision_VAF_Study_Gene_Overlap_top30_associated_noBG.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(c.top.30, meta)

pdf(file.path("data/output", "Revision_Study_Gene_Overlap_top30_associated_noBG.pdf"), width = 13, height = 6.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()



# REVISION: SINGLE SAMPLE - GENE OVERLAP BY COHORT ------------------------


## 1. BREAST
df.all <- rbind(save.it$pathogenic$breast, save.it$potential$breast)#, save.it$associated$full)
df.all <- cleanHomHet(df.all, by.vaf = TRUE, by.dp = 10)

df.all.noBG <- df.all %>%
  dplyr::filter(!(varID %in% unique(df.ctrl.full$varID)))

c.top.30.noBG <- getGeneOverlap(df.all.noBG, n.top = 30, clean.variants = T)
table(c.top.30.noBG$Gene.refGene)[order(table(c.top.30.noBG$Gene.refGene), decreasing = TRUE)]
## make heatmap
plt.heat.ct30.noBG <- plotGeneOverlapHeatmap(c.top.30.noBG, meta)
draw(plt.heat.ct30.noBG, heatmap_legend_side = "bottom")

## generate outputs
plt.vaf <- plotVAF(c.top.30.noBG)
ggsave(file.path("data/output", "BREAST_Revision_VAF_Study_Gene_Overlap_top30_noBG.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")
openxlsx::write.xlsx(c.top.30.noBG, file.path("data/output", "BREAST_Revision_VAF_Study_Gene_Overlap_top30_noBG.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(c.top.30.noBG, meta)

pdf(file.path("data/output", "BREAST_Revision_Study_Gene_Overlap_top30_noBG.pdf"), width = 13, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()


## 2. COLON
df.all <- rbind(save.it$pathogenic$colon, save.it$potential$colon)#, save.it$associated$full)
df.all <- cleanHomHet(df.all, by.vaf = TRUE, by.dp = 10)

df.all.noBG <- df.all %>%
  dplyr::filter(!(varID %in% unique(df.ctrl.full$varID)))

c.top.30.noBG <- getGeneOverlap(df.all.noBG, n.top = 30, clean.variants = T)
table(c.top.30.noBG$Gene.refGene)[order(table(c.top.30.noBG$Gene.refGene), decreasing = TRUE)]
## make heatmap
plt.heat.ct30.noBG <- plotGeneOverlapHeatmap(c.top.30.noBG, meta)
draw(plt.heat.ct30.noBG, heatmap_legend_side = "bottom")

## generate outputs
plt.vaf <- plotVAF(c.top.30.noBG)
ggsave(file.path("data/output", "COLON_Revision_VAF_Study_Gene_Overlap_top30_noBG.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")
openxlsx::write.xlsx(c.top.30.noBG, file.path("data/output", "COLON_Revision_VAF_Study_Gene_Overlap_top30_noBG.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(c.top.30.noBG, meta)

pdf(file.path("data/output", "COLON_Revision_Study_Gene_Overlap_top30_noBG.pdf"), width = 13, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()


## 3. LEUKEMIA
df.all <- rbind(save.it$pathogenic$leukemia, save.it$potential$leukemia)#, save.it$associated$full)
df.all <- cleanHomHet(df.all, by.vaf = TRUE, by.dp = 10)

df.all.noBG <- df.all %>%
  dplyr::filter(!(varID %in% unique(df.ctrl.full$varID)))

c.top.30.noBG <- getGeneOverlap(df.all.noBG, n.top = 30, clean.variants = T)
table(c.top.30.noBG$Gene.refGene)[order(table(c.top.30.noBG$Gene.refGene), decreasing = TRUE)]
## make heatmap
plt.heat.ct30.noBG <- plotGeneOverlapHeatmap(c.top.30.noBG, meta)
draw(plt.heat.ct30.noBG, heatmap_legend_side = "bottom")

## generate outputs
plt.vaf <- plotVAF(c.top.30.noBG)
ggsave(file.path("data/output", "LEUKEMIA_Revision_VAF_Study_Gene_Overlap_top30_noBG.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")
openxlsx::write.xlsx(c.top.30.noBG, file.path("data/output", "LEUKEMIA_Revision_VAF_Study_Gene_Overlap_top30_noBG.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(c.top.30.noBG, meta)

pdf(file.path("data/output", "LEUKEMIA_Revision_Study_Gene_Overlap_top30_noBG.pdf"), width = 13, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()




































## 3. get all matching varIDs from the controls 
df.ctrl.match <- df.ctrl.full %>%
  dplyr::filter(varID %in% unique(c.top.30$varID))

## 4. check out the variants distributions
table(c.top.30$Gene.refGene)[order(table(c.top.30$Gene.refGene), decreasing = TRUE)]
table(df.ctrl.match$Gene.refGene)[order(table(df.ctrl.match$Gene.refGene), decreasing = TRUE)]

length(unique(c.top.30$sampleID))
length(unique(df.ctrl.match$sampleID))


plt.vaf <- plotVAF(test2)
ggsave(file.path("data/output", "Revision_VAF_Study_Gene_Overlap_top20.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")
openxlsx::write.xlsx(test2, file.path("data/output", "Revision_VAF_Study_Gene_Overlap_top20.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(hmm, meta, plt.colors)

pdf(file.path("data/output", "Revision_Study_Gene_Overlap_top20.pdf"), width = 13, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()















c.top.30.plot <- c.top.30 %>%
  dplyr::filter(varID %in% muh$varID)

plt.heat_study <- plotGeneOverlapHeatmap(c.top.30.plot, meta, plt.colors)
draw(plt.heat_study, heatmap_legend_side = "bottom")


c.top.30.plot <- c.top.30 %>%
  dplyr::filter(!(varID %in% df.ctrl.full$varID))

plt.heat_study <- plotGeneOverlapHeatmap(c.top.30.plot, meta, plt.colors)
draw(plt.heat_study, heatmap_legend_side = "bottom")




## 2. get the top 30 genes with overlapping variants in 
c.top.30 <- getGeneOverlap(df.all, n.top = 30, clean.variants = TRUE)

## 3. get all matching varIDs from the controls 
df.ctrl.match <- df.ctrl.full %>%
  dplyr::filter(varID %in% unique(c.top.30$varID))

## 4. check out the variants distributions
table(c.top.30$Gene.refGene)[order(table(c.top.30$Gene.refGene), decreasing = TRUE)]
table(df.ctrl.match$Gene.refGene)[order(table(df.ctrl.match$Gene.refGene), decreasing = TRUE)]

length(unique(c.top.30$sampleID))
length(unique(df.ctrl.match$sampleID))

compareVariantCohortDist <- function(df.cases, df.ctrls) {
  ## subset ctrls to match varIDs in cases
  df.ctrls <- df.ctrls %>%
    dplyr::filter(varID %in% unique(df.cases$varID)) %>%
    dplyr::mutate(comp = "ctrls") %>%
    dplyr::select(sampleID, comp, varID)
  
  df.cases <- df.cases %>%
    dplyr::mutate(comp = "cases") %>%
    dplyr::select(sampleID, comp, varID)
  
  df.comp <- rbind(df.ctrls, df.cases)
  
  
  contingency_table <- table(df.comp$comp, df.comp$varID)
  
  
  # Determine group sizes
  group_sizes <- df.comp %>%
    dplyr::group_by(comp) %>%
    dplyr::summarise(n_samples = n_distinct(sampleID)) %>%
    tidyr::pivot_wider(names_from = comp, values_from = n_samples)
  
  # Create contingency table for each variant
  contingency_results <- list()
  
  unique_variants <- unique(df.comp$varID)
  
  for (variant in unique_variants) {
    # Count variant occurrences in cases and controls
    cases_count <- sum(df.comp$varID == variant & df.comp$comp == "cases")
    ctrls_count <- sum(df.comp$varID == variant & df.comp$comp == "ctrls")
    
    # Non-variant counts derived from group sizes
    non_cases_count <- group_sizes$cases - cases_count
    non_ctrls_count <- group_sizes$ctrls - ctrls_count
    
    # Construct contingency table
    contingency_table <- matrix(
      c(cases_count, non_cases_count, ctrls_count, non_ctrls_count),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(Group = c("Cases", "Controls"), Variant = c("Present", "Absent"))
    )
    
    # Perform Fisher's exact test or Chi-squared test
    if (any(contingency_table < 5)) {
      test_result <- fisher.test(contingency_table)
      test_type <- "Fisher"
    } else {
      test_result <- chisq.test(contingency_table)
      test_type <- "Chi-squared"
    }
    
    # Store results
    contingency_results[[variant]] <- list(
      variant = variant,
      contingency_table = contingency_table,
      p_value = test_result$p.value,
      test_type = test_type
    )
  }
  
  # Summarize results
  results_df <- do.call(rbind, lapply(contingency_results, function(res) {
    data.frame(
      varID = res$variant,
      p_value = res$p_value,
      test_type = res$test_type
    )
  }))
  
  # Return results
  list(
    results = results_df,
    group_sizes = group_sizes,
    contingency_details = contingency_results
  )
  
  
}






plt.vaf <- plotVAF(test2)
ggsave(file.path("data/output", "Revision_VAF_Study_Gene_Overlap_top20.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")
openxlsx::write.xlsx(test2, file.path("data/output", "Revision_VAF_Study_Gene_Overlap_top20.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(hmm, meta, plt.colors)

pdf(file.path("data/output", "Revision_Study_Gene_Overlap_top20.pdf"), width = 13, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()

unique(test2$Gene.refGene)

df.ctrl <- rbind(out.list$breast, out.list$colon, out.list$leukemia)
df.ctrl <- cleanHomHet(df.ctrl, by.vaf = TRUE, by.dp = 10)


df.case.var <- test2 %>%
  dplyr::select(sampleID, varID, Gene.refGene) %>%
  dplyr::filter(Gene.refGene == "CTBP2")


df.dummy <- df.ctrl %>%
  #dplyr::select(sampleID, varID, Gene.refGene) %>%
  dplyr::filter(varID %in% df.case.var$varID) %>%
  distinct(sampleID, varID, .keep_all = TRUE) %>%
  dplyr::filter(Gene.refGene == "CTBP2")




length(unique(df.case.var$varID))
length(unique(df.dummy$varID))


intersect(df.case.var$varID, df.dummy$varID)
setdiff(df.case.var$varID, df.dummy$varID)

hmm <- test2 %>%
  #dplyr::select(sampleID, varID, Gene.refGene) %>%
  dplyr::filter(varID %in% setdiff(df.case.var$varID, df.dummy$varID))

table(df.dummy$varID)

table(df.case.var$varID)


plotGeneOverlapHeatmap(df.dummy, meta, plt.colors)


### top. genes with high relevance
gene.select <- c("PMS1","CTBP2","COL18A1", "MAP3K1", "BCR", "PAX6","NOTCH4", "ZNF717", "ARID1B", "KAT6B", "RECQL4","KMT2C","CTBP2", "ATM","MSH6")



df.all.rev <- rbind(save.it$pathogenic$full, save.it$revel$full)
df.all.rev <- cleanHomHet(df.all.rev, by.vaf = TRUE, by.dp = 10)
test <- getGeneOverlap(df.all.rev, by.gene = gene.select, clean.variants = TRUE)

plt.vaf <- plotVAF(test)
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Study_Gene_Overlap_selected.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")

openxlsx::write.xlsx(test, file.path("data/output", "Table_XX_VAF_Study_Gene_Overlap_selected.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(test, meta, plt.colors)

pdf(file.path("data/output", "Figure_4_LEGACY_Study_Gene_Overlap_selected.pdf"), width = 12, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()

## REVISION STYLE:
gene.select <- c("CTBP2", "MUC4","KMT2C","MUC16","MAP3K1","ZNF717","COL18A1",
                 "ATG3","PAX6","LFNG","PMS1","RBMX","BCR","HRC","PSORS1C1",
                 "CDCP2","CHST15","TMEM254","OR6C76")

## from manuscript
gene.select <- c("BARD1", "BRCA1", "DICER1", "ERCC3", "FANCE", "FAT1", "MLH3", "MSH2", "NSD1", "NTRK1", "PALB2", "PDGFRA", "RECQL", "RECQL4", "SETBP1", "SLX4" )

## just top relevant over the three main cohorts
gene.select <- c("DNMT1","HLA-DRB1","CTBP2", "CTBP2", "POTED","MAP3K1","CPS1")


### top. genes with high relevance
gene.select <- c("PMS1","CTBP2","COL18A1", "MAP3K1", "BCR", "PAX6","NOTCH4", "ZNF717", "ARID1B", "KAT6B", "RECQL4","KMT2C","CTBP2", "ATM","MSH6")

df.all.rev <- rbind(save.it$pathogenic$full,save.it$associated$full,  save.it$revel$full)
df.all.rev <- cleanHomHet(df.all.rev, by.vaf = TRUE, by.dp = 10)
test <- getGeneOverlap(df.all.rev, by.gene = gene.select, clean.variants = TRUE)


test.sum <- test %>%
  dplyr::filter(Gene.refGene == "CTBP2") %>%
  #dplyr::filter(AAChange.refGene == "PMS1:NM_001321049:exon4:c.478dupA:p.L164Vfs*4") %>%
  #dplyr::filter(ExonicFunc.refGene == "nonsynonymous SNV") %>%
  ungroup()

openxlsx::write.xlsx(test.sum, file.path("data/output", "CTBP2_dummy.xlsx"))


table(test.sum$AAChange.refGene)
table(test.sum$cohort)
table(test.sum$Gene.refGene)

test.sum <- dplyr::left_join(test.sum, meta, by = "sampleID")
table(test.sum$Cancer.TypeClean)


plt.vaf <- plotVAF(test)
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Study_Gene_Overlap_selected.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")

openxlsx::write.xlsx(test, file.path("data/output", "Table_XX_VAF_Study_Gene_Overlap_selected.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(test, meta, plt.colors)

pdf(file.path("data/output", "Figure_4_ALT_Study_Gene_Overlap_selected.pdf"), width = 12, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()


test2 <- getGeneOverlap(df.all.rev, n.top = 30, clean.variants = TRUE)
table(test2$Gene.refGene)[order(table(test2$Gene.refGene), decreasing = TRUE)]
unique(test2$Gene.refGene)
plt.heat_study <- plotGeneOverlapHeatmap(test2, meta, plt.colors)


plt.vaf <- plotVAF(test2)
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Study_Gene_Overlap_top10.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")

openxlsx::write.xlsx(test2, file.path("data/output", "Table_XX_VAF_Study_Gene_Overlap_top10.xlsx"))


pdf(file.path("data/output", "Figure_4_Study_Gene_Overlap_top10.pdf"), width = 13, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()

#
# LEGACY: SINGLE SAMPLE - GENE OVERLAP ------------------------------------

meta <- readRDS(file.path("data", "RDS", "VRI_Pilot_meta_v3.RDS")) 
colnames(meta)[1] <- "sampleID"

## get all associated genes
a.genes <- unique(c(unique(tbl.mutations$colon$Gene), unique(tbl.mutations$breast$Gene),unique(tbl.mutations$leukemia$Gene)))
## get merged tables for gene overlap
save.it <- readRDS(file.path("data/RDS", "VRI_Pilot_SingleSample_Variant_Gene_Overlap.28122024.RDS"))
## get curated 'manuscript-list' for gene overlap
df.curated <- openxlsx::read.xlsx(file.path("manuscript", "Supplementary_Table.revision.xlsx"), sheet = 2, startRow = 2)

cols <- pals::tableau20(20)
p.cols <- c(cols[1],cols[3],cols[5],cols[7],cols[9],cols[11],cols[13],cols[15],cols[17],cols[19],cols[6],cols[6],cols[6])

all.labels <- c("frameshift_insertion","frameshift_deletion","nonsynonymous_SNV",
                "nonframeshift_deletion","stopgain","frameshift_deletion,nonsynonymous_SNV",
                "frameshift_deletion,stopgain","frameshift_deletion,frameshift_insertion",
                "startloss","nonframeshift_insertion","frameshift_insertion,nonsynonymous_SNV",
                "frameshift_deletion,frameshift_insertion,stopgain","stoploss",
                "nonframeshift_deletion,nonframeshift_insertion")

sorted_labels <- sapply(all.labels, function(x) {
  sorted_elements <- sort(unlist(strsplit(x, ","))) # Split by comma, sort, and recombine
  paste(sorted_elements, collapse = ",")
})

# View the result
names(sorted_labels) <- NULL
sorted_labels <- unique(sorted_labels)

plt.colors = structure(cols, names = sorted_labels)

## 1. legacy gene selection, but new data
gene.legacy <- c("PMS1", "KMT2C", "ABL1", "MKI67", "MSH6", "WRN", "RECQL", "POLQ", "ALK", "POLH","ATM")

df.all <- rbind(save.it$pathogenic$full, save.it$potential$full) #, save.it$associated$full)
df.all <- cleanHomHet(df.all, by.vaf = TRUE, by.dp = 10)
df.all.legacy <- getGeneOverlap(df.all, by.gene = gene.legacy, clean.variants = TRUE) 

plt.vaf <- plotVAF(df.all.legacy)
ggsave(file.path("data/output", "LEGACY_VAF_Study_Gene_Overlap_selected.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")
openxlsx::write.xlsx(df.all.legacy, file.path("data/output", "LEGACY_VAF_Study_Gene_Overlap_selected.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(df.all.legacy, meta, plt.colors)

pdf(file.path("data/output", "LEGACY_Study_Gene_Overlap_selected.pdf"), width = 12, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()


## 2. curated gene selection, but new data
gene.curated <- df.curated$Gene

df.all <- rbind(save.it$pathogenic$full, save.it$potential$full) #, save.it$associated$full)
df.all <- cleanHomHet(df.all, by.vaf = TRUE, by.dp = 10)
df.all.curated <- getGeneOverlap(df.all, by.gene = gene.curated, clean.variants = TRUE) 

plt.vaf <- plotVAF(df.all.curated)
ggsave(file.path("data/output", "CURATED_VAF_Study_Gene_Overlap_selected.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")
openxlsx::write.xlsx(df.all.curated, file.path("data/output", "CURATED_VAF_Study_Gene_Overlap_selected.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(df.all.curated, meta, plt.colors)

pdf(file.path("data/output", "CURATED_Study_Gene_Overlap_selected.pdf"), width = 12, height = 10.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()














# SINGLE SAMPLE STATS -----------------------------------------------------

meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS"))

## Load single sample stats creted by:
## awk 'BEGIN { printf("#!/bin/bash\n") } { printf("bcftools stats %s > %s_stats.txt\n", $0, $0) }' <(ls *DP10.vcf.gz) > run_bcfstats.sh
bla <- list.files(path = file.path("data", "RAW","annotation_single_samples"), 
                  pattern = ".gz_stats.txt")

df <- NULL
for( leSample in 1:length(bla) ) {
  tmp.typeID <- gsub(".DP10.vcf.gz_stats.txt", "", bla[leSample], fixed = TRUE)
  tmp.fID <- strsplit(tmp.typeID, "_(?!.*_)", perl=TRUE)[[1]][2]
  
  file_path <- file.path("data", "RAW", "annotation_single_samples", bla[leSample])
  
  
  tmp.df <- data.table::fread(cmd = paste("grep", "number of records:",
                                          file_path))
  
  tmp.df <- data.frame("typeID"=tmp.typeID, "fID"=tmp.fID, "variant"=tmp.df$V3, 
                       "count"=tmp.df$V4) %>%
    dplyr::mutate(variant = gsub(":","",variant, fixed = TRUE)) %>%
    dplyr::mutate(variant = gsub("number of ","",variant, fixed = TRUE)) %>%
    dplyr::filter(variant %in% c("records", "SNPs", "indels")) %>%
    tidyr::pivot_wider(names_from = variant, values_from = count)
  
  df <- rbind.data.frame(df, tmp.df)
}


meta <- dplyr::mutate(meta, fID = S5)

dummy <- dplyr::left_join(df, meta[, c("Sample", "Cancer.TypeClean", "fID")], by = "fID")


## and save as RDS- because I care
saveRDS(dummy, file = file.path("data", "RDS", "cancer_samples_varStats.RDS"), 
        compress = TRUE)

openxlsx::write.xlsx(dummy, file = file.path("data", "RAW", "cancer_samples_variant_stats.xlsx"))



# PLOT --------------------------------------------------------------------

meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS"))
dummy <- readRDS(file.path("data", "RDS", "cancer_samples_varStats.RDS"))

### prepare le plot
pdf <- dummy %>%
  dplyr::mutate(dummy, Cancer.TypeRedux = factor(
    ifelse(Cancer.TypeClean %in% c("breast", "colon","leukemia"), Cancer.TypeClean, "other"), 
    levels = c("breast","colon","leukemia","other")))

pdf_long <- tidyr::gather(pdf[, c("typeID","fID","Sample",
                                  "Cancer.TypeClean", "Cancer.TypeRedux",
                                  "records","SNPs","indels")], Metric, 
                          value, records:indels) %>%
  dplyr::mutate(Metric = factor(Metric, levels = c("records","SNPs","indels")))

variable_names <- list(
  "records" = "Total variants" ,
  "SNPs" = "SNPs",
  "indels" = "InDels"
)

variable_labeller <- function(variable,value){
  return(variable_names[value])
}

le_plot <- ggboxplot(pdf_long, x = "Cancer.TypeRedux", y = "value",
                     color = "Cancer.TypeRedux", palette = "jco",
                     add = "jitter") + 
  facet_wrap(~Metric, scales="fixed", ncol=3,  labeller= variable_labeller) +
  theme_bw() +
  guides(color=guide_legend(title="Cancer Type")) +
  theme(legend.position="bottom", 
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.line=element_blank(),
        strip.text.x = element_text(size = 22),
        panel.border=element_blank(),
        legend.title =  element_text(size =22),
        legend.text =  element_text(size =22))


ggsave(file.path("data/output", "Cancer_variants_boxplot.png"), plot = le_plot, width = 15, height = 5, dpi = "retina")
ggsave(file.path("data/output", "Cancer_variants_boxplot.pdf"), plot = le_plot, width = 15, height = 5)



# COVERAGE STATS ----------------------------------------------------------

data.dir <- "data"
meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS"))

## Load single sample stats creted by:
## awk 'BEGIN { printf("#!/bin/bash\n") } { printf("bcftools stats %s > %s_stats.txt\n", $0, $0) }' <(ls VRI_*.vcf) > run_bcfstats.sh
bla <- list.files(path = file.path("data", "RAW","alignment_metrics"), 
                  pattern = "coverage.txt")

for( leSample in 1:length(bla) ) {
  
  if( leSample == 1 ) {
    tmp.sID <- gsub("_hg38.recal.bam_coverage.txt", "", bla[leSample], fixed = TRUE)
    df <- data.table::fread(file.path("data","RAW","alignment_metrics", bla[leSample]))
    colnames(df) <- c("chrom", eval(tmp.sID))
  } else {
    tmp.sID <- gsub("_hg38.recal.bam_coverage.txt", "", bla[leSample], fixed = TRUE)
    tmp.df <- data.table::fread(file.path("data","RAW","alignment_metrics", bla[leSample]))
    colnames(tmp.df) <- c("chrom", eval(tmp.sID))
    
    df <- dplyr::left_join(df, tmp.df, by = "chrom")
  }
}

keeper <- setdiff(intersect(meta$Sample, colnames(df)), c())# 
c("DXB_TWR_21_22942335_00687", "DXB_DDC_21_4093291_00796"))

keeper <- setdiff(intersect(meta$Sample, colnames(df)), 
                  c(colnames(as.data.frame(df)[,which(df[1,] < 30)])))

df.sub <- df %>%
  as.data.frame(.) %>%
  dplyr::select(c("chrom",all_of(keeper)))

rownames(df.sub) <- df.sub$chrom
df.sub <- df.sub[c(paste0("chr",1:22), "chrX","chrY","chrM"),]


mean(rowMedians(as.matrix(df.sub[,-1])))
mean(rowMeans(as.matrix(df.sub[,-1])))
rowMax(as.matrix(df.sub[,-1]))
rowMin(as.matrix(df.sub[,-1]))

mean(colSums(df.sub[,-1]) / 25)
max(colSums(df.sub[,-1]) / 25)
min(colSums(df.sub[,-1]) / 25)

length(intersect(meta$Sample, colnames(df)))
setdiff(meta$Sample, colnames(df))

#
bad.ones <- colnames(as.data.frame(df)[,which(df[1,] < 30)])

blubb <- meta[meta$Sample %in% bad.ones,]

c("DXB_TWR_21_22942335_00687", "DXB_DDC_21_4093291_00796")
df[]



# ALIGNMENT STATS ---------------------------------------------------------

data.dir <- "data"
meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS"))

## Load single sample stats creted by:
## awk 'BEGIN { printf("#!/bin/bash\n") } { printf("bcftools stats %s > %s_stats.txt\n", $0, $0) }' <(ls VRI_*.vcf) > run_bcfstats.sh
blubb <- list.files(path = file.path("data", "RAW","alignment_metrics"), 
                    pattern = "metrics.txt")

for( leSample in 1:length(blubb) ) {
  
  if( leSample == 1 ) {
    tmp.sID <- gsub("_hg38.recal.bam_alignment_metrics.txt", "", blubb[leSample], fixed = TRUE)
    df <- data.table::fread(file.path("data","RAW","alignment_metrics", blubb[leSample]),
                            skip = "CATEGORY", nrows = 3, header = TRUE)
    df <- df %>%
      dplyr::mutate(sID = eval(tmp.sID))
  } else {
    tmp.sID <- gsub("_hg38.recal.bam_alignment_metrics.txt", "", blubb[leSample], fixed = TRUE)
    tmp.df <- data.table::fread(file.path("data","RAW","alignment_metrics", blubb[leSample]))
    tmp.df <- tmp.df %>%
      dplyr::mutate(sID = eval(tmp.sID))
    
    df <- rbind.data.frame(df, tmp.df)
  }
}

keeper <- intersect(meta$Sample, df$sID)

dummy <- df %>%
  dplyr::filter(sID %in% keeper)

dfull <- dummy %>%
  dplyr::filter(CATEGORY == "PAIR")

dfirst <- dummy %>%
  dplyr::filter(CATEGORY == "FIRST_OF_PAIR")

dfsecond <- dummy %>%
  dplyr::filter(CATEGORY == "SECOND_OF_PAIR")

table.obj <- list()
table.obj[["Paired_Reads"]] <- dfull
table.obj[["Reads_First"]] <- dfirst
table.obj[["Reads_Second"]] <- dfsecond

openxlsx::write.xlsx(table.obj, file = file.path("data", "output", "Alignment_statistic.xlsx"))





# VARIANT GENE MATCHING ---------------------------------------------------

df <- data.table::fread(file.path("data","RAW","cancer_subtypes", 
                                  "VRI_study_gene_match.txt"),header = FALSE)

length(unique(df$V4))


df <- NULL
for( leSample in 1:length(bla) ) {
  tmp.typeID <- gsub(".vcf_stats.txt", "", bla[leSample], fixed = TRUE)
  tmp.fID <- strsplit(tmp.typeID, "_(?!.*_)", perl=TRUE)[[1]][2]
  tmp.df <- data.table::fread(cmd = paste("grep", "number of records:",
                                          file.path("data","RAW",
                                                    "single_sample_annotation", 
                                                    bla[leSample])))
  tmp.df <- data.frame("typeID"=tmp.typeID, "fID"=tmp.fID, "variant"=tmp.df$V3, 
                       "count"=tmp.df$V4) %>%
    dplyr::mutate(variant = gsub(":","",variant, fixed = TRUE)) %>%
    dplyr::mutate(variant = gsub("number of ","",variant, fixed = TRUE)) %>%
    dplyr::filter(variant %in% c("records", "SNPs", "indels")) %>%
    tidyr::pivot_wider(names_from = variant, values_from = count)
  
  df <- rbind.data.frame(df, tmp.df)
}


# PANEL - MAIN FIGURE 1 ---------------------------------------------------
## halfway decent colors to use

meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS"))

## remove overhead of unspecified Cancers
meta.sub <- meta %>%
  dplyr::filter(!(Cancer == "yes" & is.na(Cancer.TypeClean))) %>%   # remove the unspecified types
  dplyr::filter(!(Cancer.TypeClean == "sicklecell") | Cancer == "no") %>%  # remove the not cancers
  dplyr::filter(!(Cancer.TypeClean == "acuteidiopatheicthrompocytopeniapurpura") | Cancer == "no")%>%
  dplyr::mutate(Cancer.TypeClean = str_to_title(Cancer.TypeClean)) %>%    # make nice labels & and adjust order with factors
  dplyr::mutate(Cancer.TypeRedux = factor(
    ifelse(Cancer.TypeClean %in% c("Breast", "Colon","Leukemia"), Cancer.TypeClean, "Other"), 
    levels = c("Breast","Leukemia","Colon","Other"))) %>%
  dplyr::mutate(Cancer.TypeClean = 
                  factor(Cancer.TypeClean, levels = c("Breast","Leukemia","Colon",
                                                      "Lung","Prostate","Uterine","Pancreatic","Thyroid",
                                                      "Renalcell","Sarcoma","Multiplemyeloma","Lymphoma","Bladder"))) 
## Sanity check
table(meta.sub$Cancer)

#### MAKE CANCER PIE-CHART
## select Cancer.TypeRedux OR Cancer.TypeClean
data <- data.frame(table(meta.sub$Cancer.TypeRedux[meta.sub$Cancer == "yes"]))
colnames(data) <- c("Cancer.Type","value")

# Compute the position of labels
data <- data %>% 
  arrange(desc(Cancer.Type)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop ) %>%
  mutate(., alt.labels = ifelse(value > 4, value, ""))


plt.pie <- ggplot(data, aes(x="", y=prop, fill=Cancer.Type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  scale_fill_manual(values=plt.cols) +
  geom_text(aes(y = ypos, label = alt.labels), color = "white", size=6) +
  labs(x = "Ratio of missing values", y = "Number of covariates", fill = "Cancer Type") +
  guides(color=guide_legend(title="Cancer Type")) +
  theme_bw() +
  theme_minimal() +
  theme(legend.position="none", 
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.line=element_blank(),
        strip.text.x = element_text(size = 22),
        panel.border=element_blank(),
        legend.title =  element_text(size =22),
        legend.text =  element_text(size =22)) + 
  
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.ticks = element_blank(), plot.title=element_text(size=14, face="bold"))

#### MAKE COVERAGE-BOX-PLOT
dummy <- readRDS(file.path("data", "RDS", "cancer_samples_varStats.RDS"))

pdf <- dummy %>%
  dplyr::filter(!(is.na(Cancer.TypeClean))) %>%   # remove the unspecified types
  dplyr::filter(!(Cancer.TypeClean == "sicklecell")) %>%  # remove the not cancers
  dplyr::filter(!(Cancer.TypeClean == "acuteidiopatheicthrompocytopeniapurpura"))%>%
  dplyr::mutate(Cancer.TypeClean = str_to_title(Cancer.TypeClean)) %>%
  dplyr::mutate(Cancer.TypeRedux = factor(
    ifelse(Cancer.TypeClean %in% c("Breast", "Colon","Leukemia"), Cancer.TypeClean, "Other"), 
    levels = c("Breast","Leukemia","Colon","Other"))) 

## make long format
pdf_long <- tidyr::gather(pdf[, c("typeID","fID","Sample",
                                  "Cancer.TypeClean", "Cancer.TypeRedux",
                                  "records","SNPs","indels")], Metric, 
                          value, records:indels) %>%
  dplyr::mutate(Metric = factor(Metric, levels = c("records","SNPs","indels")))

variable_names <- list(
  "records" = "Total variants" ,
  "SNPs" = "SNPs",
  "indels" = "InDels"
)

plt.box <- ggboxplot(pdf_long, x = "Cancer.TypeRedux", y = "value",
                     color = "Cancer.TypeRedux", palette = "jco",
                     add = "jitter") + 
  facet_wrap(~Metric, scales="free_y", ncol=3,  labeller= variable_labeller) +
  scale_colour_manual(values = plt.cols) +
  theme_bw() +
  theme_minimal() +
  guides(color=guide_legend(title="Cancer Type")) +
  theme(legend.position="top", 
        #axis.text.x=element_blank(),
        #axis.text.y=element_blank(),
        #axis.ticks=element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        #axis.line=element_blank(),
        strip.text.x = element_text(size = 22),
        panel.border=element_blank(),
        legend.title =  element_text(size =22),
        legend.text =  element_text(size =22))


#### MAKE TABLE PLOT

df.table <- data.frame("\t" = c("Study (n=298)", "Control (n=236)", "Case (n=62)", "Breast (n=18)", "Leukemia (n=18)", "Colon (n=11)", "Other (n=15)"),
                       "Total Variants" = c(71152976,70425585,727391,58271,58698,56722,553650),
                       "SNPs" = c(61549809,60886911,662898,52897,53416,51482,505103),
                       "InDels" = c(8680637,8620336,60301,5374,5282,5290,44355))

df.table <- df.table %>%
  dplyr::mutate(Total.Variants = formatC(Total.Variants, digits = 0,big.mark = ",", format = "f")) %>%
  dplyr::mutate(SNPs = formatC(SNPs, digits = 0,big.mark = ",", format = "f")) %>%
  dplyr::mutate(InDels = formatC(InDels, digits = 0,big.mark = ",", format = "f"))

colnames(df.table) <- c("","Total Variants","SNPs","InDels")

plt.table <- ggpubr::ggtexttable(df.table, rows = NULL) 
#%>%
#  tab_add_title(text = "subtitle", face = "plain", size = 10) %>%
#  tab_add_title(text = "main.title", face = "bold", padding = unit(0.1, "line")) %>%
#  tab_add_footnote(text = "*Table created using ggpubr", size = 10, face = "italic")


#### MAKE THE PANEL
mf1 <- ggarrange(ggarrange(plt.pie, plt.table, ncol = 2, widths = c(1.2,1.8), labels = c("A)", "B)"), 
                           font.label = list(size = 24, color = "black", 
                                             face = "bold", family = NULL)), plt.box,
                 nrow = 2, labels = c("", "C)"), 
                 font.label = list(size = 24, color = "black", face = "bold", family = NULL))

ggsave(file.path("data/output", "MainFigure_1_panel.png"), plot = mf1, width = 12, height = 8, dpi = "retina")
ggsave(file.path("data/output", "MainFigure_1_panel.pdf"), plot = mf1, width = 12, height = 8)



# PANEL - MAIN FIGURE 2 ---------------------------------------------------
library(ComplexHeatmap)
df <- openxlsx::read.xlsx(file.path("data/RAW", "Variant Figure.xlsx")
                          , sheet = 3)

##### DO ONLY MAIN COHORTS
pdf <- df %>%
  dplyr::mutate(mia = apply(df, 1, function(x) sum(is.na(x)))) %>%
  dplyr::filter(Cancer %in% c("Breast","Colon","Leukemia")) %>%
  dplyr::arrange(Cancer, mia) %>%
  dplyr::mutate(cohort.ID = c(paste("B",1:15, sep=""),paste("C",1:9, sep=""),paste("L",1:16, sep="")), .before = Sample.ID) %>%
  dplyr::select(c("Cancer","cohort.ID","Sample.ID","PMS1","KMT2C","ABL1",
                  "MKI67","MSH6","WRN","RECQL","POLQ","ALK",
                  "POLH","ATM"))
## fix colors
labels <- unique(c(pdf[,4],pdf[,5],pdf[,6],pdf[,7],pdf[,8],pdf[,9],pdf[,10]))
labels <- labels[!is.na(labels)]
colors = structure(plt.cols[c(1,2,3,6,11)], names = labels)
# fix labels
heat_legend_labels <- c("Frameshift Ins. (FI)", "FI/SNV", "Framshift Del. (FD)", "SNV", "non-FD")


mat <- t(as.matrix(pdf))

colnames(mat) <- mat[2,]

h.map <- ComplexHeatmap::Heatmap(
  (mat[-c(1,2,3),]), na_col = "lightgrey",column_split = mat[1,], column_title=NULL, row_km=TRUE, col = colors,
  row_title=NULL , column_names_rot = 45,
  top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "lightgrey"),
                                                      labels = c("Breast", "Colon", "Leukemia"), 
                                                      labels_gp = gpar(col = "black", fontsize = 12))),
  name = "Variant types", 
  rect_gp = gpar(col = "white", lwd = 2), 
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  heatmap_legend_param = list(labels = heat_legend_labels, at = labels, direction = "horizontal", nrow = 1)
)


pdf(file.path("data/output", "MainFigure_2_panel.pdf"), width = 12, height = 4.5)
draw(h.map, heatmap_legend_side = "bottom")
dev.off()



##### DO ALL SAMPLES MAIN COHORTS
pdf <- df %>%
  dplyr::mutate(mia = apply(df, 1, function(x) sum(is.na(x)))) %>%
  #dplyr::filter(Cancer %in% c("Breast","Colon","Leukemia")) %>%
  dplyr::mutate(Cancer = ifelse(Cancer %in% c("Breast", "Colon", "Leukemia"), Cancer, "Other")) %>%
  dplyr::arrange(Cancer, mia) %>%
  dplyr::mutate(cohort.ID = c(paste("B",1:15, sep=""),paste("C",1:9, sep=""),paste("L",1:16, sep=""),paste("O",1:12, sep="")), .before = Sample.ID) %>%
  dplyr::select(c("Cancer","cohort.ID","Sample.ID","PMS1","KMT2C","ABL1",
                  "MKI67","MSH6","WRN","RECQL","POLQ","ALK",
                  "POLH","ATM"))
## fix colors
labels <- unique(c(pdf[,4],pdf[,5],pdf[,6],pdf[,7],pdf[,8],pdf[,9],pdf[,10]))
labels <- labels[!is.na(labels)]
colors = structure(plt.cols[c(1,2,3,6,11)], names = labels)
# fix labels
heat_legend_labels <- c("Frameshift Ins. (FI)", "FI/SNV", "Framshift Del. (FD)", "SNV", "non-FD")


mat <- t(as.matrix(pdf))

colnames(mat) <- mat[2,]

h.map <- ComplexHeatmap::Heatmap(
  (mat[-c(1,2,3),]), na_col = "lightgrey",column_split = mat[1,], column_title=NULL, row_km=TRUE, col = colors,
  row_title=NULL , column_names_rot = 45,
  top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "lightgrey"),
                                                      labels = c("Breast", "Colon", "Leukemia", "Other"), 
                                                      labels_gp = gpar(col = "black", fontsize = 12))),
  name = "Variant types", 
  rect_gp = gpar(col = "white", lwd = 2), 
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  heatmap_legend_param = list(labels = heat_legend_labels, at = labels, direction = "horizontal", nrow = 1)
)


pdf(file.path("data/output", "MainFigure_2_panel_full.pdf"), width = 15, height = 4.5)
draw(h.map, heatmap_legend_side = "bottom")
dev.off()



# PCA - PREP --------------------------------------------------------------

meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS")) 

meta <- meta %>%
  dplyr::filter(!(Cancer == "yes" & is.na(Cancer.TypeClean))) %>%   # remove the unspecified types
  dplyr::filter(!(Cancer.TypeClean == "sicklecell") | Cancer == "no") %>%  # remove the not cancers
  dplyr::filter(!(Cancer.TypeClean == "acuteidiopatheicthrompocytopeniapurpura") | Cancer == "no")%>%
  dplyr::mutate(Cancer.TypeClean = str_to_title(Cancer.TypeClean)) %>%    # make nice labels & and adjust order with factors
  dplyr::mutate(Cancer.TypeRedux = factor(
    ifelse(is.na(Cancer.TypeClean), "Healthy", ifelse(Cancer.TypeClean  %in% c("Breast", "Colon","Leukemia"), Cancer.TypeClean, "Other")), 
    levels = c("Breast","Leukemia","Colon","Other", "Healthy"))) %>%
  dplyr::mutate(Cancer.TypeClean = 
                  factor(Cancer.TypeClean, levels = c("Breast","Leukemia","Colon",
                                                      "Lung","Prostate","Uterine","Pancreatic","Thyroid",
                                                      "Renalcell","Sarcoma","Multiplemyeloma","Lymphoma","Bladder"))) 

## normal filtering
## load the output
df.vec <- data.table::fread( file.path("data/RAW/plink", "VRI_all_study_PCA.eigenvec")) %>%
  dplyr::select(-V1)
df.val <- data.table::fread( file.path("data/RAW/plink", "VRI_all_study_PCA.eigenval")) 
total.var = sum(df.val$V1)
metric.df <- getPlinkPCAMetrics(df.val, df.vec, total.var)

## make plot df
plot.df <- df.vec[,1:41] %>% data.frame(stringsAsFactors = FALSE) %>%
  dplyr::rename_at(seq_len(dim(df.vec[,2:41])[2])+1,~paste("PC",seq_len(dim(df.vec[,2:41])[2]), sep = "")) %>%
  dplyr::rename_at(1, ~"Sample") %>%
  dplyr::left_join(., meta, by = "Sample")

## get metrics

## save-it
saveRDS(list(plot.df, metric.df), file = file.path("data", "RDS", "VRI_all_study_PCA.RDS"), compress = TRUE)


## less stringent filtering
## load the output


df.vec <- data.table::fread( file.path("data/RAW/plink", "VRI_all_study_alternate_PCA.eigenvec")) %>%
  dplyr::select(-V1)
df.val <- data.table::fread( file.path("data/RAW/plink", "VRI_all_study_alternate_PCA.eigenval")) 
total.var = sum(df.val$V1)
metric.df <- getPlinkPCAMetrics(df.val, df.vec, total.var)

## make plot df
plot.df <- df.vec[,1:41] %>% data.frame(stringsAsFactors = FALSE) %>%
  dplyr::rename_at(seq_len(dim(df.vec[,2:41])[2])+1,~paste("PC",seq_len(dim(df.vec[,2:41])[2]), sep = "")) %>%
  dplyr::rename_at(1, ~"Sample") %>%
  dplyr::left_join(., meta, by = "Sample")

## get metrics

## save-it
saveRDS(list(plot.df, metric.df), file = file.path("data", "RDS", "VRI_all_study_alternate_PCA.RDS"), compress = TRUE)


## with LD
## load the output
df.vec <- data.table::fread( file.path("data/RAW/plink", "VRI_all_study_alternate_pruned_PCA.eigenvec")) %>%
  dplyr::select(-V1)
df.val <- data.table::fread( file.path("data/RAW/plink", "VRI_all_study_alternate_pruned_PCA.eigenval")) 
total.var = sum(df.val$V1)
metric.df <- getPlinkPCAMetrics(df.val, df.vec, total.var)

## make plot df
plot.df <- df.vec[,1:41] %>% data.frame(stringsAsFactors = FALSE) %>%
  dplyr::rename_at(seq_len(dim(df.vec[,2:41])[2])+1,~paste("PC",seq_len(dim(df.vec[,2:41])[2]), sep = "")) %>%
  dplyr::rename_at(1, ~"Sample") %>%
  dplyr::left_join(., meta, by = "Sample")

## get metrics

## save-it
saveRDS(list(plot.df, metric.df), file = file.path("data", "RDS", "VRI_all_study_alternate_pruned_PCA.RDS"), compress = TRUE)


# PCA FIGURE --------------------------------------------------------------

meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS")) 

meta <- meta %>%
  dplyr::filter(!(Cancer == "yes" & is.na(Cancer.TypeClean))) %>%   # remove the unspecified types
  dplyr::filter(!(Cancer.TypeClean == "sicklecell") | Cancer == "no") %>%  # remove the not cancers
  dplyr::filter(!(Cancer.TypeClean == "acuteidiopatheicthrompocytopeniapurpura") | Cancer == "no")%>%
  dplyr::mutate(Cancer.TypeClean = str_to_title(Cancer.TypeClean)) %>%    # make nice labels & and adjust order with factors
  dplyr::mutate(Cancer.TypeRedux = factor(
    ifelse(is.na(Cancer.TypeClean), "Healthy", ifelse(Cancer.TypeClean  %in% c("Breast", "Colon","Leukemia"), Cancer.TypeClean, "Other")), 
    levels = c("Breast","Leukemia","Colon","Other", "Healthy"))) %>%
  dplyr::mutate(Cancer.TypeClean = 
                  factor(Cancer.TypeClean, levels = c("Breast","Leukemia","Colon",
                                                      "Lung","Prostate","Uterine","Pancreatic","Thyroid",
                                                      "Renalcell","Sarcoma","Multiplemyeloma","Lymphoma","Bladder"))) 

pca.df <- readRDS(file.path("data", "RDS", "VRI_all_study_PCA.RDS"))

plot.df <- pca.df[[1]] %>%
  dplyr::filter(Sample %in% meta$Sample) %>%
  dplyr::select(Sample:PC20) %>%
  dplyr::left_join(meta, by = "Sample") %>%
  dplyr::mutate(p.size = factor(ifelse(Cancer == "yes", 3, 2)))
metric.df <- pca.df[[2]]

p.axes <- c("PC1","PC3")
cov.col <- "Cancer.TypeRedux"

names(plt.cols) <- c("Breast","Leukemia","Colon","Other","d1","d2","Healthy")
plt.cols[7] <- "#AED6F1"

pca12 <- panelPCA(plot.df, metric.df, c("PC1","PC2"), "none")
pca23 <- panelPCA(plot.df, metric.df, c("PC2","PC3"), "bottom")
pca34 <- panelPCA(plot.df, metric.df, c("PC3","PC4"), "none")

## generate panel
options(repr.plot.width=24, repr.plot.height=12)
plt.pcaPanel_hor <- ggarrange(pca12, pca23, pca34, ncol = 3, labels = c("A)", "B)","C)"), 
                              font.label = list(size = 24, color = "black", 
                                                face = "bold", family = NULL), common.legend = TRUE, legend="bottom")

ggsave(file.path("data/output", "PCA_horizontal_panel.png"), plot = plt.pcaPanel_hor, width = 16, height = 6, dpi = "retina")
ggsave(file.path("data/output", "PCA_horizontal_panel.pdf"), plot = plt.pcaPanel_hor, width = 16, height = 6)




# GAGE - GO PATHWAY PLOT --------------------------------------------------

meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS")) 

meta <- meta %>%
  dplyr::filter(!(Cancer == "yes" & is.na(Cancer.TypeClean))) %>%   # remove the unspecified types
  dplyr::filter(!(Cancer.TypeClean == "sicklecell") | Cancer == "no") %>%  # remove the not cancers
  dplyr::filter(!(Cancer.TypeClean == "acuteidiopatheicthrompocytopeniapurpura") | Cancer == "no")%>%
  dplyr::mutate(Cancer.TypeClean = str_to_title(Cancer.TypeClean)) %>%    # make nice labels & and adjust order with factors
  dplyr::mutate(Cancer.TypeRedux = factor(
    ifelse(is.na(Cancer.TypeClean), "Healthy", ifelse(Cancer.TypeClean  %in% c("Breast", "Colon","Leukemia"), Cancer.TypeClean, "Other")), 
    levels = c("Breast","Leukemia","Colon","Other", "Healthy"))) %>%
  dplyr::mutate(Cancer.TypeClean = 
                  factor(Cancer.TypeClean, levels = c("Breast","Leukemia","Colon",
                                                      "Lung","Prostate","Uterine","Pancreatic","Thyroid",
                                                      "Renalcell","Sarcoma","Multiplemyeloma","Lymphoma","Bladder"))) 

#go.df <- readRDS(file.path("data", "RDS", "VRI_single_sample_cancer_GAGE_GO.RDS")) 

## HMMM
go.df <- readRDS(file = file.path("data", "RDS", "backgroundRemoved_go.df.RDS"))
go.df <- readRDS(file = file.path("data", "RDS", "singlesample_withBG_go.df.RDS"))
## HMMM
singlesample_withBG_go.df.RDS

## GO
## subset to reasonable numbers



input.df <- go.df$greater


go.up <- cleanPW(go.df$greater, cutoff = 1, filter.by = "total")

 #%>%
  dplyr::filter(pathway %in% c("DNA_REPAIR","DNA_METABOLIC_PROCESS", "DOUBLE-STRAND_BREAK_REPAIR",
                               "DNA_CONFORMATION_CHANGE","CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS",
                               "PROTEIN-DNA_COMPLEX_ASSEMBLY","TELOMERE_ORGANIZATION","TELOMERE_MAINTENANCE",
                               "PROTEIN-DNA_COMPLEX_SUBUNIT_ORGANIZATION","RDNA_HETEROCHROMATIN_ASSEMBLY"))


go.down <- go.df$less %>%
  dplyr::mutate(Cancer.TypeClean = str_to_title(cohort)) %>% 
  dplyr::mutate(Cancer.TypeRedux = factor(
    ifelse(is.na(Cancer.TypeClean), "Healthy", ifelse(Cancer.TypeClean  %in% c("Breast", "Colon","Leukemia"), Cancer.TypeClean, "Other")), 
    levels = c("Breast","Leukemia","Colon","Other", "Healthy"))) %>%
  dplyr::filter(Cancer.TypeRedux != "Other") %>%
  filter(pathway %in% names(table(go.df$less$pathway))[which(table(go.df$less$pathway) > 10)]) %>%
  dplyr::filter(pathway %in% c("MITOCHONDRIAL_TRANSLATION","MITOCHONDRIAL_GENE_EXPRESSION","TRANSLATIONAL_ELONGATION",
                               "MITOCHONDRIAL_TRANSLATIONAL_ELONGATION",
                               "TRANSLATIONAL_TERMINATION","TRANSLATION","PEPTIDE_BIOSYNTHETIC_PROCESS",
                               "AMIDE_BIOSYNTHETIC_PROCESS","PEPTIDE_METABOLIC_PROCESS",
                               "CELLULAR_AMIDE_METABOLIC_PROCESS"))

plt.pwUP <- plotPathwayBar(go.up, xlabel="Number of UP-regulated pathways", ylabel="GO pathways", w.label = 20 )
plt.pwDOWN <- plotPathwayBar(go.down, xlabel="Number of DOWN-regulated pathways", ylabel="GO pathways", w.label = 20 )


test <- subsetPWCohort(go.df$greater, n.top = 50, by.cohort = 0)

test <- test %>%
  arrange(desc(n.total), desc(n.cohort))

data.table::fwrite(list(unique(test$pathway)), file.path("data/output", "top50.GO_pathways_greater.txt"), sep = "\t")

test <- subsetPWCohort(go.df$less, n.top = 50, by.cohort = 0)

test <- test %>%
  arrange(desc(n.total), desc(n.cohort))

data.table::fwrite(list(unique(test$pathway)), file.path("data/output", "top50.GO_pathways_less.txt"), sep = "\t")




plotPathwayBarCohort(test, xlabel="Number of UP-regulated pathways", ylabel="GO pathways", w.label = 20 )

pw.select <- c("DNA_REPAIR","DNA_METABOLIC_PROCESS", "DOUBLE-STRAND_BREAK_REPAIR",
               "DNA_CONFORMATION_CHANGE","CELLULAR_RESPONSE_TO_DNA_DAMAGE_STIMULUS",
               "PROTEIN-DNA_COMPLEX_ASSEMBLY","TELOMERE_ORGANIZATION","TELOMERE_MAINTENANCE",
               "PROTEIN-DNA_COMPLEX_SUBUNIT_ORGANIZATION","RDNA_HETEROCHROMATIN_ASSEMBLY")







mf1 <- ggarrange(plt.pwUP, plt.pwDOWN, nrow = 2, labels = c("A)", "B)"), 
                 font.label = list(size = 24, color = "black", 
                                   face = "bold", family = NULL))

ggsave(file.path("data/output", "Pathway_GO_panel.png"), plot = mf1, width = 10, height = 15, dpi = "retina")
ggsave(file.path("data/output", "Pathway_GO_panel.pdf"), plot = mf1, width = 10, height = 15)



8.26 inch x 11.69


plotPathwayBar <- function( p.df, xlabel="Number of UP-regulated pathways", ylabel="REACTOME pathways", w.label=35) {
  
  ggplot(p.df, aes(y = forcats::fct_inorder(pathway, ordered=TRUE), 
                   fill = Cancer.TypeRedux)) +
    geom_bar(width = .85, color= "black") +
    scale_fill_manual(values=plt.cols) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
    scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "_" ," "),
                                                   width = w.label)) +
    facet_grid(cols = vars(Cancer.TypeRedux), scales = "free_x", labeller = label_parsed) +
    xlab(eval(xlabel)) + ylab(eval(ylabel)) + 
    this_theme() +
    theme(legend.position="none", 
          #axis.text.x=element_blank(),
          axis.text.y=element_text(size =12),
          axis.ticks.y=element_blank(),
          axis.title.x = element_text(size =16),
          axis.title.y = element_blank(),
          #axis.line=element_blank(),
          strip.text.x = element_text(size = 16),
          panel.border=element_blank())
  
}



theme_bw() + 
  scale_size_manual(values=c(3,4)) +
  scale_color_manual(values=plt.cols, na.value = "lightgrey") +
  scale_fill_manual(values=plt.cols, na.value = "lightgrey") +
  xlab(metric.df$var.label[which(metric.df$pc.id == p.axes[1])]) +
  ylab(metric.df$var.label[which(metric.df$pc.id == p.axes[2])]) +
  labs(fill = "Cancer Status") +
  guides(fill = guide_legend(override.aes = list(size = 10)), colour="none", size="none") +
  theme(legend.position=eval(leg.pos), 
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x = element_text(size =26),
        axis.title.y = element_text(size =26),
        #axis.line=element_blank(),
        panel.border=element_blank(),
        legend.title =  element_text(size =26),
        legend.text =  element_text(size =26))




# UMAP FIGURE -------------------------------------------------------------

library(plotly) 
library(umap) 

pca.df <- readRDS(file.path("data", "RDS", "VRI_all_study_PCA.RDS"))

plot.df <- pca.df[[1]] %>%
  dplyr::filter(Sample %in% meta$Sample) %>%
  dplyr::select(Sample:PC40) %>%
  dplyr::left_join(meta, by = "Sample") %>%
  dplyr::mutate(p.size = factor(ifelse(Cancer == "yes", 3, 2)))
metric.df <- pca.df[[2]]


u.input <- plot.df[,2:29]

uae.umap = umap(u.input, n_components = 2, random_state = 15) 


layout <- uae.umap[["layout"]] %>% data.frame()

final <- cbind(layout, plot.df) 

fig <- plot_ly(final, x = ~X1, y = ~X2, color = ~final$Cancer.TypeRedux, 
               #colors = c('#636EFA','#EF553B','#00CC96'), 
               type = 'scatter', mode = 'markers')%>%  
  
  layout(plot_bgcolor = "#e5ecf6",
         legend=list(title=list(text='species')), 
         xaxis = list( 
           title = "0"),  
         yaxis = list( 
           title = "1")) 

fig


popCols <- c("#1F77B4", "#FFBFD4", "#82CBFF", "#2CA02C", "#930000", "#000000")
names(popCols) <- c('African','American','East Asian','European','South Asian','Middle East')

plot.df <- final %>%
  dplyr::mutate(shape = ifelse(Cancer == "yes", "circle", "square"))

fig2 <- plot_ly(plot.df, x = ~X1, y = ~X2, z = ~X3, color = ~final$Cancer.TypeRedux, size = 0.5, #colors = plt.cols, 
                symbol = ~final$Cancer.TypeRedux, symbols = c("square", "circle"))
fig2 <- fig2 %>% add_markers() 
fig2 <- fig2 %>% layout(scene = list(xaxis = list(title = '0'), 
                                     yaxis = list(title = '1'), 
                                     zaxis = list(title = '2'))) 


fig2

htmlwidgets::saveWidget(
  widget = fig2, #the plotly object
  file = file.path("data/output","UAE_1KG_UMAP.html"), #the path & file name
  selfcontained = TRUE #creates a single html file
)










ggplot(data = plot.df, aes(x=get(p.axes[1]), y=get(p.axes[2]), fill=get(cov.col))) +
  geom_point(data = df.healthy, aes(x=get(p.axes[1]), y=get(p.axes[2]), size = p.size, colour=get(cov.col)), alpha = 0.95, shape = 19) +
  geom_point(data = df.case, aes(x=get(p.axes[1]), y=get(p.axes[2]), size = p.size), alpha = 0.85, shape = 21) + 
  
  theme_bw() + 
  scale_size_manual(values=c(3,4)) +
  scale_color_manual(values=plt.cols, na.value = "lightgrey") +
  scale_fill_manual(values=plt.cols, na.value = "lightgrey") +
  xlab(metric.df$var.label[which(metric.df$pc.id == p.axes[1])]) +
  ylab(metric.df$var.label[which(metric.df$pc.id == p.axes[2])]) +
  labs(colour = "Cancer Status")





panelPCA <- function(plot.df, metric.df, p.axes, leg.pos = "bottom") {
  
  df.healthy <- plot.df %>%
    dplyr::filter(Cancer == "no")
  
  df.case <- plot.df %>%
    dplyr::filter(Cancer == "yes")
  
  tmp.return <- ggplot(data = plot.df, aes(x=get(p.axes[1]), y=get(p.axes[2]), fill=get(cov.col))) +
    geom_point(data = df.healthy, aes(x=get(p.axes[1]), y=get(p.axes[2]), size = p.size, colour=get(cov.col)), alpha = 0.95, shape = 19) +
    geom_point(data = df.case, aes(x=get(p.axes[1]), y=get(p.axes[2]), size = p.size), alpha = 0.85, shape = 21) + 
    theme_bw() + 
    scale_size_manual(values=c(3,4)) +
    scale_color_manual(values=plt.cols, na.value = "lightgrey") +
    scale_fill_manual(values=plt.cols, na.value = "lightgrey") +
    xlab(metric.df$var.label[which(metric.df$pc.id == p.axes[1])]) +
    ylab(metric.df$var.label[which(metric.df$pc.id == p.axes[2])]) +
    labs(fill = "Cancer Status") +
    guides(fill = guide_legend(override.aes = list(size = 10)), colour="none", size="none") +
    theme(legend.position=eval(leg.pos), 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x = element_text(size =26),
          axis.title.y = element_text(size =26),
          #axis.line=element_blank(),
          panel.border=element_blank(),
          legend.title =  element_text(size =26),
          legend.text =  element_text(size =26))
  
  return(tmp.return)
  
}



#### OLD CODE
mat <- (as.matrix(pdf))

#rownames(mat) <- NULL#mat[,1]

ComplexHeatmap::Heatmap(
  (mat[,-1]), na_col = "lightgrey",row_split = mat[,1],
  #col = c("None" = "white", "SNP" = "blue", "Insertion" = "green", "Deletion" = "red"),
  name = "Categorical Heatmap", rect_gp = gpar(col = "white", lwd = 2),
  cluster_rows = FALSE,  # Optional: Set to TRUE if you want to cluster rows
  cluster_columns = FALSE  # Optional: Set to TRUE if you want to cluster columns
)


panelPCA <- function(plot.df, metric.df, p.axes) {
  tmp.return <- ggplot(data=plot.df, aes(x=get(p.axes[1]), y=get(p.axes[2]), fill=Cancer, colour=Cancer.Type)) +
    #geom_point(alpha=plot.df$alpha, size=plot.df$p.size, shape=plot.df$shape) + 
    theme_bw() +
    #scale_fill_manual(values = popCols, labels = c('African','American','East Asian','European','South Asian')) +
    #scale_colour_manual(values = popCols, name = "Admixture",
    #                    labels = c("Component 1","Component 2","Component 3","Component 4","Component 5"), 
    #                    breaks = c("European","South Asian","African","American","East Asian")) +
    #xlim(metric.df[eval(p.axes[1]), "axis.min"], metric.df[eval(p.axes[1]), "axis.max"]) +
    #ylim(metric.df[eval(p.axes[2]), "axis.min"], metric.df[eval(p.axes[2]), "axis.max"]) +
    #xlab(metric.df[eval(p.axes[1]), "var.label"]) +
    #ylab(metric.df[eval(p.axes[2]), "var.label"]) + 
    guides(fill="none", colour = guide_legend(override.aes = list(size = 10))) +
    #xlab(eval(p.axes[1])) +
    #ylab(eval(p.axes[2])) + 
    theme(legend.position="none", 
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x = element_text(size =26),
          axis.title.y = element_text(size =26),
          #axis.line=element_blank(),
          panel.border=element_blank(),
          legend.title =  element_text(size =26),
          legend.text =  element_text(size =26))
  
  return(tmp.return)
  
}




# SUPPLEMENTARY FIGURE 1 - IMPUTATION - MATCHING --------------------------

### get meta and subset to final data
meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS"))

## remove overhead of unspecified Cancers
meta.sub <- meta %>%
  dplyr::filter(!(Cancer == "yes" & is.na(Cancer.TypeClean))) %>%   # remove the unspecified types
  dplyr::filter(!(Cancer.TypeClean == "sicklecell") | Cancer == "no") %>%  # remove the not cancers
  dplyr::filter(!(Cancer.TypeClean == "acuteidiopatheicthrompocytopeniapurpura") | Cancer == "no")%>%
  dplyr::mutate(Cancer.TypeClean = str_to_title(Cancer.TypeClean)) %>%    # make nice labels & and adjust order with factors
  dplyr::mutate(Cancer.TypeRedux = factor(
    ifelse(Cancer.TypeClean %in% c("Breast", "Colon","Leukemia"), Cancer.TypeClean, "Other"), 
    levels = c("Breast","Leukemia","Colon","Other"))) %>%
  dplyr::mutate(Cancer.TypeClean = 
                  factor(Cancer.TypeClean, levels = c("Breast","Leukemia","Colon",
                                                      "Lung","Prostate","Uterine","Pancreatic","Thyroid",
                                                      "Renalcell","Sarcoma","Multiplemyeloma","Lymphoma","Bladder"))) 

### I. prepare BMI plot
p.cols <- c(cols[2],cols[4],viridis::viridis(7))#colfunc(7))
names(p.cols) <- c("female","male","SUW", "UW", "NW", "OW", "OCI", "OCII", "OCIII")
mu <- data.frame("Gender"=c("female","male"), 
                 "grp.mean"=c(mean(meta.sub$BMI.imp[meta.sub$Gender == "female"], na.rm = TRUE),
                              mean(meta.sub$BMI.imp[meta.sub$Gender == "male"], na.rm = TRUE)))
#########
# Plot with density area and line coloured but legend not right
p1 <- ggplot(data=data, aes(x=value)) + 
  geom_density(aes(fill=type, colour=type), alpha=0.3 ) + 
  geom_vline(data=vlines, aes(xintercept=mean_median, colour=labels), 
             linetype="dashed", size=1.5, show_guide=TRUE ) 

g1 <- ggplotGrob(p1)

# Plot with density line not coloured but legend is ok
p2 <- ggplot(data=data, aes(x=value)) + 
  geom_density(aes(fill=type), alpha=0.3 ) + 
  geom_vline(data=vlines, aes(xintercept=mean_median, colour=labels), 
             linetype="dashed", size=1.5, show_guide=TRUE )  +
  guides(fill = guide_legend(override.aes = list(linetype = 0 ))) 

g2 <- ggplotGrob(p2)


# Add legend of second plot to first plot    
g1$grobs[which(g1$layout$name=="guide-box")] <- 
  g2$grobs[which(g2$layout$name=="guide-box")]  

grid::grid.newpage()    
grid::grid.draw(g1)

#########

meta.plot <- meta.sub %>%
  dplyr::mutate(bmi.frame.col = ifelse(is.na(bmi_class), "NA","BMI"))

meta.bmi <- meta.plot %>%
  dplyr::filter(!is.na(bmi_class))

meta.na <- meta.plot %>%
  dplyr::filter(is.na(bmi_class))

plt.bmi.legend <- ggplot(data=meta.plot, aes(x=BMI.imp), show.legend = F) +
  #geom_histogram(data=meta.plot, aes(fill=Gender), color="lightgrey", binwidth=1,
  #               show.legend = FALSE) + 
  geom_vline(data=mu, aes(xintercept=grp.mean), color="black",linetype="dashed") +
  geom_jitter(data=meta.bmi, aes(x=BMI.imp,y=2.5, fill=bmi_class.imp, color=bmi.frame.col), 
              height = 2.5, shape=21, size=4, alpha = 0.85, stroke = 1.5) +
  geom_jitter(data=meta.na, aes(x=BMI.imp,y=2.5, fill=bmi_class.imp, color=bmi.frame.col), 
              height = 2.5, shape=21, size=4, alpha = 0.85, stroke = 1.5) +
  
  scale_color_manual(values=c("black","red"), na.value="red") +
  scale_fill_manual(values=p.cols, na.value="red") +
  
  facet_grid(cols = vars(Gender)) + 
  guides(color = "none") + labs(fill = "BMI class") + 
  this_theme() + theme(strip.text.x = element_text(size = 12))

g.legend <- ggplotGrob(plt.bmi.legend)


plt.bmi <- ggplot(data=meta.plot, aes(x=BMI.imp), show.legend = F) +
  geom_histogram(data=meta.plot, aes(fill=Gender), color="lightgrey", binwidth=1,
                 show.legend = FALSE) + 
  geom_vline(data=mu, aes(xintercept=grp.mean), color="black",linetype="dashed") +
  geom_jitter(data=meta.bmi, aes(x=BMI.imp,y=2.5, fill=bmi_class.imp, color=bmi.frame.col), 
              height = 2.5, shape=21, size=3, alpha = 0.85, stroke = 1.5) +
  geom_jitter(data=meta.na, aes(x=BMI.imp,y=2.5, fill=bmi_class.imp, color=bmi.frame.col), 
              height = 2.5, shape=21, size=2, alpha = 0.85, stroke = 1.25) +
  
  scale_color_manual(values=c("black","red"), na.value="red") +
  scale_fill_manual(values=p.cols, na.value="red") +
  
  facet_grid(cols = vars(Gender)) + 
  guides(color = "none") + labs(fill = "BMI class") + 
  this_theme() + theme(strip.text.x = element_text(size = 12))


g.bmi <- ggplotGrob(plt.bmi)


g.bmi$grobs[which(g.bmi$layout$name=="guide-box")] <- 
  g.legend$grobs[which(g.legend$layout$name=="guide-box")]  

grid::grid.newpage()    
grid::grid.draw(g.bmi)


### II. prepare Age plot
meta.age <- meta.sub %>% 
  mutate(., Age = 2023 - Year.imp)

mu <- data.frame("Gender"=c("female","male"), 
                 "grp.mean"=c(mean(meta.age$Age[meta.age$Gender == "female"], na.rm = TRUE),
                              mean(meta.age$Age[meta.age$Gender == "male"], na.rm = TRUE)))

p.cols <- c(cols[2],cols[4],viridis::viridis(7), "#2CA02C", "#D62728")
names(p.cols) <- c("female","male","SUW", "UW", "NW", "OW", "OCI", "OCII", "OCIII", "no","yes")

plt.age <- ggpubr::ggviolin(data = meta.age, "Cancer", "Age", color = "Cancer", fill = "Gender",
                            palette = c("black", "red", "#FC4E07"),
                            add = "boxplot", add.params = list(fill = "white")) +
  #facet_wrap(~ Cancer)
  facet_grid(cols = vars(Gender))+
  scale_fill_manual(values=c("#FFBB78","#AEC7E8"), na.value="red")+ 
  this_theme() + theme(strip.text.x = element_text(size = 12))


### III. smoking bins

plt.smoking <- ggplot(data=meta.sub, aes(x=smoking_bin, fill=Gender, color = Cancer)) +
  scale_color_manual(values=p.cols) +
  scale_fill_manual(values=p.cols, na.value="red") +
  stat_count(aes(x=smoking_bin)) +
  facet_wrap(~ Gender + Cancer) +
  labs(color = "Cancer") + 
  this_theme() + theme(strip.text.x = element_text(size = 12))

library(ggpubr)

mf1 <- ggarrange(ggarrange(plt.age, g.bmi, nrow = 2, widths = c(1,1), labels = c("A)", "B)"), 
                           font.label = list(size = 24, color = "black", 
                                             face = "bold", family = NULL)), plt.smoking,
                 nrow = 1, labels = c("", "C)"), 
                 font.label = list(size = 24, color = "black", face = "bold", family = NULL))

ggsave(file.path("data/output", "SuppFigure_1_panel.png"), plot = mf1, width = 12, height = 8, dpi = "retina")
ggsave(file.path("data/output", "SuppFigure_1_panel.pdf"), plot = mf1, width = 12, height = 8)






df.healthy <- plot.df %>%
  dplyr::filter(Cancer == "no")

df.case <- plot.df %>%
  dplyr::filter(Cancer == "yes")

ggplot(data = plot.df, aes(x=get(p.axes[1]), y=get(p.axes[2]), fill=get(cov.col))) +
  geom_point(data = df.healthy, aes(x=get(p.axes[1]), y=get(p.axes[2]), size = p.size, colour=get(cov.col)), alpha = 0.95, shape = 19) +
  geom_point(data = df.case, aes(x=get(p.axes[1]), y=get(p.axes[2]), size = p.size), alpha = 0.85, shape = 21) + 
  theme_bw() + 
  scale_size_manual(values=c(3,4)) +
  scale_color_manual(values=plt.cols, na.value = "lightgrey") +
  scale_fill_manual(values=plt.cols, na.value = "lightgrey") +
  xlab(metric.df$var.label[which(metric.df$pc.id == p.axes[1])]) +
  ylab(metric.df$var.label[which(metric.df$pc.id == p.axes[2])]) +
  labs(fill = "Cancer Status") +
  guides(fill = guide_legend(override.aes = list(size = 10)), colour="none", size="none") +
  theme(legend.position=eval(leg.pos), 
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x = element_text(size =26),
        axis.title.y = element_text(size =26),
        #axis.line=element_blank(),
        panel.border=element_blank(),
        legend.title =  element_text(size =26),
        legend.text =  element_text(size =26))




















# Example data
df <- data.frame(
  x = c(1, 2, 3, 4),
  y = c(10, 5, 8, 12),
  variable1 = c("A", "B", "A", "B"),
  variable2 = c("C", "D", "C", "D")
)

# Create ggplot object with two geom layers
p <- ggplot(df, aes(x, y)) +
  geom_bar(aes(fill = variable1), stat = "identity", position = "dodge") +
  geom_bar(aes(fill = variable2), stat = "identity", position = "dodge")

# Deactivate legend for variable2
p + guides(fill = guide_legend(override.aes = list(alpha = 0)))




