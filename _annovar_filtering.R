# INITIALIZATION ----------------------------------------------------------
# This is a nice color palette.
cols <- pals::tableau20(20)
# Take a look at it.
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

# We also want to have a certain directory structure to store output and saved
# data objects. For me this is 1. data - 1.1 output; 1.2 RDS Show warnings is 
# disabled so that we are not bothered if the directories already exist.
dir.create("data", showWarnings = FALSE)
dir.create("data/output", showWarnings = FALSE)
dir.create("data/RDS", showWarnings = FALSE)
dir.create("data/RAW", showWarnings = FALSE)

# Create a variable to our working directory for produced outputs
WORKDIR <- getwd()
FILESOURCE <- file.path(WORKDIR, "data")

# LOAD LIBRARIES ----------------------------------------------------------

require(diffuStats)
require(STRINGdb)
require(igraph)
library(GSEABase)
library(msigdb)
library(Category)
library(annotate)
library(org.Hs.eg.db)
library(DESeq2)
library(gage)
library(plyr)
library(airway)
library(tidyverse)
library(igraph)
library(readr)
library(pathview)
library(ExperimentHub)


# FUNCTIONS ---------------------------------------------------------------
## This section provides you with some helper-functions.

prepATable <- function( input.df, f.freq=0.01, assoc.col=NULL ) {
  ## replace missing exac frequencies with zero
  # input.df <- input.df %>%
  #   mutate(., varID = paste(Chr, Start, End, Ref, Alt, sep="_")) %>%
  #   mutate(., Location = paste(Chr, Start, sep=":")) %>%
  #   separate_rows(Gene.refGene, sep = ";") %>%
  #   #mutate(., varID = paste(Chr, Start, End, Ref, Alt, sep="_")) %>%
  #   mutate_at(vars(starts_with("ExAC")), ~ replace_na(., 0)) %>%
  #   dplyr::filter(., ExAC_AFR < f.freq & ExAC_AMR < f.freq & ExAC_EAS < f.freq &
  #                   ExAC_FIN < f.freq & ExAC_NFE < f.freq & ExAC_SAS < f.freq)
  

  ## replace missing exac frequencies with zero
  input.df <- input.df %>%
    mutate(., varID = paste(Chr, Start, End, Ref, Alt, sep="_")) %>%
    mutate(., Location = paste(Chr, Start, sep=":")) %>%
    separate_rows(Gene.refGene, sep = ";") %>%
    #mutate(., varID = paste(Chr, Start, End, Ref, Alt, sep="_")) %>%
    mutate_at(vars(starts_with("gnomad312")), ~ replace_na(., 0)) %>%
    #dplyr::filter(., gnomad312_AF < f.freq) #%>%
    dplyr::filter(., gnomad312_AF_afr < f.freq & gnomad312_AF_ami < f.freq & gnomad312_AF_amr < f.freq &
                    gnomad312_AF_asj < f.freq & gnomad312_AF_eas < f.freq & gnomad312_AF_fin < f.freq &
                    gnomad312_AF_mid < f.freq& gnomad312_AF_nfe < f.freq& gnomad312_AF_sas < f.freq)
  
  
  ## extract all pathogenic variants
  df.path <- input.df %>%
    dplyr::filter(., CLNSIG %in% c("Pathogenic","Pathogenic/Likely_pathogenic",
                                   "Likely_pathogenic","Likely_risk_allele",
                                   "risk_factor","Pathogenic|risk_factor"))
  ## keep all benign variants
  df.benign <- input.df %>%
    dplyr::filter(., CLNSIG %in% c("Benign","Likely_benign",
                                   "Benign/Likely_benign",
                                   "Benign/Likely_benign|association",
                                   "Benign/Likely_benign|drug_response",
                                   "Benign/Likely_benign|other",
                                   "Benign|drug_response","Benign|other",
                                   "Benign|protective",
                                   "Likely_benign|drug_response|other",
                                   "Likely_benign|other"))
  
  ## remove benign and pathogenic variants
  input.df <- input.df %>%
    dplyr::filter(., !CLNSIG %in% c("Benign","Likely_benign",
                                    "Benign/Likely_benign",
                                    "Benign/Likely_benign|association",
                                    "Benign/Likely_benign|drug_response",
                                    "Benign/Likely_benign|other",
                                    "Benign|drug_response","Benign|other",
                                    "Benign|protective",
                                    "Likely_benign|drug_response|other",
                                    "Likely_benign|other","Pathogenic",
                                    "Pathogenic/Likely_pathogenic",
                                    "Likely_pathogenic","Likely_risk_allele",
                                    "risk_factor","Pathogenic|risk_factor"))
  
  ## get the potential ones by score filtering
  df.pot <- input.df %>%
    dplyr::filter(., !is.na(CADD_phred) & !is.na(DANN_score) & !is.na(SIFT4G_score) &
                    !is.na(Polyphen2_HDIV_score) & !is.na(Polyphen2_HVAR_score)) %>%
    dplyr::filter(., as.numeric(CADD_phred) >= 20 | as.numeric(DANN_score) >= 0.95 | 
                    as.numeric(SIFT4G_score) <= 0.05 | as.numeric(Polyphen2_HDIV_score) >= 0.85 |
                    as.numeric(Polyphen2_HVAR_score) >= 0.85)
  
  ## get the remaining ones without any score
  df.unscored <- input.df %>%
    dplyr::filter(., is.na(CADD_phred) & is.na(DANN_score) & is.na(SIFT4G_score) &
                    is.na(Polyphen2_HDIV_score) & is.na(Polyphen2_HVAR_score))
  
  ## everything that has a score but is likely benign, i.e., reverse signs
  df.likely_benign <- input.df %>%
    dplyr::filter(., !is.na(CADD_phred) & !is.na(DANN_score) & !is.na(SIFT4G_score) &
                    !is.na(Polyphen2_HDIV_score) & !is.na(Polyphen2_HVAR_score)) %>%
    dplyr::filter(., as.numeric(CADD_phred) < 20 | as.numeric(DANN_score) < 0.95 | 
                    as.numeric(SIFT4G_score) > 0.05 | as.numeric(Polyphen2_HDIV_score) < 0.85 |
                    as.numeric(Polyphen2_HVAR_score) < 0.85)
  
  ## associated through curated list of genes, by cancer in description, and through function
  df.associated <- df.unscored %>%
    dplyr::filter(., (Gene.refGene %in% assoc.col) | 
                    (grepl("cancer", Function_description.refGene, ignore.case = TRUE)) | 
                    (ExonicFunc.refGene %in% c("stopgain","startloss",
                                               "nonsynonymous SNV","frameshift insertion",
                                               "frameshift deletion")) )
  
  df.unscored <- df.unscored %>%
    dplyr::filter(., !((Gene.refGene %in% assoc.col) | 
                         (grepl("cancer", Function_description.refGene, ignore.case = TRUE)) | 
                         (ExonicFunc.refGene %in% c("stopgain","startloss",
                                                    "nonsynonymous SNV","frameshift insertion",
                                                    "frameshift deletion"))) )
  
  out <- list()
  out[["pathogenic"]] <- df.path
  out[["potential"]] <- df.pot
  out[["associated"]] <- df.associated
  out[["benign"]] <- df.benign
  out[["likely_benign"]] <- df.likely_benign
  out[["unscored"]] <- df.unscored
  
  return(out)
}

## GAGE functions
prepDiffSet <- function( df.list ) {
  
  tmp.pat <- df.list$pathogenic %>%
    dplyr::mutate(diff.score = 1)
  
  tmp.pot <- df.list$potential %>%
    dplyr::mutate(diff.score = 0.50)
  
  tmp.ass <- df.list$associated %>%
    dplyr::mutate(diff.score = 0.25)
  
  df <- rbind.data.frame(tmp.pat, tmp.pot, tmp.ass)
  
  summarized_df <- df %>%
    dplyr::group_by(Gene.refGene) %>%
    dplyr::summarize(gsea.rank.score = max(diff.score)) %>%
    dplyr::arrange(desc(gsea.rank.score))
  
  input_oncogenes <- data.frame("node_id"=summarized_df$Gene.refGene, "node_score"=summarized_df$gsea.rank.score, 
                                stringsAsFactors=FALSE)  
  
  return(input_oncogenes)
}



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
  ## create importance score
  # summarized_df <- df %>%
  #   dplyr::mutate(gene.importance = n() - row_number() + 1) %>%
  #   dplyr::group_by(Gene.refGene) %>%
  #   dplyr::summarize(gene.importance = sum(gene.importance)) %>%
  #   dplyr::arrange(desc(gene.importance)) %>%
  #   dplyr::mutate(gsea.rank.score = n() - row_number() + 1)
  
  summarized_df <- df %>%
    #dplyr::mutate(gene.importance = (n() - row_number() + 1) / log2(row_number()+1)) %>%
    dplyr::mutate(gene.importance = (n() - row_number() + 1)) %>%
    dplyr::group_by(Gene.refGene) %>%
    dplyr::summarize(gsea.rank.score = sum(gene.importance)) %>%
    dplyr::arrange(desc(gsea.rank.score))
  
  vec.genes <- summarized_df$gsea.rank.score
  names(vec.genes) <- summarized_df$Gene.refGene
  
  return(vec.genes)
}

createMergeDF <- function() {
  tmp <- list()
  tmp[["greater"]] <- data.frame("pathway"=NULL,"p.geomean"=NULL,"stat.mean"=NULL,
                                 "p.val"=NULL,"q.val"=NULL,"set.size"=NULL,
                                 "exp1"=NULL, "sID"=NULL, "cohort"=NULL)
  tmp[["less"]] <- data.frame("pathway"=NULL,"p.geomean"=NULL,"stat.mean"=NULL,
                              "p.val"=NULL,"q.val"=NULL,"set.size"=NULL,
                              "exp1"=NULL, "sID"=NULL, "cohort"=NULL)
  return(tmp)
}

mergeGage <- function( df.merge, df.input, s.idx, c.idx ) {
  tmp.greater <- data.frame(df.input$greater) %>%
    tibble::rownames_to_column(., "pathway") %>%
    dplyr::mutate(sID = s.idx) %>%
    dplyr::mutate(cohort = c.idx)
  
  tmp.less <- tibble::rownames_to_column(data.frame(df.input$less), "pathway") %>%
    dplyr::mutate(sID = s.idx) %>%
    dplyr::mutate(cohort = c.idx)
  
  df.merge[["greater"]] <- rbind.data.frame(df.merge[["greater"]], tmp.greater)
  df.merge[["less"]] <- rbind.data.frame(df.merge[["less"]], tmp.less)
  
  return(df.merge)
}

plotPathwayBar <- function( p.df, xlabel="Number of UP-regulated pathways", ylabel="REACTOME pathways", w.label=35) {
  
  ggplot(p.df, aes(y = forcats::fct_inorder(pathway, ordered=TRUE), 
                   fill = cohort)) +
    geom_bar(width = .85, color= "black") +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1)) +
    scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "_" ," "),
                                                   width = w.label)) +
    facet_grid(cols = vars(cohort), scales = "free_x", labeller = label_parsed) +
    xlab(eval(xlabel)) + ylab(eval(ylabel)) + 
    this_theme()
  
}

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










