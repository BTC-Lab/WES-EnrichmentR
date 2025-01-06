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
library(ggpubr)
library(plotly)


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
    mutate_at(vars(starts_with("gnomad41")), ~ replace_na(., -1)) %>%
    #dplyr::filter(., gnomad41_genome_AF < f.freq) #%>%
    dplyr::filter(., gnomad41_genome_AF_afr < f.freq & gnomad41_genome_AF_ami < f.freq & gnomad41_genome_AF_amr < f.freq &
                    gnomad41_genome_AF_asj < f.freq & gnomad41_genome_AF_eas < f.freq & gnomad41_genome_AF_fin < f.freq &
                    gnomad41_genome_AF_mid < f.freq& gnomad41_genome_AF_nfe < f.freq& gnomad41_genome_AF_sas < f.freq)
  
  
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
                    (grepl("cancer", CLNDN, ignore.case = TRUE)) | 
                    (ExonicFunc.refGene %in% c("stopgain","startloss",
                                               "nonsynonymous SNV","frameshift insertion",
                                               "frameshift deletion")) )
  
  df.unscored <- df.unscored %>%
    dplyr::filter(., !((Gene.refGene %in% assoc.col) | 
                         (grepl("cancer", CLNDN, ignore.case = TRUE)) | 
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

plotPathwayBar <- function( p.df, xlabel="Number of Enriched Pathways", ylabel="REACTOME pathways", w.label=35) {
  
  ## get top 20
  if( length(p.df$pathway) >= 20) {
    p.df <- p.df %>%
      filter(pathway %in% names(sort(table(p.df$pathway), decreasing = TRUE))[1:20])
  } else {
    p.df <- p.df %>%
      filter(pathway %in% names(sort(table(p.df$pathway), decreasing = TRUE))[1:length(p.df$pathway)])
  }

  names(sort(table(p.df$pathway), decreasing = TRUE))[1:20]
  
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


### This will clean the enrichment pathway output and remove pathways with 
### overlaps smaller than 'cutoff'.
cleanPW <- function(input.df, cutoff = 1, filter.by = "total" ) {
  
  f.col <- ifelse(filter.by == " total", "n.total", "n.cohort")
  
  
  df.out <- input.df %>%
    dplyr::mutate(Cancer.TypeClean = str_to_title(cohort)) %>% 
    dplyr::mutate(Cancer.TypeRedux = factor(
      ifelse(is.na(Cancer.TypeClean), "Healthy", ifelse(Cancer.TypeClean  %in% c("Breast", "Colon","Leukemia"), Cancer.TypeClean, "Other")), 
      levels = c("Breast","Leukemia","Colon","Other", "Healthy"))) %>%
    dplyr::filter(Cancer.TypeRedux != "Other") %>%
    group_by(pathway) %>%
    mutate(n.total = n()) %>%
    ungroup() %>%
    group_by(pathway, cohort) %>%
    mutate(n.cohort = n()) %>%
    ungroup() %>%
    dplyr::arrange(desc(n.cohort)) %>%
    filter(eval(f.col) >= eval(cutoff))
  
  return(df.out)
}

## Create a 'hom.het' column indicating alternate variant state
cleanHomHet <- function(varTable.df, gt.col = "GT", var.col = "Variant", by.vaf = FALSE, by.dp = NULL) {
  
  varTable.df <- varTable.df %>%
    mutate( GT = sub(":.*", "", get(var.col)),
            VAF = sub(".*:", "", get(var.col)), .before = cytoBand) %>%
    dplyr::filter(VAF != ".") %>%
    mutate(VAF = as.numeric(as.character(VAF))) %>%
    dplyr::mutate(phased = ifelse(grepl("|", GT, fixed = TRUE), TRUE, FALSE),
      hom.het = sapply(GT, function(x) ifelse(x == "1/1", "hom",
                                                  ifelse(x == "1|1", "hom", 
                                                         ifelse(x == "0/0", "hom", 
                                                                ifelse(x == "0|0", "hom", "het"))))), .after = GT) %>%
    mutate(AD = sub("^[^:]*:([^:]*):.*$", "\\1", get(var.col)),
           AD1 = as.numeric(as.character(sub(",.*", "", AD))),
           AD2 = as.numeric(as.character(sub(".*,", "", AD))))
  
  if( by.vaf ) {
    varTable.df <- varTable.df %>%
      dplyr::filter(
        (hom.het == "het" & VAF >= 0.3 & VAF <= 0.7) |  
          (hom.het == "hom" & phased & VAF >= 0.3) |
          (hom.het == "hom" & !phased & VAF >= 0.85)
      )
  }
  
  if( !is.null(by.dp) ) {
    varTable.df <- varTable.df %>%
      dplyr::filter(
        (AD1 + AD2 >= eval(by.dp)))  
  }
  
  return(varTable.df)
}

## subset either to top n of the cohort or to top n total
subsetPWCohort <- function(input.df, n.top = 10, by.cohort = 1, by.pw = NULL) {
  
  if( !any(grepl("n.total", colnames(input.df))) ) {
    
    input.df <- cleanPW(input.df, cutoff = 1, filter.by = "total")
  }
  
  if( !is.null(by.pw) ) {
    top.out <- input.df %>%
      dplyr::filter(pathway %in% by.pw)
  }
  else if( by.cohort == 1 ) {
    input.df <- input.df %>%
      dplyr::mutate(pw.cohort = paste0(cohort, "_", pathway))
    subset.vec <- NULL
    ## get top n for each cohort
    for( c.idx in unique(input.df$cohort) ) {
      
      tmp.p.df <- input.df %>%
        dplyr::filter(!duplicated(pw.cohort)) %>%
        dplyr::filter(cohort == eval(c.idx)) %>%
        dplyr::arrange(desc(n.cohort)) %>%
        dplyr::slice_head(n = n.top)
      
      subset.vec <- c(subset.vec, tmp.p.df$pw.cohort)
    }
    
    top.out <- input.df %>%
      dplyr::filter(pw.cohort %in% subset.vec)
    
    
  } else {
    ## return top n overall
    tmp.p.df <- input.df %>%
      dplyr::filter(!duplicated(pathway)) %>%
      dplyr::arrange(desc(n.total)) %>%
      dplyr::slice_head(n = n.top)
    
    subset.vec <- c(subset.vec, tmp.p.df$pathway)
    
    top.out <- input.df %>%
      dplyr::filter(pathway %in% subset.vec)
    
  }
  
  top.out <- top.out  %>%
    arrange(cohort, n.cohort) %>%
    dplyr::mutate(pathway = factor(pathway, levels = c(unique(pathway))))
  
  return(top.out)
  
}


## for heatmap get the overlapping Genes
getGeneOverlap <- function(input.df, n.top = 10, by.gene = NULL, clean.variants = FALSE) {
  
  ## if we want to clean up the variants, we do before subsetting
  if( clean.variants ) {
    
    input.df <- input.df %>%
      dplyr::filter(!is.na(ExonicFunc.refGene)) %>%
      dplyr::filter(ExonicFunc.refGene != "unknown") %>%
      dplyr::filter(ExonicFunc.refGene != "synonymous SNV")
    
  }
  
  if( !is.null(by.gene) ) {
    
    input.df <- dplyr::filter(input.df, Gene.refGene %in% by.gene)
    
  } else {
    
    c.names <- colnames(input.df)
    
    input.df <- input.df %>%
      group_by(Gene.refGene) %>%
      tally(name = "count") %>%
      arrange(desc(count)) %>%
      slice_head(n = n.top) %>%
      inner_join(input.df, by = "Gene.refGene") %>%
      dplyr::select(all_of(c.names))
  }
  
  
  return(input.df)
  
}

## plot VAF for all variants in given table
plotVAF <- function(varTable.df, ad.filter = 10) {
  
  ## just sanity check
  if( !any(grepl("hom.het", colnames(varTable.df))) ) {
    varTable.df <- cleanHomHet(varTable.df)
  }
  
  varTable.df <- varTable.df %>%
    dplyr::filter(
      (hom.het == "het" & AD1 >= eval(ad.filter) & AD2 >= eval(ad.filter)) |  
        (hom.het == "hom" & AD1 + AD2 >= eval(ad.filter) * 2))
  
  ggplot(varTable.df, aes(x = hom.het, y = VAF, fill = phased, color = hom.het)) +
    geom_point(position = position_jitter(width = 0.3), size = 5, alpha = 0.75, shape = 21) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.01)) +
    labs(
      title = "Variant Allele Fraction (VAF) by Genotype Class",
      x = "Genotype Class",
      y = "Variant Allele Fraction (VAF)",
      fill = "Phased" 
    ) +
    facet_grid(. ~ cohort, scales = "free", space = "free", 
               labeller = labeller(cohort = c(
                 "breast" = "Breast",
                 "colon" = "Colon",
                 "leukemia" = "Leukemia",
                 "other" = "Other"
               ))) +
    theme_minimal() +
    theme(
      panel.background = element_blank(),
      plot.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 20),
      strip.background = element_rect(fill = "lightgrey", color = "black"), 
      legend.position = "bottom", 
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.ticks.y = element_blank(),
      strip.text.x = element_text(size = 16),
      panel.border = element_blank()
    ) +
    scale_fill_manual(values = c("TRUE" = "lightblue", "FALSE" = "white")) +
    scale_color_manual(values = c("het" = "red", "hom" = "blue")) +
    scale_x_discrete(labels = c("hom" = "Hom", "het" = "Het")) + guides(color="none") 
  
}


## INIT plt.colors for consistent legend accross plots
## ToDo: find a better way --> in a package this would be a object element.
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
names(sorted_labels) <- NULL
sorted_labels <- unique(sorted_labels)
plt.colors = structure(cols, names = sorted_labels)



library(ComplexHeatmap)

plotGeneOverlapHeatmap <- function(input.df, meta, p.colors = NULL, order.by.name = FALSE, w.label = 35) {
  
  
  ### ToDo: make tiles transparent based on revel score or VAF maybe?!
  meta.labels <- meta %>%
    dplyr::select(sampleID, fID)
  
  plot.df <- input.df %>%
    dplyr::mutate(ExonicFunc.refGene = gsub(" ", "_", ExonicFunc.refGene)) %>%
    group_by(sampleID, Gene.refGene, ExonicFunc.refGene, cohort) %>% 
    tally(name = "count") %>%  # Count occurrences
    pivot_wider( names_from = ExonicFunc.refGene, values_from = count, values_fill = 0) %>%
    pivot_longer(cols = -c(sampleID, Gene.refGene, cohort),names_to = "VariantType",values_to = "Count") %>%
    dplyr::filter(Count > 0) %>%
    mutate(grouping = paste0(Gene.refGene, "_", sampleID)) %>%
    group_by(grouping, Gene.refGene, sampleID, cohort) %>% # 
    summarise(VariantType = paste(unique(VariantType), collapse = ","),Count = paste(Count, collapse = ","),.groups = "drop") %>%
    group_by(sampleID, Gene.refGene, cohort) %>%
    summarise(VariantType = paste(unique(VariantType), collapse = ";"),.groups = "drop") %>%
    pivot_wider(names_from = Gene.refGene,values_from = VariantType,values_fill = NA)
  
  
  # Apply the transformation to all but the first two columns
  plot.df <- plot.df %>%
    mutate(across(
      .cols = -c(sampleID, cohort),  # Exclude the first two columns
      .fns = ~ sapply(., function(x) {
        if (is.na(x)) return(NA_character_)  # Handle NA values
        sorted_elements <- sort(unlist(strsplit(as.character(x), ",")))  # Split, sort, and recombine
        paste(sorted_elements, collapse = ",")
      }),
      .names = "{.col}"  # Add "sorted_" prefix to transformed columns
    ))
  
  # View the transformed dataframe
  
  c.order <- colnames(plot.df)
  
  plot.df <- plot.df %>%
    dplyr::left_join(., meta.labels, by = "sampleID") %>%
    dplyr::select(c.order[1:2], fID, c.order[-c(1,2)]) %>%
    arrange(cohort, fID)
  
  labels <- plot.df %>%
    dplyr::select(-c(sampleID, cohort, fID)) %>%
    pivot_longer(cols = everything(), values_to = "value") %>% # Reshape to a long format
    distinct(value) %>%
    dplyr::filter(!is.na(value)) %>%
    pull(value)
  
  if( is.null(p.colors) ) {
    cols <- pals::tableau20()
    p.colors = structure(cols, names = labels)
  } else {
    print("Debug: else")
    p.colors <- p.colors[which(names(p.colors) %in% labels)]
  }
  # fix labels
  heat_legend_labels <- data.frame("value" = labels) %>%
    dplyr::mutate(value = gsub("_", " ", value)) %>%
    dplyr::mutate(value = gsub(",", ", ", value)) %>%
    dplyr::mutate(value = stringr::str_to_title(value)) %>%
    dplyr::mutate(value = gsub("Snv", "SNV", value)) %>%
    pull(value)
  
  
  heat_legend_labels <- sapply(heat_legend_labels, 
                               function(x) str_wrap(str_replace_all(x, "_", " "),width = w.label))
  
  n.rows <- ifelse(length(heat_legend_labels) > 9, 3, 2)
  
  mat <- t(as.matrix(plot.df))
  
  colnames(mat) <- mat[3,]
  
  ## facet labels
  f.labels <- sapply( unique(mat[2,]), function(x) {
    paste0(toupper(substring(x, 1, 1)), tolower(substring(x, 2)))
  })
  

    if( order.by.name ) {
      ## sort genes alphabetically
      row_order <- order(rownames(mat[-c(1,2,3), ]))
      sorted_mat <- mat[-c(1,2,3), ][row_order, ] 
    } else {
      row_order <- apply(mat[-c(1,2,3), ], 1, function(row) sum(is.na(row)))
      # Order rows based on the NA counts
      sorted_mat <- mat[-c(1,2,3), ][order(row_order), ]
    }

  

  ## sort sample IDs
  col_order <- order(as.numeric(gsub("S", "", colnames(sorted_mat))))
  sorted_mat <- sorted_mat[, col_order]
  
  # Generate the heatmap with sorted rows
  h.map <- ComplexHeatmap::Heatmap(
    sorted_mat, na_col = "lightgrey", column_split = mat[2,], 
    column_title = NULL, row_km = TRUE, col = p.colors,
    row_title = NULL, column_names_rot = 90,
    top_annotation = HeatmapAnnotation(
      foo = anno_block(
        gp = gpar(fill = "lightgrey"),
        labels = eval(f.labels), 
        labels_gp = gpar(col = "black", fontsize = 12)
      )
    ),
    name = "Variant types", 
    rect_gp = gpar(col = "white", lwd = 2), 
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    heatmap_legend_param = list(labels = heat_legend_labels, at = labels, direction = "horizontal", nrow = eval(n.rows))
  )
  
  ## draw(h.map, heatmap_legend_side = "bottom")

  
  return(h.map)
}


plotPathwayBarCohort <- function( p.df, xlabel="Number of Enriched Pathways", ylabel="REACTOME pathways", w.label=35) {
  
#  names(sort(table(p.df$pathway), decreasing = TRUE))[1:20]
  
  ggplot(p.df, aes(y = forcats::fct_inorder(pathway, ordered = TRUE), 
                   fill = Cancer.TypeClean)) +
    geom_bar(width = .85, color = "black") +
    scale_fill_manual(values = plt.cols) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1)) +
    scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "_", " "),
                                                   width = w.label)) +
#    facet_grid(cols = vars(cohort), scales = "free_x", labeller = label_parsed) +
    facet_grid(. ~ cohort, scales = "free", space = "free", 
               labeller = labeller(cohort = c(
                 "breast" = "Breast (n=18)",
                 "leukemia" = "Leukemia (n=18)",
                 "colon" = "Colon (n=11)"
               ))) +
    xlab(eval(xlabel)) + ylab(eval(ylabel)) + 
    this_theme() +
    theme(
      legend.position = "none", 
      axis.text.y = element_text(size = 12),
      axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_blank(),
      strip.text.x = element_text(size = 16),
      panel.border = element_blank()
    ) +
    scale_x_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 2))
  
}




## helper function - for a single sample remove background of choice
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


## just remove the healthy ctrl background
removeCtrlBackgroundList <- function(case.list, ctrl.list) {
  
  tmp.ctrl <- rbind.data.frame(ctrl.list$pathogenic, ctrl.list$potential, ctrl.list$associated, ctrl.list$benign,ctrl.list$likely_benign, ctrl.list$unscored)
  tmp.case <- case.list
  
  for( s.idx in names(case.list) ) {
    tmp.case[[s.idx]] <- removeCtrlBackground( case.list[[s.idx]], tmp.ctrl )
  }
  
  return(tmp.case)
}

## Use specifically matched background for every sample and run through whole cohort
removeCtrlBackgroundMatched <- function(case.list, ctrl.list, df.meta, match.col) {
  
  ## just dirty duplicate for output
  tmp.case <- case.list
  ctrl.list <- rbind.data.frame(ctrl.list$pathogenic, ctrl.list$potential, 
                                ctrl.list$associated, ctrl.list$benign,
                                ctrl.list$likely_benign, ctrl.list$unscored)
  ## find the matched samples
  for( s.idx in names(case.list) ) {
    ## 1. find matched ctrl
    match.idx <- df.meta$S.merge[which(df.meta$Sample == eval(s.idx))]
    match.ids <- df.meta %>% 
      dplyr::filter(get(match.col) == eval(match.idx) & Sample != eval(s.idx)) %>%
      dplyr::select(Sample)
    ## and subset the background so that all the called variants need to have ALT
    tmp.ctrl <- ctrl.list %>%
      rowwise() %>%
      filter(all(grepl("(^[^:]*1)|(^\\./\\.)", c_across(match.ids$Sample)))) %>%
      ungroup()
    
    tmp.case[[s.idx]] <- removeCtrlBackground( case.list[[s.idx]], tmp.ctrl )
    ## le progress
    print(paste0("Finished: ",s.idx))
  }
  
  return(tmp.case)
}





