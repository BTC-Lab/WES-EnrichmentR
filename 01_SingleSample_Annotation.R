# INIT --------------------------------------------------------------------

source("scripts/_network_diffusion.R")
source("scripts/_annovar_filtering.R")
set.seed(16341)

tbl.mutations <- readRDS(file = file.path("data", "RDS", "VRI_mutations_curated.RDS"))


# BiocManager::install("maftools")
# BiocManager::install("MatchIt")
# BiocManager::install("VIM")
# BiocManager::install("tximport")
# BiocManager::install("pachterlab/sleuth")
# BiocManager::install("devtools")
library(readr)
library(tidyverse)
data.dir <- file.path(getwd(),"data")
library(pals)
library("MatchIt")
# Load required library
library(VIM)
library(plyr)

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

library("org.Hs.eg.db")
library(tximport)
library(AnnotationDbi)
library(sleuth)

cols <- pals::tableau20(20)
scales::show_col(cols)

tx2gene <- function(transcriptVersion=TRUE) {
  require(biomaRt)
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  tx2gene <- biomaRt::getBM(attributes = c("ensembl_transcript_id", 
                                           "transcript_version", "ensembl_gene_id","external_gene_name"), mart = mart)
  if(transcriptVersion) {
    tx2gene <- dplyr::mutate(tx2gene, target_id = paste(ensembl_transcript_id, 
                                                        transcript_version, sep = ".")) # include version numbering 
    tx2gene <- dplyr::select(tx2gene, target_id, ens_gene = ensembl_gene_id,
                             ext_gene=external_gene_name)
  } else {
    tx2gene <- dplyr::rename(tx2gene, target_id = ensembl_transcript_id,
                             ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
    
  }
  
  return(tx2gene)
}

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

filterEtable <- function( input.df ) {
  df.ret <- input.df %>%
    dplyr::filter(., !grepl("Expression", SNP)) %>%
    dplyr::mutate("GeneId" = stringr::str_trim(Gene)) %>%
    dplyr::filter(., !grepl("hsa", GeneId)) %>%
    dplyr::filter(., !grepl("miR", GeneId)) %>%
    tidyr::separate_rows(GeneId, sep = "/") %>%
    dplyr::mutate("GeneId" = gsub("-", "",GeneId)) %>%
    dplyr::mutate("GeneId" = gsub("ß", "B",GeneId)) %>%
    dplyr::mutate("GeneId" = gsub("α", "A",GeneId)) %>%
    dplyr::mutate("GeneId" = gsub("κ", "K",GeneId)) %>%
    dplyr::mutate("GeneId" = str_to_upper(GeneId)) %>%
    distinct(GeneId) 
  
  return(df.ret)
  
}


colfunc <- colorRampPalette(c("#2CA02C", "#D62728"))
colfunc(10)

cols <- pals::tableau20(20)
p.cols <- c(cols[1],cols[3],cols[5],cols[9],cols[11],cols[13],cols[17],cols[19],cols[6])


# FUNCTIONS ---------------------------------------------------------------

prepATable <- function( input.df, f.freq=0.01, assoc.col=NULL ) {
  
  ## replace missing exac frequencies with zero
  input.df <- input.df %>%
    mutate(., varID = paste(Chr, Start, End, Ref, Alt, sep="_")) %>%
    mutate(., Location = paste(Chr, Start, sep=":")) %>%
    separate_rows(Gene.refGene, sep = ";") %>%
    mutate_at(vars(starts_with("gnomad41")), ~ replace_na(., -1)) %>%
    dplyr::filter(., gnomad41_genome_AF_afr < f.freq & gnomad41_genome_AF_ami < f.freq & gnomad41_genome_AF_amr < f.freq &
                    gnomad41_genome_AF_asj < f.freq & gnomad41_genome_AF_eas < f.freq & gnomad41_genome_AF_fin < f.freq &
                    gnomad41_genome_AF_mid < f.freq & gnomad41_genome_AF_nfe < f.freq & gnomad41_genome_AF_sas < f.freq)
  
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
  
  ## GET REVEL SCORED VARIANTS
  df.rev <- input.df %>%
    #dplyr::filter(., !is.na(REVEL_score)) %>%
    dplyr::filter(., as.numeric(REVEL_score) >= 0.75)
  
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
  out[["revel"]] <- df.rev
  
  return(out)
}

# SINGLE-SAMPLE SUBSETTING ------------------------------------------------

### We take the output from Cohort Subsetting and now just get the single samples
## txt should already exist, but for completeness sake
meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS"))

## remove overhead of unspecified Cancers
meta <- meta %>%
  dplyr::filter(!(Cancer == "yes" & is.na(Cancer.TypeClean))) %>%   # remove the unspecified types
  dplyr::filter(!(Cancer.TypeClean == "sicklecell") | Cancer == "no") %>%  # remove the not cancers
  dplyr::filter(!(Cancer.TypeClean == "acuteidiopatheicthrompocytopeniapurpura") | Cancer == "no")%>%
  dplyr::mutate(Cancer.TypeClean = str_to_title(Cancer.TypeClean)) %>%    # make nice labels & and adjust order with factors
  dplyr::mutate(Cancer.TypeClean = 
                  factor(Cancer.TypeClean, levels = c("Breast","Leukemia","Colon",
                                                      "Lung","Prostate","Uterine","Pancreatic","Thyroid",
                                                      "Renalcell","Sarcoma","Multiplemyeloma","Lymphoma","Bladder"))) 

saveRDS(meta, file.path("data/RDS", "VRI_Pilot_meta_v3.RDS"), compress = TRUE)

## make maintype_extract.txt for cancer to use with bcftools
table(meta$Cancer.TypeClean)
sID.allCancer <- meta[which(meta$Cancer == "yes"),]

data.table::fwrite(list(sID.allCancer$Sample), file = file.path("data", "RAW", "samples_all_cancer.txt"))

## subset main file again
# while read -r sample; do
#   echo "Processing sample: $sample"
#   bcftools view --min-ac=1 --threads 14 --exclude-uncalled -c1 -s "$sample" -Oz -o "${sample}.vcf.gz" VRI_truseq_rsID_norm_study.VAF.vcf.gz
# done < samples_all_cancer.txt

## filter for DP 10
# ls *.vcf.gz | awk -F'\.' '{print "bcftools filter -e \"FORMAT/DP < 10\"  "$0" --output-type z -o "$1".DP10.vcf.gz"}' > run_DP10_filter.sh


### Just treat the sum of matched samples as background because otherwise it just gets messy


# SINGLE-SAMPLE SUBSETTING (CONTROLS) -------------------------------------

meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS"))

## remove overhead of unspecified Cancers
c.meta <- meta %>%
  dplyr::filter(!(Cancer == "yes" & is.na(Cancer.TypeClean))) %>%   # remove the unspecified types
  dplyr::filter(!(Cancer.TypeClean == "sicklecell")) %>%  # remove the not cancers
  dplyr::filter(!(Cancer.TypeClean == "acuteidiopatheicthrompocytopeniapurpura") | Cancer == "no")%>%
  dplyr::mutate(Cancer.TypeClean = str_to_title(Cancer.TypeClean)) %>%    # make nice labels & and adjust order with factors
  dplyr::mutate(Cancer.TypeClean = 
                  factor(Cancer.TypeClean, levels = c("Breast","Leukemia","Colon",
                                                      "Lung","Prostate","Uterine","Pancreatic","Thyroid",
                                                      "Renalcell","Sarcoma","Multiplemyeloma","Lymphoma","Bladder"))) 

match.col="m_breast_4"

c.breast <- c.meta %>%
  dplyr::filter(Cancer == "yes") %>%
  dplyr::filter(Cancer.TypeClean == "Breast") %>%
  pull(Sample)

b.matches <- NULL
### Find my matched controls
for( s.idx in c.breast ) {
  ## 1. find matched ctrl
  match.idx <- meta$S.merge[which(meta$Sample == eval(s.idx))]
  match.ids <- meta %>% 
    dplyr::filter(get(match.col) == eval(match.idx) & Sample != eval(s.idx)) %>%
    pull(Sample)

  b.matches <- c(b.matches, match.ids)
}

## colon
match.col="m_colon_4"

c.colon <- c.meta %>%
  dplyr::filter(Cancer == "yes") %>%
  dplyr::filter(Cancer.TypeClean == "Colon") %>%
  pull(Sample)

c.matches <- NULL
### Find my matched controls
for( s.idx in c.colon ) {
  ## 1. find matched ctrl
  match.idx <- meta$S.merge[which(meta$Sample == eval(s.idx))]
  match.ids <- meta %>% 
    dplyr::filter(get(match.col) == eval(match.idx) & Sample != eval(s.idx)) %>%
    pull(Sample)
  
  c.matches <- c(c.matches, match.ids)
}

## leukemia
match.col="m_leukemia_4"

c.leukemia <- c.meta %>%
  dplyr::filter(Cancer == "yes") %>%
  dplyr::filter(Cancer.TypeClean == "Leukemia") %>%
  pull(Sample)

l.matches <- NULL
### Find my matched controls
for( s.idx in c.leukemia ) {
  ## 1. find matched ctrl
  match.idx <- meta$S.merge[which(meta$Sample == eval(s.idx))]
  match.ids <- meta %>% 
    dplyr::filter(get(match.col) == eval(match.idx) & Sample != eval(s.idx)) %>%
    pull(Sample)
  
  l.matches <- c(l.matches, match.ids)
}

out <- data.frame("V1"=unique(c(b.matches, c.matches, l.matches)))

data.table::fwrite(list(out$V1), file = file.path("data", "RAW", "samples_all_controls.txt"))

## maybe quicker to just try subsetting with the full data --> since I did no ac filter this should be fine
df.cohort <- readRDS(file.path("data/RDS", "VRI_Pilot_Cohort_annotated.filtered.26122024.RDS"))
ctrl.list <- rbind.data.frame(df.cohort$breast_ctrl$pathogenic, df.cohort$breast_ctrl$potential, 
                              df.cohort$breast_ctrl$associated)
b.ctrls <- list()
for( c.idx in b.matches ) {
  tmp.ctrl <- ctrl.list %>%
   rowwise() %>%
    filter(all(grepl("(^[^:]*1)|(^\\./\\.)", c_across(c.idx)))) %>%
   ungroup() %>%
  dplyr::mutate(Variant = get(c.idx)) %>%
    dplyr::select(Chr:FORMAT, Variant, varID, Location)
  
  b.ctrls[[c.idx]] <- tmp.ctrl
  print(c.idx)
}

ctrl.list <- rbind.data.frame(df.cohort$colon_ctrl$pathogenic, df.cohort$colon_ctrl$potential, 
                              df.cohort$colon_ctrl$associated)
c.ctrls <- list()
for( c.idx in c.matches ) {
  tmp.ctrl <- ctrl.list %>%
    rowwise() %>%
    filter(all(grepl("(^[^:]*1)|(^\\./\\.)", c_across(c.idx)))) %>%
    ungroup() %>%
    dplyr::mutate(Variant = get(c.idx)) %>%
    dplyr::select(Chr:FORMAT, Variant, varID, Location)
  
  c.ctrls[[c.idx]] <- tmp.ctrl
  print(c.idx)
  
}

ctrl.list <- rbind.data.frame(df.cohort$leukemia_ctrl$pathogenic, df.cohort$leukemia_ctrl$potential, 
                              df.cohort$leukemia_ctrl$associated)
l.ctrls <- list()
for( c.idx in l.matches ) {
  tmp.ctrl <- ctrl.list %>%
    rowwise() %>%
    filter(all(grepl("(^[^:]*1)|(^\\./\\.)", c_across(c.idx)))) %>%
    ungroup() %>%
    dplyr::mutate(Variant = get(c.idx)) %>%
    dplyr::select(Chr:FORMAT, Variant, varID, Location)
  
  l.ctrls[[c.idx]] <- tmp.ctrl
  print(c.idx)
  
}

out.list <- list()
out.list[["breast"]] <- b.ctrls
out.list[["colon"]] <- c.ctrls
out.list[["leukemia"]] <- l.ctrls

saveRDS(out.list, file.path("data/RDS", "VRI_Pilot_Control_Cohort_SS_annotated_301224.RDS"))



## subset main file again
while read -r sample; do
  echo "Processing sample: $sample"
  bcftools view --min-ac=1 --threads 14 --exclude-uncalled -c1 -s "$sample" -Oz -o "${sample}.vcf.gz" VRI_truseq_rsID_norm_study.VAF.vcf.gz
done < samples_all_controls.txt

## filter for DP 10
# ls *.vcf.gz | awk -F'\.' '{print "bcftools filter -e \"FORMAT/DP < 10\"  "$0" --output-type z -o "$1".DP10.vcf.gz"}' > run_DP10_filter.sh

# SINGLE-SAMPLE ANNOTATION ------------------------------------------------


## then run current Annovar annotations
## ls *DP10.vcf.gz | awk -F'\.' '{print "/Users/rmolbrich/Workdir/Tools/annovar/table_annovar.pl "$0" /Volumes/ERGPBackup/annovar/humandb -buildver hg38 -out "$1" -remove -protocol refGene,cytoBand,exac03,avsnp151,dbnsfp47a,clinvar_20240611,gnomad41_genome -operation g,r,f,f,f,f,f -nastring . -polish -vcfinput"}' > run_annotation.sh




# SINGLE-SAMPLE ANNOTATIONS IMPORT ----------------------------------------

## most annoying step is renaming of the columns because annovar doesn't keep track

## chr7_152152830_152152830_G_A
# meta[which(meta$fID == "S66"),] # DXB_HLA_21_4084900_00959
# 
# f.idx <- "DXB_HLA_21_4084900_00959.hg38_multianno.txt"
# 
# dummy <- df.annot.tmp[which(df.annot.tmp$Gene.refGene == "KMT2C"),]


allFiles <- list.files(file.path("data/RAW", "annotation_single_samples"), pattern = "hg38_multianno.txt")

save.list <- list()

for( f.idx in allFiles ) {
  
  sID <- gsub(".hg38_multianno.txt", "", f.idx)
  
  df.annot.tmp <- data.table::fread(file.path("data/RAW", "annotation_single_samples", eval(f.idx)),
                                    header = TRUE, na.strings = ".")
  
  df.annot.tmp <- df.annot.tmp %>%
    dplyr::rename_with(~ sID, .col = Otherinfo13) %>%
    dplyr::rename_with(~ c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"), .col = Otherinfo4:Otherinfo12)
  
  ## put into the list
  save.list[[eval(sID)]] <- df.annot.tmp
}

saveRDS(save.list, file.path("data/RDS", "VRI_Pilot_SingleSample_annotated.26122024.RDS"), compress = TRUE)



# SINGLE-SAMPLE ANNOTATIONS FILTER ----------------------------------------


meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS"))

save.list <- readRDS(file.path("data", "RDS", "VRI_Pilot_SingleSample_annotated.26122024.RDS"))
tbl.mutations <- readRDS(file = file.path("data", "RDS", "VRI_mutations_curated.RDS"))

sample.list <- list()
breast.list <- list()
colon.list <- list()
leukemia.list <- list()
other.list <- list()

f.cutoff <- 0.001

for( leSample in names(save.list) ) {
  
  c.name <- meta$Cancer.TypeClean[which(meta$Sample == eval(leSample))]
  df <- save.list[[eval(leSample)]]
  ## filter for exac frequencies smaller than 0.001
  if( c.name == "breast" ) {
    ass.table <- tbl.mutations$breast$Gene
    
    breast.list[[eval(leSample)]] <- prepATable(df, f.cutoff, ass.table)
  } else if( c.name == "colon" ) {
    ass.table <- tbl.mutations$colon$Gene
    
    colon.list[[eval(leSample)]] <- prepATable(df, f.cutoff, ass.table)
    
  } else if( c.name == "leukemia" ) {
    ass.table <- tbl.mutations$leukemia$Gene
    
    leukemia.list[[eval(leSample)]] <- prepATable(df, f.cutoff, ass.table)
    
  } else {
    ass.table = NULL
    other.list[[eval(leSample)]] <- prepATable(df, f.cutoff, ass.table)
    
  }
  
}

sample.list <- list()
sample.list[["breast"]] <- breast.list
sample.list[["colon"]] <- colon.list
sample.list[["leukemia"]] <- leukemia.list
sample.list[["other"]] <- other.list

saveRDS(sample.list, file = file.path("data", "RDS", "VRI_Pilot_SingleSample_annotated.filtered.26122024.RDS"), compress = TRUE)

## sanity
length(sample.list$breast); length(sample.list$colon); length(sample.list$leukemia); length(sample.list$other)



# SINGLE-SAMPLE - VARIANT TO GENE OVERLAP ---------------------------------


meta <- readRDS(file.path("data", "RDS", "VRI_Pilot_meta_v3.RDS"))
colnames(meta)[1] <- "sampleID"

## Take all remaining pathogenic variants and put the in a single df
sample.list <- readRDS(file.path("data", "RDS", "VRI_Pilot_SingleSample_annotated.filtered.26122024.RDS"))

out.pat.list <- list()
out.pot.list <- list()
out.ass.list <- list()
out.rev.list <- list()
for( leCohort in names(sample.list) ) {
  ## for the cohort get the list with all the samples
  cohort.list <- sample.list[[eval(leCohort)]]
  
  df.pat.out <- NULL
  df.pot.out <- NULL
  df.ass.out <- NULL
  df.rev.out <- NULL
  ## process every sample
  for( leSample in names(cohort.list) ) {
    
    df.list <- cohort.list[[eval(leSample)]]
    # pathogenic
    df.pat.tmp <- df.list$pathogenic %>%
      dplyr::mutate(sampleID = leSample,
                    cohort = eval(leCohort),
                    class.type = "pathogenic",.before = Chr) %>%
      rename(Variant = eval(leSample)) %>%
      mutate( GT = sub(":.*", "", Variant),
              VAF = sub(".*:", "", Variant), .before = cytoBand)
      
    df.pat.out <- rbind(df.pat.out, df.pat.tmp)
    # potential
    df.pot.tmp <- df.list$potential %>%
      dplyr::mutate(sampleID = leSample,
                    cohort = eval(leCohort),
                    class.type = "potential",.before = Chr) %>%
      rename(Variant = eval(leSample)) %>%
      mutate( GT = sub(":.*", "", Variant),
              VAF = sub(".*:", "", Variant), .before = cytoBand)
    
    df.pot.out <- rbind(df.pot.out, df.pot.tmp)
    # revel75
    df.ass.tmp <- df.list$associated %>%
      dplyr::mutate(sampleID = leSample,
                    cohort = eval(leCohort),
                    class.type = "associated", .before = Chr) %>%
      rename(Variant = eval(leSample)) %>%
      mutate( GT = sub(":.*", "", Variant),
              VAF = sub(".*:", "", Variant), .before = cytoBand)
    
    df.ass.out <- rbind(df.ass.out, df.ass.tmp)
    # revel75
    df.rev.tmp <- df.list$revel %>%
      dplyr::mutate(sampleID = leSample,
                    cohort = eval(leCohort),
                    class.type = "revel",.before = Chr) %>%
      rename(Variant = eval(leSample)) %>%
      mutate( GT = sub(":.*", "", Variant),
              VAF = sub(".*:", "", Variant), .before = cytoBand)
    
    df.rev.out <- rbind(df.rev.out, df.rev.tmp)
  }
  
  
  
  
  
  out.pat.list[[eval(leCohort)]] <- df.pat.out
  out.pot.list[[eval(leCohort)]] <- df.pot.out
  out.ass.list[[eval(leCohort)]] <- df.ass.out
  out.rev.list[[eval(leCohort)]] <- df.rev.out
}

## Sanity 
names(out.pat.list)
# merge all
out.pat.list[["full"]] <- rbind(out.pat.list[["breast"]], out.pat.list[["colon"]], out.pat.list[["leukemia"]], out.pat.list[["other"]])
out.pot.list[["full"]] <- rbind(out.pot.list[["breast"]], out.pot.list[["colon"]], out.pot.list[["leukemia"]], out.pot.list[["other"]])
out.ass.list[["full"]] <- rbind(out.ass.list[["breast"]], out.ass.list[["colon"]], out.ass.list[["leukemia"]], out.ass.list[["other"]])
out.rev.list[["full"]] <- rbind(out.rev.list[["breast"]], out.rev.list[["colon"]], out.rev.list[["leukemia"]], out.rev.list[["other"]])


out.pat.list[["full"]] <- cleanHomHet(out.pat.list[["full"]], by.vaf = TRUE, by.dp = 10)
out.pot.list[["full"]] <- cleanHomHet(out.pot.list[["full"]], by.vaf = TRUE, by.dp = 10)
out.ass.list[["full"]] <- cleanHomHet(out.ass.list[["full"]], by.vaf = TRUE, by.dp = 10)
out.rev.list[["full"]] <- cleanHomHet(out.rev.list[["full"]], by.vaf = TRUE, by.dp = 10)

df <- as.data.frame(table(paste0(out.pat.list[["full"]]$ExonicFunc.refGene, "_", out.pat.list[["full"]]$Func.refGene))) %>%
  mutate(Percentage = (Freq / sum(Freq)) * 100)

table(out.pat.list[["full"]]$Func.refGene)

table(out.pat.list[["full"]]$Gene.refGene)

length(unique(out.pat.list[["full"]]$Gene.refGene))
length(unique(out.pat.list[["full"]]$sampleID))

### plot VAF for the subsets
plt.pathogenic.vaf <- plotVAF(out.pat.list[["full"]])
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Cohort_Pathogenic_Variants.png"), plot = plt.pathogenic.vaf, width = 13, height = 4.5, dpi = "retina")
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Cohort_Pathogenic_Variants.pdf"), plot = plt.pathogenic.vaf, width = 13, height = 4.5, dpi = "retina")

plt.potential.vaf <- plotVAF(out.pot.list[["full"]])
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Cohort_Potential_Variants.png"), plot = plt.potential.vaf, width = 13, height = 4.5, dpi = "retina")
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Cohort_Potential_Variants.pdf"), plot = plt.potential.vaf, width = 13, height = 4.5, dpi = "retina")

plt.associated.vaf <- plotVAF(out.ass.list[["full"]])
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Cohort_Associated_Variants.png"), plot = plt.associated.vaf, width = 13, height = 4.5, dpi = "retina")
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Cohort_Associated_Variants.pdf"), plot = plt.associated.vaf, width = 13, height = 4.5, dpi = "retina")

plt.rev.vaf <- plotVAF(out.rev.list[["full"]])
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Cohort_Revel75_Variants.png"), plot = plt.rev.vaf, width = 13, height = 4.5, dpi = "retina")
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Cohort_Revel75_Variants.pdf"), plot = plt.rev.vaf, width = 13, height = 4.5, dpi = "retina")

openxlsx::write.xlsx(out.pat.list, file.path("data/output", "VRI_Cancers_Pathogenic_Variants.xlsx"))
openxlsx::write.xlsx(out.pot.list, file.path("data/output", "VRI_Cancers_Potential_Variants.xlsx"))
openxlsx::write.xlsx(out.pot.list, file.path("data/output", "VRI_Cancers_Associated_Variants.xlsx"))
openxlsx::write.xlsx(out.rev.list, file.path("data/output", "VRI_Cancers_Revel75_Variants.xlsx"))


save.it <- list()
save.it[["pathogenic"]] <- out.pat.list
save.it[["potential"]] <- out.pot.list
save.it[["associated"]] <- out.ass.list
save.it[["revel"]] <- out.rev.list

saveRDS(save.it, file.path("data/RDS", "VRI_Pilot_SingleSample_Variant_Gene_Overlap.28122024.RDS"))


## check legacy
df.leg <- readRDS(file.path("data/RDS", "breast_cancer_variants_old.RDS"))




# ANALYSIS: SINGLE-SAMPLE - VARIANT TO GENE OVERLAP -----------------------

## Now we have pathogenic, potential (and subset of it REvel), and associated 
## variants and need to find an interesting set of genes with variants in them
## that overlap between the cohorts and the whole dataset

save.it <- readRDS(file.path("data/RDS", "VRI_Pilot_SingleSample_Variant_Gene_Overlap.28122024.RDS"))

# Play around and find good set
df.pat <- cleanHomHet(save.it$pathogenic$full, by.vaf = TRUE, by.dp = 100)
table(df.pat$Gene.refGene)[order(table(df.pat$Gene.refGene), decreasing = TRUE)]
plotVAF(df.pat)

df.pot <- cleanHomHet(save.it$potential$full, by.vaf = TRUE, by.dp = 100)
table(df.pot$Gene.refGene)[order(table(df.pot$Gene.refGene), decreasing = TRUE)]
plotVAF(df.pot)

df.ass <- cleanHomHet(save.it$associated$full, by.vaf = TRUE, by.dp = 100)
table(df.ass$Gene.refGene)[order(table(df.ass$Gene.refGene), decreasing = TRUE)]
plotVAF(df.ass)

df.rev <- cleanHomHet(save.it$revel$full, by.vaf = TRUE, by.dp = 100)
table(df.rev$Gene.refGene)[order(table(df.rev$Gene.refGene), decreasing = TRUE)] >= 20
plotVAF(df.rev)


## merge all of it
df.all <- rbind(save.it$pathogenic$full, save.it$potential$full, save.it$associated$full)
df.all <- cleanHomHet(df.all, by.vaf = TRUE, by.dp = 100)
table(df.all$Gene.refGene)[order(table(df.all$Gene.refGene), decreasing = TRUE)]
plotVAF(df.all)


df.all.rev <- rbind(save.it$pathogenic$full, save.it$associated$full, save.it$revel$full)
df.all.rev <- cleanHomHet(df.all.rev, by.vaf = TRUE, by.dp = 100)
table(df.all.rev$Gene.refGene)[order(table(df.all.rev$Gene.refGene), decreasing = TRUE)]
plotVAF(df.all.rev)
table(df.all.rev$ExonicFunc.refGene)

df.all.rev[which(df.all.rev$ExonicFunc.refGene == "stopgain"),]

# CTBP2,MUC4,KIR3DS1,ANKRD36C,ANKRD36,TEKT4,LOC100129697,HLA-DRB1,KMT2C,MUC16,
# DNMT1,NBPF1,KIR2DS3,KIR2DS4,MAP3K1,HLA-DQA2,KIR2DS1,KIR2DS2,KIR2DS5,LOC102725023,
# LOC112268355,IGFN1,CPS1,ZNF717

## get all labels and play through plotGeneOverlap to make input vector
#input.df <- getGeneOverlap(df.all.rev, n.top = 1000, clean.variants = TRUE)

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

gene.select <- c("CTBP2", "MUC4","KMT2C","MUC16","MAP3K1","ZNF717","COL18A1",
                 "ATG3","PAX6","LFNG","PMS1","RBMX","BCR","HRC","PSORS1C1",
                 "CDCP2","CHST15","TMEM254","OR6C76")

df.all.rev <- rbind(save.it$pathogenic$full, save.it$associated$full, save.it$revel$full)
df.all.rev <- cleanHomHet(df.all.rev, by.vaf = TRUE, by.dp = 10)
test <- getGeneOverlap(df.all.rev, by.gene = gene.select, clean.variants = TRUE)

plt.vaf <- plotVAF(test)
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Study_Gene_Overlap_selected.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")

openxlsx::write.xlsx(test, file.path("data/output", "Table_XX_VAF_Study_Gene_Overlap_selected.xlsx"))

plt.heat_study <- plotGeneOverlapHeatmap(test, meta, plt.colors)

pdf(file.path("data/output", "Figure_4_ALT_Study_Gene_Overlap_selected.pdf"), width = 12, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()


test2 <- getGeneOverlap(df.all.rev, n.top = 10, clean.variants = TRUE)
unique(test2$Gene.refGene)
plt.heat_study <- plotGeneOverlapHeatmap(test2, meta, plt.colors)


plt.vaf <- plotVAF(test2)
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Study_Gene_Overlap_top10.pdf"), plot = plt.vaf, width = 15, height = 7, dpi = "retina")

openxlsx::write.xlsx(test2, file.path("data/output", "Table_XX_VAF_Study_Gene_Overlap_top10.xlsx"))


pdf(file.path("data/output", "Figure_4_Study_Gene_Overlap_top10.pdf"), width = 13, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()


#### full cohort overlap stats
df <- as.data.frame(table(paste0(test$ExonicFunc.refGene, "_", test$Func.refGene))) %>%
  mutate(Percentage = (Freq / sum(Freq)) * 100)

table(test$Func.refGene)

table(test$Gene.refGene)

length(unique(test$Gene.refGene))
length(unique(test$sampleID))


### get the ones that are shared between everything
hmm <- as.data.frame(table(test$varID, test$cohort)) %>%
  dplyr::filter(Freq > 10)

unique(hmm$Var1)

muh <- dplyr::filter(test, varID %in% unique(hmm$Var1))

plt.heat_study <- plotGeneOverlapHeatmap(muh, meta, plt.colors)

pdf(file.path("data/output", "Supplementary_Figure_XX_Study_Gene_Overlap_sharedVariants.pdf"), width = 13, height = 5.5)
draw(plt.heat_study, heatmap_legend_side = "bottom")
dev.off()
openxlsx::write.xlsx(muh, file.path("data/output", "Table_XX_Study_Gene_Overlap_sharedVariants.xlsx"))


a.genes <- unique(c(unique(tbl.mutations$colon$Gene), unique(tbl.mutations$breast$Gene),unique(tbl.mutations$leukemia$Gene)))


## breast only
df.all.breast <- rbind(save.it$pathogenic$breast, save.it$revel$breast, save.it$associated$breast)
df.all.breast <- cleanHomHet(df.all.breast, by.vaf = TRUE, by.dp = 10)
table(df.all.breast$Gene.refGene)[order(table(df.all.breast$Gene.refGene), decreasing = TRUE)]


df.all.breast <- getGeneOverlap(df.all.breast, n.top = 10, clean.variants = TRUE)
plt.vaf.breast <- plotVAF(df.all.breast)
plt.heat_breast <- plotGeneOverlapHeatmap(df.all.breast, meta, plt.colors)

pdf(file.path("data/output", "Figure_2_Breast_Gene_Overlap_top10.pdf"), width = 12, height = 4.5)
draw(plt.heat_breast, heatmap_legend_side = "bottom")
dev.off()

## colon only
df.all.colon <- rbind(save.it$pathogenic$colon, save.it$revel$colon, save.it$associated$colon)
df.all.colon <- cleanHomHet(df.all.colon, by.vaf = TRUE, by.dp = 10)
table(df.all.colon$Gene.refGene)[order(table(df.all.colon$Gene.refGene), decreasing = TRUE)]


df.all.colon <- getGeneOverlap(df.all.colon, n.top = 10, clean.variants = TRUE)
plt.vaf.colon <- plotVAF(df.all.colon)
plt.heat_colon <- plotGeneOverlapHeatmap(df.all.colon, meta, plt.colors)

pdf(file.path("data/output", "Figure_2_Colon_Gene_Overlap_top10.pdf"), width = 12, height = 4.5)
draw(plt.heat_colon, heatmap_legend_side = "bottom")
dev.off()

## leukemia only
df.all.leukemia <- rbind(save.it$pathogenic$leukemia, save.it$revel$leukemia, save.it$associated$leukemia)
df.all.leukemia <- cleanHomHet(df.all.leukemia, by.vaf = TRUE, by.dp = 10)
table(df.all.leukemia$Gene.refGene)[order(table(df.all.leukemia$Gene.refGene), decreasing = TRUE)]


df.all.leukemia <- getGeneOverlap(df.all.leukemia, n.top = 10, clean.variants = TRUE)
plt.vaf.leukemia <- plotVAF(df.all.leukemia)
plt.heat_leukemia <- plotGeneOverlapHeatmap(df.all.leukemia, meta, plt.colors)

pdf(file.path("data/output", "Figure_2_Leukemia_Gene_Overlap_top10.pdf"), width = 12, height = 4.5)
draw(plt.heat_leukemia, heatmap_legend_side = "bottom")
dev.off()



g1 <- grid.grabExpr(draw(plt.heat_breast, heatmap_legend_side = "bottom"))
g2 <- grid.grabExpr(draw(plt.heat_colon, heatmap_legend_side = "bottom"))
g3 <- grid.grabExpr(draw(plt.heat_leukemia, heatmap_legend_side = "bottom"))

# Arrange the heatmaps in a panel
grid.arrange(
  g1, g2, g3, 
  ncol = 3  # Arrange in a row with 3 columns
)

### VAF Panel


# Create the arrangement with a shared legend
plt.vaf.cohort <- ggarrange(
  plt.vaf.breast + 
    theme(legend.position = "none", plot.title = element_blank()), 
  plt.vaf.colon + 
    theme(
      legend.position = "none", 
      plot.title = element_blank(), 
      axis.title.y = element_blank(),  # Remove y-axis title
      axis.text.y = element_blank(),   # Remove y-axis ticks
      axis.ticks.y = element_blank()   # Remove y-axis ticks
    ), 
  plt.vaf.leukemia + 
    theme(
      legend.position = "none", 
      plot.title = element_blank(), 
      axis.title.y = element_blank(),  # Remove y-axis title
      axis.text.y = element_blank(),   # Remove y-axis ticks
      axis.ticks.y = element_blank()   # Remove y-axis ticks
    ), 
  ncol = 3, 
  nrow = 1, 
  labels = c("A)", "B)", "C)"), 
  font.label = list(size = 24, color = "black", face = "bold"),
  common.legend = TRUE,  # Add a shared legend
  legend = "bottom"      # Position the legend at the bottom
)


# Display the combined plot
plt.vaf.cohort

ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Cohort_Gene_Overlap_top10.png"), plot = plt.vaf.cohort, width = 13, height = 4.5, dpi = "retina")
ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Cohort_Gene_Overlap_top10.pdf"), plot = plt.vaf.cohort, width = 13, height = 4.5, dpi = "retina")


openxlsx::write.xlsx(list("Breast"=df.all.breast,
                          "Colon"=df.all.colon,
                          "Leukemia"=df.all.leukemia), file.path("data/output", "Table_XX_VAF_Cohort_Gene_Overlap_top10.xlsx"))



# OLD STUFF ---------------------------------------------------------------


ggsave(file.path("data/output", "Supplementary_Figure_XX_VAF_Breast_Gene_Overlap_top10.pdf"), plot = plt.vaf.breast, width = 15, height = 7, dpi = "retina")



Fig3.plt <- ggarrange(
  # First row with 90/10 width ratio
  ggarrange(
    ad_plot.NOsa, ad_plot.sa, 
    ncol = 2, nrow = 1, 
    widths = c(95, 5), 
    labels = c("A)"), 
    font.label = list(size = 24, color = "black", face = "bold")
  ), 
  # Second row
  ggarrange(
    pca12.main, pca12.uae, ancestry_plot, 
    ncol = 3, nrow = 1, 
    labels = c("B)", "C)", "D)"), 
    font.label = list(size = 24, color = "black", face = "bold")
  ), 
  # Third row
  ggarrange(
    mt.plot, violin.plot, 
    ncol = 2, nrow = 1, 
    widths = c(1, 2), 
    labels = c("E)", "F)"), 
    font.label = list(size = 24, color = "black", face = "bold")
  ), 
  # Combine rows
  ncol = 1, nrow = 3, 
  #  labels = c("A)", "B)", "C)"), 
  font.label = list(size = 24, color = "black", face = "bold",
                    panel.spacing = unit(0.1, "cm") )
)


ggsave(file.path("data/output", "Fig3_Panel.png"), plot = Fig3.plt, width = 13, height = 12, dpi = "retina")
ggsave(file.path("data/output", "Fig3_Panel.pdf"), plot = Fig3.plt, width = 13, height = 12, dpi = "retina")


table(df.ass$ExonicFunc.refGene)



table(test2$ExonicFunc.refGene)

table(is.na(test2$ExonicFunc.refGene))


input.df <- test2










### general idea is to get the top n or some count cutoff and then see what sticks
## ideally we get the same with revel and potential --> so we can tell the reviewer 
## that we used this for selection or at least that we got the same results .. maybe?!

## cohort only
df.all.cohort <- rbind(save.it$pathogenic$leukemia, save.it$potential$leukemia, save.it$associated$leukemia)
df.all.cohort <- cleanHomHet(df.all.cohort, by.vaf = TRUE, by.dp = 100)
table(df.all.cohort$Gene.refGene)[order(table(df.all.cohort$Gene.refGene), decreasing = TRUE)]
plotVAF(df.all.cohort)




table(df.ass$ExonicFunc.refGene)






# CTBP2      KMT2C       PAK2     SCN10A


## later to get distinct variants
#  distinct(sampleID, varID, .keep_all = TRUE)

# LEGACY: SINGLE-SAMPLE - PATHOGENIC VARIANTS -----------------------------

### essentially Gnomad41 removed all the previous interesting variants with new reference frequencies
## so run without frequency filtering
meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS"))

save.list <- readRDS(file.path("data", "RDS", "VRI_Pilot_SingleSample_annotated.26122024.RDS"))
tbl.mutations <- readRDS(file = file.path("data", "RDS", "VRI_mutations_curated.RDS"))

sample.list <- list()
breast.list <- list()
colon.list <- list()
leukemia.list <- list()
other.list <- list()

for( leSample in names(save.list) ) {
  
  c.name <- meta$Cancer.TypeClean[which(meta$Sample == eval(leSample))]
  df <- save.list[[eval(leSample)]]
  ## filter for exac frequencies smaller than 0.001
  if( c.name == "breast" ) {
    ass.table <- tbl.mutations$breast$Gene
    
    breast.list[[eval(leSample)]] <- prepATable(df, 0.1, ass.table)
  } else if( c.name == "colon" ) {
    ass.table <- tbl.mutations$colon$Gene
    
    colon.list[[eval(leSample)]] <- prepATable(df, 0.1, ass.table)
    
  } else if( c.name == "leukemia" ) {
    ass.table <- tbl.mutations$leukemia$Gene
    
    leukemia.list[[eval(leSample)]] <- prepATable(df, 0.1, ass.table)
    
  } else {
    ass.table = NULL
    other.list[[eval(leSample)]] <- prepATable(df, 0.1, ass.table)
    
  }
  
  #names(df.list); dim(df.list$pathogenic); dim(df.list$potential); dim(df.list$associated); dim(df.list$unscored); dim(df.list$benign)
  #openxlsx::write.xlsx(df.list, file = file.path(FILESOURCE, "output", paste(m[1],"_",m[2],"_",m[3],"_variants.fitlered.xlsx", sep="")))
  
}

sample.list <- list()
sample.list[["breast"]] <- breast.list
sample.list[["colon"]] <- colon.list
sample.list[["leukemia"]] <- leukemia.list
sample.list[["other"]] <- other.list

out.pat.list <- list()
out.pot.list <- list()
out.rev.list <- list()
for( leCohort in names(sample.list) ) {
  ## for the cohort get the list with all the samples
  cohort.list <- sample.list[[eval(leCohort)]]
  
  df.pat.out <- NULL
  df.pot.out <- NULL
  df.rev.out <- NULL
  ## process every sample
  for( leSample in names(cohort.list) ) {
    
    df.list <- cohort.list[[eval(leSample)]]
    # pathogenic
    df.pat.tmp <- df.list$pathogenic %>%
      dplyr::mutate(sampleID = leSample,
                    cohort = eval(leCohort), .before = Chr) %>%
      rename(Variant = eval(leSample)) %>%
      mutate( GT = sub(":.*", "", Variant),
              VAF = sub(".*:", "", Variant), .before = cytoBand)
    
    df.pat.out <- rbind(df.pat.out, df.pat.tmp)
    # potential
    df.pot.tmp <- df.list$potential %>%
      dplyr::mutate(sampleID = leSample,
                    cohort = eval(leCohort), .before = Chr) %>%
      rename(Variant = eval(leSample)) %>%
      mutate( GT = sub(":.*", "", Variant),
              VAF = sub(".*:", "", Variant), .before = cytoBand)
    
    df.pot.out <- rbind(df.pot.out, df.pot.tmp)
    # revel75
    df.rev.tmp <- df.list$revel %>%
      dplyr::mutate(sampleID = leSample,
                    cohort = eval(leCohort), .before = Chr) %>%
      rename(Variant = eval(leSample)) %>%
      mutate( GT = sub(":.*", "", Variant),
              VAF = sub(".*:", "", Variant), .before = cytoBand)
    
    df.rev.out <- rbind(df.rev.out, df.rev.tmp)
  }
  
  
  
  
  
  out.pat.list[[eval(leCohort)]] <- df.pat.out
  out.pot.list[[eval(leCohort)]] <- df.pot.out
  out.rev.list[[eval(leCohort)]] <- df.rev.out
}

## Sanity 
names(out.pat.list)
# merge all
out.pat.list[["full"]] <- rbind(out.pat.list[["breast"]], out.pat.list[["colon"]], out.pat.list[["leukemia"]], out.pat.list[["other"]])
out.pot.list[["full"]] <- rbind(out.pot.list[["breast"]], out.pot.list[["colon"]], out.pot.list[["leukemia"]], out.pot.list[["other"]])
out.rev.list[["full"]] <- rbind(out.rev.list[["breast"]], out.rev.list[["colon"]], out.rev.list[["leukemia"]], out.rev.list[["other"]])


openxlsx::write.xlsx(out.pat.list, file.path("data/output", "Legacy_VRI_Cancers_Pathogenic_Variants.xlsx"))
openxlsx::write.xlsx(out.pot.list, file.path("data/output", "Legacy_VRI_Cancers_Potential_Variants.xlsx"))
openxlsx::write.xlsx(out.rev.list, file.path("data/output", "Legacy_VRI_Cancers_Revel75_Variants.xlsx"))


# Play around and find good set
dummy <- out.pot.list[["full"]]
table(dummy$Gene.refGene)[order(table(dummy$Gene.refGene), decreasing = TRUE)]


meta$Sample[which(meta$fID == "S21")]


muh <- sample.list$breast$DXB_HLA_20_4062281_00474$potential



save.list$DXB_HLA_20_4062281_00474[which(save.list$DXB_HLA_20_4062281_00474$Gene.refGene == "ADH1C"),]

fuckthis <- save.list$DXB_HLA_20_4062281_00474[which(save.list$DXB_HLA_20_4062281_00474$Chr == "chr4"),]



meh <- fuckthis[which(fuckthis$Start > 99347030),]

meta[which(meta$fID == "S66"),]


### this sample supposedly has a KMT2C variant
DXB_HLA_21_4084900_00959

lookup <- save.list$DXB_HLA_21_4084900_00959[which(save.list$DXB_HLA_21_4084900_00959$Gene.refGene == "KMT2C"),]

openxlsx::write.xlsx(lookup, file.path("data/output", "S66_lookup_KMT2C.xlsx"))



sample.list <- readRDS(file.path("data", "RDS", "VRI_Pilot_SingleSample_annotated.filtered.26122024.RDS"))

new.filt <- sample.list$breast$DXB_HLA_21_4084900_00959

names(new.filt)

le.bind <- rbind(new.filt[["pathogenic"]], new.filt[["potential"]], new.filt[["associated"]], new.filt[["benign"]], new.filt[["likely_benign"]], new.filt[["unscored"]])

le.bind[which(le.bind$Gene.refGene == "KMT2C"),]

new.filt[["benign"]][which(new.filt[["benign"]]$Gene.refGene == "KMT2C"),]




