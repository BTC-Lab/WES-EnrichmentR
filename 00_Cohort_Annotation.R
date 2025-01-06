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
    dplyr::filter(., as.numeric(REVEL_score) >= 0.25)
  
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



# COHORT SUBSETTING -------------------------------------------------------

## 0. merge samples and subset to library target bed
##bcftools_viewCommand=view -R LiftOver_HG19_to_HG38_Padded50bp_truseq-dna-exome-targeted-regions-manifest-v1-2.bed -o VRI_truseq.vcf.gz VRI_full_VQSR.vcf.gz; Date=Mon Jul 10 17:01:33 2023
##bcftools_annotateCommand=annotate -a /Users/rmolbrich/Workdir/_references/hg38/dbsnp_138.hg38.vcf.gz -c ID -o VRI_truseq_rsID.vcf.gz VRI_truseq.vcf.gz; 

## 1. add VAF tag
# bcftools +fill-tags VRI_truseq_rsID_norm.vcf.gz -Oz -o VRI_truseq_rsID_norm.VAF.vcf.gz -- -t FORMAT/VAF

#### 2. generate subsets for cancers and samples respectively
## 2.1 get all names
#bcftools query -l VRI_truseq_rsID_norm.VAF.vcf.gz > VRI_samples_full.txt
## 2.2 generate subset-lists
df.sID <- data.table::fread(file.path("data/RAW/source", "VRI_samples_full.txt"),
                                header = FALSE)

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

openxlsx::write.xlsx(meta, file.path("data/output", "letest.xlsx"))


## complete list of 298 samples
data.table::fwrite(list(meta$Sample), file = file.path("data", "RAW", "sample_subset_full.txt"))

## reduce to all matched controls plus main cancer types
dummy <- meta[(!is.na(meta$m_breast_4) | !is.na(meta$m_colon_4) | !is.na(meta$m_leukemia_4)),]
data.table::fwrite(list(dummy$Sample), file = file.path(data.dir, "RAW", "sample_subset_mainCancer.txt"))

## make maintype_extract.txt for cancer to use with bcftools
table(meta$Cancer.TypeClean)
sID.breast <- dplyr::filter(meta, Cancer.TypeClean == "Breast")
sID.colon <- dplyr::filter(meta, Cancer.TypeClean == "Colon")
sID.leukemia <- dplyr::filter(meta, Cancer.TypeClean == "Leukemia")
sID.allCancer <- dplyr::filter(meta, Cancer == "yes")

data.table::fwrite(list(sID.breast$Sample), file = file.path("data", "RAW", "samples_breast_cancer.txt"))
data.table::fwrite(list(sID.colon$Sample), file = file.path("data", "RAW", "samples_colon_cancer.txt"))
data.table::fwrite(list(sID.leukemia$Sample), file = file.path("data", "RAW", "samples_leukemia_cancer.txt"))
data.table::fwrite(list(sID.allCancer$Sample), file = file.path("data", "RAW", "samples_all_cancer.txt"))
## reduce to study set and full control set
sID.allControl <- meta[which(meta$Cancer == "no"),]
# .. and study only is basically the dataframe anyway (why didn't I do that before?!)
data.table::fwrite(list(sID.allControl$Sample), file = file.path("data", "RAW", "samples_all_control.txt"))
data.table::fwrite(list(meta$Sample), file = file.path("data", "RAW", "samples_all_study.txt"))

## add remaining/other cancers
sID.other <- meta[which(!(meta$Cancer.TypeClean %in% c("Leukemia","Breast","Colon")) & meta$Cancer == "yes"),]
data.table::fwrite(list(sID.other$Sample), file = file.path("data", "RAW", "samples_other_cancer.txt"))


#### MATCHED CONTROLS
meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS"))

sID.breast <- meta[!is.na(meta$m_breast_4),]
sID.breast <- sID.breast[is.na(sID.breast$Cancer.TypeClean),]
data.table::fwrite(list(sID.breast$Sample), file = file.path(data.dir, "RAW", "samples_breast_cancer_control.txt"))

sID.colon <- meta[!is.na(meta$m_colon_4),]
sID.colon <- sID.colon[is.na(sID.colon$Cancer.TypeClean),]
data.table::fwrite(list(sID.colon$Sample), file = file.path(data.dir, "RAW", "samples_colon_cancer_control.txt"))

sID.leukemia <- meta[!is.na(meta$m_leukemia_4),]
sID.leukemia <- sID.leukemia[is.na(sID.leukemia$Cancer.TypeClean),]
data.table::fwrite(list(sID.leukemia$Sample), file = file.path(data.dir, "RAW", "samples_leukemia_cancer_control.txt"))


### meta list for halima
length(unique(c(sID.breast$Sample, sID.colon$Sample, sID.leukemia$Sample)))

letest <- rbind.data.frame(sID.breast, sID.colon, sID.leukemia) %>%
  dplyr::filter(!duplicated(Sample))


list.out <- rbind.data.frame(letest, sID.allCancer)

openxlsx::write.xlsx(list.out, file.path("data/output", "VRI_Cohort_meta.xlsx"))


#### SUBSET for study samples
# bcftools view --samples-file sample_subset_full.txt --min-ac=1 --exclude-uncalled --output-type z --output VRI_truseq_rsID_norm_study.VAF.vcf.gz VRI_truseq_rsID_norm.VAF.vcf.gz
# 
# ## SUBSET for cohorts of interest
# bcftools view --samples-file samples_breast_cancer_control.txt --min-ac=1 --exclude-uncalled --output-type z --output VRI_breast_cancer_ctrl.vcf.gz VRI_truseq_rsID_norm_study.VAF.vcf.gz
# bcftools view --samples-file samples_breast_cancer.txt --min-ac=1 --exclude-uncalled --output-type z --output VRI_breast_cancer_case.vcf.gz VRI_truseq_rsID_norm_study.VAF.vcf.gz
# bcftools view --samples-file samples_colon_cancer_control.txt --min-ac=1 --exclude-uncalled --output-type z --output VRI_colon_cancer_ctrl.vcf.gz VRI_truseq_rsID_norm_study.VAF.vcf.gz
# bcftools view --samples-file samples_colon_cancer.txt --min-ac=1 --exclude-uncalled --output-type z --output VRI_colon_cancer_case.vcf.gz VRI_truseq_rsID_norm_study.VAF.vcf.gz
# bcftools view --samples-file samples_leukemia_cancer_control.txt --min-ac=1 --exclude-uncalled --output-type z --output VRI_leukemia_cancer_ctrl.vcf.gz VRI_truseq_rsID_norm_study.VAF.vcf.gz
# bcftools view --samples-file samples_leukemia_cancer.txt --min-ac=1 --exclude-uncalled --output-type z --output VRI_leukemia_cancer_case.vcf.gz VRI_truseq_rsID_norm_study.VAF.vcf.gz

## filter for min DP 10 again - everything that is a variant has to have DP 10
# bcftools filter -e 'FORMAT/DP < 10' VRI_breast_cancer_ctrl.vcf.gz --output-type z -o VRI_breast_cancer_ctrl.DP10.vcf.gz
# bcftools filter -e 'FORMAT/DP < 10' VRI_breast_cancer_case.vcf.gz --output-type z -o VRI_breast_cancer_case.DP10.vcf.gz
# bcftools filter -e 'FORMAT/DP < 10' VRI_colon_cancer_ctrl.vcf.gz --output-type z -o VRI_colon_cancer_ctrl.DP10.vcf.gz
# bcftools filter -e 'FORMAT/DP < 10' VRI_colon_cancer_case.vcf.gz --output-type z -o VRI_colon_cancer_case.DP10.vcf.gz
# bcftools filter -e 'FORMAT/DP < 10' VRI_leukemia_cancer_ctrl.vcf.gz --output-type z -o VRI_leukemia_cancer_ctrl.DP10.vcf.gz
# bcftools filter -e 'FORMAT/DP < 10' VRI_leukemia_cancer_case.vcf.gz --output-type z -o VRI_leukemia_cancer_case.DP10.vcf.gz

## now also make sure that at least 1/3 of samples have an alternate allele
# bcftools view -i 'COUNT(FORMAT/GT="1|1" | FORMAT/GT="0|1" | FORMAT/GT="1|0" | FORMAT/GT="1/0" | FORMAT/GT="0/1" | FORMAT/GT="1/1") >= 6' VRI_breast_cancer_case.DP10.vcf.gz -Oz -o VRI_breast_cancer_case.DP10.filtered.vcf.gz
# bcftools view -i 'COUNT(FORMAT/GT="1|1" | FORMAT/GT="0|1" | FORMAT/GT="1|0" | FORMAT/GT="1/0" | FORMAT/GT="0/1" | FORMAT/GT="1/1") >= 3' VRI_colon_cancer_case.DP10.vcf.gz -Oz -o VRI_colon_cancer_case.DP10.filtered.vcf.gz
# bcftools view -i 'COUNT(FORMAT/GT="1|1" | FORMAT/GT="0|1" | FORMAT/GT="1|0" | FORMAT/GT="1/0" | FORMAT/GT="0/1" | FORMAT/GT="1/1") >= 6' VRI_leukemia_cancer_case.DP10.vcf.gz -Oz -o VRI_leukemia_cancer_case.DP10.filtered.vcf.gz

# COHORT ANNOTATION -------------------------------------------------------

## then run current Annovar annotations
# ls *DP10.filtered.vcf.gz | awk -F'\.' '{print "/Users/rmolbrich/Workdir/Tools/annovar/table_annovar.pl "$0" /Volumes/ERGPBackup/annovar/humandb -buildver hg38 -out "$1" -remove -protocol refGene,cytoBand,exac03,avsnp151,dbnsfp47a,clinvar_20240611,gnomad41_genome -operation g,r,f,f,f,f,f -nastring . -polish -vcfinput"}' > run_annotation.sh

# COHORT ANNOTATIONS IMPORT -----------------------------------------------

## most annoying step is renaming of the columns because annovar doesn't keep track
## 1. get sampleIDs for all the subsets in the order of the files
# ls *.vcf.gz | awk -F'\.' '{print "bcftools query -l "$0" > "$1".samples.txt"}' > get_names.sh


allFiles <- list.files(file.path("data/RAW", "annotation"), pattern = "hg38_multianno.txt")

save.list <- list()

for( f.idx in allFiles ) {
  
  sID <- gsub(".hg38_multianno.txt", "", f.idx)
  
  df.annot.tmp <- data.table::fread(file.path("data/RAW", "annotation", eval(f.idx)),
                                    header = TRUE, na.strings = ".")
  
  df.sID.tmp <- data.table::fread(file.path("data/RAW", "annotation", paste0(sID, ".samples.txt")),
                                  header = FALSE)
  
  ## fix naming
  lc <- colnames(df.annot.tmp)[length(colnames(df.annot.tmp))]
  
  df.annot.tmp <- df.annot.tmp %>%
    dplyr::rename_with(~ df.sID.tmp$V1, .col = Otherinfo13:eval(lc)) %>%
    dplyr::rename_with(~ c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"), .col = Otherinfo4:Otherinfo12)
  
  ## put into the list
  save.list[[eval(sID)]] <- df.annot.tmp
}

saveRDS(save.list, file.path("data/RDS", "VRI_Pilot_Cohort_annotated.26122024.RDS"), compress = TRUE)


# COHORT ANNOTATIONS FILTER -----------------------------------------------
meta <- readRDS(file.path("data", "RDS", "meta_UAE_binned_matched_v2.RDS"))

save.list <- readRDS(file.path("data/RDS", "VRI_Pilot_Cohort_annotated.26122024.RDS"))

names(save.list)

out.list <- list()

## BREAST - cases are filtered for rarity - controls are not
out.list[["breast_cancer"]] <- prepATable(save.list[["VRI_breast_cancer_case"]], 0.001, tbl.mutations$breast$Gene)
out.list[["breast_ctrl"]] <- prepATable(save.list[["VRI_breast_cancer_ctrl"]], 1, tbl.mutations$breast$Gene)

## COLON
out.list[["colon_cancer"]] <- prepATable(save.list[["VRI_colon_cancer_case"]], 0.001, tbl.mutations$breast$Gene)
out.list[["colon_ctrl"]] <- prepATable(save.list[["VRI_colon_cancer_ctrl"]], 1, tbl.mutations$breast$Gene)

## LEUKEMIA
out.list[["leukemia_cancer"]] <- prepATable(save.list[["VRI_leukemia_cancer_case"]], 0.001, tbl.mutations$breast$Gene)
out.list[["leukemia_ctrl"]] <- prepATable(save.list[["VRI_leukemia_cancer_ctrl"]], 1, tbl.mutations$breast$Gene)

saveRDS(out.list, file.path("data/RDS", "VRI_Pilot_Cohort_annotated.filtered.26122024.RDS"), compress = TRUE)

### Make excel tables
openxlsx::write.xlsx(out.list[["breast_cancer"]], file.path("data/output", "Breast_Cancer_filtered.26122024.xlsx"))
openxlsx::write.xlsx(out.list[["colon_cancer"]], file.path("data/output", "Colon_Cancer_filtered.26122024.xlsx"))
openxlsx::write.xlsx(out.list[["leukemia_cancer"]], file.path("data/output", "Leukemia_Cancer_filtered.26122024.xlsx"))





save.list[["VRI_breast_cancer_only"]]

dummy[which(dummy$Gene.refGene == "KMT2C"),]
























