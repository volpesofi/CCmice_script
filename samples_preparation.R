# an R script to perform QTL analysis with R/qtl2
# Authors: Anna Sofia Tascini

###### load libraries #######
library(plyr)
library("qtl2")
library("qtl2convert")
library("VariantAnnotation")
library("snpStats")
library("rgr")
library("broman")
library("readxl")
library("data.table")
library(patchwork)
library(cowplot)
library("optparse")
library(stringr)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggplot2)
library(viridisLite)
library("viridis")  
library(forcats)
library(ggplot2)
library(dplyr)

###### set wd #######
setwd("/beegfs/scratch/ric.cosr/ric.bragonzi/BragonziA_1641_WGS_mouse")


###### download and import data for CC mice #######

rawdata_dir <- "data/Rawdata"
# create directory for data if it doesn't exist
if(!dir.exists(rawdata_dir)) {
  dir.create(rawdata_dir, recursive = T)
}

zenodo_url <- "https://zenodo.org/record/377036/files/"
zenodo_postpend <- "?download=1"
# download genotypes.zip if not available
genotype_file <- "genotypes.zip"
distant_file <- paste0(zenodo_url, genotype_file, zenodo_postpend)
local_file <- file.path(rawdata_dir, genotype_file)
if(!file.exists(local_file)) {
  message("downloading genotypes.zip")
  download.file(distant_file, local_file)
}

# create directory for data if it doesn't exist
genotype_dir <- file.path(rawdata_dir, "genotypes")
if(!dir.exists(genotype_dir)) {
  dir.create(genotype_dir, recursive = T)
}

# unzip if it hasn't been unzipped
file <- file.path(genotype_dir, "SEQgenotypes.csv")
if(!file.exists(file)) {
  message("unzipping genotypes")
  unzip(local_file, exdir=genotype_dir)
}

# founder genotypes from figshare
# https://doi.org/10.6084/m9.figshare.5404762.v2
fg_url <- "https://ndownloader.figshare.com/files/13623080"
local_file <- file.path(rawdata_dir, "MMnGM_processed_files.zip")
if(!file.exists(local_file)) {
  message("downloading founder genotypes")
  download.file(fg_url, local_file)
}

# unzip the founder genotype data
if(!file.exists(file.path(rawdata_dir, "MMnGM", "MMnGM_info.csv"))) {
  message("unzipping founder genotypes")
  unzip(local_file, exdir=rawdata_dir)
}

# load the allele codes
fga <- read_csv(file.path(rawdata_dir, "MMnGM", "MMnGM_allelecodes.csv"), comment = "#")
fga <- fga[fga$chr %in% c(1:19,"X"),]

# ------ chr for cicle ------
# read genotypes
message("reading genotypes")
g <- data.table::fread(file, data.table=FALSE)
g[g=="N" | g=="H"] <- NA


genotype_dir <- file.path(rawdata_dir, "genotypes")

# save file names into list for write_control_file function 
genof_list = list()
foundergenof_list = list()
pmapf_list = list()
gmapf_list = list()

p = "/beegfs/scratch/ric.cosr/ric.cosr/Bragonzi_ASTascini/bioinfo_analysis/BragonziA_1523_mouse_genomics/data/rqtl2_geno/"
p = "/beegfs/scratch/ric.cosr/ric.bragonzi/NatComm/data/rqtl2_geno/"
samples = c("S1798Nr1","S1798Nr2")

###### read and save mysamples genotype #######
for (chr in c(1:19, "X")) {
  chr_i = paste("chr", chr, sep="")
  print(paste("processing",chr_i))
  # ----- genotipyng our mouse -----
  vcf_file = paste(p,chr_i,"_Sample_Founders_annotated.vcf.gz", sep ='')
  vcf = data.table::fread(input = vcf_file,
                          sep = "\t",
                          skip = "#CHROM")
  
  # ---- Filter PASS -----
  filter =grep("PASS", vcf[["FILTER"]])
  vcf_f = vcf[filter,]
  df_vcf_cftr = as.data.frame(vcf_f)
  
  
  geno_list = list()
  for (i in c(which(colnames(df_vcf_cftr) %in% samples))) {
    sample_id =colnames(df_vcf_cftr)[i]
    geno = df_vcf_cftr[,grep(sample_id,colnames(df_vcf_cftr))]
    gene_onlyType_list = str_split(geno, pattern = ":")
    gene_onlyType = sapply(1:length(gene_onlyType_list), function(x){gene_onlyType_list[[x]][[1]]})
    df_vcf_cftr$gtmp = gene_onlyType
    df_vcf_cftr$GENO = "H"
    mydt = data.table(df_vcf_cftr)
    # print(table(mydt$gtmp))
    mydt$GENO[mydt$gtmp == "1/1"] <- mydt$ALT[mydt$gtmp == "1/1"]
    mydt$GENO[mydt$gtmp == "1|1"] <- mydt$ALT[mydt$gtmp == "1|1"]
    mydt$GENO[mydt$gtmp == "2/2"] <- mydt$ALT[mydt$gtmp == "2/2"]
    mydt$GENO[mydt$gtmp == "2|2"] <- mydt$ALT[mydt$gtmp == "2|2"]
    mydt$GENO[mydt$gtmp == "0|0"] <- mydt$REF[mydt$gtmp == "0|0"]
    mydt$GENO[mydt$gtmp == "0/0"] <- mydt$REF[mydt$gtmp == "0/0"]
    mydt$GENO[mydt$gtmp == "./."] <- mydt$REF[mydt$gtmp == "./."]
    geno_list[[sample_id]] = mydt[,c("ID","POS", "GENO")]
  }
  Geno_table = geno_list %>% reduce(left_join, by =c("ID", "POS"))
  colnames(Geno_table) <- c("marker", "position(b38)",
                            "CC037_S1_4", "CC006_S2_17")
  
  # filter Srivastava et al. (2017) on chr
  g_chr_f = g[g$chromosome == chr,]
  
  # merge my samples with Srivastava et al
  geno_i = merge(x = g_chr_f, y = Geno_table, by = c("marker","position(b38)"))
  geno_i = geno_i[order(geno_i$`position(b38)`,decreasing = F),]
  geno_i[geno_i=="N" | geno_i=="H"] <- NA
  
  # omit markers not in founder genotypes
  geno_i <- geno_i[geno_i$marker %in% fga$marker,]
  fga_chr <- fga[fga$marker %in% geno_i$marker,]
  
  # reorder rows of genotype data to match fga
  geno_i <- geno_i[match(geno_i$marker, fga_chr$marker),]
  
  # cut down to just the genotypes
  rownames(geno_i) <- geno_i$marker
  geno_i <- geno_i[,-(1:3)]
  
  # strip off individual IDs from column names
  colnames(geno_i) <- paste0(sapply(strsplit(colnames(geno_i), "Unc"), "[", 1), "Unc")
  
  # encode genotypes
  message("encoding genotypes")
  geno_i <- geno_i[,1:(length(colnames(geno_i)))]
  geno_i <- encode_geno(geno_i, fga_chr[,c("A","B")], cores=4)
  
  # omit markers with no data
  geno_i <- geno_i[rowSums(geno_i!="-") > 0, , drop=FALSE]
  
  # writing genotype
  file <- paste0(genotype_dir,"/cc_geno", chr, ".csv")
  genof_list[[chr_i]] = file
  gsub <- geno_i[rownames(geno_i) %in% fga_chr$marker[fga_chr$chr==chr], , drop=FALSE]
  write2csv(cbind(marker=rownames(gsub), gsub), file,
            overwrite=TRUE,
            comment=paste("Chromosome", chr_i, "genotypes",
                          "for Collaborative Cross (CC) lines inferred from",
                          "genotypes from Srivastava et al. (2017)",
                          "doi:10.1534/genetics.116.198838,",
                          "data at doi:10.5281/zenodo.377036",
                          "merged with Bragonzi et al. mice"))
  
  
  foundergenof_list[[chr_i]] = paste0(rawdata_dir,"/MMnGM/MMnGM_foundergeno",chr,".csv")
  pmapf_list[[chr_i]] = paste0(rawdata_dir,"/MMnGM/MMnGM_pmap",chr,".csv")
  gmapf_list[[chr_i]] = paste0(rawdata_dir,"/MMnGM/MMnGM_gmap",chr,".csv")
}

# ----- covar file -----
# download Prob36.zip if not available
zenodo_url <- "https://zenodo.org/record/377036/files/"
zenodo_postpend <- "?download=1"
prob_file <- "Prob36.zip"
distant_file <- paste0(zenodo_url, prob_file, zenodo_postpend)
local_file <- file.path(rawdata_dir, prob_file)
if(!file.exists(local_file)) {
  message("Downloading Prob36.zip")
  download.file(distant_file, local_file)
}

# create directory for data if it doesn't exist
prob_dir <- file.path(rawdata_dir, "Prob36")
if(!dir.exists(prob_dir)) {
  dir.create(prob_dir, recursive= T)
}

# unzip if it hasn't been unzipped
gzfile <- file.path(prob_dir, "CC001-Uncb38V01.csv.gz")
csvfile <- file.path(prob_dir, "CC001-Uncb38V01.csv")
if(!file.exists(gzfile) && !file.exists(csvfile)) {
  message("unzipping Prob36.zip")
  unzip(local_file, exdir=prob_dir)
}

# gunzip all of the Prob36 files
if(!file.exists(csvfile)) {
  files <- list.files(prob_dir, pattern=".csv.gz$")
  
  message("unzipping the probability files")
  
  for(file in files) {
    system(paste("gunzip", file.path(prob_dir, file)))
  }
}

# load the Prob36 files and determine X, Y, and M genotypes
message("reading the probability files")
files <- list.files(prob_dir, pattern=".csv$")
strains <- sub("\\.csv$", "", files)
probs <- setNames(vector("list", length(strains)), strains)
for(i in seq_along(files)) {
  probs[[i]] <- data.table::fread(file.path(prob_dir, files[i]), data.table=FALSE)
}

##############################
# guess the rest of the cross order
##############################
message("inferring cross order")
mprob <- t(sapply(probs, function(a) colMeans(a[a[,2]=="M", paste0(LETTERS, LETTERS)[1:8]])))
yprob <- t(sapply(probs, function(a) colMeans(a[a[,2]=="Y", paste0(LETTERS, LETTERS)[1:8]])))
xprob <- t(sapply(probs, function(a) colMeans(a[a[,2]=="X", paste0(LETTERS, LETTERS)[1:8]])))
# a bunch where we can't tell Y or M

# also need the supplementary data file
supp_file <- "SupplementalData.zip"
distant_file <- paste0(zenodo_url, supp_file, zenodo_postpend)
local_file <- file.path(rawdata_dir, supp_file)
if(!file.exists(local_file)) {
  message("downloading SupplementalData.zip")
  download.file(distant_file, local_file)
}

# extract just the CCStrains.csv file
csv_file <- "SupplmentalData/CCStrains.csv"
csv_file <- file.path(rawdata_dir, csv_file)
if(!file.exists(csv_file)) {
  message("unzipping SupplementalData.zip")
  unzip(local_file, csv_file, exdir=rawdata_dir)
}

ccstrains <- data.table::fread(csv_file, data.table=FALSE)

# check that the strain names are the same
stopifnot( all( paste0(sub("/", "-", ccstrains$Strain, fixed=TRUE), "b38V01") ==
                  strains))

# download table S2
url <- "http://www.genetics.org/highwire/filestream/438137/field_highwire_adjunct_files/10/TableS2.xlsx"
file <- basename(url)
local_file <- file.path(rawdata_dir, file)
if(!file.exists(local_file)) {
  download.file(url, local_file)
}
# TableS2.xlsx download does not always work
# save it to zenodo
tab <- as.data.frame(readxl::read_xlsx(local_file))

# check that the strain names are the same
stopifnot( all( paste0(sub("/", "-", tab$Strain, fixed=TRUE), "b38V01") ==
                  strains))

funnel <- tab[,"Funnel Code"]
mtdna <- tab[,"Mitochondria"]
ychr <- tab[,"Chromosome Y"]

# determine cross orders; use orders in Table S2 when available
cross_info <- matrix(ncol=8, nrow=length(strains))
dimnames(cross_info) <- list(ccstrains$Strain, LETTERS[1:8])
for(i in which(!is.na(funnel))) {
  cross_info[i,] <- match(unlist(strsplit(funnel[i], "")), LETTERS)
}

# when not available, use the values in table S2
use_mtdna <- mtdna; use_mtdna[mtdna=="A/D" | mtdna=="D/A"] <- "A" # use A when A/D or D/A
cross_info[is.na(cross_info[,1]),1] <- match(use_mtdna, LETTERS[1:8])[is.na(cross_info[,1])]
cross_info[is.na(cross_info[,8]),8] <- match(ychr, LETTERS[1:8])[is.na(cross_info[,8])]

# problems: CC013, CC023, CC027
# CC013 M = Y = E [genotypes say Y could be B or C, too] --- change Y to B

cross_info["CC013/GeniUnc", 8] <- 2

# check cases where mitochondrial genotype is clear
max_mprob <- apply(mprob, 1, max)
wh_mprob <- apply(mprob, 1, which.max)
stopifnot( all(cross_info[max_mprob>0.9,1] == wh_mprob[max_mprob > 0.9]) )

stopifnot( all( cross_info[,1] != cross_info[,8] ))

# are there any obvious differences?

# mitochondria
stopifnot( all( sort(max_mprob[wh_mprob != cross_info[,1]]) < 0.5))

# now the Y chr
max_yprob <- apply(yprob, 1, max)
wh_yprob <- apply(yprob, 1, which.max)
stopifnot( all( sort(max_yprob[wh_yprob != cross_info[,8]]) < 0.5))

# founders with largest X chr probabilities
wh_xprob <- t(apply(xprob, 1, order, decreasing=TRUE))

# most common X chr genotype that's not the mtDNA one in the 3rd slot
cross_info[is.na(cross_info[,3]),3] <- wh_xprob[is.na(cross_info[,3]),1]
stopifnot( all(cross_info[,1] != cross_info[,3]) )

# sort other X chr probs...put three most probable in the 2nd, 5th, and 6th slots
for(i in which(is.na(cross_info[,2]))) {
  z <- (wh_xprob[i, ] %wnin% cross_info[i, c(1,3,8)])
  cross_info[i, c(2,5,6)] <- sample(z[1:3])
  cross_info[i, c(4,7)] <- sample(z[4:5])
}

# further problems:
# CC031 Y chr is B but this is clearly on the X chromosome
# CC037 Y chr is D but this is clearly on the X chromosome
# CC056 Y chr is E but this is clearly on the X chromosome
# ... swap these with one of 2,5,6

cross_info["CC031/GeniUnc", ] <- c(1,6,5,3,7,2,4,8)
cross_info["CC037/TauUnc", ] <- c(3,4,8,7,6,1,5,2)
cross_info["CC056/GeniUnc", ] <- c(3,8,1,7,5,4,6,2)

stopifnot( all(apply(cross_info, 1, sort) == 1:8) )

message("writing cross and covariate data")

# write the cross information to a file
write2csv(cbind(id=rownames(cross_info), cross_info),
          file.path(rawdata_dir,"cc_crossinfo.csv"),
          comment=paste("Cross information for Collaborative Cross (CC) lines inferred from",
                        "genotypes from Srivastava et al. (2017)",
                        "doi:10.1534/genetics.116.198838,",
                        "data at doi:10.5281/zenodo.377036"),
          overwrite=TRUE)

Binfo = read.csv(file.path(rawdata_dir,"cross_info_Bragonzi.csv"), row.names = 1)
cross_info_Srivastava_Bragonzi = rbind(cross_info, Binfo)
write2csv(cbind(id=rownames(cross_info_Srivastava_Bragonzi), 
                cross_info_Srivastava_Bragonzi),
          file.path(rawdata_dir,"./cross_info_Srivastava_Bragonzi.csv"),
          comment=paste("Cross information for Collaborative Cross (CC) lines inferred from",
                        "genotypes from Srivastava et al. (2017)",
                        "doi:10.1534/genetics.116.198838,",
                        "data at doi:10.5281/zenodo.377036",
                        "merged with Bragonzi et al. mice"),
          overwrite=TRUE)

# write covariate info with M and Y as inferred
covar <- data.frame(id=rownames(cross_info),
                    mitochondria=tab$Mitochondria,
                    Ychr=tab[,"Chromosome Y"],
                    n_founders=tab[,"# of Founders"],
                    origin=tab[,"Origin of Strain"],
                    stringsAsFactors=FALSE)

write2csv(covar, 
          file.path(rawdata_dir,"./cc_covar.csv"),
          comment=paste("Covariate information for Collaborative Cross (CC) lines inferred from",
                        "genotypes from Srivastava et al. (2017)",
                        "doi:10.1534/genetics.116.198838,",
                        "data at doi:10.5281/zenodo.377036"),
          overwrite=TRUE)

Bcova = read.csv(file.path(rawdata_dir,"covar_Bragonzi.csv"))
covar$sex = "male"
covar_Srivastava_Bragonzi = rbind(covar, Bcova)
write2csv(covar_Srivastava_Bragonzi,
          file.path(rawdata_dir,"./covar_Srivastava_Bragonzi.csv"),
          comment=paste("Cross information for Collaborative Cross (CC) lines inferred from",
                        "genotypes from Srivastava et al. (2017)",
                        "doi:10.1534/genetics.116.198838,",
                        "data at doi:10.5281/zenodo.377036",
                        "merged with Bragonzi et al. mice"),
          overwrite=TRUE)


# ------ CONTROL FILE ------
write_control_file(description ="Backcrossing of CC037 and CC006 mice with Cftr-ΔF508 carrying mice",
                   output_file = "control.yaml",
                   crosstype = "risib8",
                   na.strings = c("-", "NA"),
                   comment.char = "#",
                   geno_file = genof_list,
                   founder_geno_file = foundergenof_list,
                   pmap_file = pmapf_list,
                   gmap_file = gmapf_list, # gmapf_list,
                   #pheno_file = "pheno.csv",
                   covar_file = file.path(rawdata_dir,"covar_Srivastava_Bragonzi.csv"),
                   crossinfo_file = file.path(rawdata_dir,"cross_info_Srivastava_Bragonzi.csv"),
                   alleles = c("A","B", "C", "D", "E", "F", "G", "H"),
                   sex_covar = "sex",
                   xchr = "X",
                   sex_codes = list(
                     "male" = "male",
                     "female" = "female"),
                   sep = ",",
                   geno_codes = list("A"=1, "H" = 2, "B"=3),
                   geno_transposed = TRUE,   # if TRUE, markers as rows
                   founder_geno_transposed = TRUE,   # if TRUE, markers as rows
                   pheno_transposed = FALSE, # if TRUE, phenotypes as rows
                   covar_transposed = FALSE, # if TRUE, covariates as rows
                   overwrite = TRUE)

