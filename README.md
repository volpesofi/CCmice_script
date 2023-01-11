# CCmice scripts

[![DOI](https://zenodo.org/badge/587369919.svg)](https://zenodo.org/badge/latestdoi/587369919)

R_qtl2 scripts used to produce the results of the manuscript "ΔF508-Cftr mutation in genetically diverse Collaborative Cross mice yields novel disease-relevant phenotypes for cystic fibrosis".

Haplotype reconstruction of Collaborative Cross mice genotypes was performed with [Rqtl2](https://kbroman.org/qtl2/) in the R environment (v. 3.6.3). 

### Samples preparation ###
Sample preparations was realised with ``samples_preparation.R`` script, as prescribed by Karl Broman in [convert_cc_data.R](https://github.com/rqtl/qtl2data/blob/main/CC/R/convert_cc_data.R) script.

VCF files were firstly filtered on the MegaMUGA and GigaMUGA (combined) SNP genomics location used for Collaborative Cross genotyping using vcftools, then imported in R with data.table (v 1.14.0), prepared with qtl2convert (v0.24) and merged with the original genome sequences of the CC lines (Srivastava et al. 2017) retrieved from [Zenodo](https://doi.org/10.5281/zenodo.377036). Founders’ genotypes for use with R/qtl2 were retrieved from [figshare](https://doi.org/10.6084/m9.figshare.5404762.v2).

### Rqtl2 analysis and Mosaic plot ###
Rqtl2 analysis and mosaic plots are performed with the ``rqtl2_analysis.R`` script. 

For each sequenced individual, the observed fraction of the autosomes assigned to each of the 36 possible diplotype states were computed with Rqtl2/calc_genoprob function. We then reduced genotype probabilities to the allele probabilities of each founder haplotype using Rqtl2/genoprob_to_alleleprob function. A haplotype mosaic plot for each sample was created with ggplot2 (v. 3.3) from the allele probability data-frame rounded to 0, 0.5 or 1. The obtained haplotype mosaic plot was compared with the original CC line, CC037 and CC006 respectively, and the small differences detected were manually inspected. The density of heterozygous variants was used to check the reliability of haplotype reconstruction by Rqtl2.

