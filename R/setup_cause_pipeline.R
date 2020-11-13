#'@title Set up analysis directory for pairs of GWAS summary statistics
#'@description A helper function that prepares a working directory for running an analysis with CAUSE.
#'@details This function downloads Snakemake pipeline and R code needed to run CAUSE on many pairs of traits.
#'@export
setup_cause_pipeline <- function(download_ldshrink=FALSE, download_eur_ld_scores=FALSE,
                                 download_plink_ref = FALSE){
  #Download LD data
  if(download_ldshrink){
    if(!dir.exists("ld/")) system("mkdir ld")
    for(chr in 1:22){
      download.file(url=paste0("https://zenodo.org/record/1464357/files/chr", chr, "_AF0.05_0.1.RDS?download=1"),
                  destfile=paste0("ld/chr", chr, "_AF0.05_0.1.RDS"))
      download.file(url=paste0("https://zenodo.org/record/1464357/files/chr", chr, "_AF0.05_snpdata.RDS?download=1"),
                  destfile=paste0("ld/chr", chr, "_AF0.05_snpdata.RDS"))
    }
  }
  if(download_eur_ld_scores){
    if(!dir.exists("ld_scores/")) system("mkdir ld_scores")
    download.file(url="https://data.broadinstitute.org/alkesgroup/LDSCORE/eur_w_ld_chr.tar.bz2",destfile="ld_scores/eur_w_ld_chr.tar.bz2")
    system("tar -xjf ld_scores/eur_w_ld_chr.tar.bz2")
    system("mv eur_w_ld_chr/ ld_socres/")
  }
  if(download_plink_ref){
    if(!dir.exists("plink_reference/")) system("mkdir plink_reference")
    download.file(url = "http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz", destfile="plink_reference/1kg.v3.tgz")
    setwd("plink_reference")
    system("tar -xzf 1kg.v3.tgz")
    setwd("..")
  }

  #Download code
  if(!dir.exists("R/")) system("mkdir R")
  files <- c("format_data.R","merge_data.R","ld_prune_one_chrom.R","ld_cat.R", "ld_prune_plink.R",
             "cause_params.R", "cause.R", "mrpresso.R", "lcv.R", "mr_package.R", "extract_results.R")
  for(f in files){
    download.file(url=paste0("https://raw.githubusercontent.com/jean997/cause/master/pipeline_code/R/", f),
                  destfile = paste0("R/", f))
  }
  download.file(url=paste0("https://raw.githubusercontent.com/jean997/cause/master/pipeline_code/Snakefile"),
                destfile = "Snakefile")
  download.file(url=paste0("https://raw.githubusercontent.com/jean997/cause/master/pipeline_code/config.yaml"),
                destfile = "config.yaml")
  download.file(url=paste0("https://raw.githubusercontent.com/jean997/cause/master/pipeline_code/cluster.yaml"),
                destfile = "cluster.yaml")
  download.file(url=paste0("https://raw.githubusercontent.com/jean997/cause/master/pipeline_code/run-snakemake.sh"),
                destfile = "run-snakemake.sh")
  download.file(url="https://raw.githubusercontent.com/jean997/cause/master/pipeline_code/gwas_pairs.csv",
                destfile="gwas_pairs.csv")
}
