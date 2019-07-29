#'@title Set up analysis directory for pairs of GWAS summary statistics
#'@description A helper function that prepares a working directory for running an analysis with CAUSE.
#'@details This function prepares a working directory to replicate the results in Section 2.3 of Morrison et al 2019
#'(https://www.biorxiv.org/content/10.1101/682237v3). For more details see https://jean997.github.io/cause/gwas_pairs.html.
#'@export
setup_gwas_pairs <- function(){
  #Download data files
  cat("Downloading data\n")
  if(!dir.exists("data/")) system("mkdir data")
  download.file(url = "https://zenodo.org/record/3351659/files/gwas_info.csv?download=1", destfile="data/gwas_info.csv")
  consortia <- c("giant", "giant", "lu",
                 "glg", "glg", "glg", "glg",
                 "ckdgen", "gefos",
                 "egg", "egg", "egg",
                 "vanderHarst", "diagram",
                 "megastroke", "magic")
  traits <- c("height", "bmi", "bfp",
              "tg", "ldl", "hdl", "tc",
              "egfrcrea",   "bone",
              "bl", "bw", "hc",
              "cad", "t2d", "as", "fg")
  tags <- paste0(consortia, "_", traits)
  for(tag in tags){
    download.file(url = paste0("https://zenodo.org/record/3351659/files/", tag, "_summary_statistics.tsv.gz?download=1"),
                  destfile=paste0("data/", tag, "_summary_statistics.tsv.gz"))
  }

  #Download LD data
  if(!dir.exists("ld/")) system("mkdir ld")
  for(chr in 1:22){
    download.file(url=paste0("https://zenodo.org/record/1464357/files/chr", chr, "_AF0.05_0.1.RDS?download=1"),
                             destfile=paste0("ld/chr", chr, "_AF0.05_0.1.RDS"))
    download.file(url=paste0("https://zenodo.org/record/1464357/files/chr", chr, "_AF0.05_snpdata.RDS?download=1"),
                  destfile=paste0("ld/chr", chr, "_AF0.05_snpdata.RDS"))
  }

  #Download code
  if(!dir.exists("ld/")) system("mkdir R")
  files <- c("cause.R", "ivw.R", "ld_cat.R", "ld_prune_one_chrom.R",
             "mregger.R", "mrpresso.R")
  for(f in files){
    download.file(url=paste0("https://raw.githubusercontent.com/jean997/cause/master/gwas_pairs_code/", f),
                  destfile = paste0("R/", f))
  }
  download.file(url=paste0("https://raw.githubusercontent.com/jean997/cause/master/gwas_pairs_code/pairs_snakemake.py"),
                destfile = "pairs_snakemake.py")
 }
