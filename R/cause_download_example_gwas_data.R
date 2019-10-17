#'@title Download GWAS data to use in examples
#'@details This function downloads Snakemake summary statistics for LDL cholesterol (Willer et al 2013 PMID 24097068),
#'coronary artery disease (van der Harst et al 2017 PMID 29212778), and asthma (Demenais et al 2018 PMID 29273806) as well as a csv file that can be
#'used with the CAUSE pipeline.
#'@export
cause_download_example_gwas_data <- function(){

  #Download code
  if(!dir.exists("raw_data/")) system("mkdir raw_data")

  #LDL
  if(!file.exists("raw_data/jointGwasMc_LDL.txt.gz")){
    download.file(url="http://csg.sph.umich.edu/abecasis/public/lipids2013/jointGwasMc_LDL.txt.gz", destfile="raw_data/jointGwasMc_LDL.txt.gz")
  }
  #CAD
  if(!file.exists("raw_data/CAD_META.gz")){
    download.file(url="ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/vanderHarstP_29212778_GCST005194/CAD_META.gz",
                destfile="raw_data/CAD_META.gz")
  }
  #Asthma
  if(!file.exists("raw_data/TAGC_meta-analyses_results_for_asthma_risk.zip")){
    download.file(url="ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/DemenaisF_29273806_GCST005212/TAGC_meta-analyses_results_for_asthma_risk.zip",
                  destfile="raw_data/TAGC_meta-analyses_results_for_asthma_risk.zip")
    system("unzip raw_data/TAGC_meta-analyses_results_for_asthma_risk.zip")
  }
  download.file(url="https://raw.githubusercontent.com/jean997/cause/master/pipeline_code/gwas_pairs.csv",
                destfile="gwas_pairs.csv")
}
