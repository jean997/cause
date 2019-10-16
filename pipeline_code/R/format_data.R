library(cause)
library(readr)

args <- commandArgs(trailingOnly=TRUE)
data_file <- args[1]
delim <- args[2]
snp <- args[3]
A1 <- args[4]
A2 <- args[5]
beta_hat <- args[6]
se <- args[7]
p_value <- args[8]
sample_size <- args[9]
output_file <- args[10]


if(p_value != "NA"){
    pstring <- paste0(", `", p_value, "`='d'")
}else{
    pstring <- ""
    p_value <- NA
}    
if(sample_size!="NA"){
    sstring <- paste0(", `", sample_size, "`='d'")
}else{
    sstring <- ""
    sample_size <- NA
}    

col_string <- paste0("cols_only(`", snp, "`='c', `",
                     A1 , "`='c', `", A2, "`='c', `",
                     beta_hat , "`='d', `", se, "`='d'",
                     pstring,  sstring, ")")


if(delim=="tab"){
    X <- read_tsv(data_file, col_types = eval(parse(text = col_string)))
}else if(delim == "space"){
    X <- read_delim(data_file, delim=" ", col_types = eval(parse(text =col_string)))
}else if(delim == ","){
    X <- read_csv(data_file, col_types = eval(parse(text = col_string)))
}else{
    X <- read_delim(data_file, delim=delim, col_types = eval(parse(text = col_string)))
}

gwas_format(X, snp, beta_hat, se, A1, A2,
             p_value = p_value, sample_size = sample_size, 
             output_file = output_file, compute_pval = TRUE)

