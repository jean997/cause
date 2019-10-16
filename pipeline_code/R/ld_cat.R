library(purrr)
library(readr)


args <- commandArgs(trailingOnly=TRUE)
out1 <- args[1]
args <- args[-1]


k <- map(args, function(x){readRDS(x)}) %>% unlist()

write_lines(k, path=out1)

