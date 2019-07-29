# Snakemake pipeline for running CAUSE on pairs of GWAS traits
#
# This snakefile performs the following steps:
#
# * Reads in all pairs of gwas files, writes out a list of SNPs in all files and with minimial LD
# * For each pair, calculate the distribution of direct effects
# * For each trait, choose a set of SNPs in minimal LD and with p < 1e-3
# * For each pair, run CAUSE
# * For each pair, create a png of posterior and with z-score in title
# * Create a plot using all results
#
# LICENSE: CC0. Do what you want with the code, but it has no guarantees.
#          https://creativecommons.org/share-your-work/public-domain/cc0/
#
# Installation with conda package manager:
#
# conda create -n cause_large python=3.6  pandas snakemake
# source activate cause_large


#To run the full pipeline, submit the following line from within the
#same directory as the Snakefile while on the head node (the paths to
#        the data files are relative to the Snakefile):
# nohup snakemake -s pairs_snakemake.py --keep-going --jobs 96 --cluster "sbatch --output={params.log}_%A.out --error={params.log}_%A.err --cpus-per-task={params.cpus} --ntasks=1 --mem-per-cpu={params.mem} --account=pi-xinhe --partition=broadwl --time=1:00:00 --job-name={params.jobname}" > pairs.out &
#


import pandas as pd

data_dir = "data/" #where the data is
ld_dir = "ld/" #where the ld data is
cause_dir = "cause/" #where CAUSE results will go
mr_dir = "mr/" #where other MR method results will go

consortia = ["giant", "giant", "lu",
            "glg", "glg", "glg", "glg",
            "icbp", "icbp", "icbp", "icbp",
            "ckdgen", "gefos",
            "egg", "egg", "egg",
            "vanderHarst", "diagram",
            "megastroke", "magic"]

traits = ["height", "bmi", "bfp",
          "tg", "ldl", "hdl", "tc",
          "dbp", "sbp", "map", "pp",
          "egfrcrea",   "bone",
          "bl", "bw", "hc",
          "cad", "t2d",
          "as", "fg"]

tags = [consortia[i] + "_" + traits[i] for i in range(len(traits))]

tag_pairs = [(tag1, tag2) for tag1 in tags for tag2 in tags if tag1!=tag2]


rule all:
    input: expand(mr_dir + '{tp[0]}__{tp[1]}_ivw.RDS', tp = tag_pairs),
           expand(mr_dir + '{tp[0]}__{tp[1]}_mrpresso.RDS', tp = tag_pairs),
           expand(mr_dir + '{tp[0]}__{tp[1]}_mregger.RDS', tp = tag_pairs),
           expand(cause_dir + '{tp[0]}__{tp[1]}_cause.RDS', tp = tag_pairs)

## Ascertained SNP list, one per trait pair

# This one writes data for trait M but only what overlaps with Y
rule data_overlap:
    input: file1 = data_dir + '{tag1}_summary_statistics.tsv.gz',
           file2 = data_dir + '{tag2}_summary_statistics.tsv.gz'
    output: out= temp(data_dir + "{tag1}__{tag2}_overlap.tsv.gz")
    params: log="tempov", mem="20G", cpus="1",
            jobname='dataov'
    shell: """
           snps() {{ gzip -cd "$@" | awk '{{ if ($7!="NA" && $7 > 0 && $6!="NA") print $3 }}' ;}}
           full() {{ gzip -cd "$@" | awk '{{ if ($7!="NA" && $7 > 0 && $6!="NA") print $0 }}' ;}}
           snps {input.file2} | awk 'NR==FNR{{F1[$0];next}}$3 in F1{{print}}' - <(full {input.file1}) | gzip > {output.out}
           """
rule ld_prune_one_chrom:
    input: data = data_dir + '{tag1}__{tag2}_overlap.tsv.gz',
           ld1 = ld_dir + 'chr{chrom}_AF0.05_0.1.RDS',
           ld2 = ld_dir + 'chr{chrom}_AF0.05_snpdata.RDS'
    output: out=temp(data_dir + "snps_{tag1}__{tag2}.pruned.{chrom}.RDS")
    params: log="ldprune", mem="10G", cpus="4",
            pval_thresh = "1e-3", r2_thresh = 0.1 , jobname='ldprune_{chrom}' 
    shell:   'Rscript R/ld_prune_one_chrom.R {input.data} {wildcards.chrom}  \
                   {params.pval_thresh} {params.r2_thresh} {input.ld1} {input.ld2} {output.out}'

rule ld_prune_combine:
    input: fls = expand( data_dir + "snps_{{tag1}}__{{tag2}}.pruned.{chr}.RDS", chr = range(1, 23))
    output: out1 = data_dir + "snps_{tag1}__{tag2}.pruned.txt"
    params: log="ld_comb", mem="2G", cpus="1",
            jobname='combine' 
    shell: "Rscript R/ld_cat.R {output.out1} {input.fls}"



#Run CAUSE
rule cause:
    input: file1 = data_dir + "{tag1}__{tag2}_overlap.tsv.gz",
           file2 = data_dir + '{tag2}__{tag1}_overlap.tsv.gz',
           snps = data_dir + 'snps_{tag1}__{tag2}.pruned.txt'
    output: params = cause_dir + '{tag1}__{tag2}_params.RDS',
            cause = cause_dir + '{tag1}__{tag2}_cause.RDS',
            data = data_dir + '{tag1}__{tag2}_data.RDS'
    params: log="cause", mem="5G", cpus="8",
            jobname='cause', seed = 100
    shell: 'Rscript R/cause.R {input.file1} {input.file2}  \
                   {input.snps} {output.params} \
                   {output.cause} {output.data} {params.seed}'


## Other MR
rule ivw:
    input: data = data_dir + '{tag1}__{tag2}_data.RDS'
    output: out = mr_dir + '{tag1}__{tag2}_ivw.RDS',
    params: log="mr", mem="1G", cpus="1",
            jobname='mr_{tag1}__{tag2}'  
    shell: 'Rscript R/ivw.R {input.data} 5e-8 {output.out} '

rule egger:
    input: data = data_dir + '{tag1}__{tag2}_data.RDS'
    output: out = mr_dir + '{tag1}__{tag2}_mregger.RDS',
    params: log="mr", mem="1G", cpus="1",
            jobname='mr_{tag1}__{tag2}'  
    shell: 'Rscript R/mregger.R {input.data} 5e-8 {output.out} '

rule mrpresso:
    input: data = data_dir + '{tag1}__{tag2}_data.RDS'
    output: out = mr_dir + '{tag1}__{tag2}_mrpresso.RDS',
    params: log="mrp", mem="5G", cpus="1",
            jobname='mrp_{tag1}__{tag2}'  
    shell: 'Rscript R/mrpresso.R {input.data} 5e-8 {output.out} '


## Plot
#rule plot:
#    input: expand(mr_dir + '{tp[0]}__{tp[1]}_ivw.RDS', tp = tag_pairs),
#           expand(mr_dir + '{tp[0]}__{tp[1]}_mrpresso.RDS', tp = tag_pairs),
#           expand(mr_dir + '{tp[0]}__{tp[1]}_mregger.RDS', tp = tag_pairs),
#           expand(cause_dir + '{tp[0]}__{tp[1]}_cause.RDS', tp = tag_pairs)
#    output: cause_dir + "gwas_fdr_cause_ivw.png", 
#            cause_dir +  "gwas_q_prop_m2.png",
#            cause_dir + "gwas_fdr_ivw_egger_mrp.png"
#    params: log="plot", mem="1G", cpus="1",
#            jobname='plot', cause_dir = cause_dir, mr_dir = mr_dir  
#    shell: 'Rscript R/plot.R {params.cause_dir} {params.mr_dir}'
 
