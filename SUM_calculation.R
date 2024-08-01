#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse)) 

option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character',
              help="ssm input file"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="MS96 output file")
)

opt = parse_args(OptionParser(option_list=option_list))
print(opt)

# load the package
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))

mut <-read_csv(opt$i)
sum = data.frame(table(mut$sample_id))
colnames(sum) = c('sample_id', 'sum')

sum_s = add_column(sum, cancer_type = 'cantype', .after = 1)

cat = mut[!duplicated(mut$sample_id),]
for (s in sum_s$sample_id){
  sum_s$cancer_type[sum_s$sample_id == s] = cat$cancer_type[cat$sample_id == s]
}

write.csv(sum_s, opt$output, row.names = FALSE)
