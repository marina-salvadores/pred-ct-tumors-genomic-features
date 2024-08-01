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
suppressPackageStartupMessages(library(MutationalPatterns))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))


## Load the corresponding reference genome.
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

mut <-read_csv(opt$i)
mut = mut[mut$mutation_type %in% c('single base substitution', 'S', 'subs', 'snvs', 'SNP'),]


if (! 'reference_allele' %in% colnames(mut) | any(is.na(mut$reference_allele)) ){
  
  gr_mut = GRanges(seqnames = mut$chr, ranges = IRanges(start = mut$start, end = mut$end), ALT = mut$mutated_to_allele,
                   sample_id = mut$sample_id, cancer_type = mut$cancer_type)
  
  gr_mut$REF = Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19, gr_mut, as.character = TRUE)

  mut = data.frame(sample_id=gr_mut$sample_id, cancer_type=gr_mut$cancer_type,
                        chr=seqnames(gr_mut), start=start(gr_mut),end=end(gr_mut), reference_allele = gr_mut$REF, mutated_to_allele = gr_mut$ALT)
  
  mut = mut[!mut$reference_allele == mut$mutated_to_allele,]
  }

mut_split = split(mut, mut$sample_id)
vcf_split = lapply(mut_split, function(x){
  gr_x <- GRanges(seqnames = x$chr, ranges = IRanges(start = x$start, end = x$end), REF = x$reference_allele, ALT = x$mutated_to_allele,
                     sample_id = x$sample_id, cancer_type = x$cancer_type)
  return(gr_x)
})

## Construct a mutation matrix from the loaded VCFs in comparison to the ref_genome
mut_mat <- data.frame(t(mut_matrix(vcf_list = vcf_split, ref_genome = ref_genome)))
mut_mat_norm = mut_mat/rowSums(mut_mat)

mut_mat_s = add_column(mut_mat_norm, sample_id = rownames(mut_mat_norm), .before = 1)
mut_mat_s = add_column(mut_mat_s, cancer_type = 'cantype', .after = 1)

cat = mut[!duplicated(mut$sample_id),]
cat$cancer_type = as.character(cat$cancer_type)
for (s in mut_mat_s$sample_id){
  mut_mat_s$cancer_type[mut_mat_s$sample_id == s] = cat$cancer_type[cat$sample_id == s]
}

write.csv(mut_mat_s, opt$output, row.names = FALSE)





