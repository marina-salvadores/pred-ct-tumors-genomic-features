#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(require(optparse)) 
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))

option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character',
              help="ssm input file"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="RMD output file")
)

opt = parse_args(OptionParser(option_list=option_list))
print(opt)

# our data
samples = read_csv(opt$input)
genes = read_csv('~/Documents/AA_PROJECT_VC/data/raw_data/Census_drivers_with_CNAs.csv')

# load CNAs data and load it into a dataframe
counter = 1
for (f in list.files("~/Documents/corrections/NUEVO/CNV/", pattern = "all_thresholded.by_genes", full.names = TRUE)){
  
  cantype = gsub(".*\\/\\/[0-9]*-([A-Z]*)-all_thresholded.by_genes.txt", "\\1", f)

  df <- read_delim(f, "\t", escape_double = FALSE, trim_ws = TRUE)
  df_red1 = df[df$`Gene Symbol` %in% genes$`Gene Symbol`,]
  
  df_red2 = data.frame(t(df_red1[, 4:ncol(df_red1)]))
  colnames(df_red2) = df_red1$`Gene Symbol`
  
  df_red3 = add_column(df_red2, sample_id = rownames(df_red2), .before = 1)
  df_red3 = add_column(df_red3, cancer_type = cantype, .after = 'sample_id')
  
  
  if (counter == 1){
    CNV = df_red3
    counter = 2
    
  } else{
    CNV = rbind(CNV,df_red3)
  }
  
}



# remove repeated datasets
CNV$sample_id = sub("-...-...-....-..$", "", CNV$sample_id)
CNV_nodups = CNV[! CNV$cancer_type %in% c('COAD', 'READ', 'GBM', 'LGG', 'KIPAN'),]
CNV_nodups = CNV_nodups[!duplicated(CNV_nodups$sample_id),]

CNV_ids = CNV_nodups[ CNV_nodups$sample_id %in% unique(samples$sample_id),]
table(CNV_ids$cancer_type)

# change the cancer types names
for (s in unique(CNV_ids$sample_id)){
  CNV_ids$cancer_type[CNV_ids$sample_id == s] <- samples$cancer_type[samples$sample_id == s]
}

# change categorical variables into loss, neutral and gain
# final_cna = CNV_ids %>% mutate_all(as.character)
# final_cna[final_cna == '0'] <-'neutral'
# final_cna[final_cna == '-1'] <-'loss'
# final_cna[final_cna == '-2'] <-'loss'
# final_cna[final_cna == '1'] <-'gain'
# final_cna[final_cna == '2'] <-'gain'

# change categorical variables into loss, neutral and gain (here -1, 0, +1)
# because sklearn machine learning algorithms does not understand categorical variables
final_cna = CNV_ids
final_cna[final_cna == '-2'] <- -1
final_cna[final_cna == '2'] <- 1

# write file
write.csv(final_cna, opt$output, row.names = FALSE)


