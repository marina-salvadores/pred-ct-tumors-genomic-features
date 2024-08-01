#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse)) # don't say "Loading required package: optparse"
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf
# vignette: http://www.icesi.edu.co/CRAN/web/packages/optparse/vignettes/optparse.pdf

option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character',
              help="ssm input file"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="RMD output file"),
  make_option(c("-f", "--norm_by_factor"), action="store_true", default=FALSE,
              help="normalize by factor"),
  make_option(c("-c", "--cov_file"), action="store", default=NA, type='character',
              help="file containing the normalization factors for the ssm data type"),
  make_option(c("-t", "--cov_threshold"), action="store", default=NA, type='integer',
              help="threshold for removing windows lower than a determine bp coverage"),
  make_option(c("-r", "--norm_by_row"), action="store", default='none', type='character',
              help="normalize by row. Options: 'mean' or 'sum' ")
)

opt = parse_args(OptionParser(option_list=option_list))
print(opt)



# main point of program is here, do this whether or not "verbose" is set
if(!is.na(opt$i) & !is.na(opt$o)) {

  
    #load the windows 
    suppressPackageStartupMessages(library(readr))
    suppressPackageStartupMessages(library(dplyr))
    windows_1Mb <- read_csv("~/Documents/AA_PROJECT_VC/data/raw_data/windows_1Mb_ranges_for_RMD.bed")
    mut <-read_csv(opt$i)
    
    #make them genomic ranges
    suppressPackageStartupMessages(require(GenomicRanges))
    gr_windows <- GRanges(seqnames = windows_1Mb$chr, ranges = IRanges(start = windows_1Mb$Start, end = windows_1Mb$End, names = windows_1Mb$window), strand = '*')
    gr_muts <- GRanges(seqnames = mut$chr, ranges = IRanges(start = mut$start, end = mut$end), sample_id = mut$sample_id, cancer_type = mut$cancer_type)
    
    # 1) calculate RMD
    # divide the dataframe by sample
    mut_split = split(gr_muts,gr_muts$sample_id)
      
    # calculate RMD per sample
    RMD_split = lapply(mut_split, function(x){
        
        count <- countOverlaps(gr_windows, x)
        count = data.frame('sample_id' = as.character(unique(x$sample_id)), 'cancer_type'= as.character(unique(x$cancer_type)), t(count))
        return(count)
        
    })
      
    # merge the split
    RMD_merge = dplyr::bind_rows(RMD_split)
    
    
    # 2) Norm by factor (if required)
    if (opt$norm_by_factor == TRUE){
      
      wind = read_csv(opt$cov_file)
      
      wind$NF = 1000000/wind$cov_bp
      head(wind)
      
      hicov = wind[wind$cov_bp > opt$cov_threshold,]
      RMDhicov = RMD_merge[, c('sample_id', 'cancer_type', hicov$window)]


      RMDnf = lapply(colnames(RMDhicov), function(x){

        if (! x %in% c('sample_id', 'cancer_type')){
          nf = hicov$NF[hicov$window == x]
          y = as.data.frame(RMDhicov[, x]*nf)
          colnames(y) = x
          return(y)    
        }else{
          sol = as.data.frame(RMDhicov[, x])
          colnames(sol) = x
          return(sol)
        }
      })
      
      
      RMDNF = dplyr::bind_cols(RMDnf)
      RMD_merge = RMDNF
    }
      

    # 3) norm by row (if required)
    
    if (opt$norm_by_row == 'mean'){
    
      RMD_merge[,3:ncol(RMD_merge)] = RMD_merge[,3:ncol(RMD_merge)]/rowMeans(RMD_merge[,3:ncol(RMD_merge)])
      print('Normalization by row applied (rowMeans)')
      
    } else if (opt$norm_by_row == 'sum'){
      
      RMD_merge[,3:ncol(RMD_merge)] = RMD_merge[,3:ncol(RMD_merge)]/rowSums(RMD_merge[,3:ncol(RMD_merge)])
      print('Normalization by row applied (rowSums)')
      
    } else {
      print('Any normalization by row applied.')
    }
    
    
    # 4) save the file 
    write.csv(RMD_merge, opt$o, row.names = FALSE)
    

    
} else {
  print("Input and Output file parameters are mandatory") 
}










