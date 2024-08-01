#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(require(optparse)) 
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(GenomicRanges))

option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character',
              help="ssm input file"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="RMD output file")
)

opt = parse_args(OptionParser(option_list=option_list))
print(opt)

opt$input = 'data/realWXS/ssm_realWXS_tcga_muts20_7523samples_mergecantypes.csv'

# our data
samples = read_csv(opt$input)
#genes = read_csv('~/Documents/AA_PROJECT_VC/data/raw_data/Census_drivers_with_CNAs.csv')

# load CNAs data and load it into a dataframe
all = lapply(list.files("~/Documents/corrections/NUEVO/CNV/", pattern = "CNASNPHg19", full.names = TRUE), function(f){
  
  df <- read_delim(f, "\t", escape_double = FALSE, trim_ws = TRUE)
  df$cancer_type = gsub(".*\\/\\/[0-9]*-([A-Z]*)-CNASNPHg19.txt", "\\1", f)
  return(df)
})

CNV = dplyr::bind_rows(all)
colnames(CNV)[c(1,2,3,4)] = c('sample_id', 'chr', 'start', 'end')
head(CNV)

# remove repeated datasets
CNV$sample_id = sub("-...-...-....-..$", "", CNV$sample_id)
CNV_nodups = CNV[! CNV$cancer_type %in% c('COAD', 'READ', 'GBM', 'LGG', 'KIPAN'),]

CNV_ids = CNV_nodups[ CNV_nodups$sample_id %in% unique(samples$sample_id),]
table(CNV_ids$cancer_type)
length(unique(CNV_ids$sample_id))

head(CNV_ids)

# add threshold values
CNV_ids$thres = 0
CNV_ids$thres[CNV_ids$Segment_Mean > 0.1] = 1
CNV_ids$thres[CNV_ids$Segment_Mean < -0.1] = -1
table(CNV_ids$thres)

# check CNA length
CNV_ids$seg_len = CNV_ids$end - CNV_ids$start
table(is.na(CNV_ids$seg_len))
CNV_ids = CNV_ids[!is.na(CNV_ids$seg_len),]
min(CNV_ids$seg_len)
max(CNV_ids$seg_len)

library(ggplot2)
# plot their distributions
ggplot(data=CNV_ids, aes(x=CNV_ids$seg_len, fill=factor(CNV_ids$thres))) + 
  geom_histogram(alpha=.6, colour = 'black') + theme_classic() #+ xlim(-5 , 100)

ggplot(data=CNV_ids, aes(x=CNV_ids$seg_len, fill=factor(CNV_ids$thres))) + 
  geom_density(alpha=.3, colour = 'black') + theme_classic() + xlim(-5 , 100000)


write.csv(CNV_ids, 'CNAs_allcantypes_perseg_thresholded.csv', row.names = FALSE)


table(CNV_ids$seg_len > 5000000)

CNV_ids = CNV_ids[CNV_ids$seg_len < 5000000,]

# map to RMD 1Mb windows
# we will map separately +1 and -1 and then check per window which is the predominant one
head(CNV_ids)
CNV_ids$chr = sub('^', 'chr', CNV_ids$chr)
table(CNV_ids$chr)
CNA_del = CNV_ids[CNV_ids$thres == -1,]
CNA_amp = CNV_ids[CNV_ids$thres == 1,]

#load the windows 
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))
windows_1Mb <- read_csv("~/Documents/AA_PROJECT_VC/data/raw_data/windows_1Mb_ranges_for_RMD.bed")

#make them genomic ranges
suppressPackageStartupMessages(require(GenomicRanges))
gr_windows <- GRanges(seqnames = windows_1Mb$chr, ranges = IRanges(start = windows_1Mb$Start, end = windows_1Mb$End, names = windows_1Mb$window), strand = '*')
gr_del <- GRanges(seqnames = CNA_del$chr, ranges = IRanges(start = CNA_del$start, end = CNA_del$end), sample_id = CNA_del$sample_id, cancer_type = CNA_del$cancer_type)
gr_amp <- GRanges(seqnames = CNA_amp$chr, ranges = IRanges(start = CNA_amp$start, end = CNA_amp$end), sample_id = CNA_amp$sample_id, cancer_type = CNA_amp$cancer_type)
head(gr_del)
head(gr_amp)

# 1) calculate RMD
# divide the dataframe by sample
mut_split_del = split(gr_del,gr_del$sample_id)
mut_split_amp = split(gr_amp,gr_amp$sample_id)

# calculate RMD per sample for loss
RMD_split_del = lapply(mut_split_del, function(x){

  count <- countOverlaps(gr_windows, x)
  count = data.frame('sample_id' = as.character(unique(x$sample_id)), 'cancer_type'= as.character(unique(x$cancer_type)), t(count))
  return(count)
  
})

# calculate RMD per sample for gain
RMD_split_amp = lapply(mut_split_amp, function(x){
  
  count <- countOverlaps(gr_windows, x)
  count = data.frame('sample_id' = as.character(unique(x$sample_id)), 'cancer_type'= as.character(unique(x$cancer_type)), t(count))
  return(count)
  
})

# merge the split
RMD_merge_del = dplyr::bind_rows(RMD_split_del)
RMD_merge_amp = dplyr::bind_rows(RMD_split_amp)

RMD_merge_amp[1:5,1:15]
RMD_merge_del[1:5,1:15]

rowSums(RMD_merge_del[,3:ncol(RMD_merge_del)])

RMD_merge_amp_backup = RMD_merge_amp
RMD_merge_del_backup = RMD_merge_del

# we make sure that the samples are in the same order
table(RMD_merge_amp$sample_id == RMD_merge_del$sample_id)

new = lapply(colnames(RMD_merge_amp), function(y){

  if (y %in% c('sample_id', 'cancer_type')){
    print(y)
    return(RMD_merge_amp[, y])
  } else {
    
    col = RMD_merge_amp[, y] - RMD_merge_del[,y]
    col[col < 0] = -1
    col[col > 0] = 1

    #df = data.frame(y = col)
    #colnames(df) = y
    return(col)
  }
})

RMD_CNAs = dplyr::bind_cols(new)
colnames(RMD_CNAs) = colnames(RMD_merge_amp)
head(RMD_CNAs)

# change the cancer types names
for (s in unique(RMD_CNAs$sample_id)){
  RMD_CNAs$cancer_type[RMD_CNAs$sample_id == s] <- samples$cancer_type[samples$sample_id == s]
}

RMD_CNAs[1:20, 1:20]

a = RMD_CNAs[duplicated(RMD_CNAs$sample_id),]
head(a)


RMD_CNAs_fil = RMD_CNAs[!duplicated(RMD_CNAs$sample_id),]
length(unique(RMD_CNAs$sample_id))
write.csv(RMD_CNAs_fil, 'results/CNAs/RMD_CNA_maxlen5Mb_7405samples.csv', row.names = FALSE)

#RMD_CNAs = read_csv('results/CNAs/RMD_CNA_7405samples.csv')

# remove the windows that contains driver genes + 1Mb window in each side 
# windows
gr_windows

# og TCGA and 64 amp del
og = read_csv('data/raw_data/og_exons_reduce.csv')
og_cnas = read_csv('data/raw_data/Census_drivers_with_CNAs.csv')
head(og)
head(og_cnas)
og_cnas1 = og_cnas[,c('Gene Symbol', 'Genome Location')]
og_cnas1 = tidyr::separate(og_cnas1, col = 'Genome Location', into = c('chr', 'pos'), sep = ':')
og_cnas1 = tidyr::separate(og_cnas1, col = 'pos', into = c('start', 'end'), sep = '-')
og_cnas1$chr = sub('^', 'chr', og_cnas1$chr)
colnames(og_cnas1)[1] = 'gene'
head(og_cnas1)
allog = rbind(og, og_cnas1)
allog = allog[!duplicated(allog$gene),]
allog$start = as.numeric(allog$start)
allog$end = as.numeric(allog$end)
allog = allog[!is.na(allog$start),]

gr_og = GRanges(seqnames = as.character(allog$chr), ranges = IRanges(start = allog$start, end = allog$end), gene = allog$gene)
gr_og = gr_og[! seqnames(gr_og) == 'chrX', ]
head(gr_og)
table(seqnames(gr_og))
count <- countOverlaps(gr_windows, gr_og)

dc = data.frame(count)
dc$windows = rownames(dc)
head(dc)
table(dc$count == 1)

# function to extract windows with count>=1 and their colindant ones
extract.with.context <- function(x, rows, after = 0, before = 0) {
  
  match.idx  <- which(rownames(x) %in% rows)
  span       <- seq(from = -before, to = after)
  extend.idx <- c(outer(match.idx, span, `+`))
  extend.idx <- Filter(function(i) i > 0 & i <= nrow(x), extend.idx)
  extend.idx <- sort(unique(extend.idx))
  
  return(x[extend.idx, , drop = FALSE])
}

windws = extract.with.context(dc, rownames(dc[dc$count>0,]), after = 1, before = 1)
windws[1:10,]

# windows that have OGs
wind_og = dc[dc$count >0,]

RMD_CNAs_fil = RMD_CNAs[!duplicated(RMD_CNAs$sample_id),]
# RMD passenger CNAs
RMDfil = RMD_CNAs_fil[, ! colnames(RMD_CNAs_fil) %in% windws$windows]
# write file
write.csv(RMDfil, 'results/CNAs/plotfinal/RMD_CNAs_passenger_1eachside_only64.csv', row.names = FALSE)


# RMD driver CNAs
RMDdriv = RMD_CNAs_fil[, colnames(RMD_CNAs_fil) %in% c('sample_id', 'cancer_type', wind_og$windows)]
head(RMDdriv)
# write file
write.csv(RMDdriv, 'results/CNAs/plotfinal/RMD_CNAs_drivers_only64.csv', row.names = FALSE)


# RMD driver CNA 64 Census
og_cnas = read_csv('data/raw_data/Census_drivers_with_CNAs.csv')
head(og_cnas)

og_cnas1 = og_cnas[,c('Gene Symbol', 'Genome Location')]
og_cnas1 = tidyr::separate(og_cnas1, col = 'Genome Location', into = c('chr', 'pos'), sep = ':')
og_cnas1 = tidyr::separate(og_cnas1, col = 'pos', into = c('start', 'end'), sep = '-')
og_cnas1$chr = sub('^', 'chr', og_cnas1$chr)
colnames(og_cnas1)[1] = 'gene'
head(og_cnas1)
og_cnas1$start = as.numeric(og_cnas1$start)
og_cnas1$end = as.numeric(og_cnas1$end)
og_cnas1 = og_cnas1[!is.na(og_cnas1$start),]

gr_og = GRanges(seqnames = as.character(og_cnas1$chr), ranges = IRanges(start = og_cnas1$start, end = og_cnas1$end), gene = og_cnas1$gene)
gr_og = gr_og[! seqnames(gr_og) == 'chrX', ]
head(gr_og)
table(seqnames(gr_og))
count <- countOverlaps(gr_windows, gr_og)

dc = data.frame(count)
dc$windows = rownames(dc)
head(dc)
table(dc$count == 1)
