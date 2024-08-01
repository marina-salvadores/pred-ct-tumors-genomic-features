#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse)) 

option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character',
              help="ssm input file"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="RMD output file"),
  make_option(c("-p", "--pathway"), action="store_true", default=FALSE,
              help="add pathways as features"),
  make_option(c("-s", "--hotspots"), action="store_true", default=FALSE,
              help="add hotspots as features")
)

opt = parse_args(OptionParser(option_list=option_list))
print(opt)


########## function 1 - calculate binary matrix of mutations in OG/TSG ##########
calculate_og <- function(gr_muts){
og = read_csv('~/Documents/AA_PROJECT_VC/data/raw_data/og_exons_reduce.csv')
# add tert promoter to the list
og = rbind(og, data.frame(gene = 'TERT_promoter', chr = 'chr5', start = 1293469, end = 1296248 ))
gr_og = GenomicRanges::GRanges(seqnames=og$chr,ranges = IRanges(start = og$start, end = og$end),name = og$gene)
names(gr_og) = gr_og$name

# it is already reduce the version that I load - if not uncoment this commands
#red = unlist(reduce(split(gr_og, elementMetadata(gr_og)$name)))
#red$gene = names(red)

mut_split = split(gr_muts,gr_muts$names)

# calculate OG counts per sample
OG_split = lapply(mut_split, function(x){

  count <- countOverlaps(gr_og, x)
  indx2 <- duplicated(names(count))|duplicated(names(count),fromLast=TRUE)
  count2 <- count[indx2]
  lista = split(count2, names(count2)) 
  lista_sum = lapply(lista, function(x){
    sum(x)
  })
  
  count_nodup = count[!indx2]
  res = dplyr::bind_cols(lista_sum)
  countfinal = cbind(data.frame(t(count_nodup)), data.frame(res))
  
  rowfinal = data.frame('sample_id' = unique(x$names), 'cancer_type'= unique(x$cancer_type), countfinal)
  return(rowfinal)
  
})

# bind rows
OG_merge = dplyr::bind_rows(OG_split)

# transform to binary
OG_merge_binary = lapply(OG_merge, function(x){
  if (is.integer(x)){
    newx = as.numeric(gsub("[1-9]+[0-9]*", 1, x))
    return(newx)
  }else{
    return(x)
  }
})

# bind cols
df_binary = dplyr::bind_cols(OG_merge_binary)
return(df_binary)
}
########## end function 1 ##########


########## function 2 - calculate pathways features ##########
calculate_pathways <- function(og_mat){
  
  pathways = read_csv('~/Documents/AA_PROJECT_VC/data/raw_data/pathways.csv')
  head(pathways)
  pathways = pathways[pathways$Pathway != 'Other',]
  p_split = split(pathways,pathways$Pathway)
  
  temp = og_mat
  temp[unique(pathways$Pathway)] <- NA
  
  for (y in p_split){
    p = unique(y$Pathway)
    gene = sub("-", "\\.",y$Gene)
    sub = temp[,gene]
    
    temp[rowSums(sub) == 0,p] <- 0
    temp[rowSums(sub) > 0,p] <- 1
    
  }
  
  return(temp)
  
}
########## end function 2 ##########


########## function 3 - calculate hotspots features ##########
calculate_hotspots <- function(gr_muts, og_mat){
  
  # OG 3d hotspots parsed by COSMIC to obtain genomic coordinates
  hotspots_positions_COSMIC = read_csv("../../data/raw_data/hotspots_COSMIC_mutin100.csv")
  gr_hotspot = GenomicRanges::GRanges(seqnames=hotspots_positions_COSMIC$chr,
                                      ranges = IRanges(start = hotspots_positions_COSMIC$start, end = (hotspots_positions_COSMIC$end)),
                                      name = hotspots_positions_COSMIC$mut_cds)
  names(gr_hotspot) = gr_hotspot$name
  
  # divide the dataframe by sample
  mut_split = split(gr_muts,gr_muts$names)
  
  # calculate hotspots mutations
  RMD_split = lapply(mut_split, function(x){

    count <- countOverlaps(gr_hotspot, x)
    rowfinal = data.frame('sample_id' = unique(x$names), 'cancer_type'= unique(x$cancer_type), t(count))
    return(rowfinal)
    
  })
  
  og_hotspots = dplyr::bind_rows(RMD_split)
  colSums(og_hotspots[,3:ncol(og_hotspots)])  
  
  if (all(og_hotspots$sample_id == og_mat$sample_id)){
    og_matrix2_hot = cbind(og_mat, og_hotspots[,3:ncol(og_hotspots)])
    return(og_matrix2_hot)
  }
  
  
}
########## end function 3 ##########






# run the data
# load the packages
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(GenomicRanges))

mut <-read_csv(opt$i)
gr_mutations = GenomicRanges::GRanges(seqnames=mut$chr, ranges = IRanges(start = mut$start, end = mut$end), names = mut$sample_id, cancer_type = mut$cancer_type)


# OG
og_matrix = calculate_og(gr_mutations)
  
if (opt$pathway){
  og_pathway = calculate_pathways(og_matrix)
}
  
if(exists('og_pathway')){
  og_matrix2 = og_pathway
} else{
  og_matrix2 = og_matrix
}
  
if (opt$hotspots){
  og_matrix3 = calculate_hotspots(gr_mutations, og_matrix2)
} else{
  og_matrix3 = og_matrix2
}

write.csv(og_matrix3, opt$output, row.names = FALSE)

