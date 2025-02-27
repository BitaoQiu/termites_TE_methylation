# shared functions for TEs.
library(GenomicRanges)
library(rtracklayer)
library(txdbmaker)
library(GenomicFeatures)
library(stringr)
library(data.table)

# Load TE annotations.
load(file.path('data/Rdata/termite_TE_annotation.Rdata'))

# Format RepeatMasker annotation into Class, Superfamily, Famliy, etc.
format_te_level = function(anno, anno_program = 'class.RM'){
  anno[[anno_program]] = gsub(' ','',anno[[anno_program]])
  anno$type.lv1 = sapply(anno[[anno_program]], FUN = function(x) strsplit(x, '/')[[1]][1])
  anno$type.lv2 = sapply(anno[[anno_program]], FUN = function(x) strsplit(strsplit(x, '/')[[1]][2],'-')[[1]][1])
  anno$type.lv2[is.na(anno$type.lv2)] = paste(anno$type.lv1[is.na(anno$type.lv2)], "Unknown", sep = ".")
  anno$type.lv2[which(anno$type.lv1 == 'Unknown')] = 'Unknown'
  type.lv2_levels = unique(anno$type.lv2[order(anno$type.lv1, anno$type.lv2,decreasing = F)])
  anno$type.lv2 = factor(anno$type.lv2, levels = type.lv2_levels)
  anno$type.lv2_full = anno[[anno_program]]
  if (grepl("RM", anno_program)){
    anno$class = unlist(RM_map[anno$type.lv1])
    anno$class.lv2 = anno$type.lv1
    anno$class.lv3 = anno$type.lv2
    anno$class.full = anno$type.lv2_full
  }
  if (grepl("MC", anno_program)){
    anno$class = anno$type.lv1
    anno$class.lv2 = anno$type.lv2
    anno$class.lv3 = anno$type.lv2
    anno$class.full = anno$type.lv2_full
  }
  return(anno)
}

# Annotate TEs.
annotate_te = function(te_data, anno = te_classified, anno_program = 'class.RM'){
  anno = format_te_level(anno, anno_program = anno_program)
  te_data$class = anno$class[match(te_data$repeat_id,anno$repeat_id)]
  te_data$class.lv2 = anno$class.lv2[match(te_data$repeat_id,anno$repeat_id)]
  te_data$class.lv3 = anno$class.lv3[match(te_data$repeat_id,anno$repeat_id)]
  te_data$class.full = anno$class.full[match(te_data$repeat_id,anno$repeat_id)]
  te_data$TE_age = anno$TE_age[match(te_data$repeat_id,anno$repeat_id)]
  te_data$TE_length = anno$Length[match(te_data$repeat_id,anno$repeat_id)]
  
  return(te_data)
}

# Read repeat annotation from RepeatMasker2.
read_TE = function(x, sp, repeat_anno = te_classified, anno_program = 'class.RM',
                   filtered = T){
  anno_repeat = read.table(x,skip = 3,header = F,sep = '',fill = T)
  anno_repeat_gd = makeGRangesFromDataFrame(anno_repeat,
                                            keep.extra.columns = T,
                                            ignore.strand = T,
                                            start.field = 'V6',end.field = 'V7',
                                            seqnames.field = 'V5', na.rm=TRUE)
  names(anno_repeat_gd@elementMetadata)[c(1:4,7,8,12)] = 
    c("SW_score",'div','del','ins',"repeat_id",'repeat_class','RM_id')
  anno_repeat_gd = annotate_te(anno_repeat_gd,anno_program = anno_program)
  if (filtered == TRUE){
    anno_repeat_gd = anno_repeat_gd[which(!anno_repeat_gd$repeat_class %in%
                                            c("Low_complexity",'Simple_repeat')),]
  }
  anno_repeat_gd$species = sp
  anno_repeat_gd$TE_class = as.character(anno_repeat_gd$class)
  return(anno_repeat_gd)
}

# Read PacBio CpG 5mC information from the output of pb-CpG-tools
# https://github.com/PacificBiosciences/pb-CpG-tools/tree/main
read_CpG = function(x, species){
  CpG = read.table(x)
  colnames(CpG)[c(4:9)] = c('m.score','haplotype','cov','m.count','um.count','m.Prob')
  CpG$sp = species
  CpG_gd = makeGRangesFromDataFrame(CpG,keep.extra.columns = T,ignore.strand = T,
                                    start.field = 'V2',end.field = 'V3',
                                    seqnames.field = 'V1')
  return(CpG_gd)
}



