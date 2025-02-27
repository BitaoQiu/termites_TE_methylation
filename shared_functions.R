# shared functions for TEs.
library(GenomicRanges)
library(rtracklayer)
library(txdbmaker)
library(GenomicFeatures)
library(stringr)
library(data.table)

# Load TE annotations.
load(file.path(work_dir,'data/Rdata/termite_TE_annotation.Rdata'))
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
    print(table(anno$type.lv1))
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

# Summarize the TE methylation level by examining the overlap of CpG sites and TEs.
summarise_repeat_CpG = function(cpg.data, repeat.data, target.var = 'repeat_id',species = NULL,
                                meth_cut = 80, cov_cut = 10, cov_up_cut = 100){
  # cpg.data: CpG site data.
  # repeat.data: TE data.
  # meth_cut: Probability threshold for a CpG site to be considered methylated.
  # cov_cut: Lower read coverage threhold.
  # cov_up_cut: upper read coverage threhold.
  # Filter CpG sites with low or extremely high coverage.
  cpg.data.filtered = cpg.data[which(cpg.data$cov > cov_cut & cpg.data$cov < cov_up_cut)]
  repeat.data.df = droplevels(data.frame(repeat.data))
  target.var.full = names(table(subset(repeat.data.df, select = target.var)))
  tmp_list = list()
  tmp_list$TEs = target.var.full
  n_max = length(tmp_list$TEs)
  tmp_list$UnM = rep(0,n_max )
  tmp_list$M = rep(0, n_max)
  for (i in 1:n_max){
    target_repeat = target.var.full[i]
    repeat.data.target = subset(repeat.data, eval(as.name(target.var)) == target_repeat)
    cpg_repeat.target = subsetByOverlaps(cpg.data.filtered,repeat.data.target)
    if (length(cpg_repeat.target) > 0){
      # Only CpG sites with > meth_cut are considered methylated.
      tmp_list$UnM[i] = sum(cpg_repeat.target$m.Prob < meth_cut)
      tmp_list$M[i] = sum(cpg_repeat.target$m.Prob >= meth_cut)
    }
    else{
      tmp_list$UnM[i] = NA
      tmp_list$M[i] = NA
    }
  }
  
  temp = data.frame(do.call(cbind,tmp_list))
  temp$UnM = as.numeric(temp$UnM)
  temp$M = as.numeric(temp$M)
  temp$UnM[is.na(temp$UnM)] = 0
  temp$M[is.na(temp$M)] = 0
  cpg_repeat.summary = temp[,c(2,3)]
  rownames(cpg_repeat.summary) = temp$TEs
  cpg_repeat.summary = rbind(cpg_repeat.summary,table(cpg.data.filtered$m.Prob > meth_cut) )
  cpg_repeat.summary$CpG = cpg_repeat.summary$M + cpg_repeat.summary$UnM
  cpg_repeat.summary$M.Prob = cpg_repeat.summary$M/cpg_repeat.summary$CpG
  tmp_nrow = dim(cpg_repeat.summary)[1]
  cpg_repeat.summary[[target.var]] = c(rownames(cpg_repeat.summary)[c(1: (tmp_nrow - 1))], "Genomic")
  target.lvl = c('class','class.lv2','class.lv3','class.full')
  cpg_repeat.summary[,target.lvl] = 
    repeat.data.df[match(cpg_repeat.summary[[target.var]],
                         repeat.data.df[[target.var]]), target.lvl]
  cpg_repeat.summary[tmp_nrow,target.lvl] = "Genomic"
  cpg_repeat.summary$M.log = cpg_repeat.summary$M.Prob/cpg_repeat.summary$M.Prob[tmp_nrow]
  cpg_repeat.summary$Species = species
  return(cpg_repeat.summary)
}

# Format the output of one-code-to-find-them-all
# (https://github.com/mptrsen/mobilome/tree/master/code/Onecodetofindthemall)

format_te_sum = function(te_sum, sp, annotatation = T, anno_program = 'class.MC'){
  te_sum = te_sum[-grep('^#',te_sum$V1),]
  te_sum = te_sum[-which(te_sum$V1 == ''),]
  te_sum = te_sum[order(te_sum$V2),]
  repeat_id = unique(sort(te_sum$V2))
  # Summarize TE abundance at the family level (repeat_id)
  te_sum.new  = data.frame(row.names = repeat_id, repeat_id = repeat_id)
  for (i in repeat_id){
    te_sum.new[i,'length'] = max(te_sum[which(te_sum$V2 == i),'V3'])
    te_sum.new[i,'fragments'] = sum(te_sum[which(te_sum$V2 == i),'V4'])
    te_sum.new[i,'number'] = sum(te_sum[which(te_sum$V2 == i),'V5'])
    te_sum.new[i,'size'] = sum(te_sum[which(te_sum$V2 == i),'V7'])
  }
  # Filter TEs with less than 4 copies.
  te_sum.new = te_sum.new[which(te_sum.new$number > 3),]
  if (annotatation == T){
    te_sum.new = annotate_te(te_sum.new,anno_program = anno_program)
  }
  te_sum.new$species = sp
  return(te_sum.new)
}

# Summarize the CpG methylation level at ?
summarize_te_CpG = function(x, anno_class = F){
  sum_vars = c("M",'UnM','Number','Size')
  summary_data = data.frame(matrix(ncol = length(sum_vars)))
  colnames(summary_data) = sum_vars
  summary_data[,sum_vars] = apply(x[,sum_vars],2, sum)
  summary_data$M.Prob_1 = median(x$M.Prob)
  summary_data$M.Prob_2 = summary_data$M/(summary_data$M + summary_data$UnM)
  summary_data$Length = median(x$Length)
  summary_data$GC = median(x$GC)
  summary_data$M.log = median(x$M.log)
  if (anno_class == T){
    summary_data$class = x$class[1]}
  return(summary_data)
}


