# Summarize TE methylation levels by examining the overlap of CpG sites and TEs.
source('shared_functions.R')
# Take M.bellicosus as an example.
# Repeat data (only scaffold Mbel_ptg000001l)
mbel_te = read_TE('data/example/Mbel_repeat_example.out','Mbel')
head(mbel_te)

# CpG data (only scaffold Mbel_ptg000001l)
mbel.cpg = read_CpG('data/example/Mbel_5mC_example.bed','MbelQ')
head(mbel.cpg)

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
  # Add the genomic CpG information to cpg_repeat.summary. 
  # Note that for this example, only the scaffold CpG information.
  cpg_repeat.summary = rbind(cpg_repeat.summary,
                             table(cpg.data.filtered$m.Prob > meth_cut)) 
  # Calculate the total number of CpG.
  cpg_repeat.summary$CpG = cpg_repeat.summary$M + cpg_repeat.summary$UnM
  # Calculate the percentage of methylated CpG.
  cpg_repeat.summary$M.Prob = cpg_repeat.summary$M/cpg_repeat.summary$CpG
  tmp_nrow = dim(cpg_repeat.summary)[1]
  cpg_repeat.summary[[target.var]] = c(rownames(cpg_repeat.summary)[c(1: (tmp_nrow - 1))], "Genomic")
  target.lvl = c('class','class.lv2','class.lv3','class.full')
  cpg_repeat.summary[,target.lvl] = 
    repeat.data.df[match(cpg_repeat.summary[[target.var]],
                         repeat.data.df[[target.var]]), target.lvl]
  # Normalized the TE methylation level by the genomic (scaffold) background, 
  # so that it's comparable between species. 
  # > 1 indidates higher than the genomic background.
  cpg_repeat.summary$M.norm = cpg_repeat.summary$M.Prob/cpg_repeat.summary$M.Prob[tmp_nrow]
  cpg_repeat.summary$Species = species
  return(cpg_repeat.summary)
}

mbel.cpg_repeat.summary.TE = summarise_repeat_CpG(mbel.cpg, mbel_te, 
                                                  target.var = 'repeat_id', species = 'Mbel')
head(mbel.cpg_repeat.summary.TE)

# For all species.
load('data/Rdata/repeat_5mC_summary_C90.Rdata')
