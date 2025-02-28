# Summarize TE abundance at the family level by formating the output of one-code-to-find-them-all
# (https://github.com/mptrsen/mobilome/tree/master/code/Onecodetofindthemall)
source('shared_functions.R')
# Use Macrotermes bellicosus as an example.
mbel_te_data = read.csv('data/repeats/Mbel_TE_summary.tsv', sep = '\t',header = F) 

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
  te_sum.filtered = te_sum.new[which(te_sum.new$number > 3),]
  if (annotatation == T){
    te_sum.filtered = annotate_te(te_sum.filtered,anno_program = anno_program)
  }
  te_sum.filtered$species = sp
  return(te_sum.filtered)
}
mbel_te_sum = format_te_sum(mbel_te_data, sp = 'Mbel', anno_program = 'class.RM')
head(mbel_te_sum)
