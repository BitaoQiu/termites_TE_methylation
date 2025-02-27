source('shared_functions.R') # This include the annotated TE information.
# Estimation of TE age based on their distribution on the termite phylogeny. ####
te_dir ='data/repeats//'
termite_species = c('Trin',"Mbel",'Mbel_S','Odon','Rgra','Csec','Znev','Mdar','Cpun')
termiteTEs_sum.list = list()
for (sp in termite_species){
  termiteTEs_sum.list[[sp]] = read.csv(file.path(te_dir, paste(sp, 'TE_summary.tsv',sep = '_')),
                                       sep = '\t',header = F)
  termiteTEs_sum.list[[sp]]  = format_te_sum(termiteTEs_sum.list[[sp]], sp = sp, anno_program = 'class.RM')
}

termite_TE.sum = rbindlist(termiteTEs_sum.list)
termite_TE.sum$species = factor(termite_TE.sum$species, 
                                levels = termite_species)
termite_TE.sum = annotate_te(termite_TE.sum,anno_program = 'class.RM')
termite_TE.sum.filtered = droplevels(termite_TE.sum)#[which(termite_TE.sum$class %in% c("CLASSI",'CLASSII')),])
target_level = 'repeat_id'
target_tes = names(table(termite_TE.sum.filtered[[target_level]]))
termite_TE.sum.matrix = data.frame(
  row.names = target_tes)

# Summarize TE abundance in a matirx.
for (sp in termite_species){
  for (te_family in rownames(termite_TE.sum.matrix)){
    repeat_ids = which(termite_TE.sum.filtered$species %in% sp & 
                         termite_TE.sum.filtered[[target_level]] %in% te_family)
    termite_TE.sum.matrix[te_family,paste(sp,'Number',sep = '_')] = sum(termite_TE.sum.filtered$number[repeat_ids])
    termite_TE.sum.matrix[te_family,paste(sp,'Size',sep = '_')] = sum(termite_TE.sum.filtered$size[repeat_ids])
  }
}
termite_TE.sum.matrix$repeat_id = rownames(termite_TE.sum.matrix)
# Read TE family length information.
termite_TE.sum.matrix =  annotate_te(termite_TE.sum.matrix,anno_program = 'class.RM')
termite_TE.sum.matrix$TE_age = NULL

# Estimate TE family age based on their distribution.
estimate_TE_age = function(x, cp_threshold = 5,
                           sp = termite_species){
  # cp_threshold: Effective copy number threshold, lower than which is considered absence in the genome.
  target_names = paste(sp,'Size', sep = '_')
  te_length = as.numeric(x['TE_length'])
  size_threshold = cp_threshold*te_length 
  # TE size threshold, calculated as effective copy number threshold x TE length.
  x = as.numeric(x[target_names])
  trin = 1
  mbel = c(2,3)
  odon = 4
  rgra = 5
  csec = 6
  znev = 7
  mdar = 8
  cpun = 9
  te_age = "Unknown"
  if (max(x[mbel]) >= size_threshold & 
      max(x[c(odon,trin,rgra,csec,znev)]) < size_threshold){
    te_age = "Mbel"
  } else if (x[odon] >= size_threshold & 
             max(x[c(mbel,trin,rgra,csec,znev)]) < size_threshold){
    te_age = "Odon"
  } else if (max(x[mbel]) >= size_threshold & x[odon] >= size_threshold & 
             max(x[c(trin,rgra,csec,znev)]) < size_threshold){
    te_age = "Macrotermitinae"
  } else if (x[trin] >= size_threshold & 
             max(x[c(mbel,odon,rgra,csec,znev)]) < size_threshold){
    te_age = "Trin"
  } else if (max(x[mbel]) >= size_threshold & min(x[c(odon,trin)]) >= size_threshold 
             & max(x[c(rgra,csec,znev)]) < size_threshold){
    te_age = "Termitidae"
  } else if (x[rgra] >= size_threshold & 
             max(x[c(mbel,odon,trin,csec,znev)]) < size_threshold){
    te_age = "Rgra"
  } else if (max(x[mbel]) >= size_threshold & 
             min(x[c(odon,trin, rgra)]) >= size_threshold & 
             max(x[c(csec,znev)]) < size_threshold){
    te_age = "GeoIsoptera"
  } else if (max(x[mbel]) >= size_threshold & 
             (min(x[c(odon,trin,rgra,csec)]) >= size_threshold| 
              min(x[c(odon,trin,rgra,znev)]) >= size_threshold)) {
    te_age = "Older"
  } else if (x[csec] >= size_threshold & 
             max(x[c(mbel,odon,trin,rgra,znev,mdar)]) < size_threshold){
    te_age = "Csec"
  } else if (x[znev] >= size_threshold & 
             max(x[c(mbel,odon,trin,rgra,csec,mdar)]) < size_threshold){
    te_age = "Znev"
  } else if (x[mdar] >= size_threshold & 
             max(x[c(mbel,odon,trin,rgra,csec,znev)]) < size_threshold){
    te_age = "Mdar"
  } else if (x[cpun] >= size_threshold & 
             max(x[c(mbel,odon,trin,rgra,csec,znev,mdar)]) < size_threshold){
    te_age = "Cpun"
  }
  return(te_age)
}

termite_TE.sum.matrix$TE_age = apply(
  termite_TE.sum.matrix,1,estimate_TE_age)

head(termite_TE.sum.matrix)
table(termite_TE.sum.matrix$TE_age)



