# Examine the association between CpG methyaltion level, TE abundance, and SVs.
source('shared_functions.R')
load('data/Rdata/Mbel_SVs_v2.Rdata')
load('data/Rdata//repeat_5mC_summary_C90.Rdata')

# SVs annotation based on the non-redundant TE library.
sv_RM = annotate_te(sv_RM)
head(sv_RM)

# Summarize TE activiies as the number of SVs that annotated in the TE family.
active_TEs = data.frame(table(sv_RM$repeat_id))
names(active_TEs)[1] = 'repeat_id'
te_classified$TE_activties = 0
te_classified$TE_activties = 
  active_TEs$Freq[match(te_classified$repeat_id,active_TEs$repeat_id)] 
te_classified$TE_activties[is.na(te_classified$TE_activties)] = 0

# Add methylation information.
te_classified[,c('M.Prob_Q','M.norm_Q')] = mbelQ.cpg_repeat.summary.TE[
  match(te_classified$repeat_id,mbelQ.cpg_repeat.summary.TE$repeat_id),c('M.Prob','M.norm')]

target_age = c('Mbel','Macrotermitinae','Termitidae',
               'GeoIsoptera','Older')
te_classified.filtered = te_classified[which(te_classified$TE_age %in% target_age),]
te_classified.filtered = droplevels(annotate_te(te_classified.filtered))

# TE spreading efficient is calculated as:
# SV copy numbers (TE_activties) normalised by TE abundance (Mbel_Number)
ggplot(te_classified.filtered,
       aes(x = (TE_activties + 1)/Mbel_Number, 
           y = M.norm_Q))+
  geom_point()+
  scale_x_log10()+
  geom_smooth(method = 'glm')+
  geom_hline(yintercept = 1, col = 'red')+
  facet_wrap(~TE_age, nrow =1)+
  scale_y_sqrt()
