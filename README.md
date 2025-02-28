# Evolutioanry arms races between transposable elements (TEs) and DNA methylation in termites
Key codes for analyzing the co-evolution between TE and TE methylation in termites.

**Outline:**

* Build a cross-species non-redundant repeat library: [shell script](scripts/build_non_redundant_TE_library.sh).
* Calling genomic CpG information from PacBio HiFi reads (with kinetic information): [shell script](scripts/call_5mC.sh).
* Summarize TE family abundance based on RepeatMasker->OneCodeToFindThemAll output: [R script](scripts/summarize_TEs.R).
* Annotate the age of TE families based on their phylogenetic distribution: [R script](scripts/annotate_TE_age.R).
* Summarize TE family CpG methylation levels: [R script](scripts/summarize_TE_CpG.R).
* Summarize TE-associated structure variants (SVs) and examine the association between TE activities and methylation levels: [R script](scripts/summarize_TE_SVs.R).

## Relevant Rdata:
* TE methylation levels in the termite genomes [repeat_5mC_summary_C90.Rdata](data/Rdata/repeat_5mC_summary_C90.Rdata).
* SVs in M.bellicosus populations and SV annotation [Mbel_SVs_v2.Rdata](data/Rdata/Mbel_SVs_v2.Rdata).
* TE familiy annotation [termite_TE_annotation.Rdata](data/Rdata/termite_TE_annotation.Rdata).
